#!/usr/bin/env python3
"""Analyze one of the three cases, or all; reads only this folder's results.

Examples:
  python3 analyze.py single_ring
  python3 analyze.py hopf_link --no-video
  python3 analyze.py all

Periodic views use a fixed [0,2*pi]^3 cube with edges split at its faces.
Outputs go to results/CASE/analysis_periodic; use --output for a fresh folder.
No GPE data or paper curves are imported. ffmpeg is required for MP4 output.
"""
from pathlib import Path
import argparse
import csv
import json
import shutil
import tomllib
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from mpl_toolkits.mplot3d.art3d import Line3DCollection

ROOT = Path(__file__).resolve().parent
CASES = ('single_ring', 'hopf_link', 'trefoil')
BOX = 2*np.pi

def periodic_segments(p, box=BOX):
    """Split each closed-curve edge at box faces and wrap the pieces.

    Use the shortest periodic edge, as for unwrapping the saved filaments.
    Midpoints select the image for each piece, so an endpoint at box belongs
    to its outgoing piece and reappears at zero on the next piece. Merely
    taking coordinates modulo box would draw spurious chords across the cube.
    """
    pieces = []
    for a, end in zip(p, np.roll(p, -1, axis=0)):
        d = (end-a+box/2) % box-box/2
        b = a+d
        cuts = [0., 1.]
        for axis in range(3):
            if abs(d[axis]) < 1e-14:
                continue
            lo, hi = sorted((a[axis], b[axis]))
            for cell in range(int(np.floor(lo/box))+1, int(np.ceil(hi/box))):
                fraction = (cell*box-a[axis])/d[axis]
                if 1e-12 < fraction < 1-1e-12:
                    cuts.append(fraction)
        cuts = np.unique(cuts)
        for left, right in zip(cuts[:-1], cuts[1:]):
            if right-left < 1e-12:
                continue
            image = box*np.floor((a+(left+right)/2*d)/box)
            piece = np.stack((a+left*d, a+right*d))-image
            if np.any(piece < -1e-10) or np.any(piece > box+1e-10):
                raise ValueError('Periodic line piece lies outside the cell')
            pieces.append(np.clip(piece, 0, box))
    return np.asarray(pieces).reshape(-1, 2, 3)

def load_case(case):
    root = ROOT/'results'/case
    config = tomllib.loads((root/'config.toml').read_text())
    records, snapshots, events, completions = {}, {}, {}, []
    for segment in sorted((root/'segments').iterdir()):
        # A resumed segment supersedes a previous tail beyond its restart time.
        a = np.atleast_1d(np.genfromtxt(segment/'diagnostics.csv', delimiter=',', names=True))
        if not len(a):
            continue
        start = float(a['t'][0])
        records = {k:v for k,v in records.items() if float(v['t']) < start-1e-12}
        snapshots = {k:v for k,v in snapshots.items() if v['t'] < start-1e-12}
        events = {k:v for k,v in events.items() if float(v['t_after']) <= start+1e-12}
        for row in a:
            records[int(row['step'])] = row
        with (segment/'snapshots.csv').open() as f:
            for row in csv.DictReader(f):
                row.update(t=float(row['t']), tstar=float(row['tstar']), path=segment/row['file'])
                snapshots[int(row['step'])] = row
        with (segment/'events.csv').open() as f:
            for row in csv.DictReader(f):
                events[float(row['t_after'])] = row
        file = segment/'completion.toml'
        if file.exists():
            completion=tomllib.loads(file.read_text())
            completion['segment']=segment.name
            completion['Julia_threads']=tomllib.loads((segment/'config.toml').read_text())['Julia_threads']
            completions.append(completion)
    a = np.array([records[k] for k in sorted(records)], dtype=a.dtype)
    return config, a, [snapshots[k] for k in sorted(snapshots)], list(events.values()), completions

def curves(path):
    a = np.atleast_2d(np.loadtxt(path, delimiter=',', skiprows=1))
    result = []
    for j in np.unique(a[:,0]):
        p = a[a[:,0]==j, 2:5]
        edges = (np.diff(p,axis=0)+BOX/2)%BOX-BOX/2
        result.append(np.vstack([p[0], p[0]+np.cumsum(edges,axis=0)]))
    return result

def resample(p, n=256):
    q = np.vstack([p,p[0]])
    ds = np.linalg.norm(np.diff(q,axis=0),axis=1)
    if np.any(ds<=1e-14):
        raise ValueError('Degenerate curve segment')
    arc = np.r_[0,np.cumsum(ds)]
    return CubicSpline(arc,q,bc_type='periodic')(np.linspace(0,arc[-1],n,endpoint=False))

def gauss(a,b,same=False):
    da=np.roll(a,-1,axis=0)-a; db=np.roll(b,-1,axis=0)-b
    r=(a+da/2)[:,None]-(b+db/2)[None,:]
    distance=np.linalg.norm(r,axis=-1)
    if same:
        np.fill_diagonal(distance,np.inf)
    return float(np.sum(np.sum(np.cross(da[:,None],db[None,:])*r,axis=-1)/distance**3)/(4*np.pi))

def geometric_history(snaps, data):
    rows=[]
    for snap, ps in zip(snaps,data):
        windings=[np.rint(np.sum((np.roll(p,-1,axis=0)-p+BOX/2)%BOX-BOX/2,axis=0)/BOX).astype(int).tolist() for p in ps]
        row={'t':snap['t'],'tstar':snap['tstar'],'loops':len(ps),'periodic_winding':windings,
             'component_coordinate_spans':[np.ptp(p,axis=0).tolist() for p in ps]}
        if any(w!=[0,0,0] for w in windings):
            row.update(writhe=None,linking=None,H_geometry=None,
                       note='Euclidean closed-loop Gauss diagnostic omitted for nonzero periodic winding')
        else:
            ns=max(48,min(256,2048//max(1,len(ps))))
            qs=[resample(p,ns) for p in ps]
            # Compact contractible loops are placed in a common periodic lift.
            for q in qs[1:]:
                q+=BOX*np.rint((qs[0].mean(axis=0)-q.mean(axis=0))/BOX)
            wr=[gauss(q,q,True) for q in qs]
            lk=[gauss(qs[i],qs[j]) for i in range(len(qs)) for j in range(i+1,len(qs))]
            row.update(samples_per_loop=ns,writhe=sum(wr),linking=sum(lk),
                       H_geometry=sum(wr)+2*sum(lk),component_writhe=wr,pair_linking=lk)
        rows.append(row)
    return rows

def time_axis(ax, gamma):
    ax.set_xlabel('Physical time t')
    ax.grid(alpha=.25)
    top=ax.secondary_xaxis('top',functions=(lambda t:gamma*t,lambda ts:ts/gamma))
    top.set_xlabel('t* = Γt/R²',fontsize=8)
    top.tick_params(labelsize=7)

def savefig(fig,out,name):
    fig.savefig(out/(name+'.png'),dpi=160)
    fig.savefig(out/(name+'.pdf'))
    plt.close(fig)

def analyze(case, out=None, video=True, frame='periodic'):
    cfg,a,snaps,events,completions=load_case(case)
    if not snaps:
        raise ValueError('No snapshots')
    out=Path(out) if out else ROOT/'results'/case/('analysis_periodic' if frame=='periodic' else 'analysis_unwrapped')
    if out.exists():
        raise FileExistsError(f'Preserving {out}; select a new --output directory')
    out.mkdir(parents=True)
    gamma=cfg['Gamma'];t=a['t'];data=[curves(s['path']) for s in snaps]
    geo=geometric_history(snaps,data)
    (out/'geometry.json').write_text(json.dumps(geo,indent=2)+'\n')
    initial_h=geo[0]['H_geometry']
    if case=='hopf_link':
        assert abs(geo[0]['linking']+1)<.005, 'Initial Hopf link orientation/geometry failed'
    if case=='single_ring':
        assert abs(initial_h)<1e-7, 'Circle has nonzero geometric helicity'
    summary={'case':case,'config':cfg,'completion_records':completions,
             'visualization_frame':frame,
             'visualization_note':'Periodic views wrap and split curves at every box face; '
                 'the fixed cube is [0,2*pi]^3. Diagnostics use the original continuous coordinates.',
             'reached_t10':bool(abs(t[-1]-10)<1e-10),
             'finite_native_diagnostics':bool(all(np.all(np.isfinite(a[k])) for k in a.dtype.names)),
             'initial_H_over_Gamma2':float(a['H_over_Gamma2'][0]),
             'final_H_over_Gamma2':float(a['H_over_Gamma2'][-1]),
             'initial_geometric_H':initial_h,'final_loops':int(a['loops'][-1]),
             'maximum_nodes':int(max(a['nodes'])), 'events':events,
             'max_relative_energy_change':float(max(abs(a['energy']/a['energy'][0]-1))),
             'max_relative_length_change':float(max(abs(a['length']/a['length'][0]-1))),
             'diagnostic_notes':'Native H=Gamma integral(v_BS.ds), divided by Gamma^2. '
               'Geometric H=sum(Wr)+2sum(Lk) uses midpoint Gauss quadrature on resampled '
               'cubic curves in a common continuous lift; no ribbon twist or periodic-image '
               'Gauss terms. This is a compact-loop geometric diagnostic, not a fully periodic '
               'helicity calculation for extended curves. Geometric values near reconnection '
               'are resolution sensitive. Energy is E_total/(rho*V); impulse is I/(rho*Gamma). '
               'Core/mesh/time match the three base configurations; no numerical convergence '
               'or paper/GPE quantitative agreement is established by these three runs.'}
    fig,axes=plt.subplots(3,3,figsize=(14,11),constrained_layout=True)
    ax=axes[0,0];ax.plot(t,a['H_over_Gamma2'],label='Native VFM')
    ax.plot([r['t'] for r in geo],[np.nan if r['H_geometry'] is None else r['H_geometry'] for r in geo],ls='--',label='Linking + writhe (approx.)')
    ax.set_ylabel('H / Γ²');ax.legend(fontsize=7)
    for ax,key,label in [(axes[0,1],'length','Length / R'),(axes[0,2],'energy','E_total / (ρV)'),
                         (axes[1,1],'loops','Loop count'),(axes[1,2],'nodes','Node count'),
                         (axes[2,1],'kappa_max','Maximum spline curvature')]:
        ax.plot(t,a[key]);ax.set_ylabel(label)
    for key in ('Ix','Iy','Iz'):
        axes[1,0].plot(t,a[key],label=key)
    axes[1,0].set_ylabel('Impulse / (ρΓ)');axes[1,0].legend(fontsize=7)
    for key in ('ds_min','ds_mean','ds_max'):
        axes[2,0].plot(t,a[key],label=key)
    axes[2,0].set_ylabel('Segment lengths');axes[2,0].legend(fontsize=7)
    for key in ('reconnection_length_loss','refinement_length_loss','removed_length'):
        axes[2,2].plot(t,a[key],label=key.replace('_',' '))
    axes[2,2].set_ylabel('Cumulative numerical length changes');axes[2,2].legend(fontsize=6)
    for ax in axes.flat:
        time_axis(ax,gamma)
        for event in events:
            ax.axvline(float(event['t_after']),color='gray',alpha=.12,lw=.7)
    fig.suptitle(f'{case} · periodic Biot–Savart 64³ mesh · a={cfg["a"]} · no mutual friction')
    savefig(fig,out,'diagnostics')
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    if case=='single_ring':
        speed=float(np.polyfit(t,a['cz'],1)[0]);prediction=gamma/(4*np.pi*cfg['R'])*(np.log(8*cfg['R']/cfg['a'])-cfg['Delta'])
        summary.update(speed_z=speed,z_initial=float(a['cz'][0]),z_final=float(a['cz'][-1]),
            displacement_z=float(a['cz'][-1]-a['cz'][0]),unbounded_speed_prediction=prediction,
            max_abs_helicity=float(max(abs(a['H_over_Gamma2']))))
        axes[0].plot(t,a['cz']-a['cz'][0],label='Measured Δz')
        axes[0].plot(t,speed*t,'--',label=f'Fit U={speed:.6f}')
        axes[0].set_ylabel('Unwrapped displacement Δz / R');axes[0].legend()
        axes[1].plot(t,np.gradient(a['cz'],t),label='Measured periodic speed')
        axes[1].axhline(prediction,ls='--',color='gray',label='Unbounded ring prediction')
        axes[1].set_ylabel('Translation speed');axes[1].legend(fontsize=8)
    else:
        for key,label in [('linking','Σ Lk'),('writhe','Σ Wr'),('H_geometry','Σ Wr + 2Σ Lk')]:
            axes[0].plot([r['t'] for r in geo],[np.nan if r[key] is None else r[key] for r in geo],label=label)
        axes[0].legend();axes[0].set_ylabel('Geometric diagnostics (approximate)')
        axes[1].step(t,a['reconnections'],where='post',label='Reconnection count')
        axes[1].step(t,a['removed_loops'],where='post',label='Removed loops')
        axes[1].set_ylabel('Cumulative count');axes[1].legend()
    for ax in axes:time_axis(ax,gamma)
    fig.suptitle(case+' · case-specific analysis');savefig(fig,out,'case_analysis')
    periodic = frame=='periodic'
    if periodic:
        low=np.zeros(3);high=np.full(3,BOX);span=high-low
        display=[ [periodic_segments(p) for p in ps] for ps in data ]
        frame_label='periodic cube [0, 2π]³'
    else:
        points=np.concatenate([p for ps in data for p in ps])
        low=points.min(axis=0);high=points.max(axis=0)
        span=np.maximum(high-low,.5);low-=.08*span;high+=.08*span;span=high-low
        frame_label='unwrapped laboratory frame'
    def draw(ax,k):
        for j,p in enumerate(data[k]):
            if periodic:
                ax.add_collection3d(Line3DCollection(display[k][j],colors=f'C{j%10}',linewidths=1.5))
            else:
                end=p[-1]+(p[0]-p[-1]+BOX/2)%BOX-BOX/2
                q=np.vstack([p,end]);ax.plot(*q.T,lw=1.5)
        if case=='single_ring' and k:
            p=data[0][0];ax.plot(*np.vstack([p,p[0]]).T,color='gray',alpha=.4,ls='--')
        ax.set(xlim=(low[0],high[0]),ylim=(low[1],high[1]),zlim=(low[2],high[2]),
               xlabel='x',ylabel='y',zlabel='z' if periodic else 'z (unwrapped)')
        ax.set_box_aspect(span);ax.view_init(elev=20,azim=-55)
        if periodic:
            for setter in (ax.set_xticks,ax.set_yticks,ax.set_zticks):
                setter([0,np.pi,BOX],labels=['0','π','2π'])
            ax.tick_params(axis='both',labelsize=8,pad=1)
        elif case=='single_ring':
            # A long laboratory trajectory makes this box narrow; sparse horizontal
            # ticks preserve equal physical scales without overlapping labels.
            ax.set_xticks([np.pi],labels=['π']);ax.set_yticks([np.pi],labels=['π'])
            ax.tick_params(axis='both',labelsize=8,pad=1)
        ax.set_title(f't = {snaps[k]["t"]:.3f}   t* = {snaps[k]["tstar"]:.3f}\n{len(data[k])} loop(s)',fontsize=10)
    fig=plt.figure(figsize=(15,6))
    fig.subplots_adjust(left=.02,right=.96,bottom=.12,top=.85,wspace=.28)
    for j,k in enumerate(np.unique(np.linspace(0,len(snaps)-1,4).astype(int))):
        draw(fig.add_subplot(1,4,j+1,projection='3d'),k)
    fig.suptitle(case+' · '+frame_label+' · equal physical coordinate scales')
    savefig(fig,out,'snapshots')
    if video:
        if shutil.which('ffmpeg') is None:
            raise RuntimeError('ffmpeg is required for MP4 output; use --no-video for plots only')
        fig=plt.figure(figsize=(11,7))
        fig.subplots_adjust(left=.02,right=.97,bottom=.10,top=.86,wspace=.15)
        ax=fig.add_subplot(121,projection='3d');track=fig.add_subplot(122)
        key='cz' if case=='single_ring' else 'H_over_Gamma2'
        values=a[key]-a[key][0] if case=='single_ring' else a[key]
        track.plot(t,values,color='gray',alpha=.6)
        marker,=track.plot([],[],'o',color='tab:blue')
        track.set(xlabel='Physical time t',ylabel='Unwrapped Δz / R' if case=='single_ring' else 'H / Γ²',xlim=(0,cfg['end_time']))
        track.grid(alpha=.25)
        def update(k):
            ax.clear();draw(ax,k)
            tt=snaps[k]['t'];marker.set_data([tt],[np.interp(tt,t,values)])
            fig.suptitle(case+' · '+frame_label+' · physical timestamps')
            return marker,
        # Uniform physical-time targets make playback independent of event-output density.
        targets=np.linspace(snaps[0]['t'],snaps[-1]['t'],201)
        times=np.array([s['t'] for s in snaps])
        selected=[int(np.argmin(abs(times-tt))) for tt in targets]
        anim=FuncAnimation(fig,update,frames=selected,interval=100)
        anim.save(out/'evolution.mp4',writer=FFMpegWriter(fps=10,codec='libx264',extra_args=['-pix_fmt','yuv420p','-threads','1']),dpi=100)
        plt.close(fig)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(f'{case}: t={t[-1]:.6f}, loops={int(a["loops"][-1])}, reconnections={int(a["reconnections"][-1])}; {out}',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case',choices=(*CASES,'all'))
    parser.add_argument('--no-video',action='store_true')
    parser.add_argument('--frame',choices=('periodic','unwrapped'),default='periodic',
                        help='Default: fixed periodic cube; unwrapped follows the laboratory coordinates')
    parser.add_argument('--output',type=Path,help='Fresh directory for one case; not valid with all')
    args=parser.parse_args()
    if args.case=='all' and args.output:
        parser.error('--output is only supported for one case')
    for case in CASES if args.case=='all' else (args.case,):
        analyze(case,args.output,not args.no_video,args.frame)
