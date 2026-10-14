"""Read-only run analysis; writes derived products only under this experiment."""
from pathlib import Path
import sys,json,csv,tomllib,argparse
import numpy as np
from scipy.spatial.distance import cdist
from scipy.interpolate import CubicSpline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,FFMpegWriter,PillowWriter
ROOT=Path(__file__).resolve().parent
G=Path('/home/yuninw/codes/simulations/ghost/run')
sys.path.insert(0,str(ROOT/'reference_tools'))
from trefoil_diagnostics import knot_check
L=2*np.pi;GAMMA=4*np.pi*2*.1171875/np.sqrt(2)

def readcurves(p):
 a=np.loadtxt(p,delimiter=',',skiprows=1,ndmin=2)
 return [unwrap(a[a[:,0]==k,1:]) for k in np.unique(a[:,0])]
def unwrap(p):
 q=np.vstack([p[0],p[0]+np.cumsum((np.diff(p,axis=0)+L/2)%L-L/2,axis=0)])
 return q

def resample(p,n=384):
 p=np.vstack([p,p[0]]);ds=np.linalg.norm(np.diff(p,axis=0),axis=1);s=np.r_[0,np.cumsum(ds)]
 return CubicSpline(s,p,bc_type='periodic')(np.linspace(0,s[-1],n,endpoint=False))
def gauss(a,b,self_=False):
 da=np.roll(a,-1,axis=0)-a;db=np.roll(b,-1,axis=0)-b
 ma=a+da/2;mb=b+db/2;r=ma[:,None]-mb[None,:];dist=np.linalg.norm(r,axis=-1)
 if self_:np.fill_diagonal(dist,np.inf)
 return float(np.sum(np.sum(np.cross(da[:,None],db[None,:])*r,axis=-1)/dist**3)/(4*np.pi))

def separation(curves):
 # Exact straight-segment distances; spline representation sampled at 2x nodes.
 p=np.concatenate(curves);d=np.concatenate([np.roll(a,-1,axis=0)-a for a in curves]);ids=np.concatenate([np.full(len(a),i) for i,a in enumerate(curves)])
 arc=np.concatenate([np.r_[0,np.cumsum(np.linalg.norm(np.roll(a,-1,axis=0)-a,axis=1))[:-1]] for a in curves])
 lens=np.array([np.linalg.norm(np.roll(a,-1,axis=0)-a,axis=1).sum() for a in curves])
 i,j=np.triu_indices(len(p),1);diff=abs(arc[i]-arc[j]);same=ids[i]==ids[j]
 nonlocal_=~same | (np.minimum(diff,lens[ids[i]]-diff)>lens[ids[i]]/6)
 # A separate strict nonadjacent diagnostic excludes two neighboring sample segments.
 local_index=np.concatenate([np.arange(len(a)) for a in curves]);ns=np.array([len(a) for a in curves]);di=abs(local_index[i]-local_index[j]);allowed=~same | (np.minimum(di,ns[ids[i]]-di)>2)
 i=i[allowed];j=j[allowed];nonlocal_=nonlocal_[allowed]
 w=(p[i]-p[j]+L/2)%L-L/2;u=d[i];v=d[j]
 aa=np.sum(u*u,axis=1);bb=np.sum(u*v,axis=1);cc=np.sum(v*v,axis=1);dd=np.sum(u*w,axis=1);ee=np.sum(v*w,axis=1)
 denom=aa*cc-bb*bb;safe=np.where(abs(denom)>1e-20,denom,np.inf)
 s=(bb*ee-cc*dd)/safe;t=(aa*ee-bb*dd)/safe
 distances=np.sum((w+s[:,None]*u-t[:,None]*v)**2,axis=1);distances[(s<0)|(s>1)|(t<0)|(t>1)|(abs(denom)<1e-20)]=np.inf
 # Four constrained minima on the parameter-square boundary.
 for sv in (0.,1.):
  tv=np.clip((ee+sv*bb)/cc,0,1);dis=np.sum((w+sv*u-tv[:,None]*v)**2,axis=1);take=dis<distances;distances[take]=dis[take];s[take]=sv;t[take]=tv[take]
 for tv in (0.,1.):
  sv=np.clip((tv*bb-dd)/aa,0,1);dis=np.sum((w+sv[:,None]*u-tv*v)**2,axis=1);take=dis<distances;distances[take]=dis[take];s[take]=sv[take];t[take]=tv
 minimum=float(np.sqrt(distances.min()))
 distances[~nonlocal_]=np.inf
 order=np.argsort(distances);sites=[]
 for k in order:
  if not np.isfinite(distances[k]):break
  # p[j] in same nearest image as p[i].
  a=p[i[k]]+s[k]*u[k];b=p[i[k]]-w[k]+t[k]*v[k];center=((a+b)/2)%L
  if all(np.linalg.norm((center-np.array(x['location'])+L/2)%L-L/2)>.35 for x in sites):
   sites.append({'distance':float(np.sqrt(distances[k])),'location':center.tolist(),'loops':[int(ids[i[k]]),int(ids[j[k]])]})
  if len(sites)==3:break
 return minimum,sites

def modes(p,nodes=96):
 p=resample(p,512);c=p.mean(axis=0);q=p-c
 imp=np.sum(np.cross(p,np.roll(p,-1,axis=0)),axis=0);imp/=np.linalg.norm(imp)
 v=np.linalg.eigh(q.T@q)[1];normal=v[:,0]
 if np.dot(normal,imp)<0:normal=-normal
 e1=v[:,2];e2=np.cross(normal,e1)
 x=q@e1;y=q@e2;z=q@normal;ang=np.unwrap(np.arctan2(y,x))
 if np.median(np.diff(ang))<0:ang=-ang
 out={'center':c.tolist(),'best_plane_rms':float(np.sqrt(np.mean(z*z))),'mean_radius':float(np.mean(np.hypot(x,y)))}
 if np.any(np.diff(ang)<=0) or (ang[-1]-ang[0])<6:out['valid']=False;return out
 theta=np.linspace(ang[0],ang[0]+2*np.pi,512,endpoint=False);ra=np.hypot(x,y)
 ax=np.r_[ang[-1]-2*np.pi,ang,ang[0]+2*np.pi]
 rf=np.fft.rfft(np.interp(theta,ax,np.r_[ra[-1],ra,ra[0]]))/512
 zf=np.fft.rfft(np.interp(theta,ax,np.r_[z[-1],z,z[0]]))/512
 out.update(valid=True,maximum_resolved_mode=min(12,int(nodes)//8),modes=[{'m':m,'radial':float(2*abs(rf[m])),'axial':float(2*abs(zf[m])),'phase':float(np.angle(rf[m]*zf[m].conjugate()))} for m in range(2,min(12,int(nodes)//8)+1)])
 return out

def loadrun(d):
 cfg=tomllib.loads((d/'config.toml').read_text());a=np.genfromtxt(d/'diagnostics.csv',delimiter=',',names=True);a=np.atleast_1d(a)
 snaps=list(csv.DictReader((d/'snapshots.csv').open()));return cfg,a,snaps

def analyze(d):
 if not (d/'complete.toml').exists():return
 cfg,a,snaps=loadrun(d);rows=[];prev=None
 for k,snap in enumerate(snaps):
  curves=readcurves(d/snap['file']);qs=[resample(p) for p in curves]
  # Curves are continuous lifts. Align compact loops to the nearest common cell.
  for j in range(1,len(qs)):qs[j]+=L*np.rint((qs[0].mean(axis=0)-qs[j].mean(axis=0))/L)
  wr=[gauss(q,q,True) for q in qs];lk=[gauss(qs[i],qs[j]) for i in range(len(qs)) for j in range(i+1,len(qs))]
  sep,sites=separation(curves)
  windings=[np.rint(np.sum((np.roll(p,-1,axis=0)-p+L/2)%L-L/2,axis=0)/L).astype(int).tolist() for p in curves]
  assert all(w==[0,0,0] for w in windings), 'Noncontractible loop: Euclidean Gauss diagnostic is not applicable'
  r={'periodic_windings':windings,'t':float(snap['t']),'tstar':float(snap['tstar']),'loops':len(curves),'writhe':sum(wr),'writhes':wr,'pair_linking':lk,'linking':sum(lk),'H_centerline':sum(wr)+2*sum(lk),'minimum_nonadjacent_segment_distance':sep,'strand_sites':sites,'modes':[modes(q,nodes=len(curves[j])//2) for j,q in enumerate(qs)]}
  if k in {0,len(snaps)-1} or prev!=len(curves):
   r['knot_checks']=[knot_check(np.vstack([q,q[0]])) for q in qs]
  if k in {0,len(snaps)-1}:
   qq=[resample(p,768) for p in curves]
   for j in range(1,len(qq)):qq[j]+=L*np.rint((qq[0].mean(axis=0)-qq[j].mean(axis=0))/L)
   r['H_centerline_768']=sum(gauss(q,q,True) for q in qq)+2*sum(gauss(qq[i],qq[j]) for i in range(len(qq)) for j in range(i+1,len(qq)))
  rows.append(r);prev=len(curves)
 out=ROOT/'analysis'/d.name;out.mkdir(exist_ok=True)
 (out/'geometry.json').write_text(json.dumps(rows,indent=2))
 events=list(csv.DictReader((d/'events.csv').open()))
 summary={'configuration':cfg,'completion':tomllib.loads((d/'complete.toml').read_text()),'initial':{n:float(a[n][0]) for n in a.dtype.names},'final':{n:float(a[n][-1]) for n in a.dtype.names},'maximum_nodes':float(max(a['nodes'])),'min_segment':float(min(a['ds_min'])),'max_segment':float(max(a['ds_max'])),'max_curvature':float(max(a['kappa_max'])),'events':events,
 'max_relative_length_change':float(np.max(abs(a['length']/a['length'][0]-1))), 'max_relative_energy_change':float(np.max(abs(a['energy']/a['energy'][0]-1))),'max_impulse_vector_relative_change':float(np.max(np.sqrt(sum((a[n]-a[n][0])**2 for n in ['Ix','Iy','Iz'])))/np.sqrt(sum(a[n][0]**2 for n in ['Ix','Iy','Iz']))),
 'max_abs_helicity':float(max(abs(a['H_over_Gamma2']))),'max_H_geometry_difference':float(max(abs(np.interp(r['t'],a['t'],a['H_over_Gamma2'])-r['H_centerline']) for r in rows))}
 if d.name.startswith('ring'):
  summary['speed_fit']=float(np.polyfit(a['t'],a['cz'],1)[0]);summary['transverse_centroid_drift']=float(max(np.hypot(a['cx']-a['cx'][0],a['cy']-a['cy'][0])))
  summary['radius_range']=[min(r['modes'][0]['mean_radius'] for r in rows),max(r['modes'][0]['mean_radius'] for r in rows)];summary['max_best_plane_rms']=max(r['modes'][0]['best_plane_rms'] for r in rows)
 (out/'summary.json').write_text(json.dumps(summary,indent=2))
 print(d.name,'t=',a['t'][-1],'rec=',a['reconnections'][-1],'loops=',a['loops'][-1],'energychange=',summary['max_relative_energy_change'],flush=True)
 return summary

def plot(d,movie=False):
 cfg,a,snaps=loadrun(d);rows=json.loads((ROOT/'analysis'/d.name/'geometry.json').read_text());out=ROOT/'figures'/d.name;out.mkdir(exist_ok=True)
 fig,axes=plt.subplots(3,3,figsize=(13,10),constrained_layout=True);axs=axes.flat
 plots=[('length','Length'),('energy','Energy / (ρV)'),('H_over_Gamma2','H / Γ² (periodic velocity)'),('loops','Loop count'),('nodes','Node count'),('kappa_max','Maximum curvature'),('reconnection_length_loss','Length removed by reconnection'),('refinement_length_loss','Length removed by refinement')]
 for ax,(key,title) in zip(axs,plots):ax.plot(a['tstar'],a[key]);ax.set(xlabel='t* = Γt/R²',ylabel=title);ax.grid(alpha=.25)
 ax=axes[2,2];ax.plot(a['tstar'],a['ds_min'],label='min');ax.plot(a['tstar'],a['ds_max'],label='max');ax.legend();ax.set(xlabel='t*',ylabel='Segment length')
 fig.suptitle(d.name+' · a='+str(cfg['a'])+' · no mutual friction');fig.savefig(out/'diagnostics.png',dpi=160);fig.savefig(out/'diagnostics.pdf');plt.close(fig)
 fig,axes=plt.subplots(1,3,figsize=(13,4),constrained_layout=True);t=[r['tstar'] for r in rows]
 for key,label in [('linking','Σ Lk'),('writhe','Σ Wr'),('H_centerline','Σ Wr + 2Σ Lk')]:axes[0].plot(t,[r[key] for r in rows],label=label)
 axes[0].legend();axes[0].set(ylabel='Geometric topology / helicity')
 axes[1].plot(t,[r['strand_sites'][0]['distance'] for r in rows]);axes[1].axhline(cfg.get('drec',.08),ls='--',color='k',label='drec');axes[1].legend();axes[1].set(ylabel='Nonlocal strand separation')
 for key in ['Ix','Iy','Iz']:axes[2].plot(a['tstar'],a[key],label=key)
 axes[2].legend();axes[2].set(ylabel='Impulse / (ρΓ)')
 for ax in axes:ax.set(xlabel='t* = Γt/R²');ax.grid(alpha=.25)
 fig.savefig(out/'topology.png',dpi=150);plt.close(fig)
 indices=np.unique(np.linspace(0,len(snaps)-1,4).astype(int));curv=[readcurves(d/s['file']) for s in snaps]
 # Display continuous lifts shifted by a common periodic translation to the domain center.
 laboratory=d.name.startswith('ring')
 centers=[np.full(3,np.pi) if laboratory else np.concatenate(p).mean(axis=0) for p in curv]
 allq=[np.concatenate(p)-centers[k]+np.pi for k,p in enumerate(curv)]
 rad=max(1.3,max(float(np.abs(q-np.pi).max()) for q in allq)*1.08)
 def draw(ax,k):
  for p in curv[k]:
   q=np.vstack([p,p[0]])-centers[k]+np.pi;ax.plot(*q.T,lw=1.5)
  ax.set(xlim=(np.pi-rad,np.pi+rad),ylim=(np.pi-rad,np.pi+rad),zlim=(np.pi-rad,np.pi+rad),xlabel='x' if laboratory else 'x (comoving)',ylabel='y',zlabel='z');ax.set_box_aspect((1,1,1));ax.view_init(25,-48)
  ax.set_title(f"t={float(snaps[k]['t']):.4f}; t*={float(snaps[k]['tstar']):.3f}\n{len(curv[k])} loop(s)")
 fig=plt.figure(figsize=(14,4));
 for j,k in enumerate(indices):draw(fig.add_subplot(1,len(indices),j+1,projection='3d'),k)
 fig.suptitle(d.name);fig.tight_layout();fig.savefig(out/'snapshots.png',dpi=170);fig.savefig(out/'snapshots.pdf');plt.close(fig)
 if movie:
  fig=plt.figure(figsize=(7,6));ax=fig.add_subplot(projection='3d');chosen=np.unique(np.linspace(0,len(snaps)-1,min(100,len(snaps))).astype(int))
  def update(k):ax.clear();draw(ax,k);fig.suptitle(d.name+(' · laboratory frame' if laboratory else ' · continuous curves, comoving view'))
  ani=FuncAnimation(fig,update,frames=chosen,interval=100);ani.save(out/'evolution.mp4',writer=FFMpegWriter(fps=10),dpi=100);ani.save(out/'evolution.gif',writer=PillowWriter(fps=10),dpi=80);plt.close(fig)
 return

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('runs',nargs='*');p.add_argument('--movie',action='store_true');p.add_argument('--plot',action='store_true');args=p.parse_args()
 for name in args.runs:
  d=ROOT/'runs'/name
  if not (ROOT/'analysis'/name/'summary.json').exists():analyze(d)
  if args.plot or args.movie:plot(d,args.movie)
