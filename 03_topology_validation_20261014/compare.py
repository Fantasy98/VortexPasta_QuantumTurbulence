from analyze import *

def comparison():
 output={}
 for case in ['ring','hopf','trefoil']:
  names=[case+s for s in ['_base','_fine','_halfdt','_d060','_d100']]
  available=[name for name in names if (ROOT/'runs'/name/'complete.toml').exists()]
  if not available:continue
  base=ROOT/'runs'/available[0];cfg,a,sn=loadrun(base);T=min(loadrun(ROOT/'runs'/x)[1]['t'][-1] for x in available)
  tx=np.linspace(0,T,301);entry={}
  fig,axes=plt.subplots(2,3,figsize=(13,8),constrained_layout=True)
  for name in available:
   d=ROOT/'runs'/name;c,b,ss=loadrun(d);ev=list(csv.DictReader((d/'events.csv').open()));earliest=min([float(e['t_after']) for e in ev]+[T]);pre=min(earliest*.9,T)
   record={'events':ev,'final_time':float(b['t'][-1]),'differences':{}}
   for key in ['length','energy','H_over_Gamma2','Ix','Iy','Iz']:
    delta=np.interp(tx,b['t'],b[key])-np.interp(tx,a['t'],a[key]);record['differences'][key]={'max_absolute':float(max(abs(delta))),'pre_event_max_absolute':float(max(abs(delta[tx<pre]))),'end_absolute':float(delta[-1])}
   # Unregistered distances use the periodic positions. Also report centroid differences separately.
   geo=[]
   for target in np.linspace(0,T,7):
    sa=min(sn,key=lambda x:abs(float(x['t'])-target));sb=min(ss,key=lambda x:abs(float(x['t'])-target))
    pa=np.concatenate([resample(p,256) for p in readcurves(base/sa['file'])]);pb=np.concatenate([resample(p,256) for p in readcurves(d/sb['file'])]);r=(pa[:,None]-pb[None,:]+L/2)%L-L/2;dist=np.linalg.norm(r,axis=-1)
    geo.append({'target_t':float(target),'t_base':float(sa['t']),'t_test':float(sb['t']),'rms_symmetric_distance':float(np.sqrt((np.mean(dist.min(axis=0)**2)+np.mean(dist.min(axis=1)**2))/2))})
   record['geometry_distances']=geo;entry[name]=record
   for ax,key in zip(axes.flat,['length','energy','H_over_Gamma2','loops','nodes','reconnection_length_loss']):ax.plot(b['tstar'],b[key],label=name);ax.set(xlabel='t* = Γt/R²',ylabel=key);ax.grid(alpha=.25)
  axes[0,0].legend(fontsize=8);fig.suptitle(case+' · fixed Γ, a, Δ; fixed drec except labeled threshold checks');fig.savefig(ROOT/'figures'/f'{case}_convergence.png',dpi=160);plt.close(fig);output[case]=entry
 (ROOT/'analysis/convergence.json').write_text(json.dumps(output,indent=2))

def ghost_comparison():
 output={}
 for case,folder in [('hopf','04_linked_rings_64_v2v3_compare'),('trefoil','06_trefoil_64_v3')]:
  g=G/folder;hist=np.loadtxt(g/'gpe/history/helicity.txt');validated=json.loads((g/'gpe/analysis/validated_snapshots.json').read_text());rows=[]
  for gp in sorted((g/'gpe/analysis').glob('loops_*.npz')):
   t=float(gp.stem.split('_')[1])
   if t>2.5:continue
   data=np.load(gp);qs=[resample(unwrap(data[k][:-1])) for k in data.files]
   for j in range(1,len(qs)):qs[j]+=L*np.rint((qs[0].mean(axis=0)-qs[j].mean(axis=0))/L)
   hc=sum(gauss(q,q,True) for q in qs)+2*sum(gauss(qs[i],qs[j]) for i in range(len(qs)) for j in range(i+1,len(qs)))
   rows.append({'t':t,'tstar':GAMMA*t,'loops':len(qs),'length':sum(np.linalg.norm(np.roll(q,-1,axis=0)-q,axis=1).sum() for q in qs),'H_geometry':hc})
  linefile=g/'analysis'/('line_helicity_check.json' if case=='hopf' else 'line_helicity.json');line=json.loads(linefile.read_text());hkey='H_line_over_Gamma2' if case=='hopf' else 'H_over_Gamma2'
  fig,axs=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
  axs[0,0].plot(hist[:,0]*GAMMA,hist[:,1]*L**3/GAMMA**2,color='gray',alpha=.5,label='GPE native regularized (not converged)')
  axs[0,0].plot([x['t']*GAMMA for x in line],[x[hkey] for x in line],'ko',ms=4,label='GPE Eq. 8 line check')
  for ax,key in [(axs[0,0],'H_geometry'),(axs[0,1],'loops'),(axs[1,0],'length')]:ax.plot([x['tstar'] for x in rows],[x[key] for x in rows],'k--',label='GPE traced geometry')
  compare=[]
  for name in [case+'_base',case+'_wide',case+'_prepared_ref',case+'_prepared_fine',case+'_prepared_halfdt',case+'_prepared_slender']:
   d=ROOT/'runs'/name
   if not (d/'complete.toml').exists():continue
   c,a,sn=loadrun(d)
   a=a[a['t'] <= (2.1 if case=='hopf' else .6)+.002]
   for ax,key in [(axs[0,0],'H_over_Gamma2'),(axs[0,1],'loops'),(axs[1,0],'length')]:ax.plot(a['tstar'],a[key],label=name)
   distances=[]
   for t in np.arange(0,min(a['t'][-1],2.5)+1e-6,.1):
    gp=g/'gpe/analysis'/f'loops_{t:.3f}.npz'
    if not gp.exists():continue
    sa=min(sn,key=lambda x:abs(float(x['t'])-t));pa=np.concatenate([resample(p,256) for p in readcurves(d/sa['file'])]);z=np.load(gp);pb=np.concatenate([resample(unwrap(z[k][:-1]),256) for k in z.files]);dist=np.linalg.norm((pa[:,None]-pb[None,:]+L/2)%L-L/2,axis=-1)
    distances.append({'gpe_t':float(t),'vfm_t':float(sa['t']),'rms':float(np.sqrt((np.mean(dist.min(axis=0)**2)+np.mean(dist.min(axis=1)**2))/2))})
   axs[1,1].plot([x['gpe_t']*GAMMA for x in distances],[x['rms'] for x in distances],label=name);compare.append({'run':name,'distance':distances})
  for ax,title in zip(axs.flat,['Helicity / Γ² (definitions labeled)','Loop count','Centerline length','RMS geometric distance to retained GPE']):ax.set(xlabel='t* = Γt/R²',ylabel=title,xlim=(0,2.1*GAMMA if case=='hopf' else .6*GAMMA));ax.grid(alpha=.25)
  axs[0,0].legend(fontsize=6);axs[1,1].legend(fontsize=6);fig.suptitle(case+' · common circulation and domain; analytic and prepared geometry distinguished');fig.savefig(ROOT/'figures'/f'{case}_ghost_comparison.png',dpi=170);plt.close(fig)
  output[case]={'gpe_traced':rows,'gpe_line_helicity':line,'comparisons':compare}
 (ROOT/'analysis/ghost_comparison.json').write_text(json.dumps(output,indent=2))

if __name__=='__main__':comparison();ghost_comparison()
