from analyze import *

def kelvin():
 results={};fig,axs=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
 for ax,case in zip(axs,['ring','hopf','trefoil']):
  p=ROOT/'analysis'/f'{case}_base/geometry.json'
  if not p.exists():continue
  rows=json.loads(p.read_text());saved=[]
  for target in ([0,1] if case=='ring' else ([.7,1.] if case=='hopf' else [.4,.6,.8])):
   r=min(rows,key=lambda x:abs(x['t']-target));valid=[x for x in r['modes'] if x.get('valid')]
   if not valid:continue
   m=max(valid,key=lambda x:x['mean_radius']);saved.append({'t':r['t'],'tstar':r['tstar'],'shape':m})
   ax.plot([x['m'] for x in m['modes']],[np.hypot(x['radial'],x['axial']) for x in m['modes']],'-o',ms=3,label=f"t={r['t']:.3f}")
  results[case]=saved;ax.set(title=case,xlabel='Azimuthal mode m',ylabel='Radial/axial amplitude',yscale='log');ax.legend(fontsize=8);ax.grid(alpha=.2)
 fig.suptitle('Centerline shape modes after removing translation and best-fit plane tilt; no dispersion/cascade inference');fig.savefig(ROOT/'figures/shape_modes.png',dpi=160);plt.close(fig)
 (ROOT/'analysis/shape_modes.json').write_text(json.dumps(results,indent=2))

def reconnection_sites():
 names=['trefoil_base','trefoil_fine','trefoil_halfdt','trefoil_d060','trefoil_d100']
 result={};fig=plt.figure(figsize=(12,4));ax=fig.add_subplot(121,projection='3d');bx=fig.add_subplot(122)
 colors=['tab:blue','tab:orange','tab:green','tab:red','tab:purple']
 for name,col in zip(names,colors):
  p=ROOT/'runs'/name/'events.csv'
  if not p.exists():continue
  ev=list(csv.DictReader(p.open()));result[name]=ev
  for i,e in enumerate(ev[:3]):
   pos=[float(e[k]) for k in ['x','y','z']];ax.scatter(*pos,color=col);ax.text(*pos,str(i+1),size=7)
  bx.plot(np.arange(1,len(ev[:3])+1),[float(e['t_after']) for e in ev[:3]],'o-',label=name,color=col)
 bx.set(xlabel='Event order within run',ylabel='Physical event time t');bx.legend(fontsize=8);bx.grid(alpha=.2);ax.set(xlabel='x',ylabel='y',zlabel='z');ax.set_title('Event order and location (not forced simultaneous)');fig.tight_layout();fig.savefig(ROOT/'figures/trefoil_event_order.png',dpi=160);plt.close(fig)
 (ROOT/'analysis/reconnection_sites.json').write_text(json.dumps(result,indent=2))

def budgets():
 result={}
 for d in (ROOT/'runs').glob('*operation_audit'):
  p=d/'operation_budget.csv'
  if not p.exists():continue
  a=np.genfromtxt(p,delimiter=',',names=True);records=[]
  diag=np.genfromtxt(d/'diagnostics.csv',delimiter=',',names=True)
  for step in np.unique(a['step']):
   q=a[a['step']==step]
   if len(q)!=3:continue
   records.append({'step':int(step),'t_step_start':float(q['t'][0]),'t_after_step':float(diag['t'][np.argmin(abs(diag['step']-(step+1)))]),'E_reconnection_change':float(q['energy'][1]-q['energy'][0]),'E_refinement_change':float(q['energy'][2]-q['energy'][1]),'L_reconnection_change':float(q['length'][1]-q['length'][0]),'L_refinement_change':float(q['length'][2]-q['length'][1]),'reconnections':int(q['reconnections'][1]-q['reconnections'][0])})
  result[d.name]={'records':records,'sum_E_reconnection_change':sum(x['E_reconnection_change'] for x in records),'sum_E_refinement_change':sum(x['E_refinement_change'] for x in records)}
 (ROOT/'analysis/operation_budgets.json').write_text(json.dumps(result,indent=2))
if __name__=='__main__':kelvin();reconnection_sites();budgets()
