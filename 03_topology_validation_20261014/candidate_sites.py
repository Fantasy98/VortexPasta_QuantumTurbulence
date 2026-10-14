from analyze import *
from scipy.optimize import linear_sum_assignment
p=ROOT/'analysis/trefoil_base/geometry.json'
if p.exists():
 rows=json.loads(p.read_text());ev=list(csv.DictReader((ROOT/'runs/trefoil_base/events.csv').open()));end=float(ev[0]['t_after']);tracked=[];last=None
 for row in rows:
  if row['t']>end:break
  sites=row['strand_sites']
  if len(sites)!=3:continue
  pos=np.array([s['location'] for s in sites])
  if last is not None:
   cost=np.linalg.norm((last[:,None]-pos[None,:]+L/2)%L-L/2,axis=-1);i,j=linear_sum_assignment(cost);sites=[sites[k] for k in j];pos=pos[j]
  tracked.append({'t':row['t'],'sites':sites});last=pos
 (ROOT/'analysis/trefoil_candidate_tracks.json').write_text(json.dumps(tracked,indent=2))
 fig,ax=plt.subplots(figsize=(7,4),constrained_layout=True)
 for k in range(3):ax.plot([r['t'] for r in tracked],[r['sites'][k]['distance'] for r in tracked],label=f'Geometric candidate {k+1}')
 ax.axhline(.08,color='k',ls='--',label='node-distance threshold')
 for e in ev[:3]:ax.axvline(float(e['t_after']),color='gray',alpha=.5)
 ax.set(xlabel='Physical time t',ylabel='Nonlocal segment distance',title='Three trefoil approach regions; spatial continuity assignment');ax.legend();ax.grid(alpha=.2);fig.savefig(ROOT/'figures/trefoil_candidate_sites.png',dpi=160);plt.close(fig)
