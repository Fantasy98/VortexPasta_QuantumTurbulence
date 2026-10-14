from analyze import *
from scipy.optimize import linear_sum_assignment
checks={}
# Controls use known analytic curves and signed Hopf-link integral, including reflection.
for name in ['ring_64','hopf_64','trefoil_128']:
 qs=readcurves(ROOT/'inputs'/f'{name}.csv');rec={}
 if name.startswith('hopf'):
  qq=[resample(p,768) for p in qs];rec['link']=gauss(*qq);assert abs(rec['link']+1)<2e-5
  qref=[q*np.array([-1,1,1]) for q in qq];rec['reflected_link']=gauss(*qref);assert abs(rec['reflected_link']-1)<2e-5
 else:
  c=knot_check(np.vstack([qs[0],qs[0][0]]));rec['knot']=c
  if name.startswith('trefoil'):assert c['trefoil']
  else:assert all(x['determinant']==1 for x in c['diagrams'])
 checks[name]=rec
# All completed runs have finite diagnostics, actual retained checkpoints and events matching counts.
checks['runs']={}
for d in (ROOT/'runs').iterdir():
 if not (d/'complete.toml').exists():continue
 c,a,sn=loadrun(d);done=tomllib.loads((d/'complete.toml').read_text());ev=list(csv.DictReader((d/'events.csv').open()))
 assert all(np.isfinite(a[n]).all() for n in a.dtype.names),d
 assert len(ev)==done['reconnections'],d
 assert (d/'checkpoint_final.vtkhdf').stat().st_size>0,d
 assert abs(float(sn[-1]['t'])-done['final_time'])<1e-10,d
 checks['runs'][d.name]={'finite':True,'event_count_consistent':True,'checkpoint_retained':True,'target_reached':done['final_time']>=c['end_time'],'final_time':done['final_time']}
# Quantify candidate site order with a spatial nearest assignment, rather than assuming event number corresponds.
if (ROOT/'runs/trefoil_base/events.csv').exists():
 base=list(csv.DictReader((ROOT/'runs/trefoil_base/events.csv').open()))[:3]
 pa=np.array([[float(e[k]) for k in ['x','y','z']] for e in base]);ass={}
 if len(base)==3:
  for name in ['trefoil_fine','trefoil_halfdt','trefoil_d060','trefoil_d100']:
   p=ROOT/'runs'/name/'events.csv'
   if not p.exists():continue
   es=list(csv.DictReader(p.open()))[:3]
   if len(es)!=3:continue
   pb=np.array([[float(e[k]) for k in ['x','y','z']] for e in es]);cost=np.linalg.norm((pa[:,None]-pb[None,:]+L/2)%L-L/2,axis=-1);ii,jj=linear_sum_assignment(cost)
   ass[name]=[{'base_event':int(i+1),'test_event':int(j+1),'distance':float(cost[i,j]),'time':float(es[j]['t_after'])} for i,j in zip(ii,jj)]
  checks['trefoil_site_matching']=ass
(ROOT/'analysis/validation_checks.json').write_text(json.dumps(checks,indent=2));print('Validated',len(checks['runs']),'completed runs, topology controls, finite data, checkpoints and event counters.')
