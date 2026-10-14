from analyze import *
from itertools import product
results={}
for name in ['ring_base','hopf_base','trefoil_base','hopf_prepared_ref','trefoil_prepared_ref']:
 d=ROOT/'runs'/name
 if not (d/'complete.toml').exists():continue
 cfg,a,sn=loadrun(d);record=[]
 for snap in [sn[0],sn[-1]]:
  qs=[resample(q,192) for q in readcurves(d/snap['file'])];mx=0.;where=None
  for shift in product([-1,0,1],repeat=3):
   if shift==(0,0,0):continue
   for i,x in enumerate(qs):
    for j,y in enumerate(qs):
     val=gauss(x,y+L*np.array(shift))
     if abs(val)>mx:mx=abs(val);where=[i,j,list(shift)]
  record.append({'t':float(snap['t']),'max_absolute_neighbor_image_link':mx,'pair_shift':where,'interpretation':'Values much below one: no integer linking with nearest periodic copies resolved.'})
  assert mx<.01,(name,mx)
 results[name]=record
(ROOT/'analysis/periodic_image_linking.json').write_text(json.dumps(results,indent=2));print(results)
