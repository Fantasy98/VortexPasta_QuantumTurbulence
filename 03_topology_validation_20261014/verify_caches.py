from pathlib import Path
import sys,json,numpy as np
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'reference_tools'))
from diagnostics import read_psi,topology
G=Path('/home/yuninw/codes/simulations/ghost/run');out={}
for name,case in [('hopf','04_linked_rings_64_v2v3_compare'),('trefoil','06_trefoil_64_v3')]:
 d=G/case/'gpe';r,loops=topology(read_psi(d/'inputs/zre.0001.out'),2);z=np.load(d/'analysis/loops_0.000.npz');cached=[z[k] for k in z.files]
 assert len(cached)==len(loops)
 errors=[float(np.max(abs(a-b))) for a,b in zip(cached,loops)]
 assert max(errors)<1e-12,(name,errors)
 out[name]={'maximum_coordinate_errors':errors,'oriented_curves_identical':True}
(ROOT/'analysis/cache_verification.json').write_text(json.dumps(out,indent=2));print(out)
