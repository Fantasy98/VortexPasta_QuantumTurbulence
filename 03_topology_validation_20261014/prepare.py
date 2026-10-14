from pathlib import Path
import sys,json,hashlib,shutil,subprocess,platform
import numpy as np
from scipy.interpolate import CubicSpline
ROOT=Path(__file__).resolve().parent
G=Path('/home/yuninw/codes/simulations/ghost/run')
sys.path.insert(0,str(ROOT/'reference_tools'))
from diagnostics import gauss_link,topology,read_psi
from trefoil_diagnostics import knot_check

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def resample(p,n):
    p=np.asarray(p); p=p if np.linalg.norm(p[-1]-p[0])<1e-8 else np.vstack([p,p[0]])
    s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(p,axis=0),axis=1))]
    q=CubicSpline(s,p,bc_type='periodic')(np.linspace(0,s[-1],max(n*16,2048)+1))
    a=np.r_[0,np.cumsum(np.linalg.norm(np.diff(q,axis=0),axis=1))]
    return np.column_stack([np.interp(np.linspace(0,a[-1],n,endpoint=False),a,q[:,k]) for k in range(3)])
def closed(p):return np.vstack([p,p[0]])
s=np.linspace(0,2*np.pi,8193);o=np.zeros_like(s)
curves={'ring':[np.pi+np.array([np.cos(s),np.sin(s),o]).T],
'hopf':[np.pi+np.array([-.5+np.cos(s),np.sin(s),o]).T,np.pi+np.array([.5+np.cos(s),o,np.sin(s)]).T],
'trefoil':[np.pi+np.array([(np.sin(s)+2*np.sin(2*s))/3,(np.cos(s)-2*np.cos(2*s))/3,-np.sin(3*s)/3]).T]}
provenance={}
for name,case in [('hopf','04_linked_rings_64_v2v3_compare'),('trefoil','06_trefoil_64_v3')]:
    stage=G/case/'gpe'; source=stage/'analysis/loops_0.000.npz'
    a=np.load(source); loops=[a[k] for k in a.files];curves[name+'_prepared']=loops
    # Independently retrace retained preparation handoff, reusing existing method.
    r,check=topology(read_psi(stage/'inputs/zre.0001.out'),2)
    provenance[name]={'cache':str(source),'sha256':sha(source),'retrace':r,'cache_lengths':[float(np.linalg.norm(np.diff(p,axis=0),axis=1).sum()) for p in loops],
      'completion':json.loads((stage/'run_complete.json').read_text()),'input_hash_matches':[]}
    for p in ['zre','zim']:
        dst=stage/f'inputs/{p}.0001.out'; expected=next(x['sha256'] for x in provenance[name]['completion']['imported_checkpoint'] if f'/{p}.' in x['destination'])
        provenance[name]['input_hash_matches'].append(sha(dst)==expected)
    assert all(provenance[name]['input_hash_matches'])
    assert r['ambiguous_cells']==r['open_components']==r['root_failures']==0
    assert len(check)==len(loops)
    for n in (64,128) if name=='hopf' else (128,256):
        pass
checks={}
for name,loops in curves.items():
    for n in ((64,128) if name.startswith(('ring','hopf')) else (128,256)):
        qs=[resample(p,n) for p in loops]
        data=np.concatenate([np.column_stack([np.full(n,i+1),q]) for i,q in enumerate(qs)])
        np.savetxt(ROOT/f'inputs/{name}_{n}.csv',data,delimiter=',',header='loop,x,y,z',comments='')
        lengths=[np.linalg.norm(np.roll(q,-1,axis=0)-q,axis=1) for q in qs]
        r={'nodes':[len(q) for q in qs],'ds_min':min(float(x.min()) for x in lengths),'ds_max':max(float(x.max()) for x in lengths),'length':sum(float(x.sum()) for x in lengths),'knots':[knot_check(closed(q)) for q in qs]}
        if len(qs)==2:r['link']=gauss_link(*[closed(resample(p,512)) for p in loops]);assert abs(abs(r['link'])-1)<1e-4
        if name.startswith('trefoil'):assert r['knots'][0]['trefoil']
        if name=='ring':assert all(x['determinant']==1 for x in r['knots'][0]['diagrams'])
        checks[f'{name}_{n}']=r
(ROOT/'analysis/initial_checks.json').write_text(json.dumps(checks,indent=2))
(ROOT/'analysis/ghost_provenance.json').write_text(json.dumps(provenance,indent=2))
# Audit retained data and actual times, rather than inferring completion from a README.
audit={}
for d in sorted(G.iterdir()):
    records=[]
    for p in d.glob('**/run_complete.json'):
        j=json.loads(p.read_text());records.append({'path':str(p),'returncode':j.get('returncode'),'duration':j.get('duration'),'final_record':j.get('final_real_checkpoint')})
    audit[d.name]={'completion_records':records,'binary_count':len(list(d.glob('**/zre.*.out'))),'parameters':[str(p) for p in d.glob('**/parameter.inp')]}
(ROOT/'analysis/environment_audit.json').write_text(json.dumps(audit,indent=2))
print(json.dumps(checks,indent=2))
