"""Check periodic display edges against synthetic crossings and retained curves."""
from pathlib import Path
import sys
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from analyze import BOX, CASES, curves, load_case, periodic_segments

def check(p):
    segments = periodic_segments(p)
    assert np.all(segments >= 0) and np.all(segments <= BOX)
    edges = (np.roll(p,-1,axis=0)-p+BOX/2)%BOX-BOX/2
    length = np.linalg.norm(edges,axis=1).sum()
    rendered_length = np.linalg.norm(segments[:,1]-segments[:,0],axis=1).sum()
    assert np.isclose(length,rendered_length,rtol=1e-12,atol=1e-12)
    assert np.max(np.linalg.norm(segments[:,1]-segments[:,0],axis=1)) <= np.max(np.linalg.norm(edges,axis=1))+1e-12
    # The same physical loop in another periodic image gives the same pieces.
    shifted = periodic_segments(p+BOX*np.array([2,-3,4]))
    assert shifted.shape == segments.shape
    assert np.allclose(shifted,segments,rtol=0,atol=1e-11)
    return float(abs(length-rendered_length))

def main():
    # Individual faces, negative coordinates, simultaneous corner crossings,
    # and endpoints exactly on a face.
    controls = [np.array([[1,1,1],[2,1,1],[1,2,1]],float)]
    for axis in range(3):
        for start,end in ((BOX-.1,BOX+.1),(-.1,.1),(0.,-.2),(BOX,BOX+.2)):
            p=np.array([[1,1,1],[1,2,1],[1,1,2]],float)
            p[:,axis]=[start,end,start]
            controls.append(p)
    controls.append(np.array([[BOX-.2]*3,[BOX+.2]*3,[BOX+.1,BOX+.3,BOX+.2]]))
    for p in controls:
        check(p)
    report={'synthetic_controls':len(controls),'cases':{}}
    for case in CASES:
        cfg,_,snaps,_,_=load_case(case)
        assert cfg['mesh']==[64,64,64] and cfg['L']==BOX
        errors=[check(p) for snap in snaps for p in curves(snap['path'])]
        report['cases'][case]={'snapshots':len(snaps),'loops_checked':len(errors),
            'max_length_error':max(errors),'mesh':cfg['mesh'],
            'plot_axis_limits':[[0,BOX]]*3}
    target=Path(__file__).with_name('wrapping_checks.json')
    with target.open('x') as f:
        json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    main()
