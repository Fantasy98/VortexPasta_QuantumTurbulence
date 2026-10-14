from analyze import *
from polygon_gauss import gauss_exact,triangle
from scipy.integrate import dblquad
# Independent quadrature of a pair of skew segments checks normalization and sign.
a=np.array([0.,0.,0.]);u=np.array([1.,.2,.1]);b=np.array([.3,.4,.7]);v=np.array([-.2,1.,.3])
f=lambda t,s:np.dot(np.cross(u,v),a+s*u-b-t*v)/np.linalg.norm(a+s*u-b-t*v)**3/(4*np.pi)
expected=dblquad(f,0,1,lambda s:0,lambda s:1,epsabs=1e-12)[0]
r00=a-b;r10=a+u-b;r11=a+u-b-v;r01=a-b-v
measured=-(triangle(r00,r10,r11)+triangle(r00,r11,r01))/(4*np.pi)
assert abs(measured-expected)<1e-12,(measured,expected)
output={'skew_segment_control':{'quadrature':expected,'solid_angle':float(measured)},'endpoints':{}}
for name in ['hopf_base','hopf_fine','trefoil_base','trefoil_fine']:
 d=ROOT/'runs'/name;c,a,sn=loadrun(d);curves=readcurves(d/sn[-1]['file']);values=[]
 for n in [384,768,1536]:
  qs=[resample(p,n) for p in curves];wr=sum(gauss_exact(q,q,True) for q in qs);lk=sum(gauss_exact(qs[i],qs[j]) for i in range(len(qs)) for j in range(i+1,len(qs)));values.append({'samples_per_loop':n,'writhe':wr,'linking':lk,'H_centerline':wr+2*lk})
 output['endpoints'][name]={'t':float(sn[-1]['t']),'H_native':float(a['H_over_Gamma2'][-1]),'exact_polygon_values':values}
(ROOT/'analysis/writhe_crosscheck.json').write_text(json.dumps(output,indent=2));print(json.dumps(output,indent=2))
