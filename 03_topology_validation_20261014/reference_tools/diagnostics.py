"""Independent Fourier diagnostics and oriented periodic plaquette vortex tracing.

All simulations are 64^3. Optional Fourier interpolation only samples the saved
64^3 trigonometric polynomial more finely for topology checks; it adds no modes.
"""
from collections import defaultdict
from pathlib import Path
from case_paths import resolve_case_path
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.optimize import root
from scipy.signal import resample
from scipy.spatial.distance import cdist

N = 64
L = 2*np.pi
XI = 7.5/N
ALPHA = 2*XI/np.sqrt(2)
BETA = 2/(XI*np.sqrt(2))
GAMMA = 4*np.pi*ALPHA
K = np.asarray(np.meshgrid(*(np.fft.fftfreq(N)*N,)*3, indexing='ij'))
K2 = np.sum(K*K, axis=0)
MASK = K2 <= (N/3)**2


def read_psi(real_path):
    real_path = resolve_case_path(real_path)
    imag_path = real_path.with_name(real_path.name.replace('zre.', 'zim.').replace('phi_re.', 'phi_im.'))
    return (np.fromfile(real_path, dtype='float64').reshape((N,)*3, order='F')
            + 1j*np.fromfile(imag_path, dtype='float64').reshape((N,)*3, order='F'))


def project(field):
    return np.fft.ifftn(np.fft.fftn(field)*MASK)


def physical(psi, acoustic=True):
    f = np.fft.fftn(psi)
    grad = np.asarray([np.fft.ifftn(1j*k*f) for k in K])
    rho = abs(psi)**2
    rho_p = project(rho).real
    ekq = float(2*ALPHA**2*np.mean(np.sum(abs(grad)**2, axis=0)))
    ep = float(ALPHA*BETA*np.mean((rho_p-1)**2))
    result = dict(mass=float(rho.mean()), kinetic_quantum=ekq, potential=ep,
                  energy=ekq+ep, rho_min=float(rho.min()), rho_max=float(rho.max()),
                  high_mode_fraction=float(np.sum(abs(f[~MASK])**2)/np.sum(abs(f)**2)))
    if acoustic:
        # sqrt(rho) v = 2 alpha Im(psi* grad psi)/|psi|; Helmholtz split.
        u = 2*ALPHA*np.imag(psi.conj()*grad)/np.maximum(abs(psi),1e-14)
        uf = np.fft.fftn(u, axes=(1,2,3))
        ku = np.sum(K*uf,axis=0)/np.where(K2==0,1,K2)
        uc = K*ku
        result['compressible_kinetic'] = float(.5*np.sum(abs(uc)**2)/N**6)
        result['incompressible_kinetic'] = float(.5*np.sum(abs(uf-uc)**2)/N**6)
    return result


def difference(a,b):
    angle = np.angle(np.vdot(a,b))
    return dict(relative_l2=float(np.linalg.norm(a-b)/np.linalg.norm(a)),
                phase_aligned_l2=float(np.linalg.norm(a-b*np.exp(-1j*angle))/np.linalg.norm(a)),
                phase_offset=float(angle),
                density_relative_l2=float(np.linalg.norm(abs(a)**2-abs(b)**2)/np.linalg.norm(abs(a)**2)))


def gauss_link(a,b):
    da=np.diff(a,axis=0); db=np.diff(b,axis=0)
    ma=(a[:-1]+a[1:])/2; mb=(b[:-1]+b[1:])/2
    r=ma[:,None]-mb[None,:]
    return float(np.sum(np.sum(np.cross(da[:,None],db[None,:])*r,axis=-1)
                        /np.linalg.norm(r,axis=-1)**3)/(4*np.pi))


def line_regularized_helicity(psi, loops, interpolation=2):
    """Independent Eq. (8)-(9) quadrature on resolved zero-line segments.

    Fourier derivatives are interpolated before evaluation. This is a secondary
    check with interpolation error, not a replacement for the native diagnostic.
    Returns H/Gamma^2. Unresolved/ambiguous loops must be excluded by the caller.
    """
    if not loops:return 0.0
    mid=np.concatenate([(p[:-1]+p[1:])/2 for p in loops])
    dl=np.concatenate([np.diff(p,axis=0) for p in loops])
    coords=(mid%L).T*(N*interpolation)/L
    def sample(f):
        field=np.fft.ifftn(f)
        if interpolation!=1:
            for axis in range(3):field=resample(field,N*interpolation,axis=axis)
        return (map_coordinates(field.real,coords,order=3,mode='grid-wrap')+
                1j*map_coordinates(field.imag,coords,order=3,mode='grid-wrap'))
    f=np.fft.fftn(psi)
    grad=np.array([sample(1j*k*f) for k in K]).T
    a=grad.real;b=grad.imag
    e=np.cross(a,b);e/=np.maximum(np.linalg.norm(e,axis=1)[:,None],1e-20)
    q=np.zeros_like(a)
    for j in range(3):
        for k in range(j,3):
            h=sample(-K[j]*K[k]*f)
            q[:,j]+=a[:,k]*h.imag-b[:,k]*h.real
            if k!=j:q[:,k]+=a[:,j]*h.imag-b[:,j]*h.real
    speed=2*ALPHA*np.sum(q*e,axis=1)/np.maximum(np.sum(abs(grad)**2,axis=1),1e-20)
    return float(np.sum(speed*np.sum(e*dl,axis=1))/GAMMA)


def topology(psi, interpolation=1):
    if interpolation != 1:
        for axis in range(3):
            psi=resample(psi, N*interpolation, axis=axis)
    n=psi.shape[0]; dx=L/n
    # A small fixed global phase changes no winding and avoids an exact zero
    # on a branch cut. No density threshold enters the plaquette count.
    phi=np.angle(psi*np.exp(.123j))
    wrap=lambda x: (x+np.pi)%(2*np.pi)-np.pi
    points=[]; signs=[]; normals=[]; cells=defaultdict(list)
    root_failures=0
    for normal,(u,v) in enumerate(((1,2),(2,0),(0,1))):
        p10=np.roll(phi,-1,axis=u);p01=np.roll(phi,-1,axis=v)
        p11=np.roll(p10,-1,axis=v)
        winding=np.rint((wrap(p10-phi)+wrap(p11-p10)+wrap(p01-p11)+wrap(phi-p01))/(2*np.pi)).astype(int)
        for ijk in np.argwhere(winding):
            i10=ijk.copy();i10[u]=(i10[u]+1)%n
            i01=ijk.copy();i01[v]=(i01[v]+1)%n
            i11=i10.copy();i11[v]=(i11[v]+1)%n
            q00=psi[tuple(ijk)];q10=psi[tuple(i10)];q01=psi[tuple(i01)];q11=psi[tuple(i11)]
            def value(t):
                a,b=t; z=q00*(1-a)*(1-b)+q10*a*(1-b)+q01*(1-a)*b+q11*a*b
                return (z.real,z.imag)
            sol=root(value,[.5,.5])
            valid=sol.success and np.max(np.abs(value(sol.x)))<1e-8 and np.all(sol.x>=-1e-6) and np.all(sol.x<=1+1e-6)
            if not valid:
                # A bilinear complex face can have a second root outside the
                # face. Locate its interior root before declaring a failure.
                for guess in ((.1,.1),(.1,.9),(.9,.1),(.9,.9)):
                    sol=root(value,guess)
                    valid=np.max(np.abs(value(sol.x)))<1e-8 and np.all(sol.x>=-1e-6) and np.all(sol.x<=1+1e-6)
                    if valid:break
            root_failures+=int(not valid)
            xy=sol.x if valid else np.array([.5,.5])
            xyz=ijk.astype(float);xyz[u]+=xy[0];xyz[v]+=xy[1]
            face=len(points);points.append(xyz*dx); signs.append(int(winding[tuple(ijk)]));normals.append(normal)
            cells[tuple(ijk)].append(face)
            below=ijk.copy();below[normal]=(below[normal]-1)%n
            cells[tuple(below)].append(face)
    graph=defaultdict(list);bad_cells=0
    for ids in cells.values():
        if len(ids)!=2:
            bad_cells+=1
            continue
        a,b=ids;graph[a].append(b);graph[b].append(a)
    points=np.asarray(points)
    visited=set();loops=[];open_components=0
    for seed in range(len(points)):
        if seed in visited:continue
        order=[];current=seed;previous=None
        while current not in visited:
            visited.add(current);order.append(current)
            neighbors=graph[current]
            if len(neighbors)!=2:break
            nxt=neighbors[0] if neighbors[0]!=previous else neighbors[1]
            previous,current=current,nxt
        if current!=seed or len(order)<3:
            open_components+=1;continue
        ids=np.array(order+[seed]);p=points[ids]
        delta=np.diff(p,axis=0);delta-=L*np.rint(delta/L)
        unwrapped=np.vstack([p[0],p[0]+np.cumsum(delta,axis=0)])
        # Signed plaquette normal fixes each loop's circulation orientation.
        tang=np.roll(delta,1,axis=0)+delta
        score=sum(tang[j,normals[idx]]*signs[idx] for j,idx in enumerate(order))
        if score<0:unwrapped=unwrapped[::-1]
        loops.append(unwrapped)
    loops.sort(key=lambda x:-np.sum(np.linalg.norm(np.diff(x,axis=0),axis=1)))
    report=dict(analysis_samples=n, winding_faces=len(points), ambiguous_cells=bad_cells,
                root_failures=root_failures, open_components=open_components, closed_loops=len(loops),
                lengths=[float(np.sum(np.linalg.norm(np.diff(p,axis=0),axis=1))) for p in loops],
                centroids=[np.average((p[:-1]+p[1:])/2,axis=0,
                                     weights=np.linalg.norm(np.diff(p,axis=0),axis=1)).tolist() for p in loops],
                oriented_area_vectors=[(.5*np.sum(np.cross(p[:-1],p[1:]),axis=0)).tolist() for p in loops],
                periodic_windings=[np.rint((p[-1]-p[0])/L).astype(int).tolist() for p in loops])
    if len(loops)==2 and bad_cells==0 and open_components==0:
        a,b=loops
        b=b+L*np.rint((a[:-1].mean(axis=0)-b[:-1].mean(axis=0))/L)
        report['gauss_link']=gauss_link(a,b)
        # Minimum point distance and local oriented tangent angle, refined by
        # optional interpolation. This is a sampled separation, not an exact minimum.
        dist=cdist(a[:-1],b[:-1]);ia,ib=np.unravel_index(np.argmin(dist),dist.shape)
        da=a[(ia+1)%(len(a)-1)]-a[(ia-1)%(len(a)-1)]
        db=b[(ib+1)%(len(b)-1)]-b[(ib-1)%(len(b)-1)]
        report['min_sampled_separation']=float(dist[ia,ib])
        report['nearest_tangent_cosine']=float(np.dot(da,db)/np.linalg.norm(da)/np.linalg.norm(db))
    return report,loops
