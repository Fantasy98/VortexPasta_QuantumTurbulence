"""Knot diagrams, nonlocal separation and core widths for traced zero lines.

Coloring matrix minor: 2*over - under_in - under_out at every crossing.
Reference: https://math.uchicago.edu/~may/VIGRE/VIGRE2011/REUPapers/Dixon.pdf
A determinant alone does not uniquely identify a knot; a three-crossing diagram
with determinant 3 identifies the trefoil. All three views are retained.
"""
import numpy as np
from scipy.spatial.distance import cdist
from scipy.ndimage import map_coordinates
from diagnostics import L, N


def analytic_trefoil(n=1000, radius=1):
    s=np.linspace(0,2*np.pi,n+1)
    return np.pi+radius*np.array([(np.sin(s)+2*np.sin(2*s))/3,
        (np.cos(s)-2*np.cos(2*s))/3,-np.sin(3*s)/3]).T


def diagram(loop, direction):
    n=np.asarray(direction,dtype=float);n/=np.linalg.norm(n)
    u=np.cross(n,[1.,0.,0.]);u/=np.linalg.norm(u);v=np.cross(n,u)
    p=loop@np.array([u,v,n]).T
    a=p[:-1];d=np.diff(p,axis=0);count=len(d)
    cross=lambda x,y:x[...,0]*y[...,1]-x[...,1]*y[...,0]
    i,j=np.triu_indices(count,2)
    use=~((i==0)&(j==count-1));i=i[use];j=j[use]
    den=cross(d[i],d[j]);use=abs(den)>1e-12;i=i[use];j=j[use];den=den[use]
    delta=a[j]-a[i]
    t=cross(delta,d[j])/den;s=cross(delta,d[i])/den
    use=(t>1e-8)&(t<1-1e-8)&(s>1e-8)&(s<1-1e-8)
    i=i[use];j=j[use];t=t[use];s=s[use]
    heights=a[i,2]+t*d[i,2]-a[j,2]-s*d[j,2]
    if len(heights) and np.min(abs(heights))<1e-7:
        return {'crossings':len(i),'degenerate':True}
    over=np.where(heights>0,i+t,j+s);under=np.where(heights>0,j+s,i+t)
    m=len(i)
    if not m:return dict(crossings=0,determinant=1,writhe=0)
    order=np.argsort(under);under=under[order];over=over[order]
    matrix=np.zeros((m,m),dtype=int)
    for k in range(m):
        oa=(np.searchsorted(under,over[k])-1)%m
        matrix[k,oa]+=2;matrix[k,k]-=1;matrix[k,(k-1)%m]-=1
    from sympy import Matrix
    det=abs(int(Matrix(matrix[:-1,:-1]).det()))
    writhe=int(np.sum(np.sign(heights*den[use])))
    return dict(crossings=m,determinant=det,writhe=writhe)


def knot_check(loop):
    views=[diagram(loop,d) for d in ([.017,.023,1],[.23,.41,1],[.7,-.2,1])]
    valid=all(v.get('determinant')==3 for v in views)
    return dict(diagrams=views,determinant_three=bool(valid),
                trefoil=bool(valid and any(v['crossings']==3 for v in views)))


def strand_separation(loop):
    d=np.linalg.norm(np.diff(loop,axis=0),axis=1);s=np.r_[0,np.cumsum(d)[:-1]];length=d.sum()
    arc=abs(s[:,None]-s[None,:]);arc=np.minimum(arc,length-arc)
    dist=cdist(loop[:-1],loop[:-1]);dist[arc<length/6]=np.inf
    i,j=np.unravel_index(dist.argmin(),dist.shape)
    return dict(nonlocal_exclusion_arclength=float(length/6),
        nonlocal_sampled_separation=float(dist[i,j]),separation_grid_cells=float(dist[i,j]/(L/N)))


def core_widths(psi,loop):
    # Normal sections through 24 arclength-spaced traced zero-line positions.
    ds=np.linalg.norm(np.diff(loop,axis=0),axis=1);arc=np.r_[0,np.cumsum(ds)]
    ids=np.searchsorted(arc,np.linspace(0,arc[-1],24,endpoint=False))
    rad=np.linspace(-.6,.6,301);width=[]
    for k in ids:
        tangent=loop[(k+1)%(len(loop)-1)]-loop[(k-1)%(len(loop)-1)]
        tangent/=np.linalg.norm(tangent)
        normal=np.cross(tangent,np.eye(3)[np.argmin(abs(tangent))]);normal/=np.linalg.norm(normal)
        for v in (normal,np.cross(tangent,normal)):
            coords=((loop[k]+rad[:,None]*v)%L).T*N/L
            z=map_coordinates(psi.real,coords,order=3,mode='grid-wrap')+1j*map_coordinates(psi.imag,coords,order=3,mode='grid-wrap')
            rho=abs(z)**2
            sides=[]
            for side in (slice(150,None),slice(150,None,-1)):
                rs=abs(rad[side]);vs=rho[side];hit=np.flatnonzero(vs>=.5)
                if len(hit) and hit[0]>0:
                    q=hit[0];sides.append(float(np.interp(.5,vs[q-1:q+1],rs[q-1:q+1])))
            if len(sides)==2:width.append(sum(sides))
    return dict(half_density_width_median=float(np.median(width)),
        half_density_width_grid_cells=float(np.median(width)/(L/N)),
        half_density_width_range=[float(min(width)),float(max(width))],successful_sections=len(width))


if __name__=='__main__':
    p=analytic_trefoil();t=knot_check(p);assert t['trefoil'],t
    s=np.linspace(0,2*np.pi,301);circle=np.array([np.cos(s),np.sin(s),np.zeros_like(s)]).T
    c=knot_check(circle);assert all(v['determinant']==1 for v in c['diagrams']),c
    print('Analytic trefoil:',t,'Circle:',c)
