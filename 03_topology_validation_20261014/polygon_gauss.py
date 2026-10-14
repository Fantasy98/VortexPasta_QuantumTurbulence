"""Exact Gauss integral for straight polygon segments via spherical solid angles.
For r(s,t)=p(s)-q(t), dOmega = -[(dp x dq) dot r]/|r|^3 ds dt.
Each parameter-square image is split into two spherical triangles.
Triangle solid angle uses the robust atan2 determinant/dot-product formula.
Shared-vertex (adjacent) straight edges have zero contribution.
"""
import numpy as np

def triangle(a,b,c):
 na=np.linalg.norm(a,axis=-1);nb=np.linalg.norm(b,axis=-1);nc=np.linalg.norm(c,axis=-1)
 numer=np.einsum('...i,...i->...',a,np.cross(b,c))
 den=na*nb*nc+np.einsum('...i,...i->...',a,b)*nc+np.einsum('...i,...i->...',b,c)*na+np.einsum('...i,...i->...',c,a)*nb
 return 2*np.arctan2(numer,den)

def gauss_exact(a,b,self_=False):
 # Chunking bounds memory for fine diagnostic sampling.
 total=0.;n=len(a);m=len(b);bn=np.roll(b,-1,axis=0)
 for start in range(0,n,128):
  inds=np.arange(start,min(n,start+128));aa=a[inds];an=a[(inds+1)%n]
  r00=aa[:,None]-b[None,:];r10=an[:,None]-b[None,:];r11=an[:,None]-bn[None,:];r01=aa[:,None]-bn[None,:]
  omega=triangle(r00,r10,r11)+triangle(r00,r11,r01)
  if self_:
   jj=np.arange(m)[None,:];dd=abs(inds[:,None]-jj);omega[np.minimum(dd,n-dd)<=1]=0
  total+=omega.sum()
 return float(-total/(4*np.pi))
