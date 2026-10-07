import numpy as np
from model import *
from crseo_data import Y
# Noll modes j=1..10 on the unit disk, evaluated at arbitrary points
def nollj(j,rr,tt):
    d={1:np.ones_like(rr),2:2*rr*np.cos(tt),3:2*rr*np.sin(tt),4:np.sqrt(3)*(2*rr**2-1),
       5:np.sqrt(6)*rr**2*np.sin(2*tt),6:np.sqrt(6)*rr**2*np.cos(2*tt),
       7:np.sqrt(8)*(3*rr**3-2*rr)*np.sin(tt),8:np.sqrt(8)*(3*rr**3-2*rr)*np.cos(tt),
       9:np.sqrt(8)*rr**3*np.sin(3*tt),10:np.sqrt(8)*rr**3*np.cos(3*tt)}
    return d[j]
def annulus_projector(eps, jmax=10, n=24, m=48):
    xr,wr=np.polynomial.legendre.leggauss(n); rr=eps+(1-eps)*0.5*(xr+1); wr=(1-eps)*0.5*wr*rr
    tt=2*np.pi*np.arange(m)/m
    RR,TT=np.meshgrid(rr,tt,indexing='ij'); ww=np.outer(wr,np.full(m,2*np.pi/m)).ravel()
    B=np.array([nollj(j,RR,TT).ravel() for j in range(1,jmax+1)])     # (J,N)
    G=(B*ww)@B.T                                                       # Gram
    Pj=np.linalg.solve(G,B*ww)                                         # LS projector (J,N)
    sel=[4,5,6,7]  # Noll 5..8 rows
    return Pj[sel], (RR*np.cos(TT)).ravel(), (RR*np.sin(TT)).ravel()
def response_centre(p,n,m,eps):
    Pp,ax,ay=annulus_projector(eps)
    Hx=X[:,None]; Hy=Yc[:,None]; Px=ax[None,:]; Py=ay[None,:]
    W=(Hx**2+Hy**2)**p*(Px**2+Py**2)**n*(Hx*Px+Hy*Py)**m
    return Pp @ W.T @ ZF.T
def fit_ann(tl,s,eps,phis):
    A=np.zeros((4,8,7,len(tl)))
    for it,(p,n,m) in enumerate(tl):
        for k,a in enumerate(phis): A[:,:,k,it]=response(p,n,m,s*np.cos(a),s*np.sin(a))
        A[:,:,6,it]=response_centre(p,n,m,eps) if eps>0 else response(p,n,m,0,0)
    M=A.reshape(-1,len(tl)); w,*_=np.linalg.lstsq(M,Y.ravel(),rcond=None)
    return w,(Y.ravel()-M@w).reshape(Y.shape)
