import numpy as np
from truth import rho, zeta, basis
OPD=np.load('truth_opd.npy')
idx=np.arange(0,rho.shape[0],4)
def fit(kmax,lmax,mmax):
    tl=basis(kmax,lmax,mmax)
    r=rho[idx]; y=OPD[:,idx].ravel()
    zx,zy=zeta[:,0][:,None],zeta[:,1][:,None]; rx,ry=r[:,0][None,:],r[:,1][None,:]
    M=np.empty((y.size,len(tl)))
    for c,(p,n,m) in enumerate(tl):
        M[:,c]=((zx**2+zy**2)**p*(rx**2+ry**2)**n*(zx*rx+zy*ry)**m).ravel()
    sc=np.abs(M).max(0); M/=sc
    w,*_=np.linalg.lstsq(M,y,rcond=None); res=y-M@w; w/=sc
    return dict(zip(tl,w)), np.sqrt((res**2).mean()), np.abs(res).max()
