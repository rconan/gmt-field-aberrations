import sys as _sys, pathlib as _pl; _P=_pl.Path(__file__).resolve().parents[1]; _sys.path[:0]=[str(_P/'model'),str(_P/'field_height'),str(_P/'rigid_body')]
import numpy as np, sys; 

from core import Model
from truth import rho, zeta, basis
OPD=np.load(str(_P/'model'/'truth_opd.npy')); OPD1=np.load(str(_P/'field_height'/'opd_m1.npy'))
idx=np.arange(0,rho.shape[0],4)
def fitall(O,tl):
    r=rho[idx]; y=O[:,idx].ravel()
    zx,zy=zeta[:,0][:,None],zeta[:,1][:,None]; rx,ry=r[:,0][None,:],r[:,1][None,:]
    M=np.empty((y.size,len(tl)))
    for c,(p,n,m) in enumerate(tl): M[:,c]=((zx**2+zy**2)**p*(rx**2+ry**2)**n*(zx*rx+zy*ry)**m).ravel()
    sc=np.abs(M).max(0); w,*_=np.linalg.lstsq(M/sc,y,rcond=None); return w/sc
TL=basis(9,12,5)
wt=fitall(OPD,TL); w1=fitall(OPD1,TL)
np.savez(str(_P/'field_height'/'full.npz'),terms=np.array(TL),w=np.array([w1,wt-w1]))
