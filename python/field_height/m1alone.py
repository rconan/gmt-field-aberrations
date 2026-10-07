"""Per-surface split of the wave aberration: W_1 = M1 alone (prime-focus wavefront, entrance-pupil
coordinates), W_2 = W - W_1. Both fitted with the polynomial basis of truth.py."""
import sys as _sys, pathlib as _pl; _P=_pl.Path(__file__).resolve().parents[1]; _sys.path[:0]=[str(_P/'model'),str(_P/'field_height'),str(_P/'rigid_body')]

import numpy as np
from raytrace import surface, M1_K, M1_C
from truth import rho, zeta, basis, P
from emulate import FV
OPD=np.load(str(_P/'model'/'truth_opd.npy'))
def trace_m1(xy, zen, azi):
    xy=np.vstack([[0.,0.],xy])
    d=np.array([np.sin(zen)*np.cos(azi),np.sin(zen)*np.sin(azi),-np.cos(zen)])
    Pp=np.c_[xy,np.zeros(len(xy))]; D=np.tile(d,(len(xy),1))
    opd=np.tan(zen)*np.cos(azi)*xy[:,0]+np.tan(zen)*np.sin(azi)*xy[:,1]
    Pp,D,L=surface(Pp,D,None,M1_K,M1_C,0.); opd+=L-L[0]
    t=(0.5/M1_C-Pp[0,2])/D[0,2]; O=Pp[0]+D[0]*t; Rs=np.linalg.norm(O-Pp[0])
    V=Pp-O; g=np.einsum('ij,ij->i',V,D); r2=np.einsum('ij,ij->i',V,V)
    s=-np.sqrt(g*g-(r2-Rs**2))-g
    opd+=s-s[0]; return opd[1:]
rows=[]
for xy in FV:
    zen=np.deg2rad(np.hypot(*xy)/60); azi=np.arctan2(xy[1],xy[0])
    rows.append(trace_m1(P,zen,azi)*1e9)
OPD1=np.array(rows); np.save('opd_m1.npy',OPD1)
idx=np.arange(0,rho.shape[0],4)
def fit(O,kmax,lmax,mmax):
    tl=basis(kmax,lmax,mmax); r=rho[idx]; y=O[:,idx].ravel()
    zx,zy=zeta[:,0][:,None],zeta[:,1][:,None]; rx,ry=r[:,0][None,:],r[:,1][None,:]
    M=np.empty((y.size,len(tl)))
    for c,(p,n,m) in enumerate(tl): M[:,c]=((zx**2+zy**2)**p*(rx**2+ry**2)**n*(zx*rx+zy*ry)**m).ravel()
    sc=np.abs(M).max(0); M/=sc
    w,*_=np.linalg.lstsq(M,y,rcond=None); res=y-M@w; w/=sc
    return dict(zip(tl,w)), np.sqrt((res**2).mean())
NAMES=[(0,2,0),(0,1,1),(0,0,2),(1,2,0),(0,1,2),(1,1,1),(0,3,0),(0,2,1),(1,3,0),(0,2,2),(1,2,1)]
lab=lambda t:'w%d%d%d'%(2*t[0]+t[2],2*t[1]+t[2],t[2])
out={}
for name,O in (('tot',OPD),('m1',OPD1)):
    for K in [(7,10,4),(9,12,5)]:
        w,rms=fit(O,*K); print(name,K,'rms %.2e'%rms, ' '.join('%s %.4f'%(lab(t),w[t]) for t in NAMES))
    out[name]=w
np.save('wsplit.npy',np.array([[out['tot'][t] for t in NAMES],[out['m1'][t] for t in NAMES]]))
