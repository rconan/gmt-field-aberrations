"""Four-part decomposition of the GMT wavefront (all orders), from ray traces with
spherical mirrors: W1s = M1 sphere alone, W1a = M1 - W1s, W2s = W(M2 sphere) - W1, W2a = W - W(M2 sphere)."""
import sys as _sys, pathlib as _pl; _P=_pl.Path(__file__).resolve().parents[1]; _sys.path[:0]=[str(_P/'model'),str(_P/'field_height'),str(_P/'rigid_body')]

import numpy as np, raytrace as RT
from raytrace import surface
import truth as T
from emulate import FV
def tel(xy,zen,azi,k1,k2,m1only=False):
    xy=np.vstack([[0.,0.],xy]); d=np.array([np.sin(zen)*np.cos(azi),np.sin(zen)*np.sin(azi),-np.cos(zen)])
    P=np.c_[xy,np.zeros(len(xy))]; D=np.tile(d,(len(xy),1))
    opd=np.tan(zen)*np.cos(azi)*xy[:,0]+np.tan(zen)*np.sin(azi)*xy[:,1]
    P,D,L=surface(P,D,None,k1,RT.M1_C,0.); opd+=L-L[0]
    if m1only:
        t=(0.5/RT.M1_C-P[0,2])/D[0,2]; O=P[0]+D[0]*t; Rs=np.linalg.norm(O-P[0])
    else:
        P,D,L=surface(P,D,None,k2,RT.M2_C,RT.M2_Z); opd+=L-L[0]
        c=P[0];dc=D[0];x,y,z=c[0],c[1],c[2]-RT.FP_Z+RT.FP_R
        g=x*dc[0]+y*dc[1]+z*dc[2];r2=x*x+y*y+z*z;s=-np.sqrt(g*g-(r2-RT.FP_R**2))-g;O=c+dc*s
        rkl=dc[0]**2+dc[1]**2
        if rkl<1e-24: Rs=23.772110269559725
        else:
            se=np.sqrt((c[0]**2+c[1]**2)/rkl);zE=c[2]+dc[2]*se;Rs=np.sqrt(O[0]**2+O[1]**2+(O[2]-zE)**2)
    V=P-O;g=np.einsum('ij,ij->i',V,D);r2=np.einsum('ij,ij->i',V,V);s=-np.sqrt(g*g-(r2-Rs**2))-g;opd+=s-s[0]
    return opd[1:]
def opds(**kw):
    return np.array([tel(T.P,np.deg2rad(np.hypot(*xy)/60),np.arctan2(xy[1],xy[0]),**kw)*1e9 for xy in FV])
TL=T.basis(9,12,5)
def fit(O):
    idx=np.arange(0,T.rho.shape[0],4); r=T.rho[idx]; y=O[:,idx].ravel()
    zx,zy=T.zeta[:,0][:,None],T.zeta[:,1][:,None]; rx,ry=r[:,0][None,:],r[:,1][None,:]
    M=np.array([((zx**2+zy**2)**p*(rx**2+ry**2)**n*(zx*rx+zy*ry)**m).ravel() for p,n,m in TL]).T
    sc=np.abs(M).max(0); w,*_=np.linalg.lstsq(M/sc,y,rcond=None); res=y-(M/sc)@w
    return w/sc, np.sqrt((res**2).mean())
if __name__=='__main__':
    k1,k2=RT.M1_K,RT.M2_K
    runs={'m1':dict(k1=k1,k2=k2,m1only=True),'m1s':dict(k1=1.,k2=k2,m1only=True),
          'tel':dict(k1=k1,k2=k2),'tel_m2s':dict(k1=k1,k2=1.)}
    W={}
    for n,kw in runs.items():
        W[n],rms=fit(opds(**kw)); print(n,'fit rms %.1e nm'%rms)
    parts=dict(w1s=W['m1s'],w1a=W['m1']-W['m1s'],w2s=W['tel_m2s']-W['m1'],w2a=W['tel']-W['tel_m2s'])
    np.savez('parts.npz',terms=np.array(TL),**parts)
    from seidel_parts import omegas
    S=omegas()
    for lab,t in (('w040',(0,2,0)),('w131',(0,1,1)),('w222',(0,0,2)),('w060',(0,3,0))):
        i=TL.index(t); print(lab,' '.join('%s %.3f'%(k,v[i]) for k,v in parts.items()))
    print('Seidel parts:',{k:{a:round(b,3) for a,b in v.items()} for k,v in S.items()})
