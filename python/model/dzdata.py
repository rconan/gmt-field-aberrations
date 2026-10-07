"""Double-Zernike data from the delrays/CEO-equivalent ray trace (entrance-pupil coordinates),
with pupil and field fits extended to Noll modes 1..J (delrays default: J=10)."""
import numpy as np
from raytrace import trace
from emulate import disc, lump, wls_projector
from inverse import noll
D_SEG=8.365; R=D_SEG/2
def nollJ(J,x,y,rmax=None):
    r=np.hypot(x,y); r=r/(rmax or r.max()); t=np.arctan2(y,x)
    return np.array([noll(j,r,t) for j in range(1,J+1)])
FV,FT=disc(20.,2.)
def segments():
    out=[]
    for sid in range(1,8):
        if sid<7:
            o=(3-2*(sid-1))*np.pi/6; org=(8.71*np.cos(o),8.71*np.sin(o))
        else: org=(0.,0.)
        V,T=disc(D_SEG,0.25,org); out.append((V,T,np.array(org)))
    return out
def dz_data(J):
    FP=wls_projector(nollJ(J,FV[:,0],FV[:,1]),lump(FV,FT))
    out=np.zeros((J,J,7))
    for k,(V,T,org) in enumerate(segments()):
        loc=V-org; P=wls_projector(nollJ(J,loc[:,0],loc[:,1]),lump(V,T))
        pup=np.array([P@trace(V,np.deg2rad(np.hypot(*f)/60),np.arctan2(f[1],f[0]))[1] for f in FV])
        out[:,:,k]=(FP@pup).T*1e9          # (pupil j, field j, segment) nm
    return out
def dz_operator(J):
    """the same mesh-based fitting operators, applied to a model wavefront (for the forward model)"""
    FP=wls_projector(nollJ(J,FV[:,0],FV[:,1]),lump(FV,FT))
    segs=[]
    for V,T,org in segments():
        loc=V-org; segs.append((wls_projector(nollJ(J,loc[:,0],loc[:,1]),lump(V,T)),V/R))
    zeta=FV/10.
    def response(p,n,m):
        out=np.zeros((J,J,7)); zx,zy=zeta[:,0][:,None],zeta[:,1][:,None]
        for k,(P,g) in enumerate(segs):
            gx,gy=g[:,0][None,:],g[:,1][None,:]
            Wv=(zx*zx+zy*zy)**p*(gx*gx+gy*gy)**n*(zx*gx+zy*gy)**m
            out[:,:,k]=(FP@(Wv@P.T)).T
        return out
    return response
