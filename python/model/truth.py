import numpy as np
from raytrace import trace
from emulate import disc, FV
R=8.365/2
pts=[]
for sid in range(1,8):
    if sid<7:
        o=(3-2*(sid-1))*np.pi/6; origin=(8.71*np.cos(o),8.71*np.sin(o))
    else: origin=(0.,0.)
    V,T=disc(8.365,0.25,origin); pts.append(V)
P=np.vstack(pts)                   # entrance coordinates, m
rho=P/R                            # global pupil coordinate in segment radii (rho_s = rho + s)
rows=[];Z=[];Rh=[]
for xy in FV:
    zen=np.deg2rad(np.hypot(*xy)/60); azi=np.arctan2(xy[1],xy[0])
    _,opd=trace(P,zen,azi); rows.append(opd*1e9)
OPD=np.array(rows)                 # (Nf, Np) nm
zeta=FV/10.
np.save('truth_opd.npy',OPD)
def basis(kmax,lmax,mmax):
    return [(p,n,m) for m in range(mmax+1) for p in range(10) for n in range(10) if 2*p+m<=kmax and 2*n+m<=lmax]
def fit(kmax,lmax,mmax):
    tl=basis(kmax,lmax,mmax)
    zx,zy=zeta[:,0][:,None],zeta[:,1][:,None]; rx,ry=rho[:,0][None,:],rho[:,1][None,:]
    cols=[((zx**2+zy**2)**p*(rx**2+ry**2)**n*(zx*rx+zy*ry)**m).ravel() for p,n,m in tl]
    M=np.array(cols).T; y=OPD.ravel()
    sc=np.abs(M).max(0); w,*_=np.linalg.lstsq(M/sc,y,rcond=None); w/=sc
    res=y-M@w
    return dict(zip(tl,w)), np.sqrt((res**2).mean()), np.abs(res).max()
