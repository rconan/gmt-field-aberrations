import sys as _sys, pathlib as _pl; _P=_pl.Path(__file__).resolve().parents[1]; _sys.path[:0]=[str(_P/'model'),str(_P/'field_height'),str(_P/'rigid_body')]
import numpy as np, sys
from emulate import disc, lump, noll10, wls_projector, FP, zeta
from noll2cplx import to_complex
D_SEG,R_OFF=8.365,8.71; R=D_SEG/2
def seg_ops(coords):
    out=[]
    for sid in range(1,8):
        o=(3-2*(sid-1))*np.pi/6 if sid<7 else 0.
        org=np.array((R_OFF*np.cos(o),R_OFF*np.sin(o))) if sid<7 else np.zeros(2)
        V,T=disc(D_SEG,0.25,tuple(org))
        if coords=='sphere': loc=-(V-org); g=loc/R-org/R
        else: loc=V-org; g=V/R
        out.append((wls_projector(noll10(loc[:,0],loc[:,1]),lump(V,T)),g))
    return out
OPS={c:seg_ops(c) for c in ('sphere','entrance')}
def response(p,n,m,coords):
    out=np.zeros((4,8,7)); zx,zy=zeta[:,0][:,None],zeta[:,1][:,None]
    for k,(P,g) in enumerate(OPS[coords]):
        gx,gy=g[:,0][None,:],g[:,1][None,:]
        W=(zx**2+zy**2)**p*(gx**2+gy**2)**n*(zx*gx+zy*gy)**m
        out[:,:,k]=(FP@(W@P.T))[:8,4:8].T
    return out
def terms(kmax,lmax,mmax):
    return [(p,n,m) for m in range(mmax+1) for p in range(10) for n in range(10) if 2*p+m<=kmax and 2*n+m<=lmax]
SETS=[(3,4,2),(3,6,2),(5,8,3)]
def fits(data,coords):
    res=[]
    for K in SETS:
        tl=terms(*K); M=np.array([response(*t,coords).ravel() for t in tl]).T
        w,*_=np.linalg.lstsq(M,data.ravel(),rcond=None); r=data.ravel()-M@w
        res.append((np.sqrt((r**2).mean()),np.abs(r).max()))
    return res
def checks(b):
    amp=lambda x,y: np.hypot(x,y)
    def ratio(jf):
        A=amp(b[0,jf-1,:6],b[1,jf-1,:6]); C=amp(b[2,jf-1,:6],b[3,jf-1,:6]); return (C/A).mean()
    a=to_complex(b); N=np.sqrt(8)*np.sqrt(4)   # Noll amplitude of a counter/co-rotating field-tilt term
    out={'defocus ratio':ratio(4),'piston ratio':ratio(1)}
    for nh,lab,Nh in ((1,'tilt',np.sqrt(8)*np.sqrt(4)),(3,'coma',np.sqrt(8)*np.sqrt(8))):
        co=a[(3,1,nh,-1)]; ctr=a[(3,1,nh,1)]
        pred=0.5*np.abs(co[:6]-co[6]); out['counter %s pred'%lab]=pred.mean()*Nh; out['counter %s meas'%lab]=np.abs(ctr[:6]).mean()*Nh
    out['seg1 astig x tilt pair']=(b[0,1,0],b[1,2,0])
    return out
if __name__=='__main__':
    from crseo_data import Y
    sets={'OLD exit (CEO)':(Y,'sphere'),'OLD entrance (replica)':(np.load(str(_P/'model'/'old_reference'/'replica_entrance.npy')),'entrance'),
          'NEW exit (replica)':(np.load('replica_sphere.npy'),'sphere'),'NEW entrance (replica)':(np.load('replica_entrance.npy'),'entrance')}
    for name,(b,c) in sets.items():
        print(name, ' fits:', '; '.join('rms %.3g max %.3g'%x for x in fits(b,c)))
        print('   ', {k:(np.round(v,4) if not isinstance(v,tuple) else tuple(np.round(v,3))) for k,v in checks(b).items()})
