"""Collimation with segment focus: Noll 4-8 at 3 field points; signature (alpha1, alpha2, delta, epsilon)
with epsilon the change of the field-independent defocus Omega_020."""
import sys as _sys, pathlib as _pl; _P=_pl.Path(__file__).resolve().parents[1]; _sys.path[:0]=[str(_P/'model'),str(_P/'field_height'),str(_P/'rigid_body')]
import numpy as np
import sens as S
from core import QX, QY, QW, SEG, R_SEG, trace, Model
from collimate import modes, HV
d=np.load(str(_P/'field_height'/'full.npz')); TL=[tuple(t) for t in d['terms']]; W=d['w']
R2=QX**2+QY**2; T=np.arctan2(QY,QX); r=np.sqrt(R2)
P48=np.array([np.sqrt(3)*(2*R2-1),np.sqrt(6)*R2*np.sin(2*T),np.sqrt(6)*R2*np.cos(2*T),
              np.sqrt(8)*(3*r**3-2*r)*np.sin(T),np.sqrt(8)*(3*r**3-2*r)*np.cos(T)])*QW
def model48(seg,Hs,alpha,dw040=0.,dw020=0.):
    sx,sy=SEG[seg]; px,py=QX+sx,QY+sy; pp=px*px+py*py; out=[]
    for h in Hs:
        Wf=0.
        for j in range(2):
            ux,uy=h[0]-alpha[j][0],h[1]-alpha[j][1]; uu=ux*ux+uy*uy; up=ux*px+uy*py
            Wf=Wf+sum(c*uu**p*pp**n*up**m for (p,n,m),c in zip(TL,W[j]))
        Wf=Wf+dw040*pp**2+dw020*pp
        out.append(P48@Wf)
    return np.array(out)
def estimate48(seg,Hs,bm,bn,iters=4):
    """unknowns: a1x a1y a2x a2y, dw040 (outer only), dw020"""
    outer=seg<6; n=6 if outer else 5
    off=bn-model48(seg,Hs,np.zeros((2,2)))
    def f(th):
        a=th[:4].reshape(2,2); d40=th[4] if outer else 0.; d20=th[-1]
        return (model48(seg,Hs,a,d40,d20)+off).ravel()
    th=np.zeros(n); y=bm.ravel()
    for _ in range(iters):
        r=y-f(th); J=np.empty((y.size,n))
        for i in range(n):
            e=np.zeros(n); e[i]=1e-6; J[:,i]=(f(th+e)-f(th-e))/2e-6
        th+=np.linalg.lstsq(J,r,rcond=None)[0]
    return th, np.sqrt(((y-f(th))**2).mean()), J
if __name__=='__main__':
    from rigid import MAS
    np.set_printoptions(precision=3,suppress=True,linewidth=150)
    out={}
    for seg in (0,6):
        bn=modes(seg,S.H)[:,:5]
        for mir in ('m1','m2'):
            cols=[]
            for i in range(5):
                e=np.zeros(5); e[i]=S.UNIT
                bm=modes(seg,S.H,**{'x1' if mir=='m1' else 'x2':e})[:,:5]
                th,rms,J=estimate48(seg,S.H,bm,bn)
                v=np.r_[th[:4]*MAS,(th[4] if seg<6 else np.nan),th[-1]]; cols.append(v)
            out[(seg,mir)]=np.array(cols).T
            print('segment %d %s: rows a1x a1y a2x a2y dw040 dw020 (mas, nm) ; cols dx dy dz rx ry\n%s'%(seg+1,mir,out[(seg,mir)]))
        print(' cond of the 15x6 (outer) / 15x5 (centre) Jacobian: %.0f'%np.linalg.cond(J))
    np.save('sens48.npy',out,allow_pickle=True)
