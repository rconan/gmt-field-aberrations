import sys
import numpy as np
import sens as S
from core import QX, QY, QW, SEG, R_SEG, trace
from rigid import MAS
G2={seg:np.load('sens.npy',allow_pickle=True).item()[(seg,'m2','G')] for seg in (0,6)}
# validation: 20 field points; Noll 4 (focus) and 9-11 as well as 5-8
rng=np.random.default_rng(7); HV=rng.uniform(-1,1,(80,2)); HV=HV[np.hypot(*HV.T)<=1][:20]
R2=QX**2+QY**2; T=np.arctan2(QY,QX); r=np.sqrt(R2)
Z={4:np.sqrt(3)*(2*R2-1),5:np.sqrt(6)*R2*np.sin(2*T),6:np.sqrt(6)*R2*np.cos(2*T),7:np.sqrt(8)*(3*r**3-2*r)*np.sin(T),
   8:np.sqrt(8)*(3*r**3-2*r)*np.cos(T),9:np.sqrt(8)*r**3*np.sin(3*T),10:np.sqrt(8)*r**3*np.cos(3*T),11:np.sqrt(5)*(6*R2**2-6*R2+1)}
PJ=np.array([Z[j]*QW for j in sorted(Z)])
def modes(seg,Hs,x1=None,x2=None):
    sx,sy=SEG[seg]; xy=np.c_[(QX+sx)*R_SEG,(QY+sy)*R_SEG]; kw={}
    if x1 is not None: kw['m1']=S.mis(x1,seg,'m1')
    if x2 is not None: kw['m2']=S.mis(x2,seg,'m2')
    out=[]
    for h in Hs:
        zen=np.deg2rad(np.hypot(*h)*10/60); azi=np.arctan2(h[1],h[0]); o,_=trace(xy,zen,azi,**kw); out.append(PJ@(o*1e9))
    return np.array(out)              # columns: Noll 4,5,...,11
def collimate(seg,x1,x2,iters=4):
    """M2 correction from Noll 5-8 at the 3 field points, Gauss-Newton with the ray-traced G2"""
    cols=[0,1,2,3,4] if seg<6 else [0,1,3,4]
    Gp=np.linalg.pinv(G2[seg][:,cols])
    bn=S.b(seg); dx=np.zeros(5); hist=[]
    for it in range(iters):
        x=x2+dx; db=(S.b(seg,x1=x1,x2=x)-bn).ravel(); hist.append(np.sqrt((db**2).mean()))
        dx[cols]-=Gp@db*S.UNIT
    return dx,hist
if __name__=='__main__':
    names=['dx','dy','dz','rx','ry']
    scale=np.array([50e-6,50e-6,50e-6,10e-6,10e-6])
    for seg in (0,6):
        mv0=modes(seg,HV)
        for case in ('M1','M2','M1+M2'):
            x1=rng.normal(size=5)*scale if 'M1' in case else np.zeros(5)
            x2=rng.normal(size=5)*scale if 'M2' in case else np.zeros(5)
            dx,hist=collimate(seg,x1,x2)
            before=modes(seg,HV,x1,x2)-mv0; after=modes(seg,HV,x1,x2+dx)-mv0
            rms=lambda a,c: np.sqrt((a[:,c]**2).mean())
            l0=S.los(seg,x1,x2); l1=S.los(seg,x1,x2+dx)
            print('seg%d %-6s M2 correction %s um/urad'%(seg+1,case,' '.join('%s %+7.2f'%(n,v*1e6) for n,v in zip(names,dx))))
            print('     fit-point rms per iteration (nm): %s'%' '.join('%.2e'%h for h in hist))
            print('     20 field points, rms change (nm): astig+coma %8.3f -> %.4f | focus %7.3f -> %7.3f | trefoil+sph %.3f -> %.4f | LOS (mas) (%+.0f,%+.0f) -> (%+.0f,%+.0f)'%(
                rms(before,[1,2,3,4]),rms(after,[1,2,3,4]),rms(before,[0]),rms(after,[0]),rms(before,[5,6,7]),rms(after,[5,6,7]),*l0,*l1))
