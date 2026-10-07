"""Ray-trace check that an M1-segment figure error is field-constant (M1 is the stop)."""
import numpy as np, sens as S
import raytrace as RT, core
from collimate import modes, HV
from core import SEG, R_SEG, M1_C
_sag,_dsag=RT.sag,RT.dsag; FIG={'on':False,'seg':0,'j':5,'a':0.}
def zern(j,u,v):
    r2=u*u+v*v; t=np.arctan2(v,u); r=np.sqrt(r2)
    return {4:np.sqrt(3)*(2*r2-1),5:np.sqrt(6)*r2*np.sin(2*t),6:np.sqrt(6)*r2*np.cos(2*t),
            7:np.sqrt(8)*(3*r**3-2*r)*np.sin(t),8:np.sqrt(8)*(3*r**3-2*r)*np.cos(t)}[j]
def fig(x,y):
    sx,sy=SEG[FIG['seg']]; u=x/R_SEG-sx; v=y/R_SEG-sy
    return FIG['a']*zern(FIG['j'],u,v)*(u*u+v*v<1.2)
def sag(x,y,k,c):
    s=_sag(x,y,k,c)
    return s+fig(x,y) if FIG['on'] and c==M1_C else s
def dsag(x,y,k,c):
    gx,gy=_dsag(x,y,k,c)
    if FIG['on'] and c==M1_C:
        e=1e-4; gx=gx+(fig(x+e,y)-fig(x-e,y))/(2*e); gy=gy+(fig(x,y+e)-fig(x,y-e))/(2*e)
    return gx,gy
RT.sag,RT.dsag=sag,dsag
_surf=core.surface
def surface(P,D,opl,k,c,z0):           # the chief ray (single-row call) never sees the figure error
    on=FIG['on']; FIG['on']=on and len(P)>1
    try: return _surf(P,D,opl,k,c,z0)
    finally: FIG['on']=on
core.surface=surface
np.set_printoptions(precision=3,suppress=True,linewidth=160)
for seg in (0,6):
    n0=modes(seg,HV)
    for j in (4,5,6,7,8):
        FIG.update(on=True,seg=seg,j=j,a=5e-9)       # 5 nm rms surface ~ 10 nm rms wavefront
        m=modes(seg,HV)-n0; FIG['on']=False
        mean=m.mean(0); var=m-mean
        print('seg%d Z%d surface 5 nm: mean Noll4..11 %s nm | field-dependent rms %.4f nm'%(seg+1,j,mean,np.sqrt((var**2).mean())))
