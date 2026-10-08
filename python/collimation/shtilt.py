"""Segment tip/tilt (Noll 2,3) at 3 field points for the motions that mimic M1 figure errors:
does a Shack-Hartmann's per-segment tilt see them? 1 nm rms of Noll 2/3 = 2 nm/R_seg rad = 0.099 mas."""
import numpy as np, sens as S
from core import QX, QY, QW, SEG, R_SEG, trace
from figgeom3 import H  # equilateral r=1
PJ=np.array([2*QX*QW,2*QY*QW])
def tilt(seg,Hs,x1=None,x2=None):
    sx,sy=SEG[seg]; xy=np.c_[(QX+sx)*R_SEG,(QY+sy)*R_SEG]; kw={}
    if x1 is not None: kw['m1']=S.mis(x1,seg,'m1')
    if x2 is not None: kw['m2']=S.mis(x2,seg,'m2')
    out=[]
    for h in Hs:
        zen=np.deg2rad(np.hypot(*h)*10/60); azi=np.arctan2(h[1],h[0]); o,_=trace(xy,zen,azi,**kw); out.append(PJ@(o*1e9))
    return np.array(out)*2e-9/R_SEG*206265e3/1e-9*1e-9   # mas
