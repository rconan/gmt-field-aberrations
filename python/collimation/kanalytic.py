import sys as _sys, pathlib as _pl; _P=_pl.Path(__file__).resolve().parents[1]; _sys.path[:0]=[str(_P/'model'),str(_P/'field_height'),str(_P/'rigid_body')]

import numpy as np
from coeffs import coef
from rigid import MAS, TH, LC
from seidel_parts import r1, z2
import sens as S
from sens import m1pivot
from core import m2_footprint
F_TH=0.60394279; MC=(-5.830005311184305-LC)/(r1/2-LC)
A={k:coef(k)*MAS for k in ('C1','V1','C2','V2')}          # mas per metre
# line of sight (nominal reference): p = -(m_C C1 + (1-m_C) C2)/(f theta), added to both field heights
A['C1']=A['C1']+(-MC/F_TH)*MAS; A['C2']=A['C2']+(-(1-MC)/F_TH)*MAS
def disp(x,P,zc,zv):
    T=x[:3]; tau=np.array([x[4],-x[3]])
    return T[:2]+(zc-P[2])*tau, T[:2]+(zv-P[2])*tau, T[2]-x[3]*P[1]+x[4]*P[0]   # C_perp, V_perp, V_z
def K_analytic(seg):
    P1,P2=m1pivot(seg),m2_footprint(seg); K=np.zeros((5,5))
    for i in range(5):
        x1=np.eye(5)[i]*1e-6
        C1,V1,V1z=disp(x1,P1,r1,0.)
        rhs=-(np.outer(A['C1'],C1)+np.outer(A['V1'],V1))          # (alpha_1, alpha_2) x (x,y)
        # unknowns: C2, V2 (lateral) from 2x2 per component
        M=np.c_[A['C2'],A['V2']]; CV=np.linalg.solve(M,rhs)       # rows C2, V2
        C2,V2=CV[0],CV[1]; tau2=(C2-V2)/(LC-z2)
        T2=V2-(z2-P2[2])*tau2
        rx2,ry2=-tau2[1],tau2[0]
        T2z=V1z+rx2*P2[1]-ry2*P2[0]                              # despace null: V2z = V1z
        K[:,i]=np.r_[T2,T2z,rx2,ry2]/1e-6
    return K
np.set_printoptions(precision=3,suppress=True,linewidth=140)
s=np.load('sens.npy',allow_pickle=True).item()
for seg in (6,0):
    D1,D2=s[(seg,'m1','D')],s[(seg,'m2','D')]
    rows=list(range(5)) if seg<6 else [0,1,2,3]; cols=list(range(5)) if seg<6 else [0,1,3,4]
    Kn=np.zeros((5,5)); Kn[np.ix_(cols,range(5))]=-np.linalg.solve(D2[np.ix_(rows,cols)],D1[rows])
    Ka=K_analytic(seg)
    if seg==6: Ka[2]=0.
    print('segment %d  analytic K:\n%s\n  ray-traced K:\n%s'%(seg+1,Ka,Kn))
