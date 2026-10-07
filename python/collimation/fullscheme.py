"""Collimation with focus: M2 (5 DOFs) from Noll 4-8 at 3 field points, plus pointing of the pair
with M1 tilts (M2 following); line of sight measured. Ray-trace test."""
import numpy as np
import sens as S
from collimate import modes, HV
G=np.load('G48.npy',allow_pickle=True).item()
WLOS=0.01                         # nm per mas: weight of the line of sight in the least squares
def P(seg,mir):
    k='x1' if mir=='m1' else 'x2'
    return np.array([(S.los(seg,**{k:np.eye(5)[i]*S.UNIT})-S.los(seg,**{k:-np.eye(5)[i]*S.UNIT}))/2 for i in range(5)]).T
def solve(seg,x1,x2,mode,iters=4):
    G1,G2=G[(seg,'Noll 4-8')]; P1,P2=P(seg,'m1'),P(seg,'m2')
    if mode=='M2 only':
        A=G2; dof1=[]
    else:                          # M2 (5) + M1 rx, ry
        dof1=[3,4]; A=np.r_[np.c_[G2,G1[:,dof1]],WLOS*np.c_[P2,P1[:,dof1]]]
    Ap=np.linalg.pinv(A,rcond=1e-10)
    b0=modes(seg,S.H)[:,:5]; d2=np.zeros(5); d1=np.zeros(5)
    for _ in range(iters):
        y1,y2=x1+d1,x2+d2
        r=(modes(seg,S.H,y1,y2)[:,:5]-b0).ravel()
        if mode!='M2 only': r=np.r_[r,WLOS*S.los(seg,y1,y2)]
        step=-Ap@r*S.UNIT
        d2+=step[:5]; d1[dof1]+=step[5:]
    return d1,d2
if __name__=='__main__':
    rng=np.random.default_rng(5); scale=np.array([50e-6,50e-6,50e-6,10e-6,10e-6])
    names=['dx','dy','dz','rx','ry']
    for seg in (0,6):
        mv0=modes(seg,HV)
        for case in ('M1','M1+M2'):
            x1=rng.normal(size=5)*scale; x2=rng.normal(size=5)*scale if case=='M1+M2' else np.zeros(5)
            b=modes(seg,HV,x1,x2)-mv0
            rr=lambda a,c: np.sqrt((a[:,c]**2).mean())
            print('seg%d %-6s before: astig+coma %8.3f focus %8.3f trefoil %6.3f spherical %6.3f nm  LOS (%+.0f,%+.0f) mas'%(seg+1,case,rr(b,[1,2,3,4]),rr(b,[0]),rr(b,[5,6]),rr(b,[7]),*S.los(seg,x1,x2)))
            for mode in ('M2 only','M2 + M1 tilts'):
                d1,d2=solve(seg,x1,x2,mode)
                a=modes(seg,HV,x1+d1,x2+d2)-mv0
                print('      %-14s after: astig+coma %8.4f focus %8.4f trefoil %6.4f spherical %6.4f nm  LOS (%+.1f,%+.1f) mas | M1 tilts %s urad'%(mode,rr(a,[1,2,3,4]),rr(a,[0]),rr(a,[5,6]),rr(a,[7]),*S.los(seg,x1+d1,x2+d2),np.round(d1[3:]*1e6,2)))
