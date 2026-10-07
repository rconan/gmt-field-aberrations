import numpy as np
from seidel_parts import omegas, R, TH, r1, r2, z2, K1, K2
import seidel_parts as S
from rigid import alpha, MAS, LC
W=omegas(); w1s,w2s,w2a,w1a=W[('M1','sph')],W[('M2','sph')],W[('M2','asph')],W[('M1','asph')]
M=np.array([[w1s['w222'],w2s['w222']+w2a['w222']],[w1s['w131'],w2s['w131']+w2a['w131']]]); Mi=np.linalg.inv(M)
# beta = sum_k  (w222_k, w131_k) * sigma_k ; sigma_k linear in C1, V1, C2, V2 (lateral, metres)
# sigma1s = -C1/(R1 th); D1 = V1/R; sigma2s = C2/(L th) - 2 C1/(R1 th); sigma2a = V2/(z2 th) - 2 C1/(R1 th)
def coef(which):
    s1s=s2s=s2a=d1=0.
    if which=='C1': s1s=-1/(r1*TH); s2s=s2a=-2/(r1*TH)
    if which=='V1': d1=1/R
    if which=='C2': s2s=1/(LC*TH)
    if which=='V2': s2a=1/(z2*TH)
    b222=w1s['w222']*s1s+w2s['w222']*s2s+w2a['w222']*s2a
    b131=w1s['w131']*s1s+w2s['w131']*s2s+w2a['w131']*s2a+4*w1a['w040']*d1
    return Mi@np.array([b222,b131])         # (alpha1, alpha2) per metre of displacement
print('field heights per um of lateral displacement (mas/um):')
for k in ('C1','V1','C2','V2'):
    a=coef(k)*1e-6*MAS; print('  %s: alpha_M1 %+9.4f   alpha_M2 %+9.4f'%(k,*a))
# check against rigid.alpha for a random motion
rng=np.random.default_rng(3); T1,T2=rng.normal(size=3)*1e-6,rng.normal(size=3)*1e-6
o1,o2=np.r_[rng.normal(size=2)*1e-6,0],np.r_[rng.normal(size=2)*1e-6,0]; P1=np.array([0,8.71,1.05]); P2=np.array([0,-1.1,20.1])
a=alpha(m1=(T1,o1,P1),m2=(T2,o2,P2))
mv=lambda X0,T,o,P: X0+np.cross(o,X0-P)+T
C1=mv(np.array([0,0,r1]),T1,o1,P1)[:2]; V1=mv(np.zeros(3),T1,o1,P1)[:2]
C2=mv(np.array([0,0,LC]),T2,o2,P2)[:2]; V2=mv(np.array([0,0,z2]),T2,o2,P2)[:2]
lin=sum(np.outer(coef(k),X) for k,X in (('C1',C1),('V1',V1),('C2',C2),('V2',V2)))
print('linear form reproduces alpha(): max diff %.1e'%np.abs(lin-a).max())
# despace: d Omega_040 / d (M1-M2 separation), from the Seidel sums
def w040(dz):
    h,u,hb,ub,n=R,0.,0.,TH,1.; tot=0.
    for c,k,d in ((-1/r1,K1,-(z2+dz)),(-1/r2,K2,0.)):
        np_=-n; up=(n*u-h*c*(np_-n))/np_; A=n*(u+h*c); D=up/np_-u/n
        tot+=-A*A*h*D+k*c**3*(np_-n)*h**4; u,n=up,np_; h+=d*u
    return -tot/8*1e9
e=1e-6; print('dOmega040/d(separation) = %+.4f nm per um'%((w040(e)-w040(-e))/(2*e)*1e-6))
