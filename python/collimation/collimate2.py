import numpy as np
import sens as S
from collimate import modes, HV, rng
MODES={'5-8':[1,2,3,4],'4-8':[0,1,2,3,4],'4-8,11':[0,1,2,3,4,7]}
def Gm(seg,sel):
    cols=[]
    for i in range(5):
        e=np.zeros(5); e[i]=S.UNIT
        cols.append(((modes(seg,S.H,x2=e)-modes(seg,S.H,x2=-e))/2)[:,sel].ravel())
    return np.array(cols).T
def collimate(seg,x1,x2,sel,iters=4):
    G=Gm(seg,sel); Gp=np.linalg.pinv(G); bn=modes(seg,S.H)[:,sel]; dx=np.zeros(5)
    for _ in range(iters):
        db=(modes(seg,S.H,x1,x2+dx)[:,sel]-bn).ravel(); dx-=Gp@db*S.UNIT
    return dx,np.linalg.cond(G)
scale=np.array([50e-6,50e-6,50e-6,10e-6,10e-6])
rng=np.random.default_rng(11)
for seg in (6,0):
    mv0=modes(seg,HV)
    x1=rng.normal(size=5)*scale; x2=rng.normal(size=5)*scale
    for case,a,c in (('M2',np.zeros(5),x2),('M1',x1,np.zeros(5)),('M1+M2',x1,x2)):
        for mname,sel in MODES.items():
            dx,cond=collimate(seg,a,c,sel)
            after=modes(seg,HV,a,c+dx)-mv0
            r=lambda cc: np.sqrt((after[:,cc]**2).mean())
            extra=' | x2+dx = %s'%np.round((c+dx)*1e6,4) if case=='M2' else ''
            print('seg%d %-5s modes %-7s cond %6.1f | residual rms (nm) astig+coma %.4f focus %.3f trefoil %.4f spherical %.4f%s'%(seg+1,case,mname,cond,r([1,2,3,4]),r([0]),r([5,6]),r([7]),extra))
