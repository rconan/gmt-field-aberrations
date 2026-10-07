import numpy as np; from core import *
exec(open('check2.py').read().split('for seg in')[0])
def tri(r,rot=90): return np.array([[r*np.cos(a),r*np.sin(a)] for a in np.deg2rad([rot,rot+120,rot+240])])
mas=600e3
for seg in (0,6):
    for rr in (0.4,0.6,0.8,1.0):
        Hr=tri(rr); bn=full.bfield(seg,Hr,np.zeros((2,2)))
        J=(estimate_ext(full,seg,Hr,bn,bn,iters=1)[3] if seg<6 else estimate(full,seg,Hr,bn,bn,iters=1)[2])
        C=np.linalg.inv(J.T@J)
        # transform to d = a1-a2, c = (w1_222 a1 + w2_222 a2)/Omega222 (astigmatism node shift)
        w1,w2=d['w'][0][0],d['w'][1][0]
        TL=[tuple(t) for t in d['terms']]; i=TL.index((0,0,2)); g1,g2=d['w'][0][i],d['w'][1][i]; G=g1+g2
        T=np.zeros((4,C.shape[0])); T[0,0]=T[1,1]=1; T[0,2]=T[1,3]=-1
        T[2,0]=T[3,1]=g1/G; T[2,2]=T[3,3]=g2/G
        Cd=T@C@T.T; sd=np.sqrt(np.diag(Cd))*mas
        sv=np.linalg.svd(J,compute_uv=False)
        print('seg%d r=%.1f sd(a1-a2)=(%.3f,%.3f) mas  sd(astig node)=(%.1f,%.1f) mas  corr=%.2f'%(seg+1,rr,*sd,Cd[0,2]/np.sqrt(Cd[0,0]*Cd[2,2])))
print('w222 per surface',g1,g2,G)
