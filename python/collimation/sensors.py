"""Edge sensors on M1 or on M2: which one resolves the figure / misalignment degeneracy?
Unknowns: x1 (5), x2 (5), figure f (Noll 4-8, 5). Data: Noll 4-8 at 3 field points (0.2 nm),
line of sight (1 mas), and segment-position sensors on one mirror (all 5 DOFs, or dz rx ry only).
Linear model (G48 sensitivities, H at r=0.8)."""
import numpy as np
d=np.load('figerr.npy',allow_pickle=True).item()
SN,SL=0.2,1.0
def run(seg,which,dofs,st,NT=400):
    o=d[seg]; G1,G2,P1,P2=o['G1'],o['G2'],o['P1'],o['P2']; F=np.tile(np.eye(5),(3,1))
    s2=np.r_[np.full(15,SN),np.full(2,SL)]
    A=np.r_[np.c_[G1,G2,F],np.c_[P1,P2,np.zeros((2,5))]]
    if which:
        sig=np.array([st,st,st,st/10,st/10])[dofs]          # um, urad
        S=np.zeros((len(dofs),15)); off=0 if which=='M1' else 5
        for i,k in enumerate(dofs): S[i,off+k]=1
        A=np.r_[A,S]; s2=np.r_[s2,sig]
    W=1/s2; Aw=A*W[:,None]
    # posterior covariance with a weak prior (1 mm, 100 urad, 1 um of figure... i.e. essentially flat)
    prior=np.array([50,50,50,10,10]*2+[10]*5,float)             # the misalignment / figure distributions themselves
    C=np.linalg.inv(Aw.T@Aw+np.diag(1/prior**2))
    sd=np.sqrt(np.diag(C))
    # derived quantities: collimation state of M2 relative to M1 (x2 - K x1) and figure
    K=-np.linalg.pinv(G2)@G1
    J=np.c_[-K,np.eye(5),np.zeros((5,5))]; sc=np.sqrt(np.diag(J@C@J.T))
    return sd[:5],sd[5:10],sd[10:],sc
np.set_printoptions(precision=2,suppress=True,linewidth=170)
print('posterior std (prior: 50 um / 10 urad misalignments, 10 nm figure per mode)')
for seg in (0,6):
    print('\nsegment %d'%(seg+1))
    print('  %-34s %-26s %-26s %-26s %s'%('sensors','M1 dx dy dz rx ry','M2 dx dy dz rx ry','figure Z4..Z8 (nm)','M2 - K M1 (collimation)'))
    for which,dofs,st,lab in [(None,[],0,'none'),
                              ('M1',[0,1,2,3,4],1.0,'M1, 5 DOF, 1 um / 0.1 urad'),('M2',[0,1,2,3,4],1.0,'M2, 5 DOF, 1 um / 0.1 urad'),
                              ('M1',[2,3,4],1.0,'M1, dz rx ry only'),('M2',[2,3,4],1.0,'M2, dz rx ry only'),
                              ('M1',[0,1,2,3,4],10.,'M1, 5 DOF, 10 um / 1 urad'),('M2',[0,1,2,3,4],10.,'M2, 5 DOF, 10 um / 1 urad')]:
        a,b,c,e=run(seg,which,dofs,st)
        print('  %-34s %-26s %-26s %-26s %s'%(lab,a,b,c,e))

def run2(seg,rows):
    """rows: list of (mirror, dof, sigma)"""
    o=d[seg]; G1,G2,P1,P2=o['G1'],o['G2'],o['P1'],o['P2']; F=np.tile(np.eye(5),(3,1))
    A=np.r_[np.c_[G1,G2,F],np.c_[P1,P2,np.zeros((2,5))]]; s2=list(np.r_[np.full(15,SN),np.full(2,SL)])
    for m,k,sg in rows:
        r=np.zeros(15); r[(0 if m=='M1' else 5)+k]=1; A=np.r_[A,r[None]]; s2.append(sg)
    Aw=A/np.array(s2)[:,None]; prior=np.array([50,50,50,10,10]*2+[10]*5,float)
    C=np.linalg.inv(Aw.T@Aw+np.diag(1/prior**2)); sd=np.sqrt(np.diag(C))
    K=-np.linalg.pinv(G2)@G1; J=np.c_[-K,np.eye(5),np.zeros((5,5))]
    return sd[10:],np.sqrt(np.diag(J@C@J.T))
full=lambda m,st: [(m,k,s) for k,s in zip(range(5),[st,st,st,st/10,st/10])]
print('\n--- combinations')
for seg in (0,6):
    for lab,rows in [('M1 5 DOF + M2 dz (1 um)',full('M1',1.)+[('M2',2,1.)]),
                     ('M1 5 DOF + M2 dz (0.1 um)',full('M1',1.)+[('M2',2,.1)]),
                     ('M2 5 DOF + M1 dz (0.1 um)',full('M2',1.)+[('M1',2,.1)]),
                     ('both 5 DOF (1 um)',full('M1',1.)+full('M2',1.))]:
        f,e=run2(seg,rows); print('seg%d %-28s figure %s   M2-K M1 %s'%(seg+1,lab,f,e))
