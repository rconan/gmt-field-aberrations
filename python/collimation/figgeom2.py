import numpy as np
d=np.load('figgeom.npy',allow_pickle=True).item(); GRID=d['grid']
def idx(r,az):
    if r==0: return 0
    return 1+[0.25,0.5,0.75,1.0].index(r)*12+int(round(az/30))%12
def mats(seg,pts):
    G1=np.concatenate([d[(seg,'m1')][:,p,:].T for p in pts]); G2=np.concatenate([d[(seg,'m2')][:,p,:].T for p in pts])
    F=np.tile(np.eye(5),(len(pts),1)); return G1,G2,F
def metrics(seg,pts):
    G1,G2,F=mats(seg,pts)
    s2=np.linalg.svd(np.c_[G2,F],compute_uv=False)
    G=np.c_[G1,G2]; U,s,_=np.linalg.svd(G,full_matrices=False); U=U[:,s>1e-6*s[0]]
    Q=F-U@(U.T@F); m=np.linalg.svd(Q,compute_uv=False)/np.sqrt(len(pts))   # unexplained fraction of a figure error
    return s2[-3:],m
CONF={
 'equilateral r=0.25':[(0.25,90),(0.25,210),(0.25,330)],
 'equilateral r=0.5':[(0.5,90),(0.5,210),(0.5,330)],
 'equilateral r=0.75':[(0.75,90),(0.75,210),(0.75,330)],
 'equilateral r=1':[(1,90),(1,210),(1,330)],
 'equilateral r=1 rot 30':[(1,0),(1,120),(1,240)],
 'centre + 2 at r=1, 90 deg':[(0,0),(1,0),(1,90)],
 'centre + 2 at r=1, 180 deg (collinear)':[(0,0),(1,0),(1,180)],
 'collinear r=0.5, r=1 at 0 deg + centre':[(0,0),(0.5,0),(1,0)],
 'small triangle at edge (r=1, 0/30/60 deg)':[(1,0),(1,30),(1,60)],
 'isosceles r=1 at 0/90/180':[(1,0),(1,90),(1,180)],
}
np.set_printoptions(precision=4,suppress=True,linewidth=170)
for seg in (0,6):
    print('segment',seg+1)
    print('  %-44s %-30s %s'%('configuration','3 smallest sv of [G2 F]','fraction of figure error not mimicked by M1+M2 motion (Noll-5 directions)'))
    for k,v in CONF.items():
        s,m=metrics(seg,[idx(*p) for p in v]); print('  %-44s %-30s %s'%(k,np.array2string(s,precision=4),np.array2string(m,precision=4)))
