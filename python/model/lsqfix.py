import numpy as np, truth as T, general as G
from rankcheck import design_p, SEGS
d=np.load('sec7.npy',allow_pickle=True).item(); wt=d['wt_new']
TL=T.basis(9,12,5); w=np.array([wt[t] for t in TL])
FJ,_=G.basis(5); PJ=G.jnm(3)
b=design_p(TL,SEGS,5,slice(0,10))@w            # (seg, pupil 10, field 21) flattened
b=b.reshape(7,10,21)
with open('double_zernike_lsq.txt','w') as f:
    f.write('# exact projections (nm) of the 136-term GMT wavefront, pupil Noll 1-10, field Noll 1-21, entrance layout\n')
    for ip,(j,n,m) in enumerate(PJ):
        f.write('(%d, %d, %d)\n'%(j,n,m))
        for jf,(jj,nn,mm) in enumerate(FJ):
            f.write(' (%2d, %2d, %2d): '%(jj,nn,mm)+' '.join('%+.12e'%v for v in b[:,ip,jf])+'\n')
def sel(k,l,m): return [t for t in TL if 2*t[0]+t[2]<=k and 2*t[1]+t[2]<=l and t[2]<=m]
out=[]
for K,pj in (((5,8,3),list(range(10))),((5,8,3),[3,4,5,6,7]),((5,8,3),[4,5,6,7]),((5,6,3),[4,5,6,7])):
    Tk=sel(*K); A=design_p(Tk,SEGS,5,pj); bb=design_p(TL,SEGS,5,pj)@w
    U,s,Vt=np.linalg.svd(A/np.linalg.norm(A,axis=0),full_matrices=False); r=int((s>1e-9*s[0]).sum()); lev=(Vt[:r]**2).sum(0)
    x,*_=np.linalg.lstsq(A,bb,rcond=1e-12)
    out.append((K,[p+1 for p in pj],len(Tk),r,[(t,x[i],lev[i]>1-1e-6) for i,t in enumerate(Tk)]))
    print(K,[p+1 for p in pj],len(Tk),'rank',r,'n_ident',sum(l>1-1e-6 for l in lev))
np.save('lsqfix.npy',out,allow_pickle=True)
# a design column for the Rust design test: term (0,1,1) w131, segment 1, pupil 8, field 2 / term (1,2,0), seg 7, pupil 4, field 4
A=design_p([(0,1,1),(1,2,0)],SEGS,5,slice(0,10)).reshape(7,10,21,2) if False else None
c1=design_p([(0,1,1)],SEGS,5,slice(0,10))[:,0].reshape(7,10,21); c2=design_p([(1,2,0)],SEGS,5,slice(0,10))[:,0].reshape(7,10,21)
print('D(w131; seg1, p8, f2) = %.15e'%c1[0,7,1]); print('D(w240; seg7, p4, f4) = %.15e'%c2[6,3,3]); print('D(w131; seg3, p7, f3) = %.15e'%c1[2,6,2])
