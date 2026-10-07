import numpy as np, truth as T
import general as G
from invert import NAMES, LAB
S=8.71/4.1825; PHI=[(3-2*k)*np.pi/6 for k in range(6)]
SEGS=[(S*np.cos(a),S*np.sin(a)) for a in PHI]+[(0.,0.)]
TL=T.basis(9,12,5); T31=[t for t in TL if 2*t[0]+t[2]<=5 and 2*t[1]+t[2]<=6 and t[2]<=3]
def design_p(terms,segs,nmax,pj):
    FJ,FB=G.basis(nmax); cols=[]
    for (p,n,m) in terms:
        blk=[]
        for sx,sy in segs:
            Hx=G.X[:,None];Hy=G.Yc[:,None];Px=G.X[None,:]+sx;Py=G.Yc[None,:]+sy
            Wv=(Hx**2+Hy**2)**p*(Px**2+Py**2)**n*(Hx*Px+Hy*Py)**m
            blk.append((G.PB[pj]@Wv.T@FB.T).ravel())
        cols.append(np.concatenate(blk))
    return np.array(cols).T
def lev(A):
    U,sv,Vt=np.linalg.svd(A,full_matrices=False); r=(sv>1e-9*sv[0]).sum(); return r,(Vt[:r]**2).sum(0)
for pj,name in ((slice(4,8),'pupil 5-8'),(slice(0,10),'pupil 1-10')):
    for nmax in (3,5):
        r,l=lev(design_p(T31,SEGS,nmax,pj))
        main=[l[T31.index(t)] for t in NAMES]
        print('%-10s field n<=%d: rank %2d/31; main 11 identifiable: %s (min leverage %.6f)'%(name,nmax,r,all(x>1-1e-6 for x in main),min(main)))
        if nmax==3 and pj.start==4: print('   not identifiable at n<=3:',[LAB(t) for t in NAMES if l[T31.index(t)]<1-1e-6])
