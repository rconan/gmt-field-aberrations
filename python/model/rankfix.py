import numpy as np, truth as T
from invert import NAMES, LAB
from rankcheck import design_p, SEGS
d=np.load('sec7.npy',allow_pickle=True).item(); wt=d['wt_new']
TL=T.basis(9,12,5); w=np.array([wt[t] for t in TL]); tru=np.array([wt[t] for t in NAMES])
def sel(k,l,m): return [t for t in TL if 2*t[0]+t[2]<=k and 2*t[1]+t[2]<=l and t[2]<=m]
def analyse(Tk,pj,nmax):
    A=design_p(Tk,SEGS,nmax,pj); b=design_p(TL,SEGS,nmax,pj)@w
    c=np.linalg.norm(A,axis=0); live=c>1e-10*c.max()      # structurally zero columns: unobservable
    An=A[:,live]/c[live]
    U,s,Vt=np.linalg.svd(An,full_matrices=False); r=int((s>1e-10*s[0]).sum())
    lev=np.zeros(len(Tk)); lev[live]=(Vt[:r]**2).sum(0)
    x=np.zeros(len(Tk)); x[live]=(Vt[:r].T@((U[:,:r].T@b)/s[:r]))/c[live]
    gap=s[r]/s[r-1] if r<len(s) else 0.
    return r,int(live.sum()),x,lev,1/(s[r-1]/s[0]),gap
if __name__=='__main__':
    for K,pj,nm in [((5,6,3),[4,5,6,7],3),((5,6,3),list(range(10)),3),((5,6,3),[4,5,6,7],5),((5,6,3),list(range(10)),5),
                    ((5,8,3),[4,5,6,7],5),((5,8,3),[3,4,5,6,7],5),((5,8,3),[4,5,6,7,8,9],5),((5,8,3),list(range(10)),5),
                    ((7,8,4),list(range(10)),5),((7,8,4),list(range(10)),7)]:
        Tk=sel(*K); r,nl,x,lev,cond,gap=analyse(Tk,pj,nm)
        est=np.array([x[Tk.index(t)] for t in NAMES]); ok=[lev[Tk.index(t)]>1-1e-6 for t in NAMES]
        print('%-10s %-31s n<=%d %2d terms (%2d observable) rank %2d cond %.0e gap %.0e | '%(K,[p+1 for p in pj],nm,len(Tk),nl,r,cond,gap)+' '.join('%s %s%.1e'%(LAB(t)[1:],'' if o else '*',abs((e-q)/q)) for t,e,q,o in zip(NAMES,est,tru,ok)))
