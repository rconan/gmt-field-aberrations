import numpy as np, truth as T, general as G
from invert import NAMES, LAB
from rankcheck import design_p, SEGS
d=np.load('sec7.npy',allow_pickle=True).item(); wt=d['wt_new']
TL=T.basis(9,12,5); w=np.array([wt[t] for t in TL]); tru=np.array([wt[t] for t in NAMES])
def sel(k,l,m): return [t for t in TL if 2*t[0]+t[2]<=k and 2*t[1]+t[2]<=l and t[2]<=m]
cases=[('5-8',slice(4,8)),('1-10',slice(0,10))]
res={}
for K in [(5,6,3),(5,8,3),(7,8,4)]:
    Tk=sel(*K)
    for pname,pj in cases:
        for nmax in (3,5):
            b=design_p(TL,SEGS,nmax,pj)@w
            A=design_p(Tk,SEGS,nmax,pj)
            sv=np.linalg.svd(A,compute_uv=False); r=int((sv>1e-9*sv[0]).sum())
            U,s,Vt=np.linalg.svd(A,full_matrices=False); lev=(Vt[:r]**2).sum(0)
            x,*_=np.linalg.lstsq(A,b,rcond=1e-12)
            est=np.array([x[Tk.index(t)] for t in NAMES]); ident=[lev[Tk.index(t)]>1-1e-6 for t in NAMES]
            res[(K,pname,nmax)]=(len(Tk),r,est,ident)
            err=(est-tru)/tru
            print('terms k<=%d l<=%d m<=%d (%3d) pupil %-4s n<=%d rank %3d | '%(*K,len(Tk),pname,nmax,r)+' '.join('%s %s%.1e'%(LAB(t)[1:],'' if ok else '*',abs(e)) for t,e,ok in zip(NAMES,err,ident)))
np.save('pupil10.npy',{'res':res,'tru':tru},allow_pickle=True)
