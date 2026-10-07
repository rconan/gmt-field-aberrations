import sys as _sys, pathlib as _pl; _P=_pl.Path(__file__).resolve().parents[1]; _sys.path[:0]=[str(_P/'model'),str(_P/'field_height'),str(_P/'rigid_body')]
import numpy as np, sys
import truth as T
from invert import cascade, NAMES, LAB
from general import design
from sec6 import response
S=8.71/4.1825; PHI=[(3-2*k)*np.pi/6 for k in range(6)]
SEGS=[(S*np.cos(a),S*np.sin(a)) for a in PHI]+[(0.,0.)]
TL=T.basis(9,12,5)
def truthfit(O):
    idx=np.arange(0,T.rho.shape[0],4); r=T.rho[idx]; y=O[:,idx].ravel()
    zx,zy=T.zeta[:,0][:,None],T.zeta[:,1][:,None]; rx,ry=r[:,0][None,:],r[:,1][None,:]
    M=np.empty((y.size,len(TL)))
    for c,(p,n,m) in enumerate(TL): M[:,c]=((zx**2+zy**2)**p*(rx**2+ry**2)**n*(zx*rx+zy*ry)**m).ravel()
    sc=np.abs(M).max(0); w,*_=np.linalg.lstsq(M/sc,y,rcond=None); return w/sc
def evaluate(tag,O,dz):
    w=truthfit(O); wt={t:v for t,v in zip(TL,w)}
    tru=np.array([wt[t] for t in NAMES])
    cas=cascade(dz); c=np.array([cas[LAB(t)].real for t in NAMES])
    T31=[t for t in TL if 2*t[0]+t[2]<=5 and 2*t[1]+t[2]<=6 and t[2]<=3]
    bfull=design(TL,SEGS,5)@w
    A31=design(T31,SEGS,5); f31,*_=np.linalg.lstsq(A31,bfull,rcond=None)
    f=np.array([f31[T31.index(t)] for t in NAMES])
    rank3=np.linalg.matrix_rank(design(T31,SEGS,3),tol=1e-9*np.linalg.norm(design(T31,SEGS,3),2))
    # synthetic (mesh-based) contribution of w531 to the cascade
    w531=wt[(1,1,3)] if (1,1,3) in wt else None
    print('\n==',tag,'  (rank n_h<=3: %d of %d)'%(rank3,len(T31)))
    print('%-6s %12s %12s %10s %12s %10s'%('term','direct fit','cascade','rel err','field n<=5','rel err'))
    for i,t in enumerate(NAMES):
        print('%-6s %12.4f %12.4f %10.1e %12.4f %10.1e'%(LAB(t),tru[i],c[i],(c[i]-tru[i])/tru[i],f[i],(f[i]-tru[i])/tru[i]))
    return wt,tru,c,f
old=evaluate('OLD (CEO reference)',np.load(str(_P/'model'/'old_reference'/'truth_opd.npy')),np.load(str(_P/'model'/'old_reference'/'replica_entrance.npy')))
new=evaluate('NEW (delrays reference)',np.load('truth_opd.npy'),np.load('replica_entrance.npy'))
np.save('sec7.npy',{'old':old[1:],'new':new[1:],'wt_new':new[0],'wt_old':old[0]},allow_pickle=True)
