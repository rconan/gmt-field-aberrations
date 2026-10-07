"""Misaligned M1 + M1 figure error (Noll 4-8). Linear model at the 20 validation points.
(i) scheme ignoring the figure (M2 + M1 tilts + LOS); (ii) figure modes added as unknowns and
corrected (M1 bending), degenerate directions dropped (truncated SVD)."""
import numpy as np
from fullscheme import WLOS
d=np.load('figerr.npy',allow_pickle=True).item()
rng=np.random.default_rng(11); scale=np.array([50,50,50,10,10.])
for seg in (0,6):
    o=d[seg]; F=np.tile(np.eye(5),(3,1)); FV=np.tile(np.eye(5),(20,1))
    A=np.r_[np.c_[o['G2'],o['G1'][:,[3,4]]],WLOS*np.c_[o['P2'],o['P1'][:,[3,4]]]]
    Af=np.c_[A,np.r_[F,np.zeros((2,5))]]
    for trial in range(3):
        x1=rng.normal(size=5)*scale; f=rng.normal(size=5)*10/np.sqrt(5)*np.sqrt(5)  # ~10 nm rms per mode
        meas=np.r_[o['G1']@x1+F@f, WLOS*(o['P1']@x1)]
        w0=o['V1']@x1+FV@f; los0=o['P1']@x1
        def apply(x,fc):
            w=w0+o['V2']@x[:5]+o['V1'][:,[3,4]]@x[5:7]-FV@fc
            return np.sqrt((w**2).reshape(-1,5).sum(1).mean()), los0+o['P2']@x[:5]+o['P1'][:,[3,4]]@x[5:7]
        r=[]
        x=-np.linalg.pinv(A)@meas; r.append(('ignore figure',)+apply(x,np.zeros(5))+(x,))
        U,s,Vt=np.linalg.svd(Af,full_matrices=False)
        for tol in (0.05,1e-3):
            k=s>tol*s[0]; xa=-(Vt[k].T/s[k])@(U[:,k].T@meas)
            r.append(('fig unknowns, s>%g'%tol,)+apply(xa[:7],-xa[7:])+(xa,))
        print('seg%d trial %d: before %.1f nm, figure %s nm'%(seg+1,trial,np.sqrt((w0**2).reshape(-1,5).sum(1).mean()),np.round(f,1)))
        for name,wr,lo,x in r:
            dev=x[:5]-(-np.linalg.pinv(A)@np.r_[o['G1']@x1,WLOS*(o['P1']@x1)])[:5]
            print('   %-22s residual %.2f nm  LOS (%+.1f,%+.1f)  M2 offset from collimated %s um/urad%s'%(name,wr,*lo,np.round(dev,1),
                  '  est. figure %s'%np.round(-x[7:],1) if len(x)>7 else ''))
