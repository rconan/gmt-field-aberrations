"""Misaligned M1 + M1 figure error, with measurement noise. Linear model, 20 validation points.
Unknowns: M2 (5), M1 rx ry, figure Noll 4-8 corrected by M1 bending. The near-degenerate
directions are resolved by a Tikhonov prior on either the figure correction or the M2 motion."""
import numpy as np
from fullscheme import WLOS
d=np.load('figerr.npy',allow_pickle=True).item()
rng=np.random.default_rng(3); scale=np.array([50,50,50,10,10.]); NOISE=0.2; NT=200
for seg in (0,6):
    o=d[seg]; F=np.tile(np.eye(5),(3,1)); FV=np.tile(np.eye(5),(20,1))
    A=np.r_[np.c_[o['G2'],o['G1'][:,[3,4]]],WLOS*np.c_[o['P2'],o['P1'][:,[3,4]]]]
    Af=np.c_[A,np.r_[F,np.zeros((2,5))]]
    def tik(M,b,wd):
        return -np.linalg.solve(M.T@M+np.diag(wd**2),M.T@b)
    res={}
    for t in range(NT):
        x1=rng.normal(size=5)*scale; f=rng.normal(size=5)*10
        meas=np.r_[o['G1']@x1+F@f+rng.normal(size=15)*NOISE, WLOS*(o['P1']@x1+rng.normal(size=2)*1.0)]
        w0=o['V1']@x1+FV@f; los0=o['P1']@x1
        xc=-np.linalg.pinv(A)@np.r_[o['G1']@x1,WLOS*(o['P1']@x1)]          # collimated (no figure)
        def ev(x,fc):
            w=w0+o['V2']@x[:5]+o['V1'][:,[3,4]]@x[5:7]-FV@fc
            return np.sqrt((w**2).reshape(-1,5).sum(1).mean()), np.hypot(*(los0+o['P2']@x[:5]+o['P1'][:,[3,4]]@x[5:7])), np.hypot(*(x[:2]-xc[:2])), np.sqrt(((fc-f)**2).sum())
        cases={}
        x=-np.linalg.pinv(A)@meas; cases['M2+M1 tilts, figure ignored']=ev(x,np.zeros(5))
        cases['two steps: then M1 bending']=None
        # two steps: residual field mean after step 1, corrected with bending
        r=(meas[:15]+A[:15]@x).reshape(3,5).mean(0); cases['two steps: then M1 bending']=ev(x,r)
        z=tik(Af,meas,np.r_[np.full(7,1e-3),np.full(5,0.05)]); cases['joint, prior on bending']=ev(z[:7],-z[7:])
        z=tik(Af,meas,np.r_[1.,1.,1e-3,1e-3,1e-3,1e-3,1e-3,np.full(5,1e-3)]); cases['joint, prior on M2 decentre']=ev(z[:7],-z[7:])
        for k,v in cases.items(): res.setdefault(k,[]).append(v)
    print('seg%d (%d trials; misalignment 50 um / 10 urad, figure 10 nm rms/mode, noise %.1f nm, 1 mas)'%(seg+1,NT,NOISE))
    print('   %-30s %10s %10s %14s %14s'%('','WFE nm','LOS mas','M2 decentre','figure err nm'))
    for k,v in res.items():
        v=np.array(v); print('   %-30s %10.2f %10.1f %14.1f %14.1f'%(k,*np.sqrt((v**2).mean(0))))
