"""Strategies for M1 figure errors with Shack-Hartmann centroids (48x48, 3 probes) compared with
Noll 4-8 coefficients + line of sight. Linear model; residuals at the 20 validation points.
Misalignments of M1 and M2: rms 50 um / 10 urad; figure: 10 nm rms per Noll 4-8 mode."""
import numpy as np, sys
from fullscheme import WLOS
sh=np.load('sh.npy',allow_pickle=True).item(); fe=np.load('figerr.npy',allow_pickle=True).item()
SIG_Z, SIG_L = 0.2, 1.0
SIG_C = float(sys.argv[1]) if len(sys.argv)>1 else 0.5        # mas per centroid
SIG_CAL = float(sys.argv[2]) if len(sys.argv)>2 else 0.0      # nm, static Noll 4-8 offsets per probe
SF = 4.0                                                      # figure prior (nm) in the joint solution
NT=300
def map_solve(A, y, sig, prior):
    Aw=A/sig[:,None]; return -np.linalg.solve(Aw.T@Aw+np.diag(1/prior**2), Aw.T@(y/sig))
def run(seg):
    o=fe[seg]; s=sh[seg]
    V1,V2,P1,P2=o['V1'],o['V2'],o['P1'],o['P2']; FV=np.tile(np.eye(5),(20,1))
    nC=s['G1'].shape[0]; n1=nC//3
    Fp=[np.zeros((nC,5)) for k in range(3)]
    for k in range(3): Fp[k][k*n1:(k+1)*n1]=s['F'][k*n1:(k+1)*n1]
    meas={'Zernike + LOS':dict(G1=np.r_[o['G1'],o['P1']],G2=np.r_[o['G2'],o['P2']],F=np.r_[np.tile(np.eye(5),(3,1)),np.zeros((2,5))],
                               sig=np.r_[np.full(15,SIG_Z),np.full(2,SIG_L)]),
          'Shack-Hartmann':dict(G1=s['G1'],G2=s['G2'],F=s['F'],sig=np.full(nC,SIG_C))}
    rng=np.random.default_rng(8); res={}
    for t in range(NT):
        x1=rng.normal(size=5)*[50,50,50,10,10]; x2=rng.normal(size=5)*[50,50,50,10,10]; f=rng.normal(size=5)*10
        cal=rng.normal(size=(3,5))*SIG_CAL
        xc=-np.linalg.pinv(o['G2'])@(o['G1']@x1)              # collimated M2 (Noll 4-8, no figure)
        w0=V1@x1+V2@x2+FV@f
        for mk,m in meas.items():
            G1,G2,F,sig=m['G1'],m['G2'],m['F'],m['sig']
            y=G1@x1+G2@x2+F@f+rng.normal(size=len(sig))*sig
            if mk=='Shack-Hartmann': y=y+sum(Fp[k]@cal[k] for k in range(3))
            else: y=y+np.r_[cal.ravel(),0,0]
            A7=np.c_[G2,G1[:,[3,4]]]
            def ev(name,d2,r1,fc):
                x1c=x1.copy(); x1c[3:]+=r1; xc=-np.linalg.pinv(o['G2'])@(o['G1']@x1c)   # collimated M2 for the corrected M1
                w=w0+V2@d2+V1[:,[3,4]]@r1-FV@fc
                los=P1@x1+P2@(x2+d2)+P1[:,[3,4]]@r1
                res.setdefault((mk,name),[]).append([np.sqrt((w**2).reshape(-1,5).sum(1).mean()),np.hypot(*los),
                    np.hypot(*(x2+d2-xc)[:2]),np.sqrt(((fc-f)**2).mean())])
            z=map_solve(A7,y,sig,np.full(7,1e4)); ev('figure ignored',z[:5],z[5:],np.zeros(5))
            r=y+A7@z; fb=-map_solve(F,r,sig,np.full(5,1e4)); ev('then M1 bending on the residual',z[:5],z[5:],fb)
            z=map_solve(np.c_[A7,F],y,sig,np.r_[np.full(7,1e4),np.full(5,SF)]); ev('joint, prior on bending',z[:5],z[5:7],-z[7:])
    return res
print('centroid noise %.2f mas, probe calibration offsets %.2f nm rms per mode; %d trials'%(SIG_C,SIG_CAL,NT))
for seg in (0,6):
    res=run(seg); print('segment %d'%(seg+1))
    print('   %-16s %-46s %8s %8s %12s %12s'%('measurements','strategy','WFE nm','LOS mas','M2 dec. um','bend.-fig nm'))
    for (mk,name),v in res.items():
        v=np.sqrt((np.array(v)**2).mean(0)); print('   %-16s %-46s %8.2f %8.1f %12.1f %12.1f'%(mk,name,*v))
