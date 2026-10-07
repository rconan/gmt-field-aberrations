"""Differential scheme: Noll 4-8 at 3 field points referenced to the aligned telescope, differenced
between field points (removes everything field-constant, i.e. M1 figure), M2 from the differences.
Then the field-constant residual is corrected with M1 bending (Noll 4-8)."""
import numpy as np
from fullscheme import WLOS
d=np.load('figerr.npy',allow_pickle=True).item()
np.set_printoptions(precision=4,suppress=True,linewidth=160)
PAIRS=np.array([[1,-1,0],[0,1,-1],[-1,0,1]])          # 1-2, 2-3, 3-1
D=np.kron(PAIRS,np.eye(5))                            # 15x15 on (point, mode) ordering
M=np.kron(np.eye(3)-1/3,np.eye(5))                    # mean removal
print('D^T D = 3 x mean removal:',np.allclose(D.T@D,3*M),' rank of pair differences:',np.linalg.matrix_rank(D))
names=['dx','dy','dz','rx','ry']
for seg in (0,6):
    o=d[seg]; G1,G2=o['G1'],o['G2']; F=np.tile(np.eye(5),(3,1)); FV=np.tile(np.eye(5),(20,1))
    U,s,Vt=np.linalg.svd(D@G2); print('\nseg%d singular values of the differential M2 sensitivities (nm per um/urad):'%(seg+1),s)
    for k in range(5): print('   s=%.4f  '%s[k]+' '.join('%s %+.2f'%(n,v) for n,v in zip(names,Vt[k]) if abs(v)>0.05))
    # projection of the M1 signature on the differences: M1 field-dependent signature that M2 cannot cancel
    DG1=D@G1; R=DG1-D@G2@np.linalg.pinv(D@G2,rcond=1e-6)@DG1
    print('   M1 differential signature M2 cannot cancel (nm per um/urad, per DOF):',np.linalg.norm(R,axis=0)/np.sqrt(15))
    # Monte Carlo
    rng=np.random.default_rng(2); scale=np.array([50,50,50,10,10.]); res={}
    for t in range(500):
        x1=rng.normal(size=5)*scale; x2=rng.normal(size=5)*scale; f=rng.normal(size=5)*10; n=rng.normal(size=15)*0.2
        b=G1@x1+G2@x2+F@f+n
        w0=o['V1']@x1+o['V2']@x2+FV@f
        xc=-np.linalg.pinv(G2)@(G1@x1)                 # collimated M2 position (no figure): x2+dx2 = xc
        for rc in (1e-6,0.05):
            Gp=np.linalg.pinv(D@G2,rcond=rc)
            dx=-Gp@(D@b)
            c=(b+G2@dx).reshape(3,5).mean(0)          # field-constant residual -> M1 bending
            w=w0+o['V2']@dx-FV@c
            k='rcond %g'%rc
            res.setdefault(k,[]).append([np.sqrt((w**2).reshape(-1,5).sum(1).mean()),
                np.sqrt((((w0+o['V2']@dx).reshape(-1,5)-(w0+o['V2']@dx).reshape(-1,5).mean(0))**2).sum(1).mean()),
                np.hypot(*(o['P1']@x1+o['P2']@(x2+dx))),
                *np.abs(x2+dx-xc), np.linalg.norm(c-f)])
    print('   500 trials (M1 & M2 misaligned 50 um/10 urad, figure 10 nm rms per Noll 4-8 mode, 0.2 nm noise); rms of:')
    print('   %-10s %8s %10s %9s  %s %13s'%('','WFE nm','fld-dep nm','LOS mas','|M2 - collimated| dx dy dz rx ry','figure err nm'))
    for k,v in res.items():
        v=np.sqrt((np.array(v)**2).mean(0)); print('   %-10s %8.2f %10.2f %9.0f  %s %10.1f'%(k,v[0],v[1],v[2],np.array2string(v[3:8],precision=1),v[8]))
