import numpy as np
from fullscheme import WLOS
d=np.load('figerr.npy',allow_pickle=True).item()
np.set_printoptions(precision=3,suppress=True,linewidth=160)
lab=['M2dx','M2dy','M2dz','M2rx','M2ry','M1rx','M1ry','Z4','Z5','Z6','Z7','Z8']
for seg in (0,6):
    o=d[seg]; F=np.tile(np.eye(5),(3,1))
    for name,M,l in (('M2+fig',np.c_[o['G2'],F],lab[:5]+lab[7:]),
                     ('full+fig',np.c_[np.r_[np.c_[o['G2'],o['G1'][:,[3,4]]],WLOS*np.c_[o['P2'],o['P1'][:,[3,4]]]],np.r_[F,np.zeros((2,5))]],lab)):
        U,s,Vt=np.linalg.svd(M,full_matrices=False)
        print('seg%d %s sing.values (raw units nm per um,urad,nm):'%(seg+1,name),s)
        for k in (-1,-2,-3):
            v=Vt[k]; print('   s=%.2e  '%s[k]+' '.join('%s %+.3f'%(a,b) for a,b in zip(l,v) if abs(b)>0.02))
        # noise: 1 nm rms per wavefront measurement, 1 mas LOS; drop the exact null direction
        Mp=np.linalg.pinv(M,rcond=1e-6)
        sd=np.sqrt(np.diag(Mp@Mp.T)); print('   std per 1 nm noise (pinv rcond 1e-6):',dict(zip(l,np.round(sd,2))))
