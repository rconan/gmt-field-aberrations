"""Differential scheme vs measurement noise: M2 dz excluded (unobservable), tilts and decentres from
the pair differences; field-constant residual -> M1 bending. Linear model, 20 validation points."""
import numpy as np
d=np.load('figerr.npy',allow_pickle=True).item()
D=np.kron(np.array([[1,-1,0],[0,1,-1],[-1,0,1]]),np.eye(5))
for seg in (0,6):
    o=d[seg]; G1,G2=o['G1'],o['G2']; F=np.tile(np.eye(5),(3,1)); FV=np.tile(np.eye(5),(20,1))
    cols=[0,1,3,4]; Gp=np.linalg.pinv((D@G2)[:,cols])
    print('seg%d   noise   WFE   |M2-coll| dx  dy   rx   ry  (um/urad)   bending applied   true figure   bending - figure (nm rms/mode)'%(seg+1))
    for sig in (0,0.01,0.05,0.2,1.0):
        rng=np.random.default_rng(4); out=[]
        for t in range(500):
            x1=rng.normal(size=5)*[50,50,50,10,10]; x2=rng.normal(size=5)*[50,50,50,10,10]; f=rng.normal(size=5)*10
            b=G1@x1+G2@x2+F@f+rng.normal(size=15)*sig
            xc=-np.linalg.pinv(G2)@(G1@x1)
            dx=np.zeros(5); dx[cols]=-Gp@(D@b)
            # dz: assume known from elsewhere (set to the collimated value) so the comparison isolates the differential part
            dx[2]=xc[2]-x2[2]
            c=(b+G2@dx).reshape(3,5).mean(0)
            w=o['V1']@x1+o['V2']@(x2+dx)+FV@(f-c)
            out.append([np.sqrt((w**2).reshape(-1,5).sum(1).mean()),*np.abs(x2+dx-xc)[cols],np.sqrt((c**2).mean()),np.sqrt((f**2).mean()),np.sqrt(((c-f)**2).mean())])
        v=np.sqrt((np.array(out)**2).mean(0))
        print('      %5.2f  %5.2f   %5.1f %5.1f %5.2f %5.2f          %7.1f        %7.1f       %7.1f'%(sig,*v))
