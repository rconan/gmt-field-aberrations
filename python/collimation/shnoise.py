"""Noise propagation from Shack-Hartmann centroids to segment Noll 2-8 (nm per mas of centroid noise),
with the subapertures of shmodel.py; also the mean slope (G-tilt) of coma."""
import numpy as np
from shmodel import geometry, slopes, R_SEG, zernike
for seg in (0,6):
    X,Y,sub,used,cx,cy=geometry(seg); u,v=(X-cx)/R_SEG,(Y-cy)/R_SEG
    cols=[slopes(2*u*1e-9,sub),slopes(2*v*1e-9,sub)]+[slopes(zernike(j,u,v)*1e-9,sub) for j in (4,5,6,7,8)]
    B=np.array(cols).T                      # mas per nm
    C=np.linalg.inv(B.T@B); sd=np.sqrt(np.diag(C)); corr=C/np.outer(sd,sd)
    print('segment %d (%d subapertures): nm per mas  tilt %.2f focus %.2f astig %.2f coma %.2f ; corr(tilt, coma) %.2f'%(
        seg+1,sub.sum(),sd[0],sd[2],sd[3],sd[5],corr[1,5]))
    n=sub.sum(); gt=lambda c: c[:n].mean(), 
    print('   mean slope of 10 nm of coma: %.2f mas (of 10 nm of Noll tilt: %.2f mas)'%(abs(B[n:,5].mean()*10),abs(B[n:,1].mean()*10)))
