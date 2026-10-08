"""Posterior uncertainties of the figure error and of the collimation state with Shack-Hartmann
centroids (48x48, 3 probes) and segment-position sensors; compare with Noll 4-8 + line of sight.
Priors: misalignments 50 um / 10 urad, figure 10 nm per Noll 4-8 mode."""
import numpy as np, sys
sh=np.load('sh.npy',allow_pickle=True).item(); fe=np.load('figerr.npy',allow_pickle=True).item()
SIG_C=float(sys.argv[1]) if len(sys.argv)>1 else 0.5
SIG_Z,SIG_L=0.2,1.0
PRIOR=np.array([50,50,50,10,10]*2+[10]*5,float)
def post(seg,kind,rows):
    o=fe[seg]
    if kind=='SH': s=sh[seg]; A=np.c_[s['G1'],s['G2'],s['F']]; sig=np.full(A.shape[0],SIG_C)
    else:
        A=np.r_[np.c_[o['G1'],o['G2'],np.tile(np.eye(5),(3,1))],np.c_[o['P1'],o['P2'],np.zeros((2,5))]]
        sig=np.r_[np.full(15,SIG_Z),np.full(2,SIG_L)]
    for m,k,sg in rows:
        r=np.zeros(15); r[(0 if m=='M1' else 5)+k]=1; A=np.r_[A,r[None]]; sig=np.r_[sig,sg]
    Aw=A/sig[:,None]; C=np.linalg.inv(Aw.T@Aw+np.diag(1/PRIOR**2))
    K=-np.linalg.pinv(o['G2'])@o['G1']; L=np.c_[-K,np.eye(5),np.zeros((5,5))]
    sd=np.sqrt(np.diag(C)); sc=np.sqrt(np.diag(L@C@L.T))
    return sd[10:],sc
full=lambda m,st: [(m,k,s) for k,s in zip(range(5),[st,st,st,st/10,st/10])]
SETS=[('none',[]),('M2, 5 DOFs',full('M2',1.)),('M1, 5 DOFs',full('M1',1.)),
      ('M1, piston and tip-tilt',[('M1',2,1.),('M1',3,.1),('M1',4,.1)]),('M1 and M2, 5 DOFs',full('M1',1.)+full('M2',1.))]
print('centroid noise %.2f mas; columns: seg1 astig, seg1 coma, seg1 M2 decentre | seg7 coma, seg7 M2 decentre'%SIG_C)
for kind in ('Zernike','SH'):
    print(kind)
    for lab,rows in SETS:
        f1,c1=post(0,kind,rows); f7,c7=post(6,kind,rows)
        print('   %-26s %5.1f %5.2f %5.1f | %5.2f %5.1f'%(lab,f1[1],f1[3],c1[0],f7[3],c7[0]))
print('M1 piston + tip-tilt sensors, tilt uncertainty sweep')
for kind in ('Zernike','SH'):
    print(kind)
    for st in (0.03,0.1,0.3,1,3):
        rows=[('M1',3,st),('M1',4,st),('M1',2,st*4.18)]
        f1,c1=post(0,kind,rows); f7,c7=post(6,kind,rows)
        print('   %5.2f urad  %5.1f %5.2f %5.1f | %5.2f %5.1f'%(st,f1[1],f1[3],c1[0],f7[3],c7[0]))
