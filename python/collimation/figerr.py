"""M1 segment figure error (Noll 4-8, field-constant) in the collimation scheme."""
import numpy as np
import sens as S
from collimate import modes, HV
from fullscheme import P, WLOS
G=np.load('G48.npy',allow_pickle=True).item()
np.set_printoptions(precision=3,suppress=True,linewidth=160)
NOLL=[4,5,6,7,8]
def GV(seg,mir):     # sensitivities at the 20 validation points, Noll 4-8 (nm per um/urad)
    k='x1' if mir=='m1' else 'x2'; cols=[]
    for i in range(5):
        e=np.eye(5)[i]*S.UNIT
        cols.append(((modes(seg,HV,**{k:e})-modes(seg,HV,**{k:-e}))[:,:5]/2).ravel())
    return np.array(cols).T
out={}
for seg in (0,6):
    G1,G2=G[(seg,'Noll 4-8')]; P1,P2=P(seg,'m1'),P(seg,'m2'); V1,V2=GV(seg,'m1'),GV(seg,'m2')
    F=np.tile(np.eye(5),(3,1)); FV=np.tile(np.eye(5),(len(HV),1))
    dof=[0,1,2,3,4] if seg<6 else [0,1,2,3,4]
    A=np.r_[np.c_[G2,G1[:,[3,4]]],WLOS*np.c_[P2,P1[:,[3,4]]]]
    Am=G2
    print('\n=== segment %d'%(seg+1))
    # 1) unknown figure error: response of the schemes to 10 nm rms of each Noll mode
    for name,M,ext in (('M2 only',Am,False),('M2 + M1 tilts + LOS',A,True)):
        print('--',name,' (per 10 nm rms of M1-segment figure, wavefront)')
        for j in range(5):
            f=10*F[:,j]; r=np.r_[f,np.zeros(2)] if ext else f
            x=-np.linalg.pinv(M,rcond=1e-10)@r
            res=10*FV[:,j]+V2@x[:5]+(V1[:,[3,4]]@x[5:] if ext else 0)
            los=P2@x[:5]+(P1[:,[3,4]]@x[5:] if ext else 0)
            R=res.reshape(-1,5); rms=lambda c: np.sqrt((R[:,c]**2).mean())
            print('  Z%d: M2 dx dy dz rx ry %s  M1 rx ry %s | residual focus %.2f astig %.2f coma %.2f nm (from 10) | LOS (%+.0f,%+.0f) mas'%(
                NOLL[j],x[:5],x[5:] if ext else '-',rms([0]),rms([1,2]),rms([3,4]),*los))
    # 2) figure error as 5 extra unknowns
    for name,M in (('M2 + figure',np.c_[G2,F]),('M2 + M1 tilts + LOS + figure',np.c_[A,np.r_[F,np.zeros((2,5))]])):
        n=M/np.linalg.norm(M,axis=0); s=np.linalg.svd(n,compute_uv=False)
        print('--',name,': %d unknowns, rank %d, singular values (normalised columns)'%(M.shape[1],np.linalg.matrix_rank(n,1e-8)),s)
    out[seg]=dict(G1=G1,G2=G2,P1=P1,P2=P2,V1=V1,V2=V2)
np.save('figerr.npy',out,allow_pickle=True)
