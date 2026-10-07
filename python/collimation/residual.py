import numpy as np
import sens as S
from collimate import modes
np.set_printoptions(precision=3,suppress=True,linewidth=150)
def Gm(seg,mir,sel):
    cols=[]
    for i in range(5):
        e=np.zeros(5); e[i]=S.UNIT; k='x1' if mir=='m1' else 'x2'
        cols.append(((modes(seg,S.H,**{k:e})-modes(seg,S.H,**{k:-e}))/2)[:,sel].ravel())
    return np.array(cols).T
names=['dx','dy','dz','rx','ry']
res={}
for seg in (0,6):
    for lab,sel in (('Noll 5-8',[1,2,3,4]),('Noll 4-8',[0,1,2,3,4])):
        G1,G2=Gm(seg,'m1',sel),Gm(seg,'m2',sel)
        sv2=np.linalg.svd(G2,compute_uv=False)
        Rm=G1-G2@np.linalg.pinv(G2,rcond=1e-8)@G1               # signature of M1 that M2 cannot cancel
        U,s,Vt=np.linalg.svd(Rm)
        print('segment %d, %s: singular values of G2 %s ; uncompensable M1 signature (nm rms per um/urad): %s'%(seg+1,lab,np.round(sv2,3),np.round(s/np.sqrt(len(sel)*3),4)))
        print('   leading uncompensable M1 motion:',' '.join('%s %+.3f'%(n,v) for n,v in zip(names,Vt[0])))
        res[(seg,lab)]=(G1,G2)
np.save('G48.npy',res,allow_pickle=True)
