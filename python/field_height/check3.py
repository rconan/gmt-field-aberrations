import numpy as np; from core import *
exec(open('check2.py').read().split('for seg in')[0])
for seg in (0,):
    bn=raytrace_b(seg,H); bnv=raytrace_b(seg,HV)
    for name,kw in cases:
        kw={k:dict(v,C=(m2_footprint(seg) if k=='m2' else m1pivot(seg))) for k,v in kw.items()}
        bm=raytrace_b(seg,H,**kw); bmv=raytrace_b(seg,HV,**kw)
        a,dw,rms,J=estimate_ext(full,seg,H,bm,bn)
        sx,sy=SEG[seg]; E=PROJ@(((QX+sx)**2+(QY+sy)**2)**2)
        pred=full.bfield(seg,HV,a)-full.bfield(seg,HV,np.zeros((2,2)))+bnv+dw[0]*E
        sig=np.sqrt(((bmv-bnv)**2).mean())
        print('seg%d %-13s a1=(%+.2e,%+.2e) a2=(%+.2e,%+.2e) dw040 %+.3f fit %.1e | signal %.2f, predict err %.2e  cond %.1e'%(
            seg+1,name,*a.ravel(),dw[0],rms,sig,np.sqrt(((pred-bmv)**2).mean()),np.linalg.cond(J)))
