import numpy as np; from core import *
exec(open('check2.py').read().split('for seg in')[0])
TL=[tuple(t) for t in d['terms']]
NAMES=[(0,2,0),(0,1,1),(0,0,2),(1,2,0),(0,1,2),(1,1,1),(0,3,0),(0,2,1),(1,3,0),(0,2,2),(1,2,1)]
sub=lambda names: Model(names, d['w'][:,[TL.index(t) for t in names]])
mods={'3rd':sub(NAMES[:3]),'11':sub(NAMES),'full':full}
res={}
for seg in (0,6):
    bn=raytrace_b(seg,H)
    for name,kw in cases:
        kw={k:dict(v,C=(m2_footprint(seg) if k=='m2' else m1pivot(seg))) for k,v in kw.items()}
        bm=raytrace_b(seg,H,**kw)
        out=[]
        for mn,M in mods.items():
            a,dw,rms,J=estimate_ext(M,seg,H,bm,bn) if seg<6 else (*estimate(M,seg,H,bm,bn)[:1],[0],estimate(M,seg,H,bm,bn)[1],None)
            out.append((mn,a,rms))
        af=out[-1][1]
        print('seg%d %-13s'%(seg+1,name),' '.join('%s: err %.1f%% rms %.2f'%(mn,100*np.abs(a-af).max()/np.abs(af).max(),rms) for mn,a,rms in out))
