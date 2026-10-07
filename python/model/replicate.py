import numpy as np
from raytrace import trace
from emulate import disc, lump, noll10, wls_projector, FV, FT, FP
D_SEG=8.365
def pipeline(coords='sphere'):
    """returns (4 pupil j=5..8, 8 field j=1..8, 7 segments) in nm, like double_zernike"""
    out=np.zeros((4,8,7)); raw=[]
    for sid in range(1,8):
        if sid<7:
            o=(3-2*(sid-1))*np.pi/6; origin=(8.71*np.cos(o),8.71*np.sin(o))
        else: origin=(0.,0.)
        V,T=disc(D_SEG,0.25,origin); w=lump(V,T)
        ci=np.where((V[:,0]==origin[0])&(V[:,1]==origin[1]))[0][0]
        if coords=='entrance':
            loc=V-V[ci]; P=wls_projector(noll10(loc[:,0],loc[:,1]),w)
        pup=[]
        for xy in FV:
            zen=np.deg2rad(np.hypot(*xy)/60); azi=np.arctan2(xy[1],xy[0])
            S,opd=trace(V,zen,azi)
            if coords=='sphere':
                loc=(S-S[ci])[:,:2]; P=wls_projector(noll10(loc[:,0],loc[:,1]),w)
            pup.append(P@opd)
        pup=np.array(pup)            # (Nf,10)
        dz=FP@pup                    # (10 field,10 pupil)
        out[:,:,sid-1]=dz[:8,4:8].T*1e9
    return out
