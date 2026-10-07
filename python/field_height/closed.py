import numpy as np; from core import *
def closed_form(seg,H,db,w1,w2):
    """third-order closed form. db: (3,4) measured-minus-nominal Noll 5..8; w1,w2: dicts of per-surface
    w131, w222. Returns alpha (2,2) [x,y] and dw040."""
    dA,dC=to_complex(db); eta=H[:,0]-1j*H[:,1]
    sx,sy=SEG[seg]; sig=sx-1j*sy
    V=np.c_[eta,np.ones_like(eta)]; (a1,a0),*_=np.linalg.lstsq(V,dA,rcond=None)
    Cm=dC.mean()
    if abs(sig)>0:
        delta=(-3*(a0-4*sig*Cm)/sig**2); dl=delta.real
    else: delta=0j; dl=0.
    b222=-6*a1; b131=-24*Cm+4*sig*dl
    M=np.array([[w1['w222'],w2['w222']],[w1['w131'],w2['w131']]])
    al=np.linalg.solve(M,np.array([b222,b131]))
    # conj convention alpha = ax - i ay
    return np.array([[a.real,-a.imag] for a in al]), dl, delta.imag
if __name__=='__main__':
    exec(open('check4.py').read().split('res={}')[0])
    i131,i222=TL.index((0,1,1)),TL.index((0,0,2))
    w1={'w131':d['w'][0][i131],'w222':d['w'][0][i222]}; w2={'w131':d['w'][1][i131],'w222':d['w'][1][i222]}
    for seg in (0,6):
        bn=raytrace_b(seg,H)
        for name,kw in cases[:4]:
            kw={k:dict(v,C=(m2_footprint(seg) if k=='m2' else m1pivot(seg))) for k,v in kw.items()}
            bm=raytrace_b(seg,H,**kw)
            a,dl,chk=closed_form(seg,H,bm-bn,w1,w2)
            a3=(estimate_ext(mods['3rd'],seg,H,bm,bn)[0] if seg<6 else estimate(mods['3rd'],seg,H,bm,bn)[0])
            print(seg+1,name,(a*600e3).round(2).ravel(),(a3*600e3).round(2).ravel(),'%.3f %.3f'%(dl,chk))
