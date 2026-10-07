"""Conversion of real Noll double-Zernike coefficients (pupil j=5..8, field j=1..8) into the
complex coefficients a_{(n_r,m_r)(n_h,m_h)} of the paper (unnormalized, conjugate basis)."""
import numpy as np
# field Noll j -> (n, |m|, 'c'/'s'/None)
FIELD = {1:(0,0,None),2:(1,1,'c'),3:(1,1,'s'),4:(2,0,None),5:(2,2,'s'),6:(2,2,'c'),7:(3,1,'s'),8:(3,1,'c')}
def to_complex(b):
    """b: (4 pupil [5,6,7,8], 8 field [1..8], ...) -> dict {(nr,mr,nh,mh): complex array}"""
    out = {}
    for (nr,mr,jc,js) in [(2,2,6,5),(3,1,8,7)]:
        Nr = np.sqrt(2*(nr+1))
        c = (b[jc-5] - 1j*b[js-5]) / Nr            # complex pupil, real Noll field: (8,...)
        done=set()
        for jh,(nh,mh,cs) in FIELD.items():
            if cs is None:
                out[(nr,mr,nh,0)] = c[jh-1]/np.sqrt(nh+1); continue
            if (nh,mh) in done: continue
            done.add((nh,mh))
            jcos=[j for j,v in FIELD.items() if v==(nh,mh,'c')][0]; jsin=[j for j,v in FIELD.items() if v==(nh,mh,'s')][0]
            Nh=np.sqrt(2*(nh+1))
            out[(nr,mr,nh,+mh)] = (c[jcos-1] - 1j*c[jsin-1])/Nh
            out[(nr,mr,nh,-mh)] = (c[jcos-1] + 1j*c[jsin-1])/Nh
    return out
