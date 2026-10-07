import numpy as np
from noll2cplx import to_complex
NAMES=[(0,2,0),(0,1,1),(0,0,2),(1,2,0),(0,1,2),(1,1,1),(0,3,0),(0,2,1),(1,3,0),(0,2,2),(1,2,1)]
LAB=lambda t:'w%d%d%d'%(2*t[0]+t[2],2*t[1]+t[2],t[2])
S=8.71/4.1825
PHI=np.array([(3-2*k)*np.pi/6 for k in range(6)])      # entrance-pupil segment azimuths

def cascade(b, s=S, phi=PHI):
    """Closed-form inversion. b: (4,8,7) real Noll coefficients, pupil j=5..8, field j=1..8,
    segments 1..6 outer (azimuths phi, offset s), 7 centre. Returns dict of omega estimates
    (complex: imaginary part is a consistency check)."""
    a=to_complex(b)                         # dict -> arrays over 7 segments
    sig=s*np.exp(-1j*phi)                   # sigma for the 6 outer segments
    A=lambda key: a[key][:6]; C=lambda key: a[key][6]
    av=lambda x: np.mean(x)
    w={}
    # 1. fifth/sixth-order terms reached by the counter-rotating field modes (outer segments)
    w['w351']=av(96*A((3,1,3,1))/sig**2)
    w['w262']=av(36*A((3,1,2,2))/sig**3)
    w['w151']=av(16*A((3,1,1,1))/sig**2 - 2/3*w['w351'])
    # 2. centre segment: field coma and tilt of pupil coma
    w['w331']=288*(C((3,1,3,-1)) - w['w351']/240)
    w['w131']=48*(C((3,1,1,-1)) - w['w331']/72 - w['w151']/40 - w['w351']/60)
    # 3. outer segments: field astigmatism of pupil coma, then centre field astigmatism of pupil astigmatism
    w['w242']=av(48*(A((3,1,2,-2))/np.conj(sig) - w['w262']*(1/30+s**2/12)))
    w['w222']=36*(C((2,2,2,-2)) - w['w242']/48 - w['w262']/60)
    # 4. field defocus of pupil astigmatism and pupil coma: 2x2 system in (w240, w260)
    yA=av(A((2,2,2,0))/sig**2) - w['w242']/36 - w['w262']*(1/16+s**2/12)
    yC=av(A((3,1,2,0))/sig)    - w['w242']/72 - w['w262']*(1/40+s**2/16)
    M=np.array([[1/18,1/8+s**2/6],[1/36,1/20+s**2/8]]); w['w240'],w['w260']=np.linalg.solve(M,[yA,yC])
    # 5. field piston of pupil astigmatism and pupil coma: 2x2 system in (w040, w060)
    yA=av(A((2,2,0,0))/sig**2) - w['w240']/6 - w['w242']/12 - w['w260']*(3/8+s**2/2) - w['w262']*(3/16+s**2/4)
    yC=av(A((3,1,0,0))/sig)    - w['w240']/12 - w['w242']/24 - w['w260']*(3/20+3*s**2/8) - w['w262']*(3/40+3*s**2/16)
    M=np.array([[1/3,3/4+s**2],[1/6,3/10+3*s**2/4]]); w['w040'],w['w060']=np.linalg.solve(M,[yA,yC])
    return w
