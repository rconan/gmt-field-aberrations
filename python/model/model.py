import numpy as np
from crseo_data import Y, PUPIL, FIELD

nr, nt = 10, 20     # exact for the polynomial degrees involved
xr, wr = np.polynomial.legendre.leggauss(nr); r = 0.5*(xr+1); wr = 0.5*wr*r
t = 2*np.pi*np.arange(nt)/nt
R, T = np.meshgrid(r, t, indexing='ij'); Wd = np.outer(wr, np.full(nt, 2*np.pi/nt))

def noll(j, rr, tt):
    return {1: np.ones_like(rr), 2: 2*rr*np.cos(tt), 3: 2*rr*np.sin(tt),
            4: np.sqrt(3)*(2*rr**2-1), 5: np.sqrt(6)*rr**2*np.sin(2*tt),
            6: np.sqrt(6)*rr**2*np.cos(2*tt), 7: np.sqrt(8)*(3*rr**3-2*rr)*np.sin(tt),
            8: np.sqrt(8)*(3*rr**3-2*rr)*np.cos(tt)}[j]
X, Yc, w = (R*np.cos(T)).ravel(), (R*np.sin(T)).ravel(), Wd.ravel()
ZP = np.array([noll(j,R,T).ravel() for j in PUPIL]) * w / np.pi   # (4,N) projector
ZF = np.array([noll(j,R,T).ravel() for j in FIELD]) * w / np.pi   # (8,N)

def response(p, n, m, sx, sy, ajx=0., ajy=0.):
    Hx = X[:,None]-ajx; Hy = Yc[:,None]-ajy
    Px = X[None,:]+sx;  Py = Yc[None,:]+sy
    W = (Hx**2+Hy**2)**p * (Px**2+Py**2)**n * (Hx*Px+Hy*Py)**m
    return ZP @ W.T @ ZF.T          # (4 pupil, 8 field)

def terms(kmax, lmax, mmax):
    return [(p,n,m) for m in range(mmax+1) for p in range(10) for n in range(10)
            if 2*p+m<=kmax and 2*n+m<=lmax]

def design(tl, s, phis):
    segs = [(s*np.cos(a), s*np.sin(a)) for a in phis] + [(0.,0.)]
    A = np.zeros((4,8,7,len(tl)))
    for it,(p,n,m) in enumerate(tl):
        for k,(sx,sy) in enumerate(segs):
            A[:,:,k,it] = response(p,n,m,sx,sy)
    return A

def fit(tl, s, phis, data=Y):
    A = design(tl, s, phis); M = A.reshape(-1,len(tl))
    wgt, *_ = np.linalg.lstsq(M, data.ravel(), rcond=None)
    return wgt, (data.ravel()-M@wgt).reshape(data.shape), A

def response_k(p, n, m, s, phi, kap):
    """outer-segment pupil foreshortened radially by kap: rho_g = s u + kap (rho.u) u + (rho.v) v"""
    ux, uy = np.cos(phi), np.sin(phi)
    ru = X*ux + Yc*uy; rv = -X*uy + Yc*ux
    gx = s*ux + kap*ru*ux - rv*uy;  gy = s*uy + kap*ru*uy + rv*ux
    Hx = X[:,None]; Hy = Yc[:,None]
    Px = gx[None,:]; Py = gy[None,:]
    W = (Hx**2+Hy**2)**p * (Px**2+Py**2)**n * (Hx*Px+Hy*Py)**m
    return ZP @ W.T @ ZF.T

def fit_k(tl, s, kap, phis, data=Y):
    A = np.zeros((4,8,7,len(tl)))
    for it,(p,n,m) in enumerate(tl):
        for k,a in enumerate(phis): A[:,:,k,it] = response_k(p,n,m,s,a,kap)
        A[:,:,6,it] = response(p,n,m,0.,0.)
    M = A.reshape(-1,len(tl))
    wgt, *_ = np.linalg.lstsq(M, data.ravel(), rcond=None)
    return wgt, (data.ravel()-M@wgt).reshape(data.shape), A
