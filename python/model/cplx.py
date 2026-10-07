import numpy as np
from math import factorial
nr, nt = 12, 28
xr, wr = np.polynomial.legendre.leggauss(nr); r = 0.5*(xr+1); wr = 0.5*wr*r
t = 2*np.pi*np.arange(nt)/nt
R, T = np.meshgrid(r, t, indexing='ij'); Wd = np.outer(wr, np.full(nt, 2*np.pi/nt))
Rf, Tf, w = R.ravel(), T.ravel(), Wd.ravel()
X, Yc = Rf*np.cos(Tf), Rf*np.sin(Tf)
def mu(s,n,m):
    am=abs(m); return (-1)**s*factorial(n-s)/(factorial(s)*factorial((n+am)//2-s)*factorial((n-am)//2-s))
def Zstar(n,m):
    return sum(mu(s,n,m)*Rf**(n-2*s) for s in range((n-abs(m))//2+1))*np.exp(-1j*m*Tf)
def a_direct(p,n,m, nr_,mr_, nh_,mh_, aj,phij,s,phis):
    ajx,ajy=aj*np.cos(phij),aj*np.sin(phij); sx,sy=s*np.cos(phis),s*np.sin(phis)
    Hx=X[:,None]-ajx; Hy=Yc[:,None]-ajy; Px=X[None,:]+sx; Py=Yc[None,:]+sy
    W=(Hx**2+Hy**2)**p*(Px**2+Py**2)**n*(Hx*Px+Hy*Py)**m
    return (w*Zstar(nh_,mh_)) @ W @ (w*Zstar(nr_,mr_)) / np.pi**2
