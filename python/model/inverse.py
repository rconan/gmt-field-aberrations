import numpy as np
from math import factorial
# continuous-disc quadrature, Noll real basis j=1..15 on pupil and field
nr_, nt_ = 14, 32
xr, wr = np.polynomial.legendre.leggauss(nr_); r = 0.5*(xr+1); wr = 0.5*wr*r
t = 2*np.pi*np.arange(nt_)/nt_
RR, TT = np.meshgrid(r, t, indexing='ij'); W = np.outer(wr, np.full(nt_, 2*np.pi/nt_)).ravel()
Rf, Tf = RR.ravel(), TT.ravel(); X, Y = Rf*np.cos(Tf), Rf*np.sin(Tf)
def noll_nm(j):
    n=0
    while (n+1)*(n+2)//2 < j: n+=1
    k=j-n*(n+1)//2-1; ms=[m for m in range(n%2,n+1,2) for _ in ((0,) if m==0 else (0,1))]
    return n, ms[k]
def noll(j, r, t):
    n,m=noll_nm(j)
    R=sum((-1)**s*factorial(n-s)/(factorial(s)*factorial((n+m)//2-s)*factorial((n-m)//2-s))*r**(n-2*s) for s in range((n-m)//2+1))
    if m==0: return np.sqrt(n+1)*R
    return np.sqrt(2*(n+1))*R*(np.cos(m*t) if j%2==0 else np.sin(m*t))
def basis(js): return np.array([noll(j,Rf,Tf) for j in js])*W/np.pi
def response(p,n,m,sx,sy,PJ,FJ,ajx=0.,ajy=0.):
    BP,BF=basis(PJ),basis(FJ)
    Hx=X[:,None]-ajx; Hy=Y[:,None]-ajy; Px=X[None,:]+sx; Py=Y[None,:]+sy
    Wv=(Hx**2+Hy**2)**p*(Px**2+Py**2)**n*(Hx*Px+Hy*Py)**m
    return BP @ Wv.T @ BF.T                       # (pupil, field)
def terms(kmax,lmax,mmax):
    return [(p,n,m) for m in range(mmax+1) for p in range(10) for n in range(12)
            if 2*p+m<=kmax and 2*n+m<=lmax]
def name(t): p,n,m=t; return 'w%d%d%d'%(2*p+m,2*n+m,m)
S_OFF=8.71/(8.365/2)
PHI_ENTRANCE=np.deg2rad(90-60*np.arange(6))
def design(tl,PJ,FJ,s=S_OFF,phis=PHI_ENTRANCE):
    segs=[(s*np.cos(a),s*np.sin(a)) for a in phis]+[(0.,0.)]
    A=np.zeros((len(PJ),len(FJ),7,len(tl)))
    for i,tt in enumerate(tl):
        for k,(sx,sy) in enumerate(segs): A[:,:,k,i]=response(*tt,sx,sy,PJ,FJ)
    return A
def identifiable(A, tol=1e-9):
    M=A.reshape(-1,A.shape[-1]); U,sv,Vt=np.linalg.svd(M,full_matrices=False)
    r=(sv>tol*sv[0]).sum(); V=Vt[:r]
    lev=(V**2).sum(0)        # =1 if e_t in row space (individually identifiable)
    return r, lev, sv
