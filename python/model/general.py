import numpy as np
from math import factorial
nr_, nt_ = 14, 40
xr,wr=np.polynomial.legendre.leggauss(nr_); r=0.5*(xr+1); wr=0.5*wr*r
t=2*np.pi*np.arange(nt_)/nt_
RR,TT=np.meshgrid(r,t,indexing='ij'); W=np.outer(wr,np.full(nt_,2*np.pi/nt_)).ravel()
Rf,Tf=RR.ravel(),TT.ravel(); X,Yc=Rf*np.cos(Tf),Rf*np.sin(Tf)
def jnm(nmax):
    out=[];j=0
    for n in range(nmax+1):
        for m in range(n%2,n+1,2):
            if m==0: j+=1; out.append((j,n,0))
            else: j+=1; out.append((j,n,m)); j+=1; out.append((j,n,m))
    return out
def noll(j,n,m,r,t):
    R=sum((-1)**s*factorial(n-s)/(factorial(s)*factorial((n+m)//2-s)*factorial((n-m)//2-s))*r**(n-2*s) for s in range((n-m)//2+1))
    if m==0: return np.sqrt(n+1)*R
    return np.sqrt(2*(n+1))*R*(np.cos(m*t) if j%2==0 else np.sin(m*t))
def basis(nmax):
    J=jnm(nmax); return J, np.array([noll(j,n,m,Rf,Tf)*W/np.pi for j,n,m in J])
PJ,PB=basis(3)     # pupil j<=10
def design(terms, segs, field_nmax):
    FJ,FB=basis(field_nmax)
    cols=[]
    for (p,n,m) in terms:
        blk=[]
        for sx,sy in segs:
            Hx=X[:,None];Hy=Yc[:,None];Px=X[None,:]+sx;Py=Yc[None,:]+sy
            Wv=(Hx**2+Hy**2)**p*(Px**2+Py**2)**n*(Hx*Px+Hy*Py)**m
            blk.append((PB[4:8]@Wv.T@FB.T).ravel())      # pupil j=5..8 x all field modes
        cols.append(np.concatenate(blk))
    return np.array(cols).T

def design_alpha(terms, segs, field_nmax, alpha=(0.,0.)):
    FJ,FB=basis(field_nmax); ax,ay=alpha
    cols=[]
    for (p,n,m) in terms:
        blk=[]
        for sx,sy in segs:
            Hx=X[:,None]-ax;Hy=Yc[:,None]-ay;Px=X[None,:]+sx;Py=Yc[None,:]+sy
            Wv=(Hx**2+Hy**2)**p*(Px**2+Py**2)**n*(Hx*Px+Hy*Py)**m
            blk.append((PB[4:8]@Wv.T@FB.T).ravel())
        cols.append(np.concatenate(blk))
    return np.array(cols).T
