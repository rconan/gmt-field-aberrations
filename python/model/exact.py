import sympy as sp
from math import factorial
u,ub,v,vb,A,Ab,S,Sb = sp.symbols('u ub v vb A Ab S Sb')
# alpha = Ab = a_j e^{-i phi_j}, sigma = Sb = s e^{-i phi_s}; A, S are their conjugates
def mu(s,n,m):
    am=abs(m); return sp.Rational((-1)**s*factorial(n-s), factorial(s)*factorial((n+am)//2-s)*factorial((n-am)//2-s))
def proj(a,b,n,m):
    """(1/pi) int_disk z^a zb^b Z*_{n,m} dA,  z = r e^{i th}"""
    if a-b != m: return 0
    return sum(2*mu(s,n,m)/sp.Integer(a+b+n-2*s+2) for s in range((n-abs(m))//2+1))
def W(p,n,m):
    hh=(v-A)*(vb-Ab); rr=(u+S)*(ub+Sb); hr=((v-A)*(ub+Sb)+(vb-Ab)*(u+S))/2
    return sp.expand(hh**p*rr**n*hr**m)
_cache={}
def coef(p,n,m, nr,mr, nh,mh):
    key=(p,n,m)
    if key not in _cache: _cache[key]=sp.Poly(W(p,n,m),u,ub,v,vb)
    P=_cache[key]; out=0
    for (a,b,c,d),cf in P.terms():
        if a-b!=mr or c-d!=mh: continue
        out += cf*proj(a,b,nr,mr)*proj(c,d,nh,mh)
    return sp.expand(out)
