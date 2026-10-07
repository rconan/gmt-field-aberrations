"""Pupil astigmatism / coma (complex, Eq. 31 convention) of one segment at a field point,
for each of the 11 wave-aberration terms, as polynomials in u = zeta_j and sigma (conj convention:
u = zeta_x - i zeta_y, sigma = s_x - i s_y)."""
import sympy as sp, pickle
x,y,ux,uy,sx,sy=sp.symbols('x y u_x u_y s_x s_y',real=True)
u,ub,g,gb=sp.symbols('u ubar sigma sigmabar')
NAMES=[(0,2,0),(0,1,1),(0,0,2),(1,2,0),(0,1,2),(1,1,1),(0,3,0),(0,2,1),(1,3,0),(0,2,2),(1,2,1)]
lab=lambda t:'w%d%d%d'%(2*t[0]+t[2],2*t[1]+t[2],t[2])
def mom(a,b):
    if a%2 or b%2: return sp.Integer(0)
    return sp.gamma(sp.Rational(a+1,2))*sp.gamma(sp.Rational(b+1,2))/(sp.pi*sp.gamma(sp.Rational(a+b,2)+2))
def proj(W,Z):
    P=sp.Poly(sp.expand(W*Z),x,y)
    return sp.nsimplify(sum(c*mom(a,b) for (a,b),c in P.terms()))
r2=x**2+y**2
Z5=sp.sqrt(6)*2*x*y; Z6=sp.sqrt(6)*(x**2-y**2)
Z7=sp.sqrt(8)*(3*r2-2)*y; Z8=sp.sqrt(8)*(3*r2-2)*x
# complex conversion: ux = (u+ub)/2, uy = i(u-ub)/2  (u = ux - i uy)
sub={ux:(u+ub)/2, uy:sp.I*(u-ub)/2, sx:(g+gb)/2, sy:sp.I*(g-gb)/2}
out={}
for t in NAMES:
    p,n,m=t
    px,py=x+sx,y+sy
    W=(ux**2+uy**2)**p*(px**2+py**2)**n*(ux*px+uy*py)**m
    b5,b6,b7,b8=[proj(W,Z) for Z in (Z5,Z6,Z7,Z8)]
    cA=sp.expand(((b6-sp.I*b5)/sp.sqrt(6)).subs(sub))
    cC=sp.expand(((b8-sp.I*b7)/sp.sqrt(8)).subs(sub))
    out[lab(t)]=(cA,cC)
    print(lab(t),' A:',cA,'\n      C:',cC)
pickle.dump(out,open('pAC.pkl','wb'))
