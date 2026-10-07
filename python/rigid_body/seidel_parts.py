"""Seidel sums of M1 and M2 split into base-sphere and aspheric parts (Welford conventions),
for a marginal ray of height R (segment radius) and a chief ray of angle theta (10')."""
import numpy as np
R=8.365/2; TH=np.deg2rad(10/60)
r1,r2,z2=36.,-4.1639009,20.26247614; K1,K2=-0.9982857,-0.71692784
def parts():
    h,u,hb,ub,n=R,0.,0.,TH,1.; out=[]
    for c,k,d in ((-1/r1,K1,-z2),(-1/r2,K2,0.)):
        np_=-n; up=(n*u-h*c*(np_-n))/np_; ubp=(n*ub-hb*c*(np_-n))/np_
        A=n*(u+h*c); Ab=n*(ub+hb*c); D=up/np_-u/n; a=k*c**3*(np_-n)
        sph=np.array([-A*A*h*D, -A*Ab*h*D, -Ab*Ab*h*D])
        asp=np.array([a*h**4, a*h**3*hb, a*h*h*hb*hb])
        out.append((sph,asp,dict(h=h,hb=hb)))
        u,ub,n=up,ubp,np_; h+=d*u; hb+=d*ub
    return out
# OPD-sign wave coefficients omega = -W: w040 = -S_I/8, w131 = -S_II/2, w222 = -S_III/2
def omegas():
    res={}
    for name,(sph,asp,_) in zip(('M1','M2'),parts()):
        for pn,S in (('sph',sph),('asph',asp)):
            res[(name,pn)]=dict(w040=-S[0]/8*1e9,w131=-S[1]/2*1e9,w222=-S[2]/2*1e9)
    return res
if __name__=='__main__':
    for k,v in omegas().items(): print(k,{a:round(b,3) for a,b in v.items()})
    p=parts(); print('chief-ray height at M2 per unit field: %.6f m (z2*theta=%.6f)'%(p[1][2]['hb'],z2*TH))
