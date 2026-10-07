"""Emulates delrays/src/bin/double_zernike.rs: same meshes, same lump-mass weighted
least-squares Zernike fits (Noll j=1..10) on pupil and field, applied to model wavefronts."""
import numpy as np, triangle
from crseo_data import Y

D_SEG, PITCH, R_OFF = 8.365, 0.25, 8.71
R = D_SEG/2

def disc(diameter, pitch, origin=(0.,0.)):
    x0,y0 = origin
    n_rim = int(round(np.pi*diameter/pitch))
    o = 2*np.pi*np.arange(n_rim)/n_rim
    rim = np.c_[0.5*diameter*np.cos(o)+x0, 0.5*diameter*np.sin(o)+y0]
    V = np.vstack([rim, [[x0,y0]]])
    S = np.c_[np.arange(n_rim), (np.arange(n_rim)+1)%n_rim]
    area = pitch**2*np.sqrt(3)/4
    t = triangle.triangulate(dict(vertices=V, segments=S), 'QDpqa%.17g'%area)
    return t['vertices'], t['triangles']

def lump(V, T):
    p = V[T]; a = 0.5*np.abs((p[:,1,0]-p[:,0,0])*(p[:,2,1]-p[:,0,1])-(p[:,2,0]-p[:,0,0])*(p[:,1,1]-p[:,0,1]))
    w = np.zeros(len(V)); np.add.at(w, T.ravel(), np.repeat(a/3,3)); return w

def noll10(x, y):
    r = np.hypot(x,y); r = r/r.max(); t = np.arctan2(y,x)
    return np.array([np.ones_like(r), 2*r*np.cos(t), 2*r*np.sin(t), np.sqrt(3)*(2*r**2-1),
        np.sqrt(6)*r**2*np.sin(2*t), np.sqrt(6)*r**2*np.cos(2*t),
        np.sqrt(8)*(3*r**3-2*r)*np.sin(t), np.sqrt(8)*(3*r**3-2*r)*np.cos(t),
        np.sqrt(8)*r**3*np.sin(3*t), np.sqrt(8)*r**3*np.cos(3*t)])

def wls_projector(Z, w):
    sw = np.sqrt(w); return np.linalg.pinv((Z*sw).T) * sw      # (J,N): coef = P @ data

# field mesh (arcmin) and projector
FV, FT = disc(20., 2.)
FZ = noll10(FV[:,0], FV[:,1]); FP = wls_projector(FZ, lump(FV,FT))
zeta = FV/10.                                     # field normalised to the 10' radius

SEG = []
for sid in range(1,8):
    if sid < 7:
        o = (3-2*(sid-1))*np.pi/6; origin = (R_OFF*np.cos(o), R_OFF*np.sin(o))
    else: origin = (0.,0.)
    V, T = disc(D_SEG, PITCH, origin)
    loc = -(V-np.array(origin))                 # centred on segment ray, pupil image rotated 180 deg
    Z = noll10(loc[:,0], loc[:,1]); P = wls_projector(Z, lump(V,T))
    rho_g = (loc - np.array(origin)*0 ) / R      # local, normalised
    s_vec = -np.array(origin)/R                  # segment centre in the (rotated) pupil frame
    SEG.append((P, rho_g + s_vec))

def response(p, n, m, segs=SEG):
    out = np.zeros((4,8,7))
    zx, zy = zeta[:,0][:,None], zeta[:,1][:,None]
    for k,(P,g) in enumerate(segs):
        gx, gy = g[:,0][None,:], g[:,1][None,:]
        W = (zx**2+zy**2)**p * (gx**2+gy**2)**n * (zx*gx+zy*gy)**m   # (Nf, Np)
        pup = W @ P.T                                 # (Nf, 10) pupil coefs per field point
        dz = FP @ pup                                 # (10 field, 10 pupil)
        out[:,:,k] = dz[:8, 4:8].T                    # pupil j=5..8, field j=1..8
    return out

def terms(kmax, lmax, mmax):
    return [(p,n,m) for m in range(mmax+1) for p in range(10) for n in range(10)
            if 2*p+m<=kmax and 2*n+m<=lmax]

def fit(tl, data=Y):
    A = np.stack([response(*t) for t in tl], -1); M = A.reshape(-1,len(tl))
    w,*_ = np.linalg.lstsq(M, data.ravel(), rcond=None)
    return w, (data.ravel()-M@w).reshape(data.shape), A

def make_segs(rot=-1, pitch=PITCH):
    segs=[]
    for sid in range(1,8):
        if sid<7:
            o=(3-2*(sid-1))*np.pi/6; origin=(R_OFF*np.cos(o),R_OFF*np.sin(o))
        else: origin=(0.,0.)
        V,T=disc(D_SEG,pitch,origin)
        loc=rot*(V-np.array(origin))
        P=wls_projector(noll10(loc[:,0],loc[:,1]),lump(V,T))
        segs.append((P,(loc+rot*np.array(origin))/R))
    return segs

def fit_with(tl, data, segs):
    A=np.stack([response(*t,segs=segs) for t in tl],-1); M=A.reshape(-1,len(tl))
    w,*_=np.linalg.lstsq(M,data.ravel(),rcond=None)
    return w,(data.ravel()-M@w).reshape(data.shape),A
