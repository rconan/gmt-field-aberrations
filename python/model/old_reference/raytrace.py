"""Numpy replica of the CEO/crseo ray trace used by delrays::double_zernike."""
import numpy as np
M1_K, M1_C = 1-0.9982857, 1/36.
M2_K, M2_C, M2_Z = 1-0.71692784, 1/-4.1639009, 20.26247614
FP_Z, FP_R = -5.83, 2.197173

def sag(x, y, k, c):
    r2 = x*x+y*y
    return c*r2/(1+np.sqrt(1-k*c*c*r2))
def dsag(x, y, k, c):          # partial derivatives of F = z - sag
    r2 = x*x+y*y; q = c/np.sqrt(1-k*c*c*r2)
    return -x*q, -y*q

def surface(P, D, opl, k, c, z0):
    """intersect + reflect a conic whose vertex is at height z0 (CEO: plane z=0 then Newton)"""
    P = P.copy(); P[...,2] -= z0
    kx, ly, m = D[...,0], D[...,1], D[...,2]
    s0 = -P[...,2]/m
    x1 = P[...,0]+kx*s0; y1 = P[...,1]+ly*s0
    s = np.zeros_like(s0)
    for _ in range(50):
        x = x1+kx*s; y = y1+ly*s; z = m*s
        F = z - sag(x, y, k, c); Fx, Fy = dsag(x, y, k, c)
        s = s - F/(Fx*kx + Fy*ly + m)
    x = x1+kx*s; y = y1+ly*s; z = m*s
    Nx, Ny = dsag(x, y, k, c); Nz = np.ones_like(x)
    a = -2*(kx*Nx+ly*Ny+m*Nz)/(Nx*Nx+Ny*Ny+Nz*Nz)
    Dn = np.stack([kx+a*Nx, ly+a*Ny, m+a*Nz], -1)
    Pn = np.stack([x, y, z+z0], -1)
    return Pn, Dn, s0+s

def trace(xy, zen, azi):
    """xy: (N,2) entrance coordinates at z=0. Returns sphere coordinates (N,3) and OPD (N,), metres."""
    xy = np.vstack([[0.,0.], xy])                       # row 0 = chief ray
    d = np.array([np.sin(zen)*np.cos(azi), np.sin(zen)*np.sin(azi), -np.cos(zen)])
    P = np.c_[xy, np.zeros(len(xy))]; D = np.tile(d, (len(xy),1))
    opd = np.tan(zen)*np.cos(azi)*xy[:,0] + np.tan(zen)*np.sin(azi)*xy[:,1]
    P, D, L = surface(P, D, None, M1_K, M1_C, 0.);       opd += L - L[0]
    P, D, L = surface(P, D, None, M2_K, M2_C, M2_Z);     opd += L - L[0]
    # reference sphere (CEO bundle::to_sphere)
    c = P[0]; dc = D[0]
    x, y, z = c[0], c[1], c[2]-FP_Z+FP_R
    g = x*dc[0]+y*dc[1]+z*dc[2]; rho2 = x*x+y*y+z*z
    s = -np.sqrt(g*g-(rho2-FP_R**2))-g
    O = c + dc*s                                         # chief-ray image point = sphere origin
    rkl = dc[0]**2+dc[1]**2
    if rkl < 1e-24: Rs = 23.772110269559725
    else:
        se = np.sqrt((c[0]**2+c[1]**2)/rkl); zE = c[2]+dc[2]*se
        Rs = np.sqrt(O[0]**2+O[1]**2+(O[2]-zE)**2)
    V = P - O
    g = np.einsum('ij,ij->i', V, D); rho2 = np.einsum('ij,ij->i', V, V)
    s = -np.sqrt(g*g-(rho2-Rs**2))-g
    P = P + D*s[:,None]; opd += s - s[0]
    return P[1:], opd[1:]
