"""Shack-Hartmann model: 48x48 subapertures across the 25.4 m GMT pupil, Fried geometry
(slope = mean of the two edge differences of the wavefront at the subaperture corners).
Only subapertures fully inside a segment are used. Centroid = wavefront slope, in mas.
Saves sh.npy: per segment, centroid sensitivities (mas per um/urad) to the 5 DOFs of M1 and M2
and to the Noll 4-8 figure modes (mas per nm), at the 3 field points (equilateral, r = 0.8)."""
import numpy as np, sens as S
from core import SEG, R_SEG, trace
NSUB, DPUP = 48, 25.4
d = DPUP/NSUB
MAS = 206264.806e3
def geometry(seg):
    cx, cy = SEG[seg][0]*R_SEG, SEG[seg][1]*R_SEG
    k = np.arange(NSUB+1)*d - DPUP/2                          # corner coordinates
    X, Y = np.meshgrid(k, k, indexing='xy')
    inside = np.hypot(X-cx, Y-cy) <= R_SEG
    sub = inside[:-1,:-1] & inside[1:,:-1] & inside[:-1,1:] & inside[1:,1:]   # fully illuminated
    used = np.zeros_like(inside)
    for di in (0,1):
        for dj in (0,1): used[di:NSUB+di, dj:NSUB+dj] |= sub
    return X, Y, sub, used, cx, cy
def slopes(W, sub):
    """W: wavefront (m) on the (49,49) corner grid; returns (2*nsub,) slopes in mas"""
    sx = ((W[:-1,1:]+W[1:,1:]) - (W[:-1,:-1]+W[1:,:-1]))/(2*d)
    sy = ((W[1:,:-1]+W[1:,1:]) - (W[:-1,:-1]+W[:-1,1:]))/(2*d)
    return np.r_[sx[sub], sy[sub]]*MAS
def wavefronts(seg, Hs, x1=None, x2=None):
    X, Y, sub, used, cx, cy = geometry(seg)
    xy = np.c_[X[used], Y[used]]; kw = {}
    if x1 is not None: kw['m1'] = S.mis(x1, seg, 'm1')
    if x2 is not None: kw['m2'] = S.mis(x2, seg, 'm2')
    out = []
    for h in Hs:
        zen = np.deg2rad(np.hypot(*h)*10/60); azi = np.arctan2(h[1], h[0])
        o, _ = trace(xy, zen, azi, **kw)
        W = np.zeros(X.shape); W[used] = o; out.append(slopes(W, sub))
    return np.concatenate(out)
def zernike(j, u, v):
    r2 = u*u+v*v; t = np.arctan2(v, u); r = np.sqrt(r2)
    return {4:np.sqrt(3)*(2*r2-1), 5:np.sqrt(6)*r2*np.sin(2*t), 6:np.sqrt(6)*r2*np.cos(2*t),
            7:np.sqrt(8)*(3*r**3-2*r)*np.sin(t), 8:np.sqrt(8)*(3*r**3-2*r)*np.cos(t)}[j]
def figure_columns(seg, nH):
    X, Y, sub, used, cx, cy = geometry(seg)
    u, v = (X-cx)/R_SEG, (Y-cy)/R_SEG
    return np.array([np.tile(slopes(zernike(j, u, v)*1e-9, sub), nH) for j in (4,5,6,7,8)]).T   # mas per nm
if __name__ == '__main__':
    out = {}
    for seg in (0, 6):
        X, Y, sub, used, cx, cy = geometry(seg)
        c0 = wavefronts(seg, S.H)
        cols = {}
        for mir, k in (('m1','x1'), ('m2','x2')):
            g = []
            for i in range(5):
                e = np.eye(5)[i]*S.UNIT
                g.append((wavefronts(seg, S.H, **{k:e}) - wavefronts(seg, S.H, **{k:-e}))/2)
            cols[mir] = np.array(g).T
        out[seg] = dict(G1=cols['m1'], G2=cols['m2'], F=figure_columns(seg, len(S.H)), nsub=int(sub.sum()))
        print('segment %d: %d subapertures, %d centroids at 3 field points'%(seg+1, sub.sum(), out[seg]['G1'].shape[0]))
    np.save('sh.npy', out, allow_pickle=True)
