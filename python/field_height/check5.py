import numpy as np; from core import *
exec(open('check2.py').read().split('for seg in')[0])
def tri(r,rot=90): return np.array([[r*np.cos(a),r*np.sin(a)] for a in np.deg2rad([rot,rot+120,rot+240])])
print('noise: 1 nm rms per Noll coefficient; std of alpha (units of 10\' radius -> mas) and dw040 (nm)')
for seg in (0,6):
    for rr in (0.4,0.6,0.8,1.0):
        Hr=tri(rr); bn=full.bfield(seg,Hr,np.zeros((2,2)))
        if seg<6: a,dw,rms,J=estimate_ext(full,seg,Hr,bn,bn,iters=1)
        else: a,rms,J=estimate(full,seg,Hr,bn,bn,iters=1)
        C=np.linalg.inv(J.T@J); sd=np.sqrt(np.diag(C))
        mas=600e3  # 10' in mas
        print('seg%d r=%.1f  sd(a1x,a1y,a2x,a2y)= %s mas'%(seg+1,rr,' '.join('%.2f'%(v*mas) for v in sd[:4])),
              (' sd(dw040)=%.3f nm'%sd[4]) if seg<6 else '', ' corr(a1,a2)x=%.3f'%(C[0,2]/np.sqrt(C[0,0]*C[2,2])))
# rigid body -> (alpha, dw) sensitivity
print('\nsensitivities per unit motion: alpha in mas per um or per urad; dw040 nm per um/urad')
units={'d':1e-6,'r':1e-6}
dofs=[('M1','dx'),('M1','dy'),('M1','dz'),('M1','rx'),('M1','ry'),('M2','dx'),('M2','dy'),('M2','dz'),('M2','rx'),('M2','ry')]
for seg in (0,6):
    bn=raytrace_b(seg,H)
    for m,dof in dofs:
        v=dict(T=(1e-6 if dof=='dx' else 0,1e-6 if dof=='dy' else 0,1e-6 if dof=='dz' else 0)) if dof[0]=='d' else {dof:1e-6}
        piv=m2_footprint(seg) if m=='M2' else m1pivot(seg)
        kw={m.lower():dict(v,C=piv)}
        bm=raytrace_b(seg,H,**kw)
        if seg<6: a,dw,rms,J=estimate_ext(full,seg,H,bm,bn)
        else: a,rms,J=estimate(full,seg,H,bm,bn); dw=[0.]
        print('seg%d %s %s: a1=(%+.3f,%+.3f) a2=(%+.3f,%+.3f) mas  dw040=%+.4f nm  fit rms %.1e nm'%(seg+1,m,dof,*(a.ravel()*600e3),dw[0],rms))
