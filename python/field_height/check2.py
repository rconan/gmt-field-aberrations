import numpy as np; from core import *
d=np.load('full.npz'); full=Model([tuple(t) for t in d['terms']],d['w'])
H=np.array([[0.8*np.cos(a),0.8*np.sin(a)] for a in np.deg2rad([90,210,330])])
rng=np.random.default_rng(1); HV=rng.uniform(-1,1,(40,2)); HV=HV[np.hypot(*HV.T)<=1][:20]
def m1pivot(seg):
    sx,sy=SEG[seg]; x,y=sx*R_SEG,sy*R_SEG; r2=x*x+y*y
    z=M1_C*r2/(1+np.sqrt(1-M1_K*M1_C**2*r2)); return (x,y,z)
cases=[('M2 dx 100um',dict(m2=dict(T=(1e-4,0,0)))),('M2 dy 100um',dict(m2=dict(T=(0,1e-4,0)))),
       ('M2 rx 10urad',dict(m2=dict(rx=1e-5))),('M2 ry 10urad',dict(m2=dict(ry=1e-5))),
       ('M2 dz 100um',dict(m2=dict(T=(0,0,1e-4)))),
       ('M1 dx 100um',dict(m1=dict(T=(1e-4,0,0)))),('M1 rx 10urad',dict(m1=dict(rx=1e-5))),
       ('M1 dz 100um',dict(m1=dict(T=(0,0,1e-4))))]
for seg in (0,6):
    bn=raytrace_b(seg,H); bnv=raytrace_b(seg,HV)
    for name,kw in cases:
        kw={k:dict(v,C=(m2_footprint(seg) if k=='m2' else m1pivot(seg))) for k,v in kw.items()}
        bm=raytrace_b(seg,H,**kw); bmv=raytrace_b(seg,HV,**kw)
        a,rms,J=estimate(full,seg,H,bm,bn)
        pred=full.bfield(seg,HV,a)-full.bfield(seg,HV,np.zeros((2,2)))+bnv
        sig=np.sqrt(((bmv-bnv)**2).mean())
        print('seg%d %-13s a1=(%+.2e,%+.2e) a2=(%+.2e,%+.2e) fit rms %.1e | signal %.2f nm, predict err %.2e nm'%(
            seg+1,name,*a.ravel(),rms,sig,np.sqrt(((pred-bmv)**2).mean())))
