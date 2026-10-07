import numpy as np; from core import *
d=np.load('full.npz'); full=Model([tuple(t) for t in d['terms']],d['w'])
H=np.array([[0.8*np.cos(a),0.8*np.sin(a)] for a in np.deg2rad([90,210,330])])
for seg in (0,1,6):
    bm=full.bfield(seg,H,np.zeros((2,2))); br=raytrace_b(seg,H)
    print(seg, np.abs(bm-br).max(), np.abs(br).max())
