import numpy as np; from core import *
exec(open('check4.py').read().split('res={}')[0])
with open('gmt_per_surface.txt','w') as f:
    f.write('# GMT per-surface wave-aberration coefficients (nm), entrance-pupil coordinates,\n')
    f.write('# rho in segment radii (R = 4.1825 m), zeta in units of the 10 arcmin field radius.\n')
    f.write('# M1 = M1 alone (prime-focus wavefront), M2 = total - M1. Columns: p n m w_M1 w_M2\n')
    for (p,n,m),a,b in zip(TL,d['w'][0],d['w'][1]): f.write('%d %d %d %+.12e %+.12e\n'%(p,n,m,a,b))
fix={}
def lit(a): return '['+', '.join('['+', '.join('%.12e'%v for v in row)+']' for row in a)+']'
out=[]
for seg in (0,6):
    bn=raytrace_b(seg,H)
    out.append('const B_NOM_%d: [[f64; 4]; 3] = %s;'%(seg+1,lit(bn)))
    for name,kw in [cases[0],cases[4],cases[6]]:
        kw2={k:dict(v,C=(m2_footprint(seg) if k=='m2' else m1pivot(seg))) for k,v in kw.items()}
        bm=raytrace_b(seg,H,**kw2)
        if seg<6: a,dw,rms,J=estimate_ext(full,seg,H,bm,bn)
        else: a,rms,J=estimate(full,seg,H,bm,bn); dw=[0.]
        tag=name.split()[0]+'_'+name.split()[1].upper()
        out.append('// %s: alpha = %s, dw040 = %.12e, rms = %.3e'%(name,lit(a),dw[0],rms))
        out.append('const B_%s_%d: [[f64; 4]; 3] = %s;'%(tag,seg+1,lit(bm)))
        out.append('const A_%s_%d: [[f64; 2]; 2] = %s;'%(tag,seg+1,lit(a)))
        out.append('const D_%s_%d: f64 = %.12e;'%(tag,seg+1,dw[0]))
open('fixtures.rs','w').write('\n'.join(out)+'\n')
print(open('fixtures.rs').read()[:1500])
