import sys as _sys, pathlib as _pl; _P=_pl.Path(__file__).resolve().parents[1]; _sys.path[:0]=[str(_P/'model'),str(_P/'field_height'),str(_P/'rigid_body')]
import os as _os
import re, numpy as np
p=_os.environ.get('DZ_TEX',str(_P.parents[1]/'double-zernike.tex')); s=open(p).read()
def fmt(v):
    r=round(v,3)
    return '$0.000$' if r==0 else '$%+.3f$'%r
out={}
for lab,c in (('tab:ceo-data}','sphere'),('tab:ceo-data-entrance}','entrance')):
    i=s.index('\\label{'+lab); j=s.index('\\end{longtable}',i)
    body=s[i:j]
    rows=[l for l in body.split('\n') if l.count('&')==8 and '$' in l]
    vals=np.array([[float(x.split('\\')[0].strip().strip('$')) for x in l.split('&')[2:]] for l in rows]).reshape(4,8,7)
    d=np.load(str(_P/'model'/('replica_%s.npy'%c)))-np.load(str(_P/'model'/'old_reference'/('replica_%s.npy'%c)))
    new=vals+d; out[c]=new
    nb=body
    for l,row in zip(rows,new.reshape(32,7)):
        parts=l.split('&'); tail=parts[-1]; end=tail[tail.index('\\\\'):] if '\\\\' in tail else ''
        nl='&'.join(parts[:2])+'& '+' & '.join(fmt(v) for v in row)+' '+end
        nb=nb.replace(l,nl,1)
    s=s[:i]+nb+s[j:]
    ch=(np.round(new,3)!=np.round(vals,3)).sum(); print(c,'entries changed at 3 decimals:',ch,' max change %.4f'%np.abs(d).max())
open(p,'w').write(s)
np.save(str(_P/'model'/'ceo_tables_new.npy'),out,allow_pickle=True)
