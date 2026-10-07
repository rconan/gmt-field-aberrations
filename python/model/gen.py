import sympy as sp
from exact import coef, A, Ab, S, Sb
aj2, s2 = sp.symbols('aj2 s2')
OMEGA = [(0,2,0),(0,1,1),(0,0,2),(1,2,0),(0,1,2),(1,1,1),(0,3,0),(0,2,1),(1,3,0),(0,2,2),(1,2,1)]
def wname(p,n,m): return '%d%d%d'%(2*p+m,2*n+m,m)
FIELDS = [(0,0),(1,-1),(1,1),(2,0),(2,-2),(2,2),(3,-1),(3,1)]

def groups(expr):
    """split polynomial in A,Ab,S,Sb into {(phase_a, phase_s): coeff(aj2,s2)} ;
    phase_a = power of alpha (=Ab) minus power of alphabar (=A), same for sigma"""
    out={}
    P=sp.Poly(sp.expand(expr),A,Ab,S,Sb)
    for (a,ab,s,sb),c in P.terms():
        ka=min(a,ab); ks=min(s,sb)
        key=(ab-a, sb-s)
        out[key]=out.get(key,0)+c*aj2**ka*s2**ks
    return {k:sp.expand(v) for k,v in out.items() if sp.expand(v)!=0}

def tex_coeff(c):
    c=sp.Poly(c,aj2,s2)
    terms=[]
    for (i,j),q in sorted(c.terms(),key=lambda t:(t[0][0]+t[0][1],t[0])):
        mono=''
        if i: mono+=r'\aj^{%d}'%(2*i) if i>1 else r'\aj^2'
        if j: mono+=r's^{%d}'%(2*j) if j>1 else r's^2'
        terms.append((q,mono))
    return terms

def tex_phase(pa,ps):
    out=''
    for sym,pw in ((r'\alpha',pa),(r'\sigma',ps)):
        if pw==0: continue
        base = sym if pw>0 else r'\bar'+sym
        out += base if abs(pw)==1 else base+'^{%d}'%abs(pw)
    return out

def frac(q):
    q=sp.Rational(q)
    if q.q==1: return str(abs(q.p))
    return r'\tfrac{%d}{%d}'%(abs(q.p),q.q)

def tex(expr):
    g=groups(expr)
    if not g: return '0'
    pieces=[]
    for (pa,ps) in sorted(g,key=lambda k:(abs(k[0])+abs(k[1]),k)):
        ph=tex_phase(pa,ps); ts=tex_coeff(g[(pa,ps)])
        if len(ts)==1:
            q,mono=ts[0]; sign='-' if q<0 else '+'
            body=(frac(q) if (q.q!=1 or (mono=='' and ph=='')) or abs(q)!=1 else '')+mono+ph
            if body=='' : body='1'
            pieces.append((sign,body))
        else:
            # common sign handling inside parentheses
            inner=''
            for k,(q,mono) in enumerate(ts):
                sg='-' if q<0 else ('+' if k else '')
                inner+=sg+(frac(q) if (abs(q)!=1 or mono=='') else '')+mono
            pieces.append(('+',r'\left('+inner+r'\right)'+ph))
    s=''
    for k,(sign,body) in enumerate(pieces):
        s+=('-' if sign=='-' else ('' if k==0 else '+'))+body
    return s

def tex(expr):
    g=groups(expr)
    if not g: return '0'
    pieces=[]
    for (pa,ps) in sorted(g,key=lambda k:(abs(k[0])+abs(k[1]),k)):
        ph=tex_phase(pa,ps); ts=tex_coeff(g[(pa,ps)])
        neg = all(q<0 for q,_ in ts)
        if neg: ts=[(-q,m) for q,m in ts]
        if len(ts)==1:
            q,mono=ts[0]
            lead = frac(q) if (abs(q)!=1 or (mono=='' and ph=='')) else ''
            body = lead+mono+ph
        else:
            inner=''
            for k,(q,mono) in enumerate(ts):
                sg='-' if q<0 else ('+' if k else '')
                inner+=sg+(frac(q) if (abs(q)!=1 or mono=='') else '')+mono
            body = (r'\left('+inner+r'\right)'+ph) if ph else inner
        pieces.append(('-' if neg else '+',body))
    s=''
    for k,(sign,body) in enumerate(pieces):
        s+=('-' if sign=='-' else ('' if k==0 else '+'))+body
    return s

def collimated(expr):
    return sp.expand(expr.subs({A:0,Ab:0}))
