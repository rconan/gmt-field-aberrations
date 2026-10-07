from gen import *
from exact import coef
NEW = {(0,3,0),(0,2,1),(1,3,0),(0,2,2),(1,2,1)}
LABEL = {(0,0):'piston',(1,-1):'tilt',(1,1):'tilt',(2,0):'defocus',(2,-2):'astig.',(2,2):'astig.',(3,-1):'coma',(3,1):'coma'}
def table(nr, mr, fields, label, caption):
    L=[]
    L.append(r'\begin{longtable}{@{}l l >{\raggedright\arraybackslash}p{0.55\textwidth} >{\raggedright\arraybackslash}p{0.22\textwidth}@{}}')
    L.append(r'\caption{%s}\label{%s}\\'%(caption,label))
    head=r'\toprule Field mode & & $a_{(%d,%d)i_\h}(\vaj,\vs)$ & $\vaj=0$ \\ \midrule'%(nr,mr)
    L.append(head+r' \endfirsthead')
    L.append(r'\multicolumn{4}{@{}l}{\small\emph{(continued)}}\\'+head+r' \endhead')
    L.append(r'\bottomrule \endlastfoot')
    first_group=True
    for f in fields:
        rows=[(w,coef(*w,nr,mr,*f)) for w in OMEGA]
        rows=[(w,c) for w,c in rows if c!=0]
        if not rows: continue
        if not first_group: L.append(r'\midrule')
        first_group=False
        for k,(w,c) in enumerate(rows):
            fm = (r'%s $(%d,%s)$'%(LABEL[f],f[0],('%+d'%f[1]) if f[1] else '0')) if k==0 else ''
            om = r'$\omega_{j%s}%s$'%(wname(*w), r'^{\dagger}' if w in NEW else '')
            L.append(r'%s & %s & $%s$ & $%s$ \\'%(fm,om,tex(c),tex(collimated(c))))
    L.append(r'\end{longtable}')
    return '\n'.join(L)
FULL=[(0,0),(1,-1),(1,1),(2,0),(2,-2),(2,2),(3,-1),(3,1)]
open('tab_astig.tex','w').write(table(2,2,FULL,'tab:pupil-astig',
 r'Double-Zernike coefficients $a_{(2,2)i_\h}$ of pupil astigmatism ($m\le2$, $k\le3$, $l\le6$; $^\dagger$: fifth- and sixth-order terms).'))
open('tab_coma.tex','w').write(table(3,1,FULL,'tab:pupil-coma',
 r'Double-Zernike coefficients $a_{(3,1)i_\h}$ of pupil coma ($m\le2$, $k\le3$, $l\le6$; $^\dagger$: fifth- and sixth-order terms).'))
open('tab_sph.tex','w').write(table(4,0,[(0,0),(1,-1),(2,0),(2,-2),(3,-1)],'tab:pupil-sph',
 r'Double-Zernike coefficients $a_{(4,0)i_\h}$ of pupil spherical aberration ($m\le2$, $k\le3$, $l\le6$; $^\dagger$: fifth- and sixth-order terms). Coefficients with $m_\h>0$ are the complex conjugates of those with $-m_\h$.'))
print('ok')
