"""Ray-traced OPD with the old (CEO) reference sphere: paraxial focus rounded to -5.83 m and
radius 2.197173 m. Run from this directory; writes truth_opd.npy here (used by sec6.py and sec7.py)."""
import sys, pathlib
H=pathlib.Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parent)]       # the old raytrace.py shadows ../raytrace.py
import truth
