# Python analysis scripts for the double-Zernike paper

These numpy/sympy scripts produced the numbers in *GMT Field Aberrations & Misalignments — Double Zernike
expansion* (`double-zernike.tex`). They replicate the CEO/crseo ray trace used by `delrays` with plain
numpy, so they run without a GPU.

Requirements: Python ≥ 3.10, `numpy`, and `sympy` for `model/exact.py` and `model/gen.py`.

**Run every script from its own directory** (`cd python/collimation && python3 fullscheme.py`). Scripts
read and write their data files in the current directory. Scripts that need modules from another folder
add `model/`, `field_height/` and `rigid_body/` to `sys.path` themselves.

## Layout

| Folder | Paper | Contents |
|---|---|---|
| `model/` | §§2–7, Appendix | Ray-trace replica, double-Zernike model, comparison with CEO, inversion |
| `model/old_reference/` | §7.4, Appendix | The old CEO reference sphere (focus −5.83 m, radius 2.197173 m), kept for the before/after tables |
| `field_height/` | §8.1–8.6 | Field heights of misaligned segments from pupil astigmatism and coma at 3 field points |
| `rigid_body/` | §8.7 | Field heights as functions of M1/M2 rigid-body motions (third order and four-part model), line of sight |
| `collimation/` | §8.8 | Collimating M1/M2 segment pairs with M2, focus, the full scheme, M1 figure errors |

### `model/`
- `raytrace.py`: numpy replica of the CEO ray trace. The reference sphere is centred on the paraxial focus and has the medial radius, as in `delrays`.
- `emulate.py`, `replicate.py`, `dzdata.py`: emulate `delrays/src/bin/double_zernike.rs` (meshes, lumped-mass weighted projections).
- `crseo_data.py`: CEO double-Zernike tables (old reference).
- `truth.py`: ray-traced OPD over the field. Writes `truth_opd.npy` and fits the wave-aberration coefficients ω.
- `truth2.py`: further fits of `truth_opd.npy`.
- `model.py`, `general.py`, `inverse.py`, `annulus.py`, `cplx.py`, `noll2cplx.py`: double-Zernike projections of the model and Noll ↔ complex conventions.
- `exact.py`, `gen.py`, `mktables.py`: symbolic double-Zernike coefficients and the pupil astigmatism, coma and spherical tables (sympy).
- `sec6.py`: §6 comparison with the CEO ray trace.
- `sec7.py`: §7 inversion. Writes `sec7.npy`.
- `invert.py`: closed-form cascade.
- `rankcheck.py`, `rankfix.py`, `pupil10.py`, `pupil10n.py`, `lsqfix.py`: identifiability and rank, pupil modes 1–10, and the least-squares inversion. `lsqfix.py` writes `double_zernike_lsq.txt`, the test data in `delrays/data/`.
- `apptables.py`: updates the CEO appendix tables in the paper. Set `DZ_TEX` to the paper's path.
- `old_reference/truth_old.py`: writes `old_reference/truth_opd.npy` with the old reference. `sec7.py` needs it.

### `field_height/`
- `core.py`: per-surface polynomial model `Model`, ray trace of misaligned segments `raytrace_b`, and Gauss–Newton estimators of the field heights.
- `m1alone.py`: M1-alone (prime-focus) wavefront. Writes `opd_m1.npy`.
- `models.py`: fits the per-surface coefficients. Writes `full.npz`.
- `sym.py`: symbolic pupil astigmatism and coma at a field point.
- `closed.py`: third-order closed-form solution.
- `check1.py`–`check6.py`: model checks, ray-traced misalignments, noise and field-point layout.
- `export.py`: writes `gmt_per_surface.txt` and `fixtures.rs` (the `delrays` test fixtures).

### `rigid_body/`
- `seidel_parts.py`: Seidel sums split into base-sphere and aspheric parts.
- `rigid.py`: analytic NAT field heights from rigid-body motions.
- `coeffs.py`: their coefficients.
- `parts.py`, `fourpart.py`: four-part decomposition, all orders.
- `nat_check.py`: validation against the ray trace.
- `los.py`, `los2.py`: line of sight of a misaligned segment.

### `collimation/`
- `sens.py`: sensitivity matrices G, D, P of a segment pair. Writes `sens.npy`.
- `collimate.py`, `collimate2.py`: M2 collimation from Noll 5–8 and the ray-trace test.
- `kanalytic.py`: analytic compensation matrix K.
- `focus.py`: focus signature ε. Writes `sens48.npy`.
- `residual.py`: M1 signature that M2 cannot compensate. Writes `G48.npy`.
- `fullscheme.py`: M2 plus M1 tilts with measured line of sight, Noll 4–8.
- `figtrace.py`: ray trace showing that an M1 figure error is the same at every field point.
- `figerr.py`–`figerr4.py`: M1 figure error in the collimation scheme. Covers the response, degeneracies and strategies with noise.
- `figgeom.py`–`figgeom3.py`: dependence on where the 3 field points are, and figure errors that M1+M2 motions mimic.
- `diffscheme.py`–`diffscheme3.py`: differential scheme. Differences between field points remove the figure error; M2 is derived from them and the field mean is corrected by M1 bending. Includes the noise sweep (`diffscheme2.py`) and the comparison with the field-height collimation (`diffscheme3.py`).
- `sensors.py`: edge sensors on M1 or M2. Posterior uncertainty of the misalignments, the figure error and the collimation state for different sensor sets.
- `shmodel.py`: Shack-Hartmann model (48x48 subapertures, Fried geometry). Writes `sh.npy`, the centroid sensitivities to the M1/M2 DOFs and the figure modes at the 3 probes.
- `shnoise.py`: noise propagation from centroids to segment Noll 2-8, and the mean slope of coma.
- `shtilt.py`: segment tilt at the 3 field points for the motions that mimic figure errors.
- `shstrat.py`: the figure-error strategies with centroids vs Noll 4-8 + line of sight. Arguments: centroid noise (mas), probe calibration offsets (nm).
- `shsensors.py`: posterior uncertainties with centroids and segment-position sensors. Argument: centroid noise (mas).

## Regenerating the large files

Run these in order, each from `python/`. The small `.npy` and `.npz` files are committed.

```
cd model && python3 truth.py                                     # model/truth_opd.npy
cd model/old_reference && python3 truth_old.py                   # old-reference truth_opd.npy
cd field_height && python3 m1alone.py && python3 models.py       # opd_m1.npy, full.npz
```
