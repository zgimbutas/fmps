# fmps

FMPS (Fast Multi-Particle Scattering): Matlab/Octave front end plus Fortran 90 kernels (via mwrap-generated MEX) for electromagnetic scattering from many spheres/particles in R^3.

## Build

- Run `make` (Linux: Octave + Matlab, OpenMP) or `make macosx` from `fmps/`. All builds must be run from `fmps/`; `makefile.mwrap` does the work and the valid `TARGET=` values live in `fmps/mwrap.inc`.
- Requires gfortran > 4.4.0; the mwrap tool source is bundled in `contrib/mwrap-0.33.3/`.
- Tests are `fmps/test*.m`, run from Matlab/Octave with `fmps/` on the path. Example drivers are in `examples/`.

## Solver pipeline

`fmps/fmps_solver.m` is the canonical workflow; keep this function chain intact:

1. Choose the expansion order `nterms` and build the spherical grid with `e3fgrid`.
2. Define geometry: sphere centers, radii and types. Type codes: PEC = 1, dielectric = 2, file-based scattering matrix = 3 (matrix produced by the Muller solver).
3. Evaluate the incoming plane wave on the sphere boundaries with `em3d_planearb_targeval`.
4. Convert E, H fields to vector spherical harmonic (local) expansions with `em3ehformta`.
5. Apply the reflection/scattering matrices with `em3d_multa_r`.
6. Solve with GMRES using `em3d_multa_fmps` as the matrix-vector product.
7. Post-process: evaluate fields at targets with `em3d_mpole_targeval`.

The solver struct `A` carries all problem state (geometry, material parameters, grid data, expansion order, reflection matrices) and is threaded through every Matlab-level solver call; add new state to `A` rather than to globals.

## Matrix-vector product tiers

`em3d_multa_fmps` selects one of three multipole-to-local translation routines via its internal `imulta_type`:

- `em3d_multa_mpta` — Matlab with partial F90; very slow, debugging only.
- `em3d_multa_mptaf90` — pure Fortran 90 direct translation; small problems.
- `em3d_multa_mptafmm` — Fortran 90 with FMM acceleration; the default and the only choice for large problems (takes `A.iprec`).

Any change to one tier must be checked against the others, since they are supposed to agree to the requested tolerance.

## Reflection coefficients and rotations

- `rcoefs_pec`, `rcoefs_pmc`, `rcoefs_die_arb`, `rcoefs_die_arb_coated` compute sphere scattering matrices analytically and back types 1 and 2.
- `rmatr_file` reads a precomputed reflection matrix for type 3; those files come from the Muller solver in `muller/` (prebuilt executables and drivers; Fortran source in `src/muller/`) and let arbitrary-shape particles be embedded in the multi-sphere solve.
- `emabrotaf` / `emabrotab` apply forward/inverse rotations to multipole expansions so type 3 particles can be oriented by Euler angles (`A.rot`).
