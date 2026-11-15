# Marlics 2.0

Marlics 2.0 extends the Landau–de Gennes solver of version 1.x in four directions: one code
for the CPU and the GPU, two-dimensional simulations, hydrodynamics (backflow) and active
nematics. Input files written for 1.x run unchanged and, with the default options, give the
same results: the refactored code was compared with 1.x in 11 configurations and agrees to
about 1e-15.

## Highlights

- **CPU and GPU from one source.** The functional derivatives, energies and integrators are
  written once and run on the CPU (OpenMP) or on an NVIDIA GPU (CUDA).
- **Two-dimensional simulations** (`dimension 2`) with the full five-component Q tensor.
- **Backflow** (`flow stokes`): Beris–Edwards nematodynamics coupled to Stokes flow, in
  periodic boxes, between flat walls and in any confined geometry.
- **Active nematics** (`activity`, `chiral_activity`): achiral and chiral active stresses.
- **New geometries**: ellipse, regular polygon and a disk with two anchorings.
- **Variational surface dynamics**: surface sites follow the exact gradient of the discrete
  energy.
- **Tests**: `make test` and `make test-gpu`.

## Building and running

```
make                              # marlics (CPU, OpenMP)
make test                         # consistency tests
make gpu                          # marlics_gpu (CPU or GPU, needs CUDA)
make gpu CUDA_HOST=g++-6          # e.g. CUDA 9.1 with an old host compiler
make gpu CUDA_ARCH=-arch=sm_75    # CUDA >= 10: the card's architecture
make test-gpu                     # GPU against CPU
./marlics input.txt
./marlics_gpu input.txt
```

`marlics_gpu` runs on the CPU unless the input file says `device gpu`. Both executables read
the same input files.

## New features

### GPU

| Input line | Values | Effect |
|---|---|---|
| `device` | `cpu` (default), `gpu` | where the simulation runs (`gpu` needs `marlics_gpu`) |
| `gpu_precision` | `double` (default), `mixed` | `mixed` evaluates the functional derivatives in single precision: about twice as fast per step, with a force error of about 1e-6 |

- FIRE, DP5 (`DP5m`, `DP5s`), RK2 and Euler run on the GPU with every geometry except the
  experimental wedges.
- Memory on the GPU: FIRE needs about 3 arrays of 5·Nx·Ny·Nz doubles (301³: about 4 GB),
  DP5 about 9 (301³: about 10 GB).
- Measured on a 301³ droplet with FIRE: 0.16 s per step on a GTX 1660 Ti against 0.56 s per
  step on the CPU.

### Two-dimensional simulations

`dimension 2` simulates a single layer (Nz = 1, set automatically) with no variation along z.
Q keeps its five components, so chirality and an out-of-plane director are allowed: the
result is exactly the z-invariant 3D system. Energies refer to a slice of thickness `dz`
(default `dx`). Without `dimension`, or with `dimension 3`, everything works as in 1.x.

| `geometry` | Domain | Boundaries |
|---|---|---|
| `disk` (alias `cylinder`) | disk, Nx = Ny | edge: boundary 0 |
| `ellipse` | semi-axes `ellipse_a` (x), `ellipse_b` (y), in nm | edge: boundary 0 |
| `polygon` | regular polygon with `polygon_sides` sides, inscribed in the disk | edge: boundary 0 |
| `channel` | strip, periodic in x | walls at y = 0 and y = Ny−1 |
| `box` | rectangle | x walls: boundary 0; y walls: boundary 1 |
| `bulk` | fully periodic | none |
| `custom` | file of points, as in 3D | from the file |

- `hybrid_arc f` (disk): the arc |φ| < πf around +x is boundary 1, so the two parts of the
  edge take their own `anchoring_type` and `Wo1`.
- With `dimension 3` and Nz > 1 these geometries are z-invariant 3D domains (cylinder,
  elliptic cylinder, channel, square tube). `sphere`, `slab`, `hemisphere` and the wedges need
  three dimensions.
- `plane_confinement W` (10⁶ J/m³) adds the energy W·Q_zz, which keeps the director in the
  plane. Without it, planar degenerate walls let the director escape along z.
- New initial conditions: `defects` (one line `defect x y charge` per defect, in lattice
  units, on a background director at angle `phi_i`) and `cholesteric_2d` (helix axis in the
  plane at angle `phi_i`, pitch `p0_i`).
- Examples in `examples/2d/`: a defect pair in a disk, a quenched box and a cholesteric
  channel.

### Backflow (`flow stokes`)

Q is coupled to the flow of the liquid crystal in the Stokes limit, following the
Beris–Edwards equations in the form of Denniston, Orlandini and Yeomans, Phys. Rev. E 63,
056702 (2001). The stress includes the distortion stress for every elastic constant. The flow
is off by default (`flow none`) and then costs nothing.

| Input line | Units | Meaning |
|---|---|---|
| `flow` | | `none` (default) or `stokes` |
| `viscosity` | Pa·s | isotropic viscosity (required with the flow) |
| `xi` | | flow-aligning parameter (required with the flow); the start of the run reports whether it is flow aligning or tumbling |
| `wall_velocity ux uy uz` | nm/µs | velocity of the upper wall between flat walls (Couette shear) |
| `flow_solver` | | `direct` or `iterative` for confined geometries (chosen automatically) |

The solver depends on the geometry; the input does not change:

| Geometry | Solver |
|---|---|
| `bulk` (2D or 3D) | FFT |
| `channel`, `slab` | FFT along the walls and the exact discrete Stokes problem across them, mode by mode |
| `disk`, `ellipse`, `polygon`, `box`, `sphere`, `hybrid_sphere`, `custom`, ... | periodic FFT solver with Lagrange forces on two layers of surface sites; direct (precomputed matrix, up to 10 000 rows, all 2D cases) or conjugate gradients |

- No slip on every wall; the pressure and the continuity equation also hold on the surface
  sites next to the bulk, so an isotropic stress moves no fluid.
- FFTs: pocketfft on the CPU (`src/third_party/`, BSD license), cuFFT on the GPU.
- With the flow, the snapshots have three more columns, `ux,uy,uz` (nm/µs). Without it the
  output files are unchanged.
- The flow changes the dynamics, not the equilibrium states, so it needs a time integrator
  (DP5, RK2, Euler); FIRE is refused.

### Active nematics

| Input line | Units | Meaning |
|---|---|---|
| `activity` (alias `zeta`) | 10⁶ Pa | active stress −ζq, q = 2Q/3; ζ > 0 extensile, ζ < 0 contractile |
| `chiral_activity` | 10⁶ Pa | chiral active stress −ζ′(ε·q − q·ε)/2 (Hoffmann et al., Soft Matter 16, 764 (2020)) |
| `elastic_stress` | 0 or 1 (default) | 0 drops the passive stress from the flow |
| `xi_wall`, `xi_wall_width`, `xi_wall_inner` | –, nm, nm | a different flow-aligning parameter in a layer (or ring) next to the walls |
| `xi_core`, `xi_core_S` | – | a different flow-aligning parameter in the defect cores (S below `xi_core_S`·S_eq, smooth mask) |

- The active stresses act through the flow and need `flow stokes`.
- `xi_wall*` and `xi_core*` run on the CPU only.
- The active length √(L1/|ζ|) should span several lattice spacings.

### Surface dynamics

`surface_dynamics variational` makes every surface site follow the exact gradient of the
discrete energy (the bulk energy of its volume and its anchoring energy). The energy then
never rises, and the stationary states are minima of the reported energy. On coarse curved
surfaces it is markedly more accurate with planar degenerate anchoring: in a disk of radius 16
sites two +1/2 defects settle at 0.66 R, the strong-anchoring value, while the local scheme
pushes them onto the wall. The default remains `local`, the scheme of 1.x.

### Minimization (FIRE)

- FIRE is now FIRE 2.0 and uses the same functional derivative and boundary conditions as the
  time integrators, so it converges to their equilibrium states. `FIRE_el` remains as an
  alias with different velocity-mixing parameters.
- The L2 term of the FIRE force was corrected, and L3 and Ls were added to it.

### Analysis tools (`auxiliar/`)

- `defects/track_defects.py`: point defects of 2D runs, linked into trajectories, with
  annihilations, creations and defects crossing the boundary.
- `defects/track_lines.py`: disclination lines of 3D runs, with the local profile
  (+1/2 wedge, −1/2 wedge, twist) along each line.
- `codegen/gen_gradient_derivative.py`: generates the distortion stress for every elastic
  constant.

### Tests

`make test` runs the consistency tests on the CPU in about a minute:

| Test | Checks |
|---|---|
| `test_energy_gradient` | the dynamics is the gradient of the reported energy |
| `test_anchoring` | anchoring energies against their functional derivatives |
| `test_relaxation` | DP5 never raises the energy; FIRE reaches a stationary state |
| `test_initial_conditions` | the initial conditions |
| `test_2d` | 2D equals the z-invariant 3D system; 2D geometries and initial conditions |
| `test_flow` | distortion stress against finite differences; the Stokes solver; no energy rise with flow; equilibria unchanged |
| `test_flow_walls` | exact Poiseuille and Couette profiles; alignment under shear at the flow-aligning angle |
| `test_flow_confined` | the confined solver in disk, box and droplet; direct against iterative; rigid rotation |
| `test_active` | bend and splay instabilities; self-propulsion of +1/2 and not of −1/2; the chiral stress; `xi_wall` and `xi_core` |
| `test_surface_dynamics` | the variational surface rate against the exact energy gradient |

`make test-gpu` compares the GPU with the CPU (forces, energies and relaxations, with and
without flow and activity, in 2D and 3D). `auxiliar/defects/test_track_*.py` check the
trackers on synthetic fields.

## Other changes and fixes

- A single energy routine for every geometry; this fixes index collisions in the energies of
  `sphere`, `custom` and `hybrid_sphere`.
- `homeotropic_ansatz` rewritten as the bipolar cholesteric droplet.
- Lower memory use in the integrators.
- FIRE stops with a message when the state becomes NaN, also under `-Ofast`.
- Fixed buffer sizes of output file names.
- `hemisphere`, `wedge` and `wedge2` are marked experimental (see README).

## Known limitations

- `xi_wall*`, `xi_core*` and DP5 with noise run on the CPU only.
- Large 3D droplets with flow are slow: the confined solver needs tens to a hundred FFT
  solves per flow evaluation (301³ is out of reach).
- `gpu_precision mixed` is useful only with convergence tolerances looser than about 1e-6.
