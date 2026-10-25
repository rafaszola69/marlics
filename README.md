Welcome to MarLicS
------------------

 Marlics is a software written in C++ to to solve either the Beris-Edwards equation of nematodynamics without flow, or to find minimum energy states of the Landau-de Gennes free energy with the (FIRE) engine. Marlics comes with the most common setups where these simulations are performed, such as confined slab cells and spherical liquid crystal droplets. Different simulations can be performed passing a file describing the system geometry. FInally, we provide the most common anchoring energies as boundary conditions, such as Rapini-Papoular and Founier-Galatola.

The program takes as input a descriptive file giving the simulations parameters and initial conditions generating a series of different snapshots. Some examples of these input files can be found in the examples folder.

In the aux folder you can find some useful scripts with help visualize and start liquid crystal simulations. In renderize scripts you will find some scripts that renderives a view of the director in a given plane. There is also another script which generates iso surfaces of the thermodynamic order parameter in the same folder. In the boundaries folder the user can find some simple program, also written in C++, which can generate geometry files for the sphere and cylinder geometry. Both folders have a readme file explaining how the scripts/program are used.

File structure:
--------------

src/ -> folder containing MarLiCs the source code;

	src/boundary.(cpp/h) -> a superclass used to define the anchoring energies;

	src/boundary_%%%.(cpp/h) -> the class defining the anchoring type %%%;

	src/container.(cpp/h) -> Responsible for printing the output CSV file;

	src/driver.(cpp/h) -> Responsible for reading the input_file, setuping the simulation variables. Contains the definition of struct Simulation_Parameters, which contains all the require variables.
	src/energy.(cpp/h) -> a superclass used to define the functional derivative to be used along with the integrator;

	src/energy_ldg.(cpp/h) -> a class with the Landau de Gennes functional derivative;

	src/geometry.(cpp/h) -> a superclass used to define the system geometry;

	src/geometry_%%%.(cpp/h) -> the class which define the %%% geometry;

	EXPERIMENTAL geometries (unfinished, do not use for production runs):
	  hemisphere -> at the rim of the flat base, fill_ki couples surface sites to outside points (frozen Q);
	  wedge      -> periodic along y although the thickness varies with y; the time integrators diverge;
	  wedge2     -> symmetric wedge of opening angle wTheta, but normals and surface areas are copied from wedge and ignore wTheta.

	src/initial_conditions.(cpp/h) -> here we define all pre-programed initial conditions; 

	src/integrator.(cpp/h) -> a superclass used to define the any temporal integrator

	src/integrator_%%%.(cpp/h) -> the class which define the %%% integrator;

	src/fire_integrator.(cpp/h) -> a class containing the FIRE 2.0 minimizer. It uses the same functional derivative (and boundary conditions) as the time integrators, so it converges to their equilibrium states. "FIRE_el" is accepted as an alias with different velocity-mixing parameters;

	src/physics_%%%.h, src/lattice.h -> the functional derivatives, energies and per-site stencils, written once for the CPU and the GPU;

	src/exec_cpu.h, src/gpu/ -> the executors that run the integrators on the CPU (OpenMP) or on the GPU (CUDA). FIRE, DP5 (m and s), RK2 and Euler are written once (src/*_impl.h) for both;

auxiliar/ -> some auxiliar files to aid the user;

	auxiliar/defects/track_defects.py -> point defects of two-dimensional runs: finds them in every
	snapshot from the in-plane part of Q (also when the director leaves the plane), links them into
	trajectories and reports annihilations, creations and defects crossing the boundary. Times come
	from the Marlics output ("Snapshot k: t" lines). Example:
	    python3 auxiliar/defects/track_defects.py run_dir --log run.log --plot tracks.png
	The snapshots must be close enough that a defect moves less than --max-jump (5 sites) between
	them; --periodic x,y for periodic domains; --start 1 skips a random initial state.
	Test: python3 auxiliar/defects/test_track_defects.py (synthetic fields with known defects);

	auxiliar/defects/track_lines.py -> disclination lines of three-dimensional runs: crossed plaquettes
	(the director comes back reversed around them) joined into lines, with the local profile
	Omega.t along each line (+1: +1/2 wedge, -1: -1/2 wedge, 0: twist); writes the line points and,
	per snapshot, the number of lines, their length and the fractions of each profile. Example:
	    python3 auxiliar/defects/track_lines.py run_dir --log run.log --plot lines.png
	Test: python3 auxiliar/defects/test_track_lines.py (wedge and twist lines, a loop, a uniform field);

	auxiliar/renderize/ -> contain routines to renderize the directors and the isosurface of the parameter order S;

	auxiliar/boundaries/ -> contain 2 examples off c++ files that can be used to produce boundary files;

examples/ -> a series of MarLiCs usage examples;

notebook/ -> Some Mathematica nb files used to evaluate symbolically the functional derivatives.

tests/ -> consistency tests, run with "make test" (about 15 s), including the two-dimensional simulations (test_2d). They check that the dynamics (fill_ki) is the gradient of the reported energy, that the anchoring energies match their functional derivatives, that DP5 never increases the energy and that FIRE converges to a stationary state. Run them after changing energies, boundaries, geometries or integrators;
	
makefile -> a makefile used to compile the code (make), or to compile and run the tests (make test)

Two-dimensional simulations ("dimension 2" in the input file):
	A single layer (Nz=1, set automatically) with no variation along z. The Q tensor keeps its five
	components, so the physics is the 3D one (chirality, out-of-plane director): it is exactly the
	z-invariant 3D system (tests/test_2d.cpp). Energies refer to a slice of thickness dz (default dx).
	Geometries: disk (alias cylinder; Nx=Ny), channel (walls at y=0 and y=Ny-1, periodic x),
	box (walls x=0,Nx-1 -> boundary 0, y=0,Ny-1 -> boundary 1), bulk (periodic) and custom.
	These are periodic in z, so with dimension 3 and Nz>1 they are a cylinder, a channel and a
	square tube. sphere, slab, hemisphere and the wedges need three dimensions.
	plane_confinement W (10^6 J/m^3): energy W*Q_zz keeping the director in the plane. Without it,
	planar degenerate (fg) walls let the director escape along z, which is what the 3D physics does.
	It acts like a field along z with delta_epsilon<0, but independently of the field.
	Initial conditions: random, homogeneous, read_from_file, defects ("defect x y charge" lines,
	lattice units, director angle phi_i plus the defects' winding) and cholesteric_2d (helix axis in
	the plane at the angle phi_i, pitch p0_i, handedness of q0).
	examples/2d/ has input files for a defect pair in a disk, a quenched box and a cholesteric channel.

Backflow ("flow stokes" in the input file):
	Couples Q to the flow of the liquid crystal, in the Stokes limit (Beris-Edwards nematodynamics,
	Denniston, Orlandini and Yeomans, PRE 63, 056702 (2001)):
	    dQ/dt = Gamma H - u.grad Q + S(grad u, Q),    eta lap u - grad p + div sigma = 0,   div u = 0
	The stress sigma includes the distortion stress for every elastic constant (src/physics_ldg_gradient.h,
	generated by auxiliar/codegen/gen_gradient_derivative.py). The flow is solved with FFTs (pocketfft
	on the CPU, cuFFT on the GPU) in two kinds of domains:
	    periodic boxes: geometry bulk, 2D or 3D;
	    two flat walls: geometries channel (walls at y=0, Ny-1; 2D or 3D) and slab (walls at z=0, Nz-1),
	        with no-slip on the rows of surface sites (src/flow_walls.h). FFTs along the walls and, for
	        each Fourier mode, the exact discrete Stokes problem across the channel. Optional
	        "wall_velocity ux uy uz" (nm/us): the upper wall moves (Couette shear); it must be
	        parallel to the walls. The flow terms act on the bulk sites; the surface sites keep
	        their anchoring dynamics.
	Other geometries (sphere, disk, box, ...) are refused with the flow.
	Switched on by "flow stokes"; off by default ("flow none" or no flow line: no extra cost).
	Parameters, both required with the flow:
	    viscosity   isotropic viscosity, Pa s. The backflow is stronger for a smaller viscosity
	                compared with the rotational one (mu_1).
	    xi          flow-aligning parameter: the liquid crystal aligns with a shear flow for
	                xi > 3S/(2+S), at the angle theta with cos(2 theta) = 3S/((2+S) xi), and tumbles
	                below it (e.g. S=0.53: aligning for xi > 0.63, 13 degrees for xi=0.7). The start of
	                the run reports which regime the input is in.
	    activity    (optional, default 0) active stress -activity*q, q=2Q/3, in 10^6 Pa: active nematic.
	                activity>0 is extensile (bend instability of the uniform state; +1/2 defects
	                self-propel towards their head), activity<0 contractile. The active length
	                sqrt(L1/|activity|) should span several lattice spacings.
	                Checked by tests/test_active.cpp (instabilities, defect self-propulsion).
	With the flow, the snapshots have three more columns, ux,uy,uz: the velocity (nm/us) of the
	written state. Without it they are unchanged. The flow changes the dynamics, not the equilibrium states, so it needs a time
	integrator (DP5, RK2, Euler), not FIRE.
	tests/test_flow.cpp checks the Stokes solver, that the energy never increases with flow and that
	the equilibrium does not change. tests/test_flow_walls.cpp checks the solver between walls (exact
	Poiseuille and Couette profiles, the constrained-minimum property for random forces) and that the
	director aligns under shear at the flow-aligning angle. Between walls an extensile channel
	becomes unstable to a spontaneous shear flow above a threshold activity scaling as 1/width^2
	(Voituriez, Joanny and Prost, EPL 70, 404 (2005)). With backflow a +1/2/-1/2 pair annihilates sooner and closer to
	the -1/2, the +1/2 moving about twice as far (Toth, Denniston and Yeomans, PRL 88, 105504 (2002)).
	src/third_party/pocketfft_hdronly.h: pocketfft (BSD license, see its header), compiled only in
	src/flow_cpu.cpp.

GPU (NVIDIA, CUDA):
	"make gpu" builds marlics_gpu, which runs on the CPU or, with "device gpu" in the input file, on the GPU.
	FIRE, DP5, DP5s, RK2 and Euler run on the GPU, with every geometry except the experimental wedges.
	"make test-gpu" checks that the GPU and the CPU give the same forces, energies and relaxations.
	Old CUDA versions need an old host compiler and a generic architecture, e.g. CUDA 9.1:
	    make gpu CUDA_HOST=g++-6        (CUDA_ARCH defaults to compute_70 PTX, compiled by the driver)
	With CUDA >= 10 one can set the card's architecture, e.g. make gpu CUDA_ARCH=-arch=sm_75.
	gpu_precision double (default) or mixed: mixed evaluates the functional derivatives in float, about
	twice as fast per step, but its force error (~1e-6) slows down FIRE near min_force 1e-6; use it only
	with looser tolerances.
	GPU memory: FIRE needs about 3 arrays of 5*Nx*Ny*Nz doubles (301^3: ~4 GB), DP5 about 9 (301^3: ~10 GB).
	cuda/ holds the first CUDA port (2022, float only), kept as a reference.

Windows_makefile -> a makefile to compile a windows executable under a GNU environment

MarlicsVS/ ->folder containing the files required to compile  marlics in the virtual studio environment.
