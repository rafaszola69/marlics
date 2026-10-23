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

	auxiliar/renderize/ -> contain routines to renderize the directors and the isosurface of the parameter order S;

	auxiliar/boundaries/ -> contain 2 examples off c++ files that can be used to produce boundary files;

examples/ -> a series of MarLiCs usage examples;

notebook/ -> Some Mathematica nb files used to evaluate symbolically the functional derivatives.

tests/ -> consistency tests, run with "make test" (about 15 s). They check that the dynamics (fill_ki) is the gradient of the reported energy, that the anchoring energies match their functional derivatives, that DP5 never increases the energy and that FIRE converges to a stationary state. Run them after changing energies, boundaries, geometries or integrators;
	
makefile -> a makefile used to compile the code (make), or to compile and run the tests (make test)

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
