#ifndef LC_DRIVER_
#define LC_DRIVER_
#ifdef _WIN32
#define strcasecmp _stricmp
#endif
#include <iostream>
#define ERROr std::cout<<"cannot allocate memory\n"; exit(4);
#include <vector>
#include <string>
#include <map>

enum class parameter_status {unset=0, set=1};
enum class elastic_const_status {unset=0, Ls=1, Ks=2};
enum class chirality_status {unset=0, p0=1, q0=2};
enum class viscosity_status {unset=0, mu=1,gamma=2};

struct Simulation_Parameters
{

  //Domain Paramaters:
  
  int Nx;              /* grid size*/
  int Ny;              /* grid size*/
  int Nz;              /* grid size*/
  parameter_status grid_flag[3]={parameter_status::unset,
                    parameter_status::unset,
                    parameter_status::unset};
    
  double dx;           /* 10^-9 m */
  double dy;           /* 10^-9 m */
  double dz;           /* 10^-9 m */
  parameter_status grid_spacing_flag[3]={parameter_status::unset,
                    parameter_status::unset,
                    parameter_status::unset};

  
  char geometry[200];
  //dimension 2: a single layer (Nz=1) with no variation along z; Q keeps its 5 components (3D physics).
  int dimension=3;
  //Energy W*Q_zz (10^6 J/m^3) keeping the director in the xy plane for W>0, like a field along z with delta_epsilon<0.
  double plane_confinement=0;
  int Nz2=0;
  double wTheta=0;
  parameter_status geometry_flag=parameter_status::unset;
  
  /*Just for use in Sphere geometry */
  double R_in;           /* 10^-9 m */
  double R_ex;           /* 10^-9 m */
  parameter_status radius_flag[2]={parameter_status::unset,
                                   parameter_status::unset};
  /*Ellipse geometry: semi-axes along x and y (default Nx*dx/2, Ny*dy/2) */
  double ellipse_a=0;    /* 10^-9 m */
  double ellipse_b=0;    /* 10^-9 m */
  /*Disk geometry: fraction of the perimeter, an arc centred on +x, that belongs to boundary 1 (0: one boundary) */
  double hybrid_arc=0;
  /*Polygon geometry: number of sides of the regular polygon inscribed in the circle of the disk (one vertex on +x) */
  int polygon_sides=0;
  
  //Lc Parameters:

  double a=0;            /* 10^6 J/Km^3 */
  double B=0;            /* 10^6 J/m^3  */
  double C=0;            /* 10^6 J/m^3  */
  double L1=0;           /* 10^-12 N    */
  double L2=0;           /* 10^-12 N    */
  double L3=0;           /* 10^-12 N    */
  double Lq=0;           /* 10^-12 N    */
  double Ls=0;           /* 10^-12 N    */
  double k11=0, k22=0;      /* 10^-12 N    */
  double k33=0, k24=0;     /* 10^-12 N    */
  double p0=0;          /* 10^-9 m     */
  double q0=0;          /* 10^9 1/m    */
  double T=-1;            /* (T-T*)  */
  elastic_const_status elastic_flag=elastic_const_status::unset;
  chirality_status chirality_flag=chirality_status::unset;
  parameter_status thermal_flag[4]={parameter_status::unset,
                    parameter_status::unset,
                    parameter_status::unset,
                    parameter_status::unset};
  
  double mu_1=0;         /* Pa s    */
  double mu_1_s=0;       /* Pa s    */

  double gamma_1=0;
  double gamma_1_s=0;
  viscosity_status viscosity_flag[2]={viscosity_status::unset,
                      viscosity_status::unset};

  
  //Anchoring parameters:
  
  char bound_fname[200];
  std::map<int,double> theta_0;      /* degrees */
  std::map<int,double> phi_0;        /* degrees */  
  std::map<int,double> Wo1;         /* 10^-3 m */
  std::map<int,double> Wo2;         /* 10^-3 m */

  std::map<int,std::string> anchoring_type;
  std::map<int,std::string> strong_type;  // Tipo de ancoramento forte (planar, etc.)
  
  //electric field parameters:
  parameter_status field_flag[4]={parameter_status::unset,    // [0] elecfieldx flag
                                  parameter_status::unset,    // [1] elecfieldy flag
                                  parameter_status::unset,    // [2] elecfieldz flag
                                  parameter_status::unset};   // [3] deltaepslon flag
  double elecfieldx=0;
  double elecfieldy=0;
  double elecfieldz=0;
  double deltaepslon=0;
  
  //Initial Conditions:

  char   initial_conditions[200];
  char   ic_file_name[200];
  int    rng_seed;
  double theta_i;       /* degrees */
  double phi_i;         /* degrees */ 
  double p0_i;         /* nm */ 
  //"defects" initial condition: in-plane point defects at (x,y) (lattice units) with charges k.
  std::vector<double> defect_x, defect_y, defect_k;
  double R0_i;         /* nm */ 
  int    m_i;         /* nm */ 
  parameter_status ic_flag[8]={parameter_status::unset,     // [0] initinal conditions flag
                               parameter_status::unset,     // [1] ic file name  flag
                               parameter_status::unset,     // [2] theta_i flag
                               parameter_status::unset,     // [3] phi_i flag
                               parameter_status::unset,     // [4] rng_seed flag
                               parameter_status::unset,     // [5] p0_i flag
                               parameter_status::unset,     // [6] m_i flag
                               parameter_status::unset,     // [7] R0_i flag
                             };
  
  //Integrator Paramters:

  char integrator_type[200];  
  char device[20]="cpu";    //cpu or gpu
  char gpu_precision[20]="double";  //double or mixed (functional derivatives in float)
  //Backflow: "none", or "stokes" (Stokes flow in a periodic box, between flat walls or in any other
  //geometry, see flow_physics.h).
  char flow[20]="none";
  double viscosity=0.;      /* isotropic viscosity, Pa s */
  double xi=0.;             /* flow-aligning parameter (no default: aligning or tumbling depends on it) */
  bool xi_set=false;
  double activity=0.;       /* active stress -activity*q (10^6 Pa); needs the flow */
  double chiral_activity=0.; /* chiral active stress -chiral_activity*(eps.q-q.eps)/2, eps_ij=eps_ijz (10^6 Pa) */
  double elastic_stress=1.;  /* 1, or 0 to drop the passive stress from the flow (active stress only) */
  double xi_wall=0., xi_wall_width=0.; /* xi_wall at the bulk sites within xi_wall_width (nm) of a surface site */
  double xi_wall_inner=0.;  /* ... and farther than xi_wall_inner (nm): a ring instead of a layer */
  double xi_core=0., xi_core_S=0.; /* xi_core at the bulk sites with S < xi_core_S*S_eq (defect cores; S from tr Q^2) */
  double wall_velocity[3]={0.,0.,0.};  /* velocity of the upper wall (channel, slab), nm/us */
  int flow_solver=0;        /* confined geometries: 0 automatic, 1 direct, 2 iterative (flow_confined.h) */
  //Surface sites: "local" relaxes them with the local boundary condition (one-sided normal
  //derivative), "variational" with the exact gradient of the discrete energy (lattice.h).
  bool surface_variational=false;
  double Atol=0.001;
  double Rtol=0.001;
  double prefac=0.8;
  double facmin=0.4;
  double facmax=3.0;
  double noise_factor=0.01;
  double min_force=1e-9;
  parameter_status integrator_flag=parameter_status::unset;
  int integrator_parameters_flag=0;
  
  //Output name parameters:
  
  char output_folder[200]=".";
  char output_fname[200]="director_field_$$.csv";
  parameter_status output_name_status[2]={parameter_status::unset, parameter_status::unset};
  
  //Execution Parameters:
  double timeprint;    /*10^-6 s */
  double timeprint_increase_factor; /*linear=10^-6 s, logarithmic=admensional */
  char print_time_type[200];
  int first_output_file_number;  
  parameter_status timeprint_status[4]={parameter_status::unset,
                    parameter_status::unset,
                    parameter_status::unset,
                    parameter_status::unset};

  
  double ti;           /*10^-6 s */
  double tf;           /*10^-6 s */
  double dt;           /*10^-6 s */
  parameter_status time_status[3]={parameter_status::unset,
                   parameter_status::unset,
                   parameter_status::unset};
  
  double S_eq;
  double * Q_00;

  
} ;

class driver
  {

  public:

  
    struct Simulation_Parameters sim_param;
    class GEOMETRY * LcS_Geometry;
    class Integrator * LcS_Integrator;
    class container * Data_Container;
    double * Qij;
 
    void update_timeprint (void);
    void (*next_time_print)  (double *, double);
    void setup_LC ( void);
    void setup_Simulation ( void);
    int parse_input_file(char[]);
    void error_check(int error_handler, char parser[]);

    driver(void);
    ~driver(void){};

  };

  void init_container(struct Simulation_Parameters * sim_param);
  void write_state(double , const double  *, const int *);
  void finish_container( void );


  void log_next_time_print(double *, double );
  void linear_next_time_print( double *, double );







#endif
