#ifndef geometry_

#define geometry_
#include "energy.h"
#include "energy_ldg.h"
#include "lattice.h"

#include <vector>
#include <string>
#include <cstring>


class GEOMETRY
{

 public:
  
  const int Nx;
  const int Ny;
  const int Nz;
  const double dx_1;
  const double dy_1;
  const double dz_1;
  
  void read_check(int , int );
  
  //Right-hand side of the relaxation equation at every site, shared by all integrators. The default
  //uses the lattice description (lattice.h); only the experimental wedge geometries override it.
  virtual void fill_ki(double * k_i, const  double * Qij) const;
  //Energy of each site (10^-21 J): k_i[2*s] = bulk part, k_i[2*s+1] = anchoring part. Shared by all geometries.
  virtual void Energy_calc(double * k_i, const  double * Qij) const;
  virtual int backwardPTP(double *, const double * ) const =0 ;

  const int *point_type;
  int number_of_boundaries;  
  
  virtual void boundary_init( struct Simulation_Parameters * );
  virtual ~GEOMETRY() {};

 protected:

  std::string geometry_name;
  std::vector<class BOUNDARY *> bc_conditions;
  std::string boundary_needed_to_be_defined;

  const landau_de_gennes bulk_energy;

  //Energy_calc returns energies (10^-21 J): bulk densities times dV, anchoring densities times dA.
  const double dV;
  double dA;  //surface area represented by one surface site; geometries with curved surfaces override it
  virtual bool is_periodic(int axis) const { return false; }

  //For a surface site (point_type>=2): returns the boundary condition index, and fills the
  //outward normal used by fill_ki and the surface area the site represents.
  virtual int surface_site(int i, int j, int k, double v[3], double & area) const;

  //For a surface site: the sites used for its first derivatives along each axis, dQ/dx_a from
  //Q[np[a]]-Q[nm[a]] (central when both differ from the site, one-sided otherwise).
  virtual void surface_stencil(int i, int j, int k, int nm[3], int np[3]) const;

  //Lattice description shared by the CPU and GPU code, built by boundary_init (build_lattice).
  std::vector<SurfaceSite> surface;
  std::vector<BoundaryParams> boundary_params;
  void build_lattice( void );

 public:

  LatticeView view( void ) const;

 protected:

  virtual int * fill_point_type( void ) const = 0 ;
  GEOMETRY * geometry_pointer;
  GEOMETRY(const struct Simulation_Parameters * );
};

#endif
