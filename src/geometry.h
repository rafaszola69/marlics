#ifndef geometry_

#define geometry_
#include "energy.h"
#include "energy_ldg.h"

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
  
  virtual void fill_ki(double * k_i, const  double * Qij) const = 0;
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
  double bulk_energy_octant(const double * Qij, int s, int xm, int xp, int ym, int yp, int zm, int zp, const double v[3]) const;

  //Neighbour of (i,j,k) along axis (0,1,2) in direction dir (+1/-1). Outside points and
  //non-periodic domain edges return the site itself (zero one-sided difference).
  int neighbor(int i, int j, int k, int axis, int dir) const;
  virtual bool is_periodic(int axis) const { return false; }

  //For a surface site (point_type>=2): returns the boundary condition index, and fills the
  //outward normal used by fill_ki and the surface area the site represents.
  virtual int surface_site(int i, int j, int k, double v[3], double & area) const;

  virtual int * fill_point_type( void ) const = 0 ;
  GEOMETRY * geometry_pointer;
  GEOMETRY(const struct Simulation_Parameters * );
};

#endif
