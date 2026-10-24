#ifndef GEOMETRY_2D_
#define GEOMETRY_2D_
// Geometries bounded in the xy plane and periodic along z. With "dimension 2" (Nz=1) they are
// two-dimensional; with Nz>1 they are their z-invariant 3D versions (cylinder, channel, square tube).

#include "driver.h"
#include "geometry.h"

//Disk of radius R (2D droplet; infinite cylinder along z in 3D). One boundary (0), radial normal.
class Geometry_Disk : public GEOMETRY
{
 public:

  Geometry_Disk(const struct Simulation_Parameters * );
  ~Geometry_Disk(void) {};
  virtual int surface_site(int, int, int, double [3], double &) const ;
  virtual void surface_stencil(int, int, int, int [3], int [3]) const ;
  virtual bool is_periodic(int axis) const { return axis==2; }
  virtual int backwardPTP(double *, const double * ) const { return 0; }

 private:

  double R_in, R_ex;
  const double HNx, HNy;
  virtual int * fill_point_type( void ) const ;
};

//Channel with walls at y=0 (boundary 0, normal -y) and y=Ny-1 (boundary 1, normal +y), periodic in x.
class Geometry_Channel : public GEOMETRY
{
 public:

  Geometry_Channel(const struct Simulation_Parameters * );
  ~Geometry_Channel(void) {};
  virtual int surface_site(int, int, int, double [3], double &) const ;
  virtual void surface_stencil(int, int, int, int [3], int [3]) const ;
  virtual bool is_periodic(int axis) const { return axis!=1; }
  virtual int backwardPTP(double *, const double * ) const { return 0; }

 private:

  virtual int * fill_point_type( void ) const ;
};

//Rectangle with walls at x=0, Nx-1 (boundary 0, normals -x, +x) and y=0, Ny-1 (boundary 1, normals
//-y, +y). The corners belong to boundary 0, with diagonal normals.
class Geometry_Box : public GEOMETRY
{
 public:

  Geometry_Box(const struct Simulation_Parameters * );
  ~Geometry_Box(void) {};
  virtual int surface_site(int, int, int, double [3], double &) const ;
  virtual void surface_stencil(int, int, int, int [3], int [3]) const ;
  virtual bool is_periodic(int axis) const { return axis==2; }
  virtual int backwardPTP(double *, const double * ) const { return 0; }

 private:

  virtual int * fill_point_type( void ) const ;
};

#endif
