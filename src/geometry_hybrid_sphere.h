#ifndef HYBRID_SPHERE_
#define HYBRID_SPHERE_

#include "driver.h"
#include "geometry.h"

class Geometry_Hybrid_Sphere : public GEOMETRY
{
  

 private:

   double R_ex;
   double R_in;

  const double HNx;
  const double HNy;
  const double HNz;


  const double dx;
  const double dy;
  const double dz;

  
 public:
  virtual void surface_stencil(int, int, int, int [3], int [3]) const ;
  virtual int surface_site(int, int, int, double [3], double &) const ;
  virtual int backwardPTP(double *, const double * ) const ;
     

  Geometry_Hybrid_Sphere(const struct Simulation_Parameters * );
  ~Geometry_Hybrid_Sphere(void);

   private:

  virtual int * fill_point_type( void ) const ;  
  
  
};

#endif
