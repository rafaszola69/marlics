#ifndef HEMISPHERE_
#define HEMISPHERE_

#include "driver.h"
#include "geometry.h"

class Geometry_Hemisphere : public GEOMETRY
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
    

  Geometry_Hemisphere(const struct Simulation_Parameters * );
  ~Geometry_Hemisphere(void);

   private:

  virtual int * fill_point_type( void ) const ;  
  
  
};

#endif
