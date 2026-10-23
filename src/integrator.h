#ifndef INTEGRATOR_
#define INTEGRATOR_

#include <cstdlib>          
#include "geometry.h"


class Integrator
{

 protected:

  const int Nx, Ny, Nz;
  double dt;
  
  GEOMETRY * sample_geometry;
  virtual ~Integrator() {};
  
 public:
  
  virtual bool evolve(double *,  double *, double) = 0;
  
  //Work array for the Runge-Kutta stages; allocated only by the integrators that use it.
  double *Qtij;

 Integrator( GEOMETRY  * lc_pointer ) :  sample_geometry(lc_pointer),

    Nx(lc_pointer->Nx),
    Ny(lc_pointer->Ny),
    Nz(lc_pointer->Nz),
    Qtij(NULL)
    {
    };
    
};
 
#endif
