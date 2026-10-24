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

  //With backflow: the flow velocity for the state Qij (host array, 3 per site, nm/us), for the
  //snapshots. NULL without flow.
  virtual const double * velocity(double * Qij) { return NULL; }
  
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
