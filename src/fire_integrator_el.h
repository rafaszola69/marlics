#ifndef FIRE_EL_
#define FIRE_EL_

#include "integrator.h"
#include <stdio.h>           
#include <stdlib.h>

class FIRE_el : public Integrator
{

 private:

 public:

  

  double *force, *vQij, v2, f2, *Qijt, *energy, min_force;
  double dtmin, potence, alphainit,scale, alphaMin, dt_init;
  double scaleAlpha, scaleP, dtMax, scaleM, alpha;
  int Np, Nmin, proceed;
  
  virtual bool evolve( double *, double *, double);    
  virtual ~FIRE_el(void) {};
  FIRE_el( GEOMETRY  * , const struct Simulation_Parameters * );


};


#endif
