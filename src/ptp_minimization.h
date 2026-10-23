#ifndef PTP_
#define PTP_

#include "integrator.h"
#include <stdio.h>           
#include <stdlib.h>

class PTP_min : public Integrator
{

 private:

 public:

  double *force, *vQij, v2, f2, *Qijt, *energy, min_force;
  double dtmin, potence, alphainit,scale, alphaMin, dt_init;
  double scaleAlpha, scaleP, dtMax, scaleM, alpha;
  int Np, Nmin, proceed;
  
  virtual bool evolve( double *, double *, double);    
  virtual ~PTP_min(void) {};
  PTP_min( GEOMETRY  * , const struct Simulation_Parameters * );


};


#endif
