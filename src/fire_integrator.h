#ifndef FIRE_
#define FIRE_

#include "integrator.h"
#include <stdio.h>
#include <stdlib.h>

//FIRE 2.0 minimizer (Guénolé et al., Comput. Mater. Sci. 175, 109584 (2020)).
//Exec runs the element-wise work: CpuExec (OpenMP) or the CUDA executor.
template<class Exec>
class FIRE : public Integrator
{
 private:

  Exec exec;
  double *force, *vQij;
  double min_force;

  //Time step and velocity mixing control. They persist between calls of evolve,
  //so printing a snapshot does not restart the minimization.
  double dtmin, dtMax, scaleP, scaleM;
  double alpha, alphainit, alphaMin, scaleAlpha;
  int Np, Nmin, steps;

 public:

  virtual bool evolve( double *, double *, double);
  virtual ~FIRE(void);
  FIRE( GEOMETRY  * , const struct Simulation_Parameters *, double alpha_init=0.3, double alpha_min=0.1 );
};

#endif
