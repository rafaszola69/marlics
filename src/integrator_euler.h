#ifndef EULER_
#define EULER_

#include "integrator.h"
#include <stdio.h>
#include <stdlib.h>

//Euler integrator with a fixed step. Exec runs the element-wise work: CpuExec (OpenMP) or the
//CUDA executor.
template<class Exec>
class Euler : public Integrator
{
 private:

  Exec exec;
  double *k_1;

 public:

  virtual bool evolve( double *, double *, double);
  virtual const double * velocity(double * Qij) { return exec.host_velocity(Qij); }
  virtual ~Euler(void);
  Euler( GEOMETRY  * , const struct Simulation_Parameters * );
};

#endif
