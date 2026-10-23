#ifndef RK2_
#define RK2_

#include "integrator.h"
#include <stdio.h>
#include <stdlib.h>

//Runge-Kutta2 (midpoint) integrator with a fixed step. Exec runs the element-wise work: CpuExec (OpenMP) or the
//CUDA executor.
template<class Exec>
class RK2 : public Integrator
{
 private:

  Exec exec;
  double *k_1, *k_2;

 public:

  virtual bool evolve( double *, double *, double);
  virtual ~RK2(void);
  RK2( GEOMETRY  * , const struct Simulation_Parameters * );
};

#endif
