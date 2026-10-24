#ifndef DP5_
#define DP5_

#include "integrator.h"
#include <stdio.h>
#include <stdlib.h>

//Embedded 7-stage FSAL Runge-Kutta pair: stage s uses a[s-1][0..s-1] (s=1..6, k_1 is stage 0);
//the last stage is the 5th-order solution, and e gives its error estimate from k_1..k_7.
struct ButcherTableau
{
  const char * name;
  double a[6][6];
  double e[7];
};

extern const ButcherTableau dormand_prince_m;  //DP5(4)7M, Dormand & Prince (1980)
extern const ButcherTableau dormand_prince_s;  //DP5(4)7S

//Adaptive Dormand-Prince integrator. Exec runs the element-wise work: CpuExec (OpenMP) or the
//CUDA executor.
template<class Exec>
class DP5 : public Integrator
{
 private:

  Exec exec;
  const ButcherTableau & tableau;
  const double Atol;
  const double Rtol;
  const double prefac;
  const double facmin;
  const double facmax;
  double *Qt, *k[7];

 public:

  virtual bool evolve( double *, double *, double);
  virtual const double * velocity(double * Qij) { return exec.host_velocity(Qij); }
  virtual ~DP5(void);
  DP5( GEOMETRY  * , const struct Simulation_Parameters *, const ButcherTableau & =dormand_prince_m );
};

#endif
