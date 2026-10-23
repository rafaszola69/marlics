// FIRE member definitions, included by the files that instantiate it for an executor.
#include "driver.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include "geometry.h"
#include "integrator.h"
#include "fire_integrator.h"

//NaN/inf test on the bit pattern: with -Ofast (finite-math-only) std::isfinite is optimized away.
static bool not_finite(double x)
{
  uint64_t bits;
  memcpy(&bits,&x,sizeof(bits));
  return ((bits>>52)&0x7ff)==0x7ff;
}

//FIRE minimizer. The force is the same right-hand side used by the time integrators (fill_ki),
//so FIRE converges to the same equilibrium as DP5 & co., boundary conditions included.
template<class Exec>
FIRE<Exec>::FIRE( GEOMETRY  * lc_pointer, const struct Simulation_Parameters *sim_param, double alpha_init, double alpha_min ) : Integrator( lc_pointer), exec(lc_pointer)
{
  dt=sim_param->dt;
  force=exec.alloc(5L*Nx*Ny*Nz);
  vQij=exec.alloc(5L*Nx*Ny*Nz);

  dtmin = dt*0.01;
  dtMax = dt*100;
  scaleP = 1.1;      //dt increase factor after Nmin consecutive downhill steps
  scaleM = 0.5;      //dt decrease factor on an uphill step
  Nmin = 20;         //also the initial delay before dt can be decreased
  alphainit = alpha_init;
  alphaMin = alpha_min;
  scaleAlpha = 0.99;
  alpha = alphainit;
  Np = 0;
  steps = 0;
  min_force=sim_param->min_force;
  printf("Initializing simulation with FIRE(Fast Inertial Relaxation Engine) method on the %s.\n", Exec::name());
  std::cout << "dt=" << dt << ", dtmin=" << dtmin << ", dtmax=" << dtMax << ", alpha=" << alphainit << ", alpha_min=" << alphaMin << ", min_force=" << min_force << " \n";
};

template<class Exec>
FIRE<Exec>::~FIRE(void)
{
  exec.release(force);
  exec.release(vQij);
}

template<class Exec>
bool FIRE<Exec>::evolve( double * Qij_host, double *time, double tf )
{
  const long N=5L*Nx*Ny*Nz;
  double * const Q=exec.attach(Qij_host);
  double * const F=force;
  double * const V=vQij;
  bool converged=false;

  exec.fill_ki(F,Q);

  while(*time<tf && !converged)
    {
      const SumSumMax r=exec.sum_sum_max(N, MARLICS_LAMBDA (long ll, double& p, double& f2, double& fmax)
                                          {
                                            p=V[ll]*F[ll];
                                            f2=F[ll]*F[ll];
                                            fmax=fabs(F[ll]);
                                          });
      const double potence=r.sum0, f2=r.sum1, f2max=r.max;

      if(not_finite(f2))
        {
          printf("NaN found in the system at time %g.\nAborting the program.\n",*time);
          exit(10);
        }
      if(f2max<min_force)
        {
          converged=true;
          printf("Stable condition achivied at time %g\n",*time);
          break;
        }
      else if(potence>0)
        {
          Np++;
          if(Np>Nmin)
            {
              dt=fmin(dt*scaleP,dtMax);
              alpha=fmax(alphaMin,scaleAlpha*alpha);
            }
        }
      else
        {
          Np=0;
          if(steps>=Nmin)
            {
              dt=fmax(dt*scaleM,dtmin);
              alpha=alphainit;
            }
        }

      const double h=dt, a=alpha;

      //Uphill: step half a step back and stop.
      if(potence<=0)
        exec.for_each(N, MARLICS_LAMBDA (long ll)
                      {
                        Q[ll]-=0.5*h*V[ll];
                        V[ll]=0;
                      });

      //Semi-implicit Euler with FIRE velocity mixing.
      const double v2=exec.sum(N, MARLICS_LAMBDA (long ll)
                               {
                                 V[ll]+=h*F[ll];
                                 return V[ll]*V[ll];
                               });
      const double mix=a*sqrt(v2/f2);
      exec.for_each(N, MARLICS_LAMBDA (long ll)
                    {
                      V[ll]=(1-a)*V[ll]+mix*F[ll];
                      Q[ll]+=h*V[ll];
                    });

      exec.fill_ki(F,Q);

      *time+=dt;
      steps++;
      if( steps%50==0 )
        {
          printf("time=%g, dt=%g, force=%g, max_force=%g, velocity=%g, P=%g, alpha=%g, Np=%d\n",*time, dt, sqrt(f2), f2max, sqrt(v2), potence, alpha, Np);
          fflush(stdout);
        }
    }

  std::cout << "time=" << *time << ", Energy=" << exec.energy(Q) << std::endl;
  exec.detach(Qij_host,Q);
  return !converged;
};
