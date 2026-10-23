#include "driver.h"
#include "omp.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdint>
#include "geometry.h"
#include "integrator.h"
#include "fire_integrator.h"

#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)

//NaN/inf test on the bit pattern: with -Ofast (finite-math-only) std::isfinite is optimized away.
static bool not_finite(double x)
{
  uint64_t bits;
  memcpy(&bits,&x,sizeof(bits));
  return ((bits>>52)&0x7ff)==0x7ff;
}

//FIRE minimizer. The force is the same right-hand side used by the time integrators (fill_ki),
//so FIRE converges to the same equilibrium as DP5 & co., boundary conditions included.
FIRE::FIRE( GEOMETRY  * lc_pointer, const struct Simulation_Parameters *sim_param, double alpha_init, double alpha_min ) : Integrator( lc_pointer)
{
  dt=sim_param->dt;
  if((force= (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  if((vQij = (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  if((energy = (double *)calloc(2*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}

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
  printf("Initializing simulation with FIRE(Fast Inertial Relaxation Engine) method.\n");
  std::cout << "dt=" << dt << ", dtmin=" << dtmin << ", dtmax=" << dtMax << ", alpha=" << alphainit << ", alpha_min=" << alphaMin << ", min_force=" << min_force << " \n";
};

bool FIRE::evolve( double * Qij, double *time, double tf )
{
  int ll;
  const int N=5*Nx*Ny*Nz;
  double Total_Energy;
  double potence, v2, f2, f2max;
  bool converged=false;

  #pragma omp parallel default(shared) private(ll)
  {
    sample_geometry->fill_ki(force,Qij);
    #pragma omp barrier

    while(*time<tf && !converged)
      {
        #pragma omp single
          { potence=0; f2=0; f2max=0; }

        #pragma omp for schedule(static) reduction(+: potence, f2) reduction(max: f2max)
          for( ll=0; ll<N;ll++)
            {
              potence+=vQij[ll]*force[ll];
              f2+=force[ll]*force[ll];
              f2max=MAX(f2max,fabs(force[ll]));
            }

        #pragma omp single
          {
            if(not_finite(f2))
              {
                printf("NaN found in the system at time %g.\nAborting the program.\n",*time);
                exit(10);
              }
            if(f2max<min_force)
              {
                converged=true;
                printf("Stable condition achivied at time %g\n",*time);
              }
            else if(potence>0)
              {
                Np++;
                if(Np>Nmin)
                  {
                    dt=min(dt*scaleP,dtMax);
                    alpha=MAX(alphaMin,scaleAlpha*alpha);
                  }
              }
            else
              {
                Np=0;
                if(steps>=Nmin)
                  {
                    dt=MAX(dt*scaleM,dtmin);
                    alpha=alphainit;
                  }
              }
          }

        if(!converged)
          {
            //Uphill: step half a step back and stop.
            if(potence<=0)
              {
                #pragma omp for schedule(static)
                  for( ll=0; ll<N;ll++)
                    {
                      Qij[ll]-=0.5*dt*vQij[ll];
                      vQij[ll]=0;
                    }
              }

            //Semi-implicit Euler with FIRE velocity mixing.
            #pragma omp single
              v2=0;
            #pragma omp for schedule(static) reduction(+: v2)
              for( ll=0; ll<N;ll++)
                {
                  vQij[ll]+=dt*force[ll];
                  v2+=vQij[ll]*vQij[ll];
                }

            #pragma omp for schedule(static)
              for( ll=0; ll<N;ll++)
                {
                  vQij[ll]=(1-alpha)*vQij[ll]+alpha*sqrt(v2/f2)*force[ll];
                  Qij[ll]+=dt*vQij[ll];
                }

            sample_geometry->fill_ki(force,Qij);

            #pragma omp single
              {
                *time+=dt;
                steps++;
                if( steps%50==0 )
                  {
                    printf("time=%g, dt=%g, force=%g, max_force=%g, velocity=%g, P=%g, alpha=%g, Np=%d\n",*time, dt, sqrt(f2), f2max, sqrt(v2), potence, alpha, Np);
                    fflush(stdout);
                  }
              }
          }
      }

    Total_Energy=0;
    sample_geometry->Energy_calc(energy,Qij);
    #pragma omp for schedule(static) reduction(+: Total_Energy)
      for( ll=0; ll<2*Nx*Ny*Nz;ll++) Total_Energy+=energy[ll];

    #pragma omp single
      std::cout << "time=" << *time << ", Energy=" << Total_Energy << std::endl;
  }
  return !converged;
};
