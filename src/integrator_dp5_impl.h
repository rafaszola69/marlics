// DP5 member definitions, included by the files that instantiate it for an executor.
#include "driver.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "geometry.h"
#include "integrator.h"
#include "integrator_dp5.h"

template<class Exec>
DP5<Exec>::DP5( GEOMETRY  * lc_pointer, const struct Simulation_Parameters *sim_param, const ButcherTableau & t ) : Integrator( lc_pointer),
                                   exec(lc_pointer, sim_param),
                                   tableau(t),
                                   Atol(sim_param->Atol),
                                   Rtol(sim_param->Rtol),
                                   prefac(sim_param->prefac),
                                   facmin(sim_param->facmin),
                                   facmax(sim_param->facmax)
{
  dt=sim_param->dt;
  Qt=exec.alloc(5L*Nx*Ny*Nz);
  for(int s=0; s<7; s++) k[s]=exec.alloc(5L*Nx*Ny*Nz);

  std::cout << "Dormand-Prince " << tableau.name << " on the " << Exec::name() << ".\n";
  std::cout << "dt=" << dt << " \n";

  if(sim_param->integrator_parameters_flag<5)
    {
      std::cout << "You defined " <<  sim_param->integrator_parameters_flag << " of the 5 optional paramters of the DOrmand-Prince integrator.\nThe remaining values will be set to their standardvalues.\n";
    }

  std::cout << "Atol=" << Atol << " \n";
  std::cout << "Rtol=" << Rtol << " \n";
  std::cout << "facmin=" << facmin << "\n";
  std::cout << "facmax=" << facmax << "\n";
  std::cout << "prefac=" << prefac << "\n\n";
};

template<class Exec>
DP5<Exec>::~DP5(void)
{
  exec.release(Qt);
  for(int s=0; s<7; s++) exec.release(k[s]);
}

//Stage values passed by value to the element-wise lambdas (no pointers to host objects).
struct DP5Stage
{
  const double * k[7];
  double c[7];
  int n;
};

//The step size persists between calls of evolve, so printing a snapshot does not restart it.
//A step shortened to land on tf does not reduce the step size of the next call.
template<class Exec>
bool DP5<Exec>::evolve( double * Qij_host, double *time, double tf )
{
  const long N=5L*Nx*Ny*Nz;
  double * const Q=exec.attach(Qij_host);
  double * const T=Qt;
  int information_step=1;
  double hfactor=1.0;

  exec.fill_ki(k[0],Q);

  while(*time<tf)
    {
      const bool clipped= (tf-*time) < dt;
      const double h= clipped ? tf-*time : dt;

      //Stages 2 to 7: Qt = Q + h*sum_j a_sj k_j.
      for(int s=1; s<7; s++)
        {
          DP5Stage st;
          st.n=s;
          for(int j=0; j<s; j++) { st.k[j]=k[j]; st.c[j]=h*tableau.a[s-1][j]; }
          exec.for_each(N, MARLICS_LAMBDA (long ll)
                        {
                          double sum=0;
                          for(int j=0; j<st.n; j++) sum+=st.c[j]*st.k[j][ll];
                          T[ll]=Q[ll]+sum;
                        });
          exec.fill_ki(k[s],T);
        }

      //Error estimate, scaled by the tolerances.
      DP5Stage er;
      er.n=7;
      for(int j=0; j<7; j++) { er.k[j]=k[j]; er.c[j]=tableau.e[j]; }
      const double atol=Atol, rtol=Rtol;
      const double global_error=exec.max(N, MARLICS_LAMBDA (long ll)
                                         {
                                           double local_error=0;
                                           for(int j=0; j<er.n; j++) local_error+=er.c[j]*er.k[j][ll];
                                           return fabs(local_error)/(atol+fabs(T[ll])*rtol);
                                         });

      const bool accepted= global_error<1.0;
      if(accepted)
        {
          *time+=h;

          if( information_step%10==0 )
            {
              std::cout << "time=" << *time << ", dt=" << h << ", global_error=" << global_error <<
              ", hfactor=" << hfactor << std::endl;
              information_step=0;
            }
          information_step++;
        }

      hfactor=fmin(facmax,fmax(facmin,prefac*pow(global_error,-0.2000000000)));
      if(!(accepted && clipped)) dt=h*hfactor;
      if(dt<1e-14)
        {
          printf("Convergence problem!!!\n Aborting execution!!!\n");
          exit(5);
        }

      if(accepted)
        {
          exec.for_each(N, MARLICS_LAMBDA (long ll) { Q[ll]=T[ll]; });
          //First same as last: the last stage is the derivative at the new Q.
          double * k1=k[0]; k[0]=k[6]; k[6]=k1;
        }
    }

  std::cout << "time=" << *time << ", Energy=" << exec.energy(Q) << std::endl;
  exec.detach(Qij_host,Q);
  return true;
};
