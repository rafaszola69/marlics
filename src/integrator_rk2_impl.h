// RK2 member definitions, included by the files that instantiate it for an executor.
#include "driver.h"
#include <cstdio>
#include "geometry.h"
#include "integrator.h"
#include "integrator_rk2.h"

template<class Exec>
RK2<Exec>::RK2( GEOMETRY  * lc_pointer, const struct Simulation_Parameters *sim_param ) : Integrator( lc_pointer), exec(lc_pointer)
{
  dt=sim_param->dt;
  Qtij=exec.alloc(5L*Nx*Ny*Nz);
  k_1=exec.alloc(5L*Nx*Ny*Nz);
  k_2=exec.alloc(5L*Nx*Ny*Nz);
  std::cout << "Runge-Kutta2 on the " << Exec::name() << ".\n";
  std::cout << "dt=" << dt << " \n";
};

template<class Exec>
RK2<Exec>::~RK2(void)
{
  exec.release(Qtij);
  exec.release(k_1);
  exec.release(k_2);
}

template<class Exec>
bool RK2<Exec>::evolve( double * Qij_host, double *time, double tf )
{
  const long N=5L*Nx*Ny*Nz;
  double * const Q=exec.attach(Qij_host);
  double * const T=Qtij;
  double * const K1=k_1;
  double * const K2=k_2;
  const double h=dt;
  int information_step=1;

  while(*time<tf)
    {
      //1st Stage:
      exec.fill_ki(K1,Q);

      //2nd Stage:
      exec.for_each(N, MARLICS_LAMBDA (long ll) { T[ll]=Q[ll]+0.5*h*K1[ll]; });
      exec.fill_ki(K2,T);

      exec.for_each(N, MARLICS_LAMBDA (long ll) { Q[ll]+=h*K2[ll]; });

      *time+=dt;
      if( information_step%25==0 )
        {
          std::cout << "time=" << *time << ", dt=" << dt << std::endl;
          information_step=0;
        }
      information_step++;
    }

  std::cout << "time=" << *time << ", Energy=" << exec.energy(Q) << std::endl;
  exec.detach(Qij_host,Q);
  return true;
};
