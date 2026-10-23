// Euler member definitions, included by the files that instantiate it for an executor.
#include "driver.h"
#include <cstdio>
#include "geometry.h"
#include "integrator.h"
#include "integrator_euler.h"

template<class Exec>
Euler<Exec>::Euler( GEOMETRY  * lc_pointer, const struct Simulation_Parameters *sim_param ) : Integrator( lc_pointer), exec(lc_pointer)
{
  dt=sim_param->dt;
  k_1=exec.alloc(5L*Nx*Ny*Nz);
  std::cout << "Euler on the " << Exec::name() << ".\n";
  std::cout << "dt=" << dt << " \n";
};

template<class Exec>
Euler<Exec>::~Euler(void)
{
  exec.release(k_1);
}

template<class Exec>
bool Euler<Exec>::evolve( double * Qij_host, double *time, double tf )
{
  const long N=5L*Nx*Ny*Nz;
  double * const Q=exec.attach(Qij_host);
  double * const K1=k_1;
  const double h=dt;
  int information_step=1;

  while(*time<tf)
    {
      exec.fill_ki(K1,Q);
      exec.for_each(N, MARLICS_LAMBDA (long ll) { Q[ll]+=h*K1[ll]; });

      *time+=dt;
      if( information_step%10==0 )
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
