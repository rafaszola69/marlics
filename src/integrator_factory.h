#ifndef INTEGRATOR_FACTORY_H
#define INTEGRATOR_FACTORY_H
// Integrators available for an executor (CpuExec or the CUDA executor), by input file name.
// Returns NULL for a name that is not one of them.
#include <cstring>
#include <iostream>
#include "driver.h"
#include "integrator.h"
#include "fire_integrator.h"
#include "integrator_dp5.h"
#include "integrator_rk2.h"
#include "integrator_euler.h"

template<class Exec>
Integrator * make_integrator(GEOMETRY * geometry, const Simulation_Parameters * sim_param)
{
  const char * name=sim_param->integrator_type;

  if ( strcasecmp(name,"DP5m") == 0 || strcasecmp(name,"DP5") == 0 || strcasecmp(name,"Dormand-Prince_m") == 0
       || strcasecmp(name,"Rk54m") == 0 || strcasecmp(name,"Rk5m") == 0 )
    {
      std::cout << "Loading the Dormand-Prince5(4)m integrator.\n\n";
      return new DP5<Exec>(geometry, sim_param, dormand_prince_m);
    }
  if ( strcasecmp(name,"DP5s") == 0 || strcasecmp(name,"Dormand-Prince_s") == 0
       || strcasecmp(name,"Rk54s") == 0 || strcasecmp(name,"Rk5s") == 0 )
    {
      std::cout << "Loading the Dormand-Prince5(4)s integrator.\n\n";
      return new DP5<Exec>(geometry, sim_param, dormand_prince_s);
    }
  if ( strcasecmp(name,"RK2") == 0 || strcasecmp(name,"Runge-Kutta2") == 0 )
    {
      std::cout << "Loading the Runge-Kutta2 integrator.\n\n";
      return new RK2<Exec>(geometry, sim_param);
    }
  if ( strcasecmp(name,"RK1") == 0 || strcasecmp(name,"Euler") == 0 )
    {
      std::cout << "Loading the Euler integrator.\n\n";
      return new Euler<Exec>(geometry, sim_param);
    }
  if ( strcasecmp(name,"FIRE") == 0 )
    {
      std::cout << "Loading the FIRE integrator.\n\n";
      return new FIRE<Exec>(geometry, sim_param);
    }
  if ( strcasecmp(name,"FIRE_el") == 0 )
    {
      //Kept for old input files: same as FIRE, with the alpha values the former FIRE_el used.
      std::cout << "Loading the FIRE integrator (FIRE_el alpha parameters).\n\n";
      return new FIRE<Exec>(geometry, sim_param, 0.4, 0.3);
    }
  return NULL;
}

#endif
