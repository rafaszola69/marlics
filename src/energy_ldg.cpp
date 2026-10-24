#include "driver.h"
#include "physics_ldg.h"
#include "energy.h"
#include "energy_ldg.h"
#include <stdio.h>
#include <math.h>

landau_de_gennes::landau_de_gennes(const struct Simulation_Parameters * lc) : ENERGY(lc)
{
  std::cout << "L1=" << L1 << " \n";
  std::cout << "L2=" << L2 << " \n";
  std::cout << "L3=" << L3 << " \n";
  std::cout << "Lq=" << Lq << " \n";
  std::cout << "Ls=" << Ls << " \n\n";
  std::cout << "p0=" << p0 << " \n";
  std::cout << "dx=" << lc->dx << " \n";
  std::cout << "dy=" << lc->dy << " \n";
  std::cout << "dz=" << lc->dz << " \n";
  std::cout << "S_eq=" << S_eq << " \n\n";
        
  std::cout << "a=" << a << " \n";
  std::cout << "b=" << bb << " \n";
  std::cout << "c=" << cc << " \n";
  std::cout << "T=" << lc->T << " \n\n";
  if(lc->plane_confinement!=0) std::cout << "plane_confinement=" << lc->plane_confinement << " \n\n";

  std::cout << "Lambda=" << 1/Lambda << " \n";
  std::cout << "Lambda_s=" << 1/Lambda_s << " \n\n";

};

                
                

double landau_de_gennes::functional_derivative_00(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const
{
  return ldg_functional_derivative_00(bulk_params, QN, dQ, ddQ, v);
}

double landau_de_gennes::functional_derivative_01(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const
{
  return ldg_functional_derivative_01(bulk_params, QN, dQ, ddQ, v);
}

double landau_de_gennes::functional_derivative_02(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const
{
  return ldg_functional_derivative_02(bulk_params, QN, dQ, ddQ, v);
}

double landau_de_gennes::functional_derivative_11(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const
{
  return ldg_functional_derivative_11(bulk_params, QN, dQ, ddQ, v);
}

double landau_de_gennes::functional_derivative_12(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const
{
  return ldg_functional_derivative_12(bulk_params, QN, dQ, ddQ, v);
}
  
double landau_de_gennes::energy_calculation(const double  QN[5],const double  dQ[15], const double v[3]) const
{
  return ldg_energy_calculation(bulk_params, QN, dQ, v);
}
