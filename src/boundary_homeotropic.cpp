#include "driver.h"
#include "physics_homeotropic.h"
#include "geometry.h"
#include "boundary.h"
#include "boundary_homeotropic.h"
#include <math.h>
#include <iostream>
#include <vector>

Boundary_Homeotropic::Boundary_Homeotropic(const Simulation_Parameters * sim_param, int boundary_number) : BOUNDARY ( sim_param, boundary_number)
{                      
  
  condition_name="Homeotropic";

  //Asserting anchoring energy is set and get its value:
  try
    {
      Wo1= sim_param->Wo1.at(boundary_id);
      std::cout << "Wo1= " << Wo1 <<".\n";
      bc_params.kind=BC_HOMEOTROPIC;
      bc_params.Wo1=Wo1;
      
    }
  catch(std::out_of_range dummy_var )
    {

      assert_parameter_is_set(false, "Wo1");

    }

};

                
                

//Evaluate dQ00 at boundary
double Boundary_Homeotropic::functional_derivative_00(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const
{
  return homeotropic_functional_derivative_00(bulk_params, bc_params, QN, dQ, v);
}
//Evaluate dQ01 at boundary
double Boundary_Homeotropic::functional_derivative_01(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double  v[3]) const
{
  return homeotropic_functional_derivative_01(bulk_params, bc_params, QN, dQ, v);
}
//Evaluate dQ02 at boundary
double Boundary_Homeotropic::functional_derivative_02(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double  v[3]) const
{
  return homeotropic_functional_derivative_02(bulk_params, bc_params, QN, dQ, v);
}

//Evaluate dQ11 at boundary
double Boundary_Homeotropic::functional_derivative_11(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double  v[3]) const
{
  return homeotropic_functional_derivative_11(bulk_params, bc_params, QN, dQ, v);
}
//Evaluate dQ12 at boundary
double Boundary_Homeotropic::functional_derivative_12(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double  v[3]) const
{
  return homeotropic_functional_derivative_12(bulk_params, bc_params, QN, dQ, v);
}

double Boundary_Homeotropic::energy_calculation(const double  QN[5],const double  dQ[15],const double  v[3]) const
{
  return homeotropic_energy_calculation(bulk_params, bc_params, QN, dQ, v);
}
