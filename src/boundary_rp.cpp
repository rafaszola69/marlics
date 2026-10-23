#include "driver.h"
#include "physics_rp.h"
#include "geometry.h"
#include "boundary.h"
#include "boundary_rp.h"
#include <math.h>
#include <iostream>
#include <vector>

Boundary_Rp::Boundary_Rp(const Simulation_Parameters * sim_param, int boundary_number) : BOUNDARY ( sim_param, boundary_number)
{                      
  double phi_0;
  double theta_0;

  
  condition_name="Rapini-Papoular(rp)";

  //Asserting anchoring energy is set and get its value:
  try
    {
      Wo1= sim_param->Wo1.at(boundary_id);
      std::cout << "Wo1= " << Wo1 <<".\n";
      
    }
  catch(std::out_of_range dummy_var )
    {

      assert_parameter_is_set(false, "Wo1");

    }

  //Asserting easy axis angle phi_0 is set and get its value:
  try
    {
      phi_0= sim_param->phi_0.at(boundary_id);
      std::cout << "phi_0= " << phi_0 <<".\n";
      
    }
  catch(std::out_of_range dummy_var )
    {

      assert_parameter_is_set(false, "phi_0");

    }

    //Asserting easy axis angle theta_0 is set and get its value:
  try
    {
      theta_0= sim_param->theta_0.at(boundary_id);
      std::cout << "theta_0= " << theta_0 <<".\n\n";
      
    }
  catch(std::out_of_range dummy_var )
    {

      assert_parameter_is_set(false, "theta_0");

    }

  double n[3];
  double Pi=M_PI;

  
  n[0]=cos(phi_0*Pi/180)*sin(theta_0*Pi/180);
  n[1]=sin(phi_0*Pi/180)*sin(theta_0*Pi/180);
  n[2]=cos(theta_0*Pi/180);

  Q0_00= (0.5*S_eq*(3.0*n[0]*n[0]-1.0));
  Q0_01= (0.5*S_eq*(3.0*n[0]*n[1]));
  Q0_02= (0.5*S_eq*(3.0*n[0]*n[2]));
  Q0_11= (0.5*S_eq*(3.0*n[1]*n[1]-1.0));
  Q0_12= (0.5*S_eq*(3.0*n[1]*n[2]));

  bc_params.kind=BC_RP;
  bc_params.Wo1=Wo1;
  bc_params.Q0[0]=Q0_00; bc_params.Q0[1]=Q0_01; bc_params.Q0[2]=Q0_02; bc_params.Q0[3]=Q0_11; bc_params.Q0[4]=Q0_12;
  
};

                
                

//Evaluate dQ00 at boundary
double Boundary_Rp::functional_derivative_00(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double v[3]) const
{
  return rp_functional_derivative_00(bulk_params, bc_params, QN, dQ, v);
}
//Evaluate dQ01 at boundary
double Boundary_Rp::functional_derivative_01(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double  v[3]) const
{
  return rp_functional_derivative_01(bulk_params, bc_params, QN, dQ, v);
}
//Evaluate dQ02 at boundary
double Boundary_Rp::functional_derivative_02(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double  v[3]) const
{
  return rp_functional_derivative_02(bulk_params, bc_params, QN, dQ, v);
}
//Evaluate dQ11 at boundary
double Boundary_Rp::functional_derivative_11(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double  v[3]) const
{
  return rp_functional_derivative_11(bulk_params, bc_params, QN, dQ, v);
}
//Evaluate dQ12 at boundary
double Boundary_Rp::functional_derivative_12(const double  QN[5],const double  dQ[15],const double  ddQ[30], const double  v[3]) const
{
  return rp_functional_derivative_12(bulk_params, bc_params, QN, dQ, v);
}

double Boundary_Rp::energy_calculation(const double  QN[5],const double  dQ[15],const double  v[3]) const
{
  return rp_energy_calculation(bulk_params, bc_params, QN, dQ, v);
}
