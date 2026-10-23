#include "driver.h"
#include "omp.h"
#include <cstdio>           
#include <cstdlib>          
#include <cstring>
#include <cmath>            
#include "geometry.h"
#include "integrator.h"
#include "ptp_minimization.h"
#pragma warning(disable : 3002)
#define cl 2
#define cll 2
#define fixed_chunk_size 1
#define chunk_size 1
#define new_chunk_size 4*Nx*Ny

#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)

//Dormand-Prince integrator
PTP_min::PTP_min( GEOMETRY  * lc_pointer, const struct Simulation_Parameters *sim_param ) : Integrator( lc_pointer)
{
  dt=sim_param->dt;
  //allocate the ith-stage array:
  if((force= (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  if((vQij = (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  if((Qijt = (double *)calloc(5*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  if((energy = (double *)calloc(2*Nx*Ny*Nz, sizeof(double)))==NULL){ERROr}
  
  dtmin = dt*0.01;
  dtMax = dt * 100;
  alphainit = 0.3;
  scaleAlpha = 0.99;
  alphaMin = 0.1;
  scaleP =  1.1 ;
  scaleM =  0.9 ;
  Nmin = 10;
  min_force=sim_param->min_force;
  dt_init = dt;
  alpha = alphainit;
  printf("Initializing simulation with FIRE(Fast Inertial Relaxation Engine) method.\n");
  std::cout << "dt=" << dt << " \n";
};

bool PTP_min::evolve( double * Qij, double *time, double tf )
{
  bool continue_condition=false;
  int ll;
  printf("first run\n");
  while (sample_geometry->backwardPTP(vQij,Qij)>2){
  #pragma omp single 
    printf("HERE %g,%g\n",Qij[5*(134+Nx*(56+Ny*68))],vQij[5*(134+Nx*(56+Ny*68))]);
  #pragma omp for simd schedule(simd:dynamic,new_chunk_size)          
    for( ll=0; ll<5*Nx*Ny*Nz;ll++) Qij[ll]=Qij[ll]+vQij[ll];
  #pragma omp single 
    printf("HERE %g,%g\n",Qij[5*(134+Nx*(56+Ny*68))],vQij[5*(134+Nx*(56+Ny*68))]);
  
  }

  
  return continue_condition;
};
