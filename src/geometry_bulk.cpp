#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> 
#include <ctime>
#include <cstring>

#include "driver.h"
#include "geometry.h"
#include "geometry_bulk.h"
#include "energy_ldg.h"
#include <gsl/gsl_randist.h>
#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)

#define chunk_size 4*Nx*Ny


Geometry_Bulk::Geometry_Bulk(const struct Simulation_Parameters * sim_param) : GEOMETRY (sim_param)
{
  point_type=fill_point_type( );
  geometry_pointer=&(*this);

  geometry_name="Bulk";
  number_of_boundaries=0;
  bc_conditions=std::vector<class BOUNDARY *>(number_of_boundaries);

};




int * Geometry_Bulk::fill_point_type( void )  const 
{
  int i,j,k;
  int * point_kind;
  
  if((point_kind= (int *)calloc(Nx*Ny*Nz, sizeof(int)))==NULL){ERROr}

  for( i= 0; i< Nx; i++)
    {
      for( j= 0; j< Ny; j++)
        {
          for( k= 0; k< Nz; k++)
            {     
                
              point_kind[(k*Ny+j)*Nx+i]=1;
                                                              
            }
        }    
    }     
    
  return point_kind;
}


int Geometry_Bulk::backwardPTP(double *, const double * ) const {}




