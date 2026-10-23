#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> 
#include <ctime>
#include <cstring>

#include "driver.h"
#include "geometry.h"
#include "geometry_slab.h"
#include "boundary.h"
#include "energy_ldg.h"
#include <gsl/gsl_randist.h>
#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)


slab::slab(const struct Simulation_Parameters * sim_param) : GEOMETRY (sim_param)
{
  point_type=fill_point_type( );
  geometry_pointer=&(*this);

  geometry_name="Slab";
  number_of_boundaries=2;
  bc_conditions=std::vector<class BOUNDARY *>(number_of_boundaries);
  boundary_needed_to_be_defined="0 and 1";
  
  
};


void slab::surface_stencil(int i, int j, int k, int nm[3], int np[3]) const
{
  const int s=(k*Ny+j)*Nx+i;
  if(point_type[s]==2)
    {
      //Central differences in the (periodic) plane, one-sided upwards in z.
      nm[0]=(k*Ny+j)*Nx+(i+Nx-1)%Nx;  np[0]=(k*Ny+j)*Nx+(i+1)%Nx;
      nm[1]=(k*Ny+(j+Ny-1)%Ny)*Nx+i;  np[1]=(k*Ny+(j+1)%Ny)*Nx+i;
      nm[2]=s;                         np[2]=s+Nx*Ny;
    }
  else
    {
      //Central differences in the (periodic) plane, one-sided downwards in z.
      nm[0]=(k*Ny+j)*Nx+(i+Nx-1)%Nx;  np[0]=(k*Ny+j)*Nx+(i+1)%Nx;
      nm[1]=(k*Ny+(j+Ny-1)%Ny)*Nx+i;  np[1]=(k*Ny+(j+1)%Ny)*Nx+i;
      nm[2]=s-Nx*Ny;                   np[2]=s;
    }
}

    
int * slab::fill_point_type( void )  const 
{
  int i,j,k;
  int * point_kind;

  
  if((point_kind= (int *)calloc(Nx*Ny*Nz, sizeof(int)))==NULL){ERROr}

  for( i= 0; i< Nx; i++)
    {
      for( j= 0; j< Ny; j++)
        {

          k=0;
          point_kind[(k*Ny+j)*Nx+i]=2;

          for( k= 1; k< Nz-1; k++)
            {       
              
              point_kind[(k*Ny+j)*Nx+i]=1;
                                                              
            }
          
          k=Nz-1;
          point_kind[(k*Ny+j)*Nx+i]=3;
        }
    }     
    
  return point_kind;
}


int slab::backwardPTP(double *, const double * ) const {}

int slab::surface_site(int i, int j, int k, double v[3], double & area) const
{
  v[0]=0.0;
  v[1]=0.0;
  v[2]= point_type[(k*Ny+j)*Nx+i]==2 ? -1.0 : 1.0;
  area=dA;
  return point_type[(k*Ny+j)*Nx+i]==2 ? 0 : 1;
}
