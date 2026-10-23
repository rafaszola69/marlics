#include "driver.h"
#include "geometry.h"
#include "geometry_hemisphere.h"
#include "boundary.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> 
#include <ctime>
#include <cstring>

#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)


Geometry_Hemisphere::Geometry_Hemisphere(const struct Simulation_Parameters * sim_param) : GEOMETRY (sim_param),
                                                                                   HNx ((Nx-1.)/2),
                                                                                   HNy ((Ny-1.)/2.),
                                                                                   HNz (0.0),
                                                                                   dx(sim_param->dx),
                                                                                   dy(sim_param->dy),
                                                                                   dz(sim_param->dz)
{
  geometry_pointer=&(*this);
  geometry_name="Hemisphere";
  std::cout << "WARNING: the Hemisphere geometry is experimental and has known problems (see README).\n";
  number_of_boundaries=2;
  bc_conditions=std::vector<class BOUNDARY *>(number_of_boundaries);
  boundary_needed_to_be_defined="0 and 1";
  
  switch( sim_param->radius_flag[0] )
    {
     case parameter_status::unset:
         
       std::cout <<"Sphere internal radius not seted.\n"
                 <<"Using standard value.\n\n";

       R_in=( (Nx-1.)/2. - 0.4)*dx;
     break;

     case parameter_status::set:

       R_in=sim_param->R_in;
       break;
    }

  
    switch (sim_param->radius_flag[1])
     {
       case parameter_status::unset:
   
         std::cout <<"Sphere external radius not seted.\n"
                   <<"Using standard value.\n\n";
         R_ex=( (Nx-1)/2.+0.95)*dx;

         break;
   
       case parameter_status::set:
   
         R_ex=sim_param->R_ex;
         break;
     }

   std::cout <<"R_in=" << R_in <<std::endl
             <<"R_ex=" << R_ex <<std::endl;

   point_type=fill_point_type( );

   //Each curved-surface site represents an equal share of the cap area.
   int n_surface=0;
   for(int s=0; s<Nx*Ny*Nz; s++) if(point_type[s]==3) n_surface++;
   dA=2*M_PI*pow(0.5*(R_in+R_ex),2)/n_surface;
   std::cout << "Surface area per curved-surface site dA=" << dA << std::endl;
}

Geometry_Hemisphere::~Geometry_Hemisphere()
{
  
  
}

void Geometry_Hemisphere::surface_stencil(int i, int j, int k, int nm[3], int np[3]) const
{
  const int s=(k*Ny+j)*Nx+i;
  if(point_type[s]==2)
    {
      //Central differences in the (periodic) plane, one-sided upwards in z.
      nm[0]=(k*Ny+j)*Nx+(i+Nx-1)%Nx;  np[0]=(k*Ny+j)*Nx+(i+1)%Nx;
      nm[1]=(k*Ny+(j+Ny-1)%Ny)*Nx+i;  np[1]=(k*Ny+(j+1)%Ny)*Nx+i;
      nm[2]=s;                         np[2]=s+Nx*Ny;
      return;
    }

  //One-sided differences towards the centre of the droplet.
  const int c[3]={i,j,k}, step[3]={1,Nx,Nx*Ny};
  const double half[3]={HNx, HNy, HNz};
  for(int a=0; a<3; a++)
    {
      nm[a]= c[a]>=half[a] ? s-step[a] : s;
      np[a]= c[a]>=half[a] ? s         : s+step[a];
    }
}

int * Geometry_Hemisphere::fill_point_type( void )  const 
{
  int i,j,k;
  double delta_x,delta_y,delta_z;
  double rr;
  int * point_kind;

  
  if((point_kind= (int *)calloc(Nx*Ny*Nz, sizeof(int)))==NULL){ERROr}

  for( i= 0; i< Nx; i++)
  {
    for( j= 0; j< Ny; j++)
    {
      for( k= 0; k< Nz; k++)
      {       point_kind[(k*Ny+j)*Nx+i]=2;//Setting surface points
            
        delta_x=(i-HNx)*dx;
        delta_y=(j-HNy)*dy;
        delta_z=(k-HNz)*dz;
        rr=sqrt( delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);
        
    
        if(k ==0 && rr <=R_ex)
        {
          point_kind[(k*Ny+j)*Nx+i]=2;//Setting bottom points
        }
        else if(rr < R_in)
        {                                 
          point_kind[(k*Ny+j)*Nx+i]=1;//Setting bulk points
        }
        else if(rr >=R_in && rr <=R_ex)    
        {
          point_kind[(k*Ny+j)*Nx+i]=3;//Setting surface points
        }
        else
        {
          point_kind[(k*Ny+j)*Nx+i]=0;//Setting outside points
        }
      }  
    }         
  }

  return point_kind;
}

int Geometry_Hemisphere::backwardPTP(double *, const double * ) const {}



  

int Geometry_Hemisphere::surface_site(int i, int j, int k, double v[3], double & area) const
{
  //Flat bottom (point type 2) and curved cap (point type 3), as in fill_ki.
  if(point_type[(k*Ny+j)*Nx+i]==2)
    {
      v[0]=0.0;
      v[1]=0.0;
      v[2]=-1.0;
      area=dx*dy;
      return 0;
    }
  const double delta_x=(i-HNx)*dx;
  const double delta_y=(j-HNy)*dy;
  const double delta_z=(k-HNz)*dz;
  const double rr=sqrt(delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);

  v[0]=delta_x/rr;
  v[1]=delta_y/rr;
  v[2]=delta_z/rr;
  area=dA;
  return 1;
}
