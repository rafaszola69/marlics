#include "driver.h"
#include "geometry.h"
#include "geometry_hybrid_sphere.h"
#include "boundary.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> 
#include <ctime>
#include <cstring>

#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)


Geometry_Hybrid_Sphere::Geometry_Hybrid_Sphere(const struct Simulation_Parameters * sim_param) : GEOMETRY (sim_param),
                                                                                   HNx ((Nx-1.)/2),
                                                                                   HNy ((Ny-1.)/2.),
                                                                                   HNz ((Nz-1.)/2.),
                                                                                   dx(sim_param->dx),
                                                                                   dy(sim_param->dy),
                                                                                   dz(sim_param->dz)
{

  geometry_pointer=&(*this);
  geometry_name="Sphere";
  number_of_boundaries=2;
  bc_conditions=std::vector<class BOUNDARY *>(number_of_boundaries);
  boundary_needed_to_be_defined="2";

  
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

   //Each surface site represents an equal share of the droplet surface area.
   int n_surface=0;
   for(int s=0; s<Nx*Ny*Nz; s++) if(point_type[s]>=2) n_surface++;
   dA=4*M_PI*pow(0.5*(R_in+R_ex),2)/n_surface;
   std::cout << "Surface area per surface site dA=" << dA << std::endl;
};




Geometry_Hybrid_Sphere::~Geometry_Hybrid_Sphere()
{
  
  
}

 

void Geometry_Hybrid_Sphere::surface_stencil(int i, int j, int k, int nm[3], int np[3]) const
{
  const int s=(k*Ny+j)*Nx+i;
  //One-sided differences towards the centre of the droplet.
  const int c[3]={i,j,k}, step[3]={1,Nx,Nx*Ny};
  const double half[3]={Nx/2, Ny/2, Nz/2};
  for(int a=0; a<3; a++)
    {
      nm[a]= c[a]>=half[a] ? s-step[a] : s;
      np[a]= c[a]>=half[a] ? s         : s+step[a];
    }
}

int * Geometry_Hybrid_Sphere::fill_point_type( void )  const 
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
            {       
              delta_x=(i-HNx)*dx;
              delta_y=(j-HNy)*dy;
              delta_z=(k-HNz)*dz;
              rr=sqrt( delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);
              

              if(rr < R_in)
                {                                 

                  point_kind[(k*Ny+j)*Nx+i]=1;//Setting bulk points
                          
                }
              else if(rr >=R_in && rr <=R_ex)    
                {
                  if (k>HNz) point_kind[(k*Ny+j)*Nx+i]=2;//Setting top surface points
                  else point_kind[(k*Ny+j)*Nx+i]=3;       //Setting bottom surface points

                          
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

int Geometry_Hybrid_Sphere::backwardPTP(double *, const double * ) const {}

int Geometry_Hybrid_Sphere::surface_site(int i, int j, int k, double v[3], double & area) const
{
  const double delta_x=(i-HNx)*dx;
  const double delta_y=(j-HNy)*dy;
  const double delta_z=(k-HNz)*dz;
  const double rr=sqrt(delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);

  v[0]=delta_x/rr;
  v[1]=delta_y/rr;
  v[2]=delta_z/rr;
  area=dA;
  return point_type[(k*Ny+j)*Nx+i]==2 ? 0 : 1;
}
