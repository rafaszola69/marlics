#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> 
#include <ctime>
#include <cstring>

#include "driver.h"
#include "geometry.h"
#include "geometry_wedge2.h"
#include "boundary.h"
#include "energy_ldg.h"
#include <gsl/gsl_randist.h>
#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)


Geometry_Wedge2::Geometry_Wedge2(const struct Simulation_Parameters * sim_param) : GEOMETRY (sim_param), Nz2 (sim_param->Nz2), wTheta (sim_param->wTheta)
{
  point_type=fill_point_type( );
  geometry_pointer=&(*this);

  geometry_name="Wedge2";
  std::cout << "WARNING: the Wedge2 geometry is experimental and has known problems (see README).\n";
  number_of_boundaries=2;
  bc_conditions=std::vector<class BOUNDARY *>(number_of_boundaries);
  boundary_needed_to_be_defined="0 and 1";
  
  
};


void  Geometry_Wedge2::fill_ki(double * k_i,
                    const double * Qij)  const 
{
  const double theta=atan2(Nz-Nz2,Ny);
	  
#pragma omp for schedule(dynamic,1) collapse(2) 
  for( int k= 0; k< Nz; k++)
    {
      for( int j= 0; j< Ny; j++)
        {
          for( int i= 0; i< Nx; i++)
            {
							if (point_type[(Nx*(Ny*k+j)+i)]==1)
								{
									double dQ[15];
									double ddQ[30];
									double QN[5];
									double v[3];
					
									//check_bulk_limits( i,  j,  k);
									
									int ip1= (i+1)%Nx;
									int jp1= (j+1)%Ny;
									int kp1= (k+1);
									int im1= (i+Nx-1)%Nx;
									int jm1= (j+Ny-1)%Ny;
									int km1= (k-1);
					
									for(int ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
					
									//Calcule first derivatives of Qij:
									for(int ll=0; ll<=4;ll++) dQ[ll]= 0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
									if(j==0)     for(int ll=0; ll<=4;ll++) dQ[5+ll]= (Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+j)*Nx+i)+ll])*dy_1;
					else		if(j==Ny-1)  for(int ll=0; ll<=4;ll++) dQ[5+ll]= (Qij[5*((k*Ny+j)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
									else         for(int ll=0; ll<=4;ll++) dQ[5+ll]= 0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
								
									for(int ll=0; ll<=4;ll++) dQ[10+ll]= 0.5*(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1;
					
									//Calculate second derivatives of Qij:
									for(int ll=0; ll<=4;ll++) ddQ[ll]= (Qij[5*((k*Ny+j)*Nx+ip1)+ll]-2.0*QN[ll]+Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1*dx_1;
									for(int ll=0; ll<=4;ll++) ddQ[10+ll]= 0.25*(Qij[5*((kp1*Ny+j)*Nx+ip1)+ll]+Qij[5*((km1*Ny+j)*Nx+im1)+ll]-Qij[5*((kp1*Ny+j)*Nx+im1)+ll]-Qij[5*((km1*Ny+j)*Nx+ip1)+ll])*dx_1*dz_1;
									for(int ll=0; ll<=4;ll++) ddQ[25+ll]= (Qij[5*((kp1*Ny+j)*Nx+i)+ll]-2.0*QN[ll]+Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1*dz_1;
									if(j==0){	
										for(int ll=0; ll<=4;ll++) ddQ[ll+5]= 0.5*(Qij[5*((k*Ny+jp1)*Nx+ip1)+ll]+Qij[5*((k*Ny+j)*Nx+im1)+ll]-Qij[5*((k*Ny+jp1)*Nx+im1)+ll]-Qij[5*((k*Ny+j)*Nx+ip1)+ll])*dx_1*dy_1;
										for(int ll=0; ll<=4;ll++) ddQ[15+ll]= (2*QN[ll]-5*Qij[5*((k*Ny+j+1)*Nx+i)+ll]+4*Qij[5*((k*Ny+j+2)*Nx+i)+ll]-Qij[5*((k*Ny+j+3)*Nx+i)+ll])*dy_1*dy_1*dy_1;
										for(int ll=0; ll<=4;ll++) ddQ[20+ll]= 0.5*(Qij[5*((kp1*Ny+jp1)*Nx+i)+ll]+Qij[5*((km1*Ny+j)*Nx+i)+ll]-Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+jp1)*Nx+i)+ll])*dy_1*dz_1;
									}else if(j==Ny-1){	
										for(int ll=0; ll<=4;ll++) ddQ[ll+5]= 0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]+Qij[5*((k*Ny+jm1)*Nx+im1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll]-Qij[5*((k*Ny+jm1)*Nx+ip1)+ll])*dx_1*dy_1;
										for(int ll=0; ll<=4;ll++) ddQ[15+ll]= (2*QN[ll]-5*Qij[5*((k*Ny+j-1)*Nx+i)+ll]+4*Qij[5*((k*Ny+j-2)*Nx+i)+ll]-Qij[5*((k*Ny+j-3)*Nx+i)+ll])*dy_1*dy_1*dy_1;
										for(int ll=0; ll<=4;ll++) ddQ[20+ll]= 0.5*(Qij[5*((kp1*Ny+j)*Nx+i)+ll]+Qij[5*((km1*Ny+jm1)*Nx+i)+ll]-Qij[5*((kp1*Ny+jm1)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dy_1*dz_1;
									}else {	
										for(int ll=0; ll<=4;ll++) ddQ[ll+5]= 0.25*(Qij[5*((k*Ny+jp1)*Nx+ip1)+ll]+Qij[5*((k*Ny+jm1)*Nx+im1)+ll]-Qij[5*((k*Ny+jp1)*Nx+im1)+ll]-Qij[5*((k*Ny+jm1)*Nx+ip1)+ll])*dx_1*dy_1;
										for(int ll=0; ll<=4;ll++) ddQ[15+ll]= (Qij[5*((k*Ny+jp1)*Nx+i)+ll]-2.0*QN[ll]+Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1*dy_1;
										for(int ll=0; ll<=4;ll++) ddQ[20+ll]= 0.25*(Qij[5*((kp1*Ny+jp1)*Nx+i)+ll]+Qij[5*((km1*Ny+jm1)*Nx+i)+ll]-Qij[5*((kp1*Ny+jm1)*Nx+i)+ll]-Qij[5*((km1*Ny+jp1)*Nx+i)+ll])*dy_1*dz_1;
									}
								  //~ if ((i==0)*(j==Ny-1)*(k==Nz/2))	for(int ll=0; ll<=4;ll++) printf("%d %g %g\n",(-2*Qij[5*((k*Ny+jm1)*Nx+i)+ll]+QN[ll]+Qij[5*((k*Ny+j-2)*Nx+i)+ll])*dy_1*dy_1,(-QN[ll]+Qij[5*((k*Ny+j-2)*Nx+i)+ll])*dy_1*dy_1);
								
									k_i[5*(Nx*(Ny*k+j)+i)+0]=bulk_energy.functional_derivative_00(QN,dQ,ddQ,v); 
									k_i[5*(Nx*(Ny*k+j)+i)+1]=bulk_energy.functional_derivative_01(QN,dQ,ddQ,v); 
									k_i[5*(Nx*(Ny*k+j)+i)+2]=bulk_energy.functional_derivative_02(QN,dQ,ddQ,v); 
									k_i[5*(Nx*(Ny*k+j)+i)+3]=bulk_energy.functional_derivative_11(QN,dQ,ddQ,v); 
									k_i[5*(Nx*(Ny*k+j)+i)+4]=bulk_energy.functional_derivative_12(QN,dQ,ddQ,v);
								}else if(point_type[(Nx*(Ny*k+j)+i)]==2)
								{
									double dQ[15];
									double ddQ[30];
									double QN[5];
									double v[3];
							
									int ip1= (i+1)%Nx;
									int jp1= (j+1)%Ny;
    				      int kp1= k+1;
									int im1= (i+Nx-1)%Nx;
									int jm1= (j+Ny-1)%Ny;
									//~ int km1= (k-1);
									v[0]=0.0;
									v[1]=0.0;
									v[2]=-1.0;
										
									for(int ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
							
									//Calcule first derivatives of Qij:
									for(int ll=0; ll<=4;ll++) dQ[ll]=0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
										
									for(int ll=0; ll<=4;ll++) dQ[5+ll]=0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
										
									for(int ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((kp1*Ny+j)*Nx+i)+ll]-Qij[5*((k*Ny+j)*Nx+i)+ll])*dz_1;
										
									k_i[5*(Nx*(Ny*k+j)+i)+0]= bc_conditions[0]->functional_derivative_00(QN,dQ,ddQ,v);
									k_i[5*(Nx*(Ny*k+j)+i)+1]= bc_conditions[0]->functional_derivative_01(QN,dQ,ddQ,v);
									k_i[5*(Nx*(Ny*k+j)+i)+2]= bc_conditions[0]->functional_derivative_02(QN,dQ,ddQ,v);
									k_i[5*(Nx*(Ny*k+j)+i)+3]= bc_conditions[0]->functional_derivative_11(QN,dQ,ddQ,v);
									k_i[5*(Nx*(Ny*k+j)+i)+4]= bc_conditions[0]->functional_derivative_12(QN,dQ,ddQ,v);
								}else if(point_type[(Nx*(Ny*k+j)+i)]==3)
								{
									double dQ[15];
									double ddQ[30];
									double QN[5];
									double v[3];
							
									int ip1= (i+1)%Nx;
									int jp1= (j+1)%Ny;
//~     				      int kp1= k;
									int im1= (i+Nx-1)%Nx;
									int jm1= (j+Ny-1)%Ny;
									int km1= (k-1);

									v[0]=0.0;
									v[1]=-sin(theta);
									v[2]=cos(theta);
										
									for(int ll=0; ll<=4;ll++) QN[ll]=Qij[5*(Nx*(Ny*k+j)+i)+ll];
							
									//Calcule first derivatives of Qij:
									for(int ll=0; ll<=4;ll++) dQ[ll]=0.5*(Qij[5*((k*Ny+j)*Nx+ip1)+ll]-Qij[5*((k*Ny+j)*Nx+im1)+ll])*dx_1;
										
									for(int ll=0; ll<=4;ll++) dQ[5+ll]=0.5*(Qij[5*((k*Ny+jp1)*Nx+i)+ll]-Qij[5*((k*Ny+jm1)*Nx+i)+ll])*dy_1;
										
									for(int ll=0; ll<=4;ll++) dQ[10+ll]=(Qij[5*((k*Ny+j)*Nx+i)+ll]-Qij[5*((km1*Ny+j)*Nx+i)+ll])*dz_1;
										
									k_i[5*(Nx*(Ny*k+j)+i)+0]= bc_conditions[1]->functional_derivative_00(QN,dQ,ddQ,v);
									k_i[5*(Nx*(Ny*k+j)+i)+1]= bc_conditions[1]->functional_derivative_01(QN,dQ,ddQ,v);
									k_i[5*(Nx*(Ny*k+j)+i)+2]= bc_conditions[1]->functional_derivative_02(QN,dQ,ddQ,v);
									k_i[5*(Nx*(Ny*k+j)+i)+3]= bc_conditions[1]->functional_derivative_11(QN,dQ,ddQ,v);
									k_i[5*(Nx*(Ny*k+j)+i)+4]= bc_conditions[1]->functional_derivative_12(QN,dQ,ddQ,v);
								}else
								{	
									k_i[5*(Nx*(Ny*k+j)+i)+0]= 0;
									k_i[5*(Nx*(Ny*k+j)+i)+1]= 0;
									k_i[5*(Nx*(Ny*k+j)+i)+2]= 0;
									k_i[5*(Nx*(Ny*k+j)+i)+3]= 0;
									k_i[5*(Nx*(Ny*k+j)+i)+4]= 0;
								}
								
            }
        }
    }
	
}

    
int * Geometry_Wedge2::fill_point_type( void )  const 
{
  int i,j,k;
  int * point_kind;
  double height;
  double midZ=Nz/2.0;
  
  if((point_kind= (int *)calloc(Nx*Ny*Nz, sizeof(int)))==NULL){ERROr}

  for( i= 0; i< Nx; i++)
    {
      for( j= 0; j< Ny; j++)
        {
					height=tan(M_PI*wTheta/360)*j;
          for( k= 0;  k< midZ; k++) point_kind[(k*Ny+j)*Nx+i]=2;
          for( k= midZ; k< Nz; k++) point_kind[(k*Ny+j)*Nx+i]=3;
          for( k= midZ-height+1; k< midZ+height-1; k++) if (k>-1 && k<Nz) point_kind[(k*Ny+j)*Nx+i]=1;
          for( k= 0; k< midZ-height-2 ; k++) if (k>-1 && k<Nz) point_kind[(k*Ny+j)*Nx+i]=0;
          for( k= midZ+height+2; k< Nz; k++) if (k>-1 && k<Nz) point_kind[(k*Ny+j)*Nx+i]=0;
        }
    }     
    
  return point_kind;
}


int Geometry_Wedge2::backwardPTP(double *, const double * ) const {}

int Geometry_Wedge2::surface_site(int i, int j, int k, double v[3], double & area) const
{
  //Same normals as fill_ki: bottom (0,0,-1); top (0,-sin(theta),cos(theta)), with area dx*dy/cos(theta).
  const double theta=atan2(Nz-Nz2,Ny);

  v[0]=0.0;
  if(point_type[(k*Ny+j)*Nx+i]==2)
    {
      v[1]=0.0;
      v[2]=-1.0;
      area=dA;
      return 0;
    }
  v[1]=-sin(theta);
  v[2]=cos(theta);
  area=dA/cos(theta);
  return 1;
}
