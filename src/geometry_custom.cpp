#include "driver.h"
#include "geometry.h"
#include "geometry_custom.h"
#include "boundary.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> 
#include <ctime>
#include <cstring>
#include <algorithm>   
#pragma warning(disable : 4996)
#define MAX(a,b) ((a) > (b) ? a : b)
#define min(a,b) ((a) < (b) ? a : b)


Geometry_Custom::Geometry_Custom(const struct Simulation_Parameters * sim_param) : GEOMETRY (sim_param),                                                                                                                HNx ((Nx-1.)/2),
                                                                                   HNy ((Ny-1.)/2.),
                                                                                   HNz ((Nz-1.)/2.),
                                                                                   dx(sim_param->dx),
                                                                                   dy(sim_param->dy),
                                                                                   dz(sim_param->dz)
{

  geometry_pointer=&(*this);
  geometry_name="Custom";
  v=(double**)calloc(sim_param->Nx*sim_param->Ny*sim_param->Nz,sizeof(double*));
  for (int i=0;i< (Nx*Ny*Nz); i++)v[i]=(double*)calloc(3,sizeof(double));
  point_type=fill_pt_and_normals(v, sim_param->bound_fname,number_of_boundaries);
  test_derivatives(sim_param->integrator_type);
  
  boundary_needed_to_be_defined=std::to_string(number_of_boundaries);
  bc_conditions=std::vector<class BOUNDARY *>(number_of_boundaries);
};

Geometry_Custom::~Geometry_Custom()
{
  
}

void Geometry_Custom::surface_stencil(int i, int j, int k, int nm[3], int np[3]) const
{
  const int s=(k*Ny+j)*Nx+i;
  //One-sided differences opposite to the normal (periodic box); no derivative along zero components.
  const int c[3]={i,j,k}, n[3]={Nx,Ny,Nz}, step[3]={1,Nx,Nx*Ny};
  for(int a=0; a<3; a++)
    {
      np[a]= v[s][a]>=0 ? s : s+(((c[a]+1)%n[a])-c[a])*step[a];
      nm[a]= v[s][a]>0  ? s+(((c[a]-1+n[a])%n[a])-c[a])*step[a] : s;
    }
}

void  Geometry_Custom::test_derivatives( const char integrator_type[])  const 
{
  int counter=0;

//  #pragma omp master
  for( int k= 0; k< Nz; k++)
  {
    for( int j= 0; j< Ny; j++)
    {
      for( int i= 0; i< Nx; i++)
      {
        
        int ip,jp,kp, im, jm, km;
        
        if(point_type[(k*Ny+j)*Nx+i] == 1)
        {

          //check_bulk_limits( i,  j,  k);  
          ip= (i+1)%Nx;
          jp= (j+1)%Ny;
          kp= (k+1)%Nz;
          im= (i-1+Nx)%Nx;
          jm= (j-1+Ny)%Ny;
          km= (k-1+Nz)%Nz;
          
          if(point_type[(k*Ny+j)*Nx+ip]*point_type[(k*Ny+j)*Nx+im]*
             point_type[(k*Ny+jp)*Nx+i]*point_type[(k*Ny+jm)*Nx+i]*
             point_type[(kp*Ny+j)*Nx+i]*point_type[(km*Ny+j)*Nx+i]*
             point_type[(k*Ny+jp)*Nx+ip]*point_type[(k*Ny+jp)*Nx+im]*
             point_type[(k*Ny+jm)*Nx+ip]*point_type[(k*Ny+jm)*Nx+im]*
             point_type[(kp*Ny+j)*Nx+ip]*point_type[(kp*Ny+j)*Nx+im]*
             point_type[(km*Ny+j)*Nx+ip]*point_type[(km*Ny+j)*Nx+im]*
             point_type[(kp*Ny+jp)*Nx+i]*point_type[(kp*Ny+jm)*Nx+i]*
             point_type[(km*Ny+jp)*Nx+i]*point_type[(km*Ny+jm)*Nx+i]==0)
             {
               printf("Warning!!! The bulk cell (%d,%d,%d) is trying to interact with a empty cell (type 0)!\n",i,j,k);
               counter++;
             if(!point_type[(k*Ny+j)*Nx+im] )printf("%d %d %d \n",im,j,k);
             if(!point_type[(k*Ny+jm)*Nx+i] )printf("%d %d %d \n",i,jm,k);
             if(!point_type[(km*Ny+j)*Nx+i] )printf("%d %d %d \n",i,j,km);
             if(!point_type[(k*Ny+jp)*Nx+im])printf("%d %d %d \n",im,jp,k);
             if(!point_type[(k*Ny+jm)*Nx+im])printf("%d %d %d \n",im,jm,k);
             if(!point_type[(kp*Ny+j)*Nx+im])printf("%d %d %d \n",im,j,kp);
             if(!point_type[(km*Ny+j)*Nx+im])printf("%d %d %d \n",im,j,km);
             if(!point_type[(kp*Ny+jm)*Nx+i])printf("%d %d %d \n",i,jm,kp);
             if(!point_type[(km*Ny+jm)*Nx+i])printf("%d %d %d \n",i,jm,km);
             if(!point_type[(k*Ny+j)*Nx+ip] )printf("%d %d %d \n",ip,j,k);
             if(!point_type[(k*Ny+jp)*Nx+i] )printf("%d %d %d \n",i,jp,k);
             if(!point_type[(kp*Ny+j)*Nx+i] )printf("%d %d %d \n",i,j,kp);
             if(!point_type[(k*Ny+jp)*Nx+ip])printf("%d %d %d \n",ip,jp,k);
             if(!point_type[(k*Ny+jm)*Nx+ip])printf("%d %d %d \n",ip,jm,k);
             if(!point_type[(kp*Ny+j)*Nx+ip])printf("%d %d %d \n",ip,j,kp);
             if(!point_type[(km*Ny+j)*Nx+ip])printf("%d %d %d \n",ip,j,km);
             if(!point_type[(kp*Ny+jp)*Nx+i])printf("%d %d %d \n",i,jp,kp);
             if(!point_type[(km*Ny+jp)*Nx+i])printf("%d %d %d \n",i,jp,km);
             }
             
        }          
        else if (point_type[(k*Ny+j)*Nx+i] == 0)
        {
          
        }
        else 
        {
          //Boundary condtions:
          
          ip= (v[(k*Ny+j)*Nx+i][0] >=0) ? i : (i+1)%Nx ;
          jp= (v[(k*Ny+j)*Nx+i][1] >=0) ? j : (j+1)%Ny ;
          kp= (v[(k*Ny+j)*Nx+i][2] >=0) ? k : (k+1)%Nz ;
          
          im= (v[(k*Ny+j)*Nx+i][0] > 0) ? (i-1+Nx)%Nx : i ;
          jm= (v[(k*Ny+j)*Nx+i][1] > 0) ? (j-1+Ny)%Ny : j ;
          km= (v[(k*Ny+j)*Nx+i][2] > 0) ? (k-1+Nz)%Nz : k ;
          if(strcasecmp(integrator_type,"fire"))
          {
            if(point_type[(k*Ny+j)*Nx+ip]*point_type[(k*Ny+j)*Nx+im]*
               point_type[(k*Ny+jp)*Nx+i]*point_type[(k*Ny+jm)*Nx+i]*
               point_type[(kp*Ny+j)*Nx+i]*point_type[(km*Ny+j)*Nx+i]*
               point_type[(k*Ny+jp)*Nx+ip]*point_type[(k*Ny+jp)*Nx+im]*
               point_type[(k*Ny+jm)*Nx+ip]*point_type[(k*Ny+jm)*Nx+im]*
               point_type[(kp*Ny+j)*Nx+ip]*point_type[(kp*Ny+j)*Nx+im]*
               point_type[(km*Ny+j)*Nx+ip]*point_type[(km*Ny+j)*Nx+im]*
               point_type[(kp*Ny+jp)*Nx+i]*point_type[(kp*Ny+jm)*Nx+i]*
               point_type[(km*Ny+jp)*Nx+i]*point_type[(km*Ny+jm)*Nx+i]==0)
              {
                printf("Warning!!! The surface cell (%d,%d,%d) is trying to interact with a empty cell (type 0)!\n",i,j,k);
                counter++;
              }
          }else{
            if(point_type[(k*Ny+j)*Nx+ip]*point_type[(k*Ny+j)*Nx+im]*
               point_type[(k*Ny+jp)*Nx+i]*point_type[(k*Ny+jm)*Nx+i]*
               point_type[(kp*Ny+j)*Nx+i]*point_type[(km*Ny+j)*Nx+i]==0)
              {
                printf("Warning!!! The surface cell (%d,%d,%d) is trying to interact with a empty cell (type 0)!\n",i,j,k);
                counter++;
              }
          }
        }
	    }
    }
  }
  #pragma omp master
  {
    if(counter!=0)
    {
    printf("%d cells are trying to interact with an empty cell!!\n Please, check your boundary file.\n",counter);
    exit(10);
    }
  }
}

int * Geometry_Custom::fill_pt_and_normals( double **v , const char bound_file_name[], int &number_of_boundaries)  const 
{
  int i, nn;
  int * point_kind, max_point_kind=0, line_number, read_tester=1;
  int NoV,pos_nx, pos_ny, pos_nz, pos_pt;
  FILE *bound_input = fopen(bound_file_name,"r");
  if(bound_input==0){perror(bound_file_name);exit(2);}
  printf("Reading boundaries from %s\n",bound_file_name);
  char *line=(char*)malloc(500);
  std::string testline;
  line=fgets(line,500,bound_input);
  testline=line;
  NoV=1+std::count(testline.begin(),testline.end(),',');
  double *var= new double[NoV];
  char strvar[500];
  rewind(bound_input);
  for (i=0; i<NoV;i++) 
  {
    char names[500];
    int dummy=fscanf(bound_input,"%[^,\n ],",names);
    if(strcasecmp(names,"nx")==0)pos_nx=i;
    if(strcasecmp(names,"ny")==0)pos_ny=i;
    if(strcasecmp(names,"nz")==0)pos_nz=i;
    if(strcasecmp(names,"pt")==0)pos_pt=i;
  }
  printf("Using the positions of nx %d ny %d nz %d and pt %d.\n ",pos_nx,pos_ny,pos_nz,pos_pt);
  if((point_kind= (int *)calloc(Nx*Ny*Nz, sizeof(int)))==NULL){ERROr}
  line_number=0;
  line=fgets(line,500,bound_input);
  for (i=0; i<Nx*Ny*Nz; i++){
    point_kind[i]=1;
  }  
  while (read_tester!=EOF){
    for (nn=0;nn<NoV;nn++) 
    {
      read_tester=fscanf(bound_input,"%[^,\n ],",strvar);
      if(read_tester==1)var[nn]=atof(strvar);
      else continue;
    }
    int ii=(int)var[0];
    int jj=(int)var[1];
    int kk=(int)var[2];
    line_number++;
    if (ii>=Nx|| jj>=Ny||kk>=Nz){fprintf(stderr,"Point out of the simulation box at line %d\n",line_number);exit(10);}
    if ((int)var[pos_pt]<0){fprintf(stderr,"Negative Point Type founde at line %d\n\
Please, use 0 as empty space, 1 as bulk and 2+ as surface point type.\n",line_number);exit(10);}
    line=fgets(line,500,bound_input);
             v[(Ny*kk+jj)*Nx+ii][0]=var[pos_nx];
             v[(Ny*kk+jj)*Nx+ii][1]=var[pos_ny];
             v[(Ny*kk+jj)*Nx+ii][2]=var[pos_nz];
    point_kind[(Ny*kk+jj)*Nx+ii]=(int)var[pos_pt];
    max_point_kind=MAX(max_point_kind,point_kind[(Ny*kk+jj)*Nx+ii]);
  }
  number_of_boundaries=max_point_kind-1;
  printf("%d boundary(boundaries) found in %s!!\n\n",number_of_boundaries, bound_file_name);
  fflush(stdout);
  return point_kind;
}

int * Geometry_Custom::fill_point_type( void )  const  { return 0;}

int Geometry_Custom::backwardPTP(double *, const double * ) const {}

int Geometry_Custom::surface_site(int i, int j, int k, double v[3], double & area) const
{
  const int s=(k*Ny+j)*Nx+i;

  v[0]=this->v[s][0];
  v[1]=this->v[s][1];
  v[2]=this->v[s][2];
  area=dA;
  return point_type[s]-2;
}
