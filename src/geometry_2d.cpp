#include "driver.h"
#include "geometry.h"
#include "geometry_2d.h"
#include "boundary.h"
#include <cmath>
#include <cstdlib>

//Stencil of a surface site along x and y: one-sided away from the domain edges of non-periodic
//axes, central otherwise (wrapping around periodic axes). Along z it is always periodic.
static void edge_stencil(const int N[3], const bool periodic[3], const bool low[2], const bool high[2],
                         int i, int j, int k, int nm[3], int np[3])
{
  const int c[3]={i,j,k}, step[3]={1,N[0],N[0]*N[1]};
  const int s=(k*N[1]+j)*N[0]+i;

  for(int a=0; a<3; a++)
    {
      const bool lo= a<2 && low[a], hi= a<2 && high[a];
      if(lo)      { nm[a]=s;          np[a]=s+step[a]; }
      else if(hi) { nm[a]=s-step[a];  np[a]=s; }
      else
        {
          const int cm= periodic[a] ? (c[a]-1+N[a])%N[a] : c[a]-1;
          const int cp= periodic[a] ? (c[a]+1)%N[a]      : c[a]+1;
          nm[a]=s+(cm-c[a])*step[a];
          np[a]=s+(cp-c[a])*step[a];
        }
    }
}

//################################ Disk ################################

Geometry_Disk::Geometry_Disk(const struct Simulation_Parameters * sim_param) : GEOMETRY (sim_param),
  HNx((Nx-1.)/2.), HNy((Ny-1.)/2.)
{
  geometry_pointer=&(*this);
  geometry_name="Disk";
  number_of_boundaries=1;
  bc_conditions=std::vector<class BOUNDARY *>(number_of_boundaries);
  boundary_needed_to_be_defined="0";

  if(Nx!=Ny || sim_param->dx!=sim_param->dy)
    {
      std::cout << "The disk geometry needs Nx=Ny and dx=dy.\nAborting the program.\n";
      exit(1);
    }

  //Surface shell from R_in to R_ex: thick enough that bulk sites only touch bulk or surface sites,
  //including the diagonal neighbours of the mixed derivatives.
  R_in= sim_param->radius_flag[0]==parameter_status::set ? sim_param->R_in : (HNx-0.5)*sim_param->dx;
  R_ex= sim_param->radius_flag[1]==parameter_status::set ? sim_param->R_ex : (HNx+0.95)*sim_param->dx;
  std::cout << "R_in=" << R_in << std::endl << "R_ex=" << R_ex << std::endl;

  point_type=fill_point_type( );

  //Each surface site represents an equal share of the side of the cylinder (perimeter per
  //thickness dz in two dimensions).
  int n_surface=0;
  for(int s=0; s<Nx*Ny*Nz; s++) if(point_type[s]>=2) n_surface++;
  dA=2*M_PI*0.5*(R_in+R_ex)*Nz*sim_param->dz/n_surface;
  std::cout << "Surface area per surface site dA=" << dA << std::endl;
}

int * Geometry_Disk::fill_point_type( void ) const
{
  int * point_kind;
  if((point_kind= (int *)calloc(Nx*Ny*Nz, sizeof(int)))==NULL){ERROr}

  const double d=1/dx_1;
  for(int k=0; k<Nz; k++)
    for(int j=0; j<Ny; j++)
      for(int i=0; i<Nx; i++)
        {
          const double rr=d*sqrt((i-HNx)*(i-HNx)+(j-HNy)*(j-HNy));
          point_kind[(k*Ny+j)*Nx+i]= rr<R_in ? 1 : (rr<=R_ex ? 2 : 0);
        }
  return point_kind;
}

int Geometry_Disk::surface_site(int i, int j, int k, double v[3], double & area) const
{
  const double x=i-HNx, y=j-HNy, rr=sqrt(x*x+y*y);
  v[0]=x/rr;
  v[1]=y/rr;
  v[2]=0;
  area=dA;
  return 0;
}

//One-sided differences towards the axis of the disk.
void Geometry_Disk::surface_stencil(int i, int j, int k, int nm[3], int np[3]) const
{
  const int N[3]={Nx,Ny,Nz};
  const bool periodic[3]={false,false,true};
  const bool low[2]={i<HNx, j<HNy}, high[2]={i>=HNx, j>=HNy};
  edge_stencil(N,periodic,low,high,i,j,k,nm,np);
}

//################################ Channel ################################

Geometry_Channel::Geometry_Channel(const struct Simulation_Parameters * sim_param) : GEOMETRY (sim_param)
{
  geometry_pointer=&(*this);
  geometry_name="Channel";
  number_of_boundaries=2;
  bc_conditions=std::vector<class BOUNDARY *>(number_of_boundaries);
  boundary_needed_to_be_defined="0 and 1";

  point_type=fill_point_type( );
  dA=sim_param->dx*sim_param->dz;
}

int * Geometry_Channel::fill_point_type( void ) const
{
  int * point_kind;
  if((point_kind= (int *)calloc(Nx*Ny*Nz, sizeof(int)))==NULL){ERROr}

  for(int k=0; k<Nz; k++)
    for(int j=0; j<Ny; j++)
      for(int i=0; i<Nx; i++)
        point_kind[(k*Ny+j)*Nx+i]= j==0 ? 2 : (j==Ny-1 ? 3 : 1);
  return point_kind;
}

int Geometry_Channel::surface_site(int i, int j, int k, double v[3], double & area) const
{
  v[0]=0;
  v[1]= j==0 ? -1 : 1;
  v[2]=0;
  area=dA;
  return j==0 ? 0 : 1;
}

void Geometry_Channel::surface_stencil(int i, int j, int k, int nm[3], int np[3]) const
{
  const int N[3]={Nx,Ny,Nz};
  const bool periodic[3]={true,false,true};
  const bool low[2]={false, j==0}, high[2]={false, j==Ny-1};
  edge_stencil(N,periodic,low,high,i,j,k,nm,np);
}

//################################ Box ################################

Geometry_Box::Geometry_Box(const struct Simulation_Parameters * sim_param) : GEOMETRY (sim_param)
{
  geometry_pointer=&(*this);
  geometry_name="Box";
  number_of_boundaries=2;
  bc_conditions=std::vector<class BOUNDARY *>(number_of_boundaries);
  boundary_needed_to_be_defined="0 and 1";

  point_type=fill_point_type( );
}

int * Geometry_Box::fill_point_type( void ) const
{
  int * point_kind;
  if((point_kind= (int *)calloc(Nx*Ny*Nz, sizeof(int)))==NULL){ERROr}

  for(int k=0; k<Nz; k++)
    for(int j=0; j<Ny; j++)
      for(int i=0; i<Nx; i++)
        {
          const bool xwall= i==0 || i==Nx-1, ywall= j==0 || j==Ny-1;
          point_kind[(k*Ny+j)*Nx+i]= xwall ? 2 : (ywall ? 3 : 1);
        }
  return point_kind;
}

int Geometry_Box::surface_site(int i, int j, int k, double v[3], double & area) const
{
  const double nx= i==0 ? -1 : (i==Nx-1 ? 1 : 0);
  const double ny= j==0 ? -1 : (j==Ny-1 ? 1 : 0);
  const double norm=sqrt(nx*nx+ny*ny);
  v[0]=nx/norm;
  v[1]=ny/norm;
  v[2]=0;
  //Walls normal to x represent dy*dz, walls normal to y dx*dz; a corner half of each.
  const double dx=1/dx_1, dy=1/dy_1, dz=1/dz_1;
  area= (nx!=0 && ny!=0) ? 0.5*(dx+dy)*dz : (nx!=0 ? dy*dz : dx*dz);
  return point_type[(k*Ny+j)*Nx+i]-2;
}

void Geometry_Box::surface_stencil(int i, int j, int k, int nm[3], int np[3]) const
{
  const int N[3]={Nx,Ny,Nz};
  const bool periodic[3]={false,false,true};
  const bool low[2]={i==0, j==0}, high[2]={i==Nx-1, j==Ny-1};
  edge_stencil(N,periodic,low,high,i,j,k,nm,np);
}
