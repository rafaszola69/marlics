// homeotropic_ansatz (bipolar cholesteric droplet, Krakhalev et al., Sci. Rep. 7, 14582 (2017)):
// finite everywhere (the poles are singular points of the construction), N=2(Nz-1)dz/p0_i pi-turns
// along the bipolar axis, director perpendicular to the axis there, and at the surface the points where
// the director is tangent follow the double helix of Eq. (2): azimuth = pi/2 - N*theta (mod pi).
#include "driver.h"
#include "geometry.h"
#include <cmath>
#include <cstdio>
#include <cstring>

static void director(const double* q, double n[3])
{
  const double Q[3][3]={{q[0],q[1],q[2]},{q[1],q[3],q[4]},{q[2],q[4],-q[0]-q[3]}};
  n[0]=0.3; n[1]=0.5; n[2]=0.8;
  for(int it=0; it<200; it++)
    {
      double m[3];
      for(int a=0;a<3;a++) { m[a]=n[a]; for(int b=0;b<3;b++) m[a]+=Q[a][b]*n[b]; }
      const double nn=sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
      for(int a=0;a<3;a++) n[a]=m[a]/nn;
    }
}

int main()
{
  char fname[]="tests/inputs/bipolar_droplet.txt";
  driver d;
  d.parse_input_file(fname);
  d.setup_LC();
  d.setup_Simulation();
  const GEOMETRY* g=d.LcS_Geometry;
  const int Nx=g->Nx, Ny=g->Ny, Nz=g->Nz;
  const double N=2*(Nz-1)*d.sim_param.dz/d.sim_param.p0_i, c=0.5*(Nx-1);

  bool finite=true;
  for(int s=0; s<5*Nx*Ny*Nz; s++) if(!(fabs(d.Qij[s])<10)) finite=false;

  //Director angle along the axis x=y=0, unwrapped from pole to pole.
  double turns=0, prev=0, nz_max=0;
  for(int k=0; k<Nz; k++)
    {
      const int s=(k*Ny+(int)c)*Nx+(int)c;
      if(!g->point_type[s]) continue;
      double n[3]; director(&d.Qij[5*s], n);
      nz_max=fmax(nz_max, fabs(n[2]));
      const double ang=atan2(n[1],n[0]);
      if(k>0) { double da=ang-prev; while(da> M_PI/2) da-=M_PI; while(da<-M_PI/2) da+=M_PI; turns+=da; }
      prev=ang;
    }
  turns=fabs(turns)/M_PI;

  //Tangent surface points near the equator vs Eq. (2).
  double dev_sum=0; int count=0;
  for(int k=0;k<Nz;k++) for(int j=0;j<Ny;j++) for(int i=0;i<Nx;i++)
    {
      const int s=(k*Ny+j)*Nx+i;
      if(g->point_type[s]!=2) continue;
      const double x=i-c, y=j-c, z=k-c, r=sqrt(x*x+y*y+z*z), th=acos(z/r);
      if(fabs(th-M_PI/2)>0.6) continue;
      double n[3]; director(&d.Qij[5*s], n);
      if(fabs((n[0]*x+n[1]*y+n[2]*z)/r)>0.15) continue;
      double dev=fmod(atan2(y,x)-(M_PI/2-N*th), M_PI);
      if(dev<0) dev+=M_PI;
      dev_sum+=fmin(dev, M_PI-dev); count++;
    }
  const double dev_mean=count ? dev_sum/count*180/M_PI : 90;

  const bool ok_f=finite, ok_t=fabs(turns-N)<0.05, ok_z=nz_max<1e-3, ok_d=count>100 && dev_mean<10;
  fprintf(stderr,"test_initial_conditions: homeotropic_ansatz (bipolar droplet, N=%.2f)\n", N);
  fprintf(stderr,"  finite everywhere                     %s\n", ok_f ? "ok" : "FAILED");
  fprintf(stderr,"  pi-turns along the axis: %.3f          %s\n", turns, ok_t ? "ok" : "FAILED");
  fprintf(stderr,"  max |n_z| on the axis: %.1e          %s\n", nz_max, ok_z ? "ok" : "FAILED");
  fprintf(stderr,"  tangent points vs Eq.(2): %d, mean %.1f deg  %s\n", count, dev_mean, ok_d ? "ok" : "FAILED");
  return (ok_f && ok_t && ok_z && ok_d) ? 0 : 1;
}
