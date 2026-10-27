// Two-dimensional simulations (dimension 2):
//  - a single layer gives the same forces and energies as the z-invariant 3D system (Nz=3) for the
//    disk, ellipse, channel and box geometries, every elastic term, chirality and plane confinement on;
//  - with planar degenerate (fg) walls the director escapes along z, unless plane_confinement keeps
//    it in the plane;
//  - the ellipse with equal semi-axes is the disk;
//  - the "defects" initial condition has the requested winding around each defect;
//  - "cholesteric_2d" has the handedness preferred by the chiral term.
#include "driver.h"
#include "geometry.h"
#include "integrator.h"
#include "exec_cpu.h"
#include "integrator_factory.h"
#include "initial_conditions.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <vector>

static int failed=0;

static void check(const char * what, bool ok, const char * fmt, double value)
{
  char v[100];
  sprintf(v,fmt,value);
  fprintf(stderr,"  %-52s %-12s %s\n", what, v, ok ? "ok" : "FAILED");
  failed+=!ok;
}

static driver * setup(const char * input, std::function<void(Simulation_Parameters&)> change=NULL)
{
  driver * d=new driver;
  char name[300];
  strcpy(name,input);
  d->parse_input_file(name);
  if(change) change(d->sim_param);
  d->setup_LC();
  d->setup_Simulation();
  return d;
}

//Smooth z-invariant field with an out-of-plane component.
static void smooth_field(driver * d)
{
  const int Nx=d->sim_param.Nx, Ny=d->sim_param.Ny, Nz=d->sim_param.Nz;
  for(int k=0;k<Nz;k++) for(int j=0;j<Ny;j++) for(int i=0;i<Nx;i++)
    {
      const double th=0.9*sin(0.37*i+0.2*j)+0.4, ph=0.8*cos(0.29*j)+0.5*sin(0.17*i), S=0.45;
      const double n[3]={sin(th)*cos(ph),sin(th)*sin(ph),cos(th)};
      double * q=&d->Qij[5*((k*Ny+j)*Nx+i)];
      q[0]=0.5*S*(3*n[0]*n[0]-1); q[1]=1.5*S*n[0]*n[1]; q[2]=1.5*S*n[0]*n[2]; q[3]=0.5*S*(3*n[1]*n[1]-1); q[4]=1.5*S*n[1]*n[2];
    }
}

static void z_invariance(const char * input)
{
  std::vector<double> k[2], e[2];
  int Nxy=0;
  for(int c=0; c<2; c++)
    {
      driver * d=setup(input, [c](Simulation_Parameters& p){ if(c==1) { p.dimension=3; p.Nz=3; p.grid_flag[2]=parameter_status::set; } });
      smooth_field(d);
      const long N=(long)d->sim_param.Nx*d->sim_param.Ny*d->sim_param.Nz;
      Nxy=d->sim_param.Nx*d->sim_param.Ny;
      k[c].resize(5*N); e[c].resize(2*N);
      #pragma omp parallel
      {
        d->LcS_Geometry->fill_ki(k[c].data(),d->Qij);
        #pragma omp barrier
        d->LcS_Geometry->Energy_calc(e[c].data(),d->Qij);
      }
      delete d;
    }
  //The 2D layer against the middle layer of the 3D system.
  double dk=0, mk=0, de=0, me=0;
  for(int s=0; s<Nxy; s++)
    {
      for(int l=0;l<5;l++) { dk=fmax(dk,fabs(k[0][5*s+l]-k[1][5*(Nxy+s)+l])); mk=fmax(mk,fabs(k[0][5*s+l])); }
      for(int l=0;l<2;l++) { de=fmax(de,fabs(e[0][2*s+l]-e[1][2*(Nxy+s)+l])); me=fmax(me,fabs(e[0][2*s+l])); }
    }
  char what[200];
  sprintf(what,"%s: 2D = z-invariant 3D, force",strrchr(input,'/')+1);
  check(what, dk/mk<1e-12, "%.1e", dk/mk);
  sprintf(what,"%s: 2D = z-invariant 3D, energy",strrchr(input,'/')+1);
  check(what, de/me<1e-12, "%.1e", de/me);
}

//The ellipse with a=b=N*dx/2 is the disk: same sites, forces and energies.
static void ellipse_is_disk()
{
  std::vector<double> k[2], e[2];
  std::vector<int> pt[2];
  for(int c=0; c<2; c++)
    {
      driver * d=setup("tests/inputs/2d_disk.txt", [c](Simulation_Parameters& p){ if(c==1) strcpy(p.geometry,"ellipse"); });
      smooth_field(d);
      const long N=(long)d->sim_param.Nx*d->sim_param.Ny;
      k[c].resize(5*N); e[c].resize(2*N);
      pt[c].assign(d->LcS_Geometry->point_type,d->LcS_Geometry->point_type+N);
      #pragma omp parallel
      {
        d->LcS_Geometry->fill_ki(k[c].data(),d->Qij);
        #pragma omp barrier
        d->LcS_Geometry->Energy_calc(e[c].data(),d->Qij);
      }
      delete d;
    }
  double dk=0, mk=0, de=0, me=0;
  for(size_t s=0; s<k[0].size(); s++) { dk=fmax(dk,fabs(k[0][s]-k[1][s])); mk=fmax(mk,fabs(k[0][s])); }
  for(size_t s=0; s<e[0].size(); s++) { de=fmax(de,fabs(e[0][s]-e[1][s])); me=fmax(me,fabs(e[0][s])); }
  check("ellipse a=b: same sites as the disk", pt[0]==pt[1], "%g", 0.);
  check("ellipse a=b: same force as the disk", dk/mk<1e-12, "%.1e", dk/mk);
  check("ellipse a=b: same energy as the disk", de/me<1e-12, "%.1e", de/me);
}

static void relax(driver * d)
{
  Integrator * it=make_integrator<CpuExec>(d->LcS_Geometry,&d->sim_param);
  double t=0;
  it->evolve(d->Qij,&t,d->sim_param.tf);
}

//Mean n_z^2 over the liquid crystal sites, from Q_zz=(S/2)(3 n_z^2-1) with S=S_eq.
static double mean_nz2(driver * d)
{
  const long N=(long)d->sim_param.Nx*d->sim_param.Ny;
  const double S=d->sim_param.S_eq;
  double m=0; int n=0;
  for(long s=0; s<N; s++)
    if(d->LcS_Geometry->point_type[s]!=0)
      {
        const double Qzz=-(d->Qij[5*s]+d->Qij[5*s+3]);
        m+=(2*Qzz/S+1)/3; n++;
      }
  return m/n;
}

static void planar_escape(double W)
{
  driver * d=setup("tests/inputs/2d_box.txt", [W](Simulation_Parameters& p){
      p.L2=p.L3=p.Ls=p.Lq=0;
      p.plane_confinement=W;
      p.anchoring_type[1]="fg"; p.Wo1[1]=1.0;
      p.tf=20000;
    });
  relax(d);
  const double nz2=mean_nz2(d);
  char what[200];
  if(W==0)
    {
      sprintf(what,"box, fg walls, no confinement: escapes, <n_z^2>");
      check(what, nz2>0.9, "%.3f", nz2);
    }
  else
    {
      sprintf(what,"box, fg walls, plane_confinement %g: planar, <n_z^2>",W);
      check(what, nz2<0.01, "%.1e", nz2);
    }
  delete d;
}

//Winding of the in-plane director around a square loop of half side h centred at (x,y), in units of 2 pi.
static double winding(driver * d, double x, double y, int h)
{
  const int Nx=d->sim_param.Nx;
  std::vector<int> I, J;
  const int i0=(int)floor(x)-h+1, i1=(int)floor(x)+h, j0=(int)floor(y)-h+1, j1=(int)floor(y)+h;
  for(int i=i0;i<i1;i++) { I.push_back(i); J.push_back(j0); }
  for(int j=j0;j<j1;j++) { I.push_back(i1); J.push_back(j); }
  for(int i=i1;i>i0;i--) { I.push_back(i); J.push_back(j1); }
  for(int j=j1;j>j0;j--) { I.push_back(i0); J.push_back(j); }
  double w=0, prev=0;
  for(size_t m=0; m<=I.size(); m++)
    {
      const double * q=&d->Qij[5*(J[m%I.size()]*Nx+I[m%I.size()])];
      const double a=atan2(2*q[1], q[0]-q[3]);   //twice the director angle
      if(m>0) { double da=a-prev; while(da>M_PI) da-=2*M_PI; while(da<-M_PI) da+=2*M_PI; w+=da; }
      prev=a;
    }
  return w/(4*M_PI);
}

static void defects()
{
  const double X[3]={6.5,12.5,18.5}, K[3]={0.5,-0.5,1.0};
  driver * d=setup("tests/inputs/2d_box.txt", [&](Simulation_Parameters& p){
      strcpy(p.initial_conditions,"defects");
      for(int m=0;m<3;m++) { p.defect_x.push_back(X[m]); p.defect_y.push_back(12.5); p.defect_k.push_back(K[m]); }
    });
  for(int m=0;m<3;m++)
    {
      const double w=winding(d,X[m],12.5,2);
      char what[200];
      sprintf(what,"defects IC: winding around a %+.1f defect",K[m]);
      check(what, fabs(w-K[m])<1e-9, "%+.3f", w);
    }
  delete d;
}

//Energy of the in-plane helix with the handedness of q0 and with the opposite one.
static void cholesteric_handedness()
{
  driver * d=setup("tests/inputs/2d_box.txt", [](Simulation_Parameters& p){
      strcpy(p.geometry,"bulk");
      p.Nx=24; p.Ny=4;
      p.L2=p.L3=p.Ls=0; p.Lq=p.L1;
      p.plane_confinement=0;
      strcpy(p.initial_conditions,"cholesteric_2d");
      p.p0_i=60; p.ic_flag[5]=parameter_status::set; p.phi_i=0;
    });
  CpuExec exec(d->LcS_Geometry);
  const double E_same=exec.energy(d->Qij);
  Simulation_Parameters opposite=d->sim_param;
  opposite.q0=-opposite.q0;
  cholesteric_2d_ic(&opposite,d->Qij,d->LcS_Geometry);
  const double E_opposite=exec.energy(d->Qij);
  check("cholesteric_2d: E(handedness of q0) < E(opposite)", E_same<E_opposite, "%.3g", E_same-E_opposite);
  delete d;
}

int main()
{
  fprintf(stderr,"test_2d: two-dimensional simulations\n");
  z_invariance("tests/inputs/2d_disk.txt");
  z_invariance("tests/inputs/2d_ellipse.txt");
  ellipse_is_disk();
  z_invariance("tests/inputs/2d_channel.txt");
  z_invariance("tests/inputs/2d_box.txt");
  planar_escape(0);
  planar_escape(0.02);
  defects();
  cholesteric_handedness();
  return failed ? 1 : 0;
}
