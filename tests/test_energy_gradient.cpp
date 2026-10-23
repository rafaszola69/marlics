// fill_ki (the right-hand side used by every integrator, FIRE included) must be the gradient of
// the discrete energy reported by Energy_calc: fill_ki = -Lambda * P[dE/dQ] / dV, where P is the
// projection onto symmetric traceless tensors written in the 5 independent components.
// Checked on a smooth periodic field in the bulk geometry, one elastic term at a time.
// The quadratic terms (L1, L2, Ls, Lq, thermotropic) are exact gradients. The cubic L3 term is
// only consistent: its mismatch must decrease at second order with the grid resolution.
#include "driver.h"
#include "geometry.h"
#include "geometry_bulk.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

struct Case { const char* name; double L1, L2, L3, Ls, Lq, q0, a, B, C; bool exact; };

static double total_energy(const GEOMETRY& g, std::vector<double>& e, const double* Q)
{
  g.Energy_calc(e.data(), Q);
  double E=0;
  for(double x : e) E+=x;
  return E;
}

static double run(const Case& c, int N)
{
  Simulation_Parameters p;
  p.Nx=p.Ny=p.Nz=N; p.dx=p.dy=p.dz=5.0;
  p.L1=c.L1; p.L2=c.L2; p.L3=c.L3; p.Ls=c.Ls; p.Lq=c.Lq; p.q0=c.q0;
  p.a=c.a; p.B=c.B; p.C=c.C; p.T=-1;
  p.mu_1=0.3; p.mu_1_s=3.0; p.S_eq=0.5;
  Geometry_Bulk g(&p);

  std::vector<double> Q(5*N*N*N), k(5*N*N*N), e(10*N*N*N);
  const double tp=2*M_PI/N, S=0.5;
  for(int z=0;z<N;z++) for(int y=0;y<N;y++) for(int x=0;x<N;x++)
    {
      double th=0.7*sin(tp*x)+0.4*cos(tp*z)+0.3, ph=0.5*cos(tp*y)+0.6*sin(tp*(x+z));
      double n[3]={sin(th)*cos(ph), sin(th)*sin(ph), cos(th)};
      double* q=&Q[5*((z*N+y)*N+x)];
      q[0]=S*(n[0]*n[0]-1/3.); q[1]=S*n[0]*n[1]; q[2]=S*n[0]*n[2]; q[3]=S*(n[1]*n[1]-1/3.); q[4]=S*n[1]*n[2];
    }
  g.fill_ki(k.data(), Q.data());

  const double Lambda=1/p.mu_1, dV=p.dx*p.dy*p.dz, h=1e-5;
  double num=0, den=0;
  for(int sample=0; sample<40; sample++)
    {
      const int s=(sample*997)%(N*N*N);
      double gr[5];
      for(int l=0;l<5;l++)
        {
          const double q0=Q[5*s+l];
          Q[5*s+l]=q0+h; const double Ep=total_energy(g,e,Q.data());
          Q[5*s+l]=q0-h; const double Em=total_energy(g,e,Q.data());
          Q[5*s+l]=q0;
          gr[l]=(Ep-Em)/(2*h);
        }
      const double expect[5]={-(2*gr[0]-gr[3])/3, -gr[1]/2, -gr[2]/2, -(2*gr[3]-gr[0])/3, -gr[4]/2};
      for(int l=0;l<5;l++)
        {
          const double ref=Lambda*expect[l]/dV;
          num+=pow(k[5*s+l]-ref,2); den+=pow(k[5*s+l],2);
        }
    }
  return sqrt(num/den);
}

int main()
{
  const Case cases[]={
    {"L1",          6,0,0,0,0,0,      0,0,0,            true },
    {"L2",          0,6,0,0,0,0,      0,0,0,            true },
    {"Ls",          0,0,0,6,0,0,      0,0,0,            true },
    {"thermotropic",0,0,0,0,0,0,      0.086,-2.12,1.74, true },
    {"Lq (chiral)", 0,0,0,0,6,0.03,   0,0,0,            true },
    {"L3",          0,0,3,0,0,0,      0,0,0,            false},
  };
  int failed=0;
  fprintf(stderr,"test_energy_gradient: fill_ki vs -Lambda*dE/dQ/dV (relative difference)\n");
  for(const Case& c : cases)
    {
      const double err=run(c,16);
      bool ok;
      if(c.exact)
        {
          ok=err<1e-7;
          fprintf(stderr,"  %-14s %.2e  %s\n", c.name, err, ok ? "ok" : "FAILED");
        }
      else
        {
          const double err2=run(c,32);
          ok=err<0.2 && err2<err/3;
          fprintf(stderr,"  %-14s %.2e (N=16) -> %.2e (N=32), second order: %s\n", c.name, err, err2, ok ? "ok" : "FAILED");
        }
      failed+=!ok;
    }
  return failed ? 1 : 0;
}
