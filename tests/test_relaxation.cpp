// Dynamics on small but complete simulations (geometry, anchoring, chirality, L1-L3 from k11-k33):
//  - DP5 relaxes as a gradient flow of the reported energy: the total energy never increases;
//  - FIRE converges, and its final state is stationary for fill_ki (max force below min_force).
#include "driver.h"
#include "geometry.h"
#include "integrator.h"
#include "integrator_dp5.h"
#include "fire_integrator.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

static double total_energy(const GEOMETRY* g, std::vector<double>& e, const double* Q)
{
  g->Energy_calc(e.data(), Q);
  double E=0;
  for(double x : e) E+=x;
  return E;
}

static bool dp5_monotone(const char* input)
{
  char fname[512];
  strncpy(fname, input, sizeof(fname)-1);
  fname[sizeof(fname)-1]='\0';
  driver d;
  d.parse_input_file(fname);
  d.setup_LC();
  d.setup_Simulation();
  GEOMETRY* g=d.LcS_Geometry;
  Integrator* dp5=new DP5(g, &d.sim_param);
  std::vector<double> e(10*g->Nx*g->Ny*g->Nz);

  double t=0, E=total_energy(g,e,d.Qij), E0=E;
  int increases=0;
  for(int it=0; it<40; it++)
    {
      dp5->evolve(d.Qij, &t, t+0.5);
      const double En=total_energy(g,e,d.Qij);
      if(En>E+1e-10*fabs(E)) increases++;
      E=En;
    }
  const bool ok=increases==0;
  fprintf(stderr,"  %-24s DP5: E %.6g -> %.6g, %d increases in 40 intervals  %s\n", input, E0, E, increases, ok ? "ok" : "FAILED");
  return ok;
}

static bool fire_converges(const char* input)
{
  char fname[512];
  strncpy(fname, input, sizeof(fname)-1);
  fname[sizeof(fname)-1]='\0';
  driver d;
  d.parse_input_file(fname);
  d.setup_LC();
  d.setup_Simulation();
  GEOMETRY* g=d.LcS_Geometry;
  Integrator* fire=new FIRE(g, &d.sim_param);
  const int N=5*g->Nx*g->Ny*g->Nz;

  double t=0;
  bool running=true;
  while(running && t<5000) running=fire->evolve(d.Qij, &t, t+100);

  std::vector<double> k(N);
  g->fill_ki(k.data(), d.Qij);
  double fmax=0;
  for(double x : k) fmax=std::max(fmax, fabs(x));
  const bool ok=!running && fmax<d.sim_param.min_force;
  fprintf(stderr,"  %-24s FIRE: converged at t=%g, max|force|=%.2e  %s\n", input, t, fmax, ok ? "ok" : "FAILED");
  return ok;
}

int main()
{
  const char* inputs[]={"tests/inputs/droplet.txt", "tests/inputs/slab.txt"};
  int failed=0;
  fprintf(stderr,"test_relaxation:\n");
  for(const char* in : inputs)
    {
      failed+=!dp5_monotone(in);
      failed+=!fire_converges(in);
    }
  return failed ? 1 : 0;
}
