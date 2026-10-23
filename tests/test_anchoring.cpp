// The anchoring energies must be consistent with the surface functional derivatives used by the
// dynamics: for the non-gradient part, functional_derivative = -Lambda_s * P[d f_s / dQ].
// The same check is done for the thermotropic bulk energy as a reference.
#include "driver.h"
#include "energy.h"
#include "energy_ldg.h"
#include "boundary.h"
#include "boundary_rp.h"
#include "boundary_fournier_galatola.h"
#include "boundary_homeotropic.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>

static double mismatch(const ENERGY& E, double Lambda)
{
  double dQ[15]={0}, ddQ[30]={0}, worst=0;
  srand(7);
  for(int trial=0; trial<50; trial++)
    {
      double Q[5], v[3];
      for(int k=0;k<5;k++) Q[k]=0.4*(rand()/(double)RAND_MAX-0.5);
      for(int k=0;k<3;k++) v[k]=rand()/(double)RAND_MAX-0.5;
      const double nv=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
      for(int k=0;k<3;k++) v[k]/=nv;

      double g[5];
      const double h=1e-6;
      for(int k=0;k<5;k++)
        {
          double Qp[5], Qm[5];
          for(int l=0;l<5;l++) { Qp[l]=Q[l]; Qm[l]=Q[l]; }
          Qp[k]+=h; Qm[k]-=h;
          g[k]=(E.energy_calculation(Qp,dQ,v)-E.energy_calculation(Qm,dQ,v))/(2*h);
        }
      const double expect[5]={-(2*g[0]-g[3])/3, -g[1]/2, -g[2]/2, -(2*g[3]-g[0])/3, -g[4]/2};
      const double fd[5]={E.functional_derivative_00(Q,dQ,ddQ,v), E.functional_derivative_01(Q,dQ,ddQ,v),
                          E.functional_derivative_02(Q,dQ,ddQ,v), E.functional_derivative_11(Q,dQ,ddQ,v),
                          E.functional_derivative_12(Q,dQ,ddQ,v)};
      double num=0, den=0;
      for(int k=0;k<5;k++) { num+=pow(fd[k]-Lambda*expect[k],2); den+=pow(fd[k],2); }
      worst=fmax(worst, sqrt(num/den));
    }
  return worst;
}

int main()
{
  Simulation_Parameters p;
  p.a=0.086; p.B=-2.12; p.C=1.74; p.T=-1; p.mu_1=0.3; p.mu_1_s=3.0; p.S_eq=0.53;
  p.dx=p.dy=p.dz=5;
  p.Wo1[0]=10; p.theta_0[0]=30; p.phi_0[0]=60;

  const landau_de_gennes bulk(&p);
  const Boundary_Rp rp(&p,0);
  const Boundary_Fg fg(&p,0);
  const Boundary_Homeotropic ho(&p,0);
  struct { const char* name; const ENERGY* E; double Lambda; } cases[]={
    {"thermotropic bulk", &bulk, 1/p.mu_1},
    {"Rapini-Papoular",   &rp,   1/p.mu_1_s},
    {"Fournier-Galatola", &fg,   1/p.mu_1_s},
    {"homeotropic",       &ho,   1/p.mu_1_s},
  };

  int failed=0;
  fprintf(stderr,"test_anchoring: functional derivative vs -Lambda*df/dQ (worst relative difference)\n");
  for(auto& c : cases)
    {
      const double err=mismatch(*c.E, c.Lambda);
      const bool ok=err<1e-7;
      failed+=!ok;
      fprintf(stderr,"  %-18s %.2e  %s\n", c.name, err, ok ? "ok" : "FAILED");
    }
  return failed ? 1 : 0;
}
