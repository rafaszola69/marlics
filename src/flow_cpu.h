#ifndef FLOW_CPU_H
#define FLOW_CPU_H
// Backflow on the CPU: Stokes flow driven by the liquid crystal stress in a periodic box, solved
// with FFTs (pocketfft, src/third_party, used only in flow_cpu.cpp so that nvcc never compiles it),
// and its terms added to fill_ki. See flow_physics.h.
#include <complex>
#include <vector>
#include <iostream>
#include <cstdlib>
#include "driver.h"
#include "geometry.h"
#include "flow_physics.h"

class FlowCpu
{
 public:

  FlowCpu(GEOMETRY * g, const Simulation_Parameters * sim_param) : L(g->view()), n((long)L.Nx*L.Ny*L.Nz)
  {
    P.eta=sim_param->viscosity;
    P.xi=sim_param->xi;
    if(!L.periodic[0] || !L.periodic[1] || !L.periodic[2] || L.n_surface>0)
      {
        std::cout << "The flow (backflow) is only available in periodic domains (geometry bulk).\nAborting the program.\n";
        exit(1);
      }
    sigma.resize(9*n); u.resize(3*n);
    for(int a=0; a<3; a++) fk[a].resize(n);
    std::cout << "Backflow: Stokes flow in a periodic box (FFT), viscosity=" << P.eta << " Pa s, xi=" << P.xi << ".\n";
  }

  //Adds the flow terms to k_i (which holds Lambda*H on entry). Called outside parallel regions.
  void add(double * k_i, const double * Qij)
  {
    const LatticeView L=this->L;
    const FlowParams P=this->P;
    double * sg=sigma.data();

    #pragma omp parallel for schedule(static)
    for(long s=0; s<n; s++) site_stress(L,P,Qij,k_i,s,&sg[9*s]);

    #pragma omp parallel for schedule(static)
    for(long s=0; s<n; s++)
      {
        double f[3];
        site_force(L,sg,s,f);
        for(int a=0; a<3; a++) fk[a][s]=f[a];
      }
    solve();

    const double * us=u.data();
    #pragma omp parallel for schedule(static)
    for(long s=0; s<n; s++)
      {
        double r[5];
        site_flow_rhs(L,P,Qij,us,s,r);
        for(int c=0; c<5; c++) k_i[5*s+c]+=r[c];
      }
  }

  //Stokes flow for a given force density (3 per site), e.g. for tests.
  void solve_force(const double * f)
  {
    for(long s=0; s<n; s++)
      for(int a=0; a<3; a++) fk[a][s]=f[3*s+a];
    solve();
  }

  //Velocity of the last solution (nm/us), 3 per site.
  const double * velocity() const { return u.data(); }

 private:

  //Velocity from the force density held in fk.
  void solve()
  {
    const LatticeView L=this->L;
    const FlowParams P=this->P;
    fft(true);
    #pragma omp parallel for schedule(static)
    for(long s=0; s<n; s++)
      {
        const int m0=s%L.Nx, m1=(s/L.Nx)%L.Ny, m2=s/((long)L.Nx*L.Ny);
        double fr[3], fi[3];
        for(int a=0; a<3; a++) { fr[a]=fk[a][s].real(); fi[a]=fk[a][s].imag(); }
        stokes_mode(L,P.eta,m0,m1,m2,fr,fi);
        for(int a=0; a<3; a++) fk[a][s]=std::complex<double>(fr[a],fi[a]);
      }
    fft(false);

    double * us=u.data();
    #pragma omp parallel for schedule(static)
    for(long s=0; s<n; s++)
      for(int a=0; a<3; a++) us[3*s+a]=fk[a][s].real();
  }

  //In-place FFT of the three components of fk (flow_cpu.cpp).
  void fft(bool forward);

  const LatticeView L;
  FlowParams P;
  const long n;
  std::vector<double> sigma, u;
  std::vector<std::complex<double>> fk[3];
};

#endif
