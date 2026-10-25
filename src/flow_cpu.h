#ifndef FLOW_CPU_H
#define FLOW_CPU_H
// Backflow on the CPU: Stokes flow driven by the liquid crystal stress in a periodic box, solved
// with FFTs (pocketfft, src/third_party, used only in flow_cpu.cpp so that nvcc never compiles it),
// or between two flat walls (flow_walls.h), and its terms added to fill_ki. See flow_physics.h.
#include <complex>
#include <vector>
#include <iostream>
#include <cstdlib>
#include "driver.h"
#include "geometry.h"
#include "flow_physics.h"
#include "flow_walls.h"

class FlowCpu
{
 public:

  FlowCpu(GEOMETRY * g, const Simulation_Parameters * sim_param) : L(g->view()), n((long)L.Nx*L.Ny*L.Nz), walls(NULL)
  {
    P.eta=sim_param->viscosity;
    P.xi=sim_param->xi;
    P.zeta=sim_param->activity;
    int axis;
    if(WallStokes::detect(L,axis))
      {
        walls=new WallStokes(L,axis,P.eta,sim_param->wall_velocity);
        n_plane=n/walls->W.M;
        scratch.resize(5*walls->W.M*n_plane);
      }
    sigma.resize(9*n); u.resize(3*n);
    for(int a=0; a<3; a++) fk[a].resize(n);
    std::cout << "Backflow: Stokes flow " << (walls ? "between walls" : "in a periodic box (FFT)") << ", viscosity=" << P.eta << " Pa s, xi=" << P.xi << ".\n";
    flow_alignment_report(P.xi,L.bulk.S_eq);
  }

  ~FlowCpu() { delete walls; }

  //Adds the flow terms to k_i (which holds Lambda*H on entry). Called outside parallel regions.
  void add(double * k_i, const double * Qij)
  {
    const LatticeView L=this->L;
    const FlowParams P=this->P;
    double * sg=sigma.data();

    #pragma omp parallel for schedule(static)
    for(long s=0; s<n; s++) if(L.point_type[s]==1) site_stress(L,P,Qij,k_i,s,&sg[9*s]);
    if(walls)
      {
        const WallFlowView W=walls->W;
        #pragma omp parallel for schedule(static)
        for(long s=0; s<n; s++) if(L.point_type[s]!=1) wall_stress(L,W,s,sg);
      }

    #pragma omp parallel for schedule(static)
    for(long s=0; s<n; s++)
      {
        double f[3]={0,0,0};
        if(L.point_type[s]==1) site_force(L,sg,s,f);
        for(int a=0; a<3; a++) fk[a][s]=f[a];
      }
    solve();

    const double * us=u.data();
    #pragma omp parallel for schedule(static)
    for(long s=0; s<n; s++)
      {
        if(L.point_type[s]!=1) continue;
        double r[5];
        site_flow_rhs(L,P,Qij,us,s,r);
        for(int c=0; c<5; c++) k_i[5*s+c]+=r[c];
      }
  }

  //Stokes flow for a given force density (3 per site), e.g. for tests. Between walls the force on
  //the wall rows is ignored.
  void solve_force(const double * f)
  {
    for(long s=0; s<n; s++)
      for(int a=0; a<3; a++) fk[a][s]= L.point_type[s]==1 ? f[3*s+a] : 0.;
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
    if(walls) { solve_walls(); return; }
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

  //Between walls: in-plane transforms, one solve per in-plane mode, wall rows set to the wall velocities.
  void solve_walls()
  {
    const LatticeView L=this->L;
    const WallFlowView W=walls->W;
    const double eta=P.eta;
    const long np=n_plane;
    double * const F[3]={(double *)fk[0].data(),(double *)fk[1].data(),(double *)fk[2].data()};
    double * sc=scratch.data();
    fft(true);
    #pragma omp parallel for schedule(static)
    for(long t=0; t<np; t++) wall_stokes_mode(L,W,eta,t,F,&sc[5*W.M*t],np);
    fft(false);

    double * us=u.data();
    const long st=axis_step(L,W.axis);
    #pragma omp parallel for schedule(static)
    for(long s=0; s<n; s++)
      {
        const int k=(int)((s/st)%W.M);
        for(int a=0; a<3; a++) us[3*s+a]= k==W.M-1 ? W.U[a] : (k==0 ? 0. : fk[a][s].real());
      }
  }

  //In-place FFT of the three components of fk (flow_cpu.cpp): over every axis in a periodic box,
  //over the two in-plane axes between walls.
  void fft(bool forward);

  const LatticeView L;
  FlowParams P;
  const long n;
  WallStokes * walls;
  long n_plane;
  std::vector<double> sigma, u, scratch;
  std::vector<std::complex<double>> fk[3];
};

#endif
