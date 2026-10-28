#ifndef FLOW_CPU_H
#define FLOW_CPU_H
// Backflow on the CPU: Stokes flow driven by the liquid crystal stress in a periodic box, solved
// with FFTs (pocketfft, src/third_party, used only in flow_cpu.cpp so that nvcc never compiles it),
// between two flat walls (flow_walls.h) or in any other geometry (flow_confined.h), and its terms
// added to fill_ki. See flow_physics.h.
#include <complex>
#include <vector>
#include <iostream>
#include <cstdlib>
#include "driver.h"
#include "geometry.h"
#include "flow_physics.h"
#include "flow_walls.h"
#include "flow_confined.h"

class FlowCpu
{
 public:

  FlowCpu(GEOMETRY * g, const Simulation_Parameters * sim_param) : L(g->view()), n((long)L.Nx*L.Ny*L.Nz), walls(NULL), confined(NULL)
  {
    P.eta=sim_param->viscosity;
    P.xi=sim_param->xi;
    P.zeta=sim_param->activity;
    P.zeta_c=sim_param->chiral_activity;
    P.el=sim_param->elastic_stress;
    int axis;
    const double * U=sim_param->wall_velocity;
    if(WallStokes::detect(L,axis))
      {
        walls=new WallStokes(L,axis,P.eta,U);
        n_plane=n/walls->W.M;
        scratch.resize(5*walls->W.M*n_plane);
      }
    else
      {
        if(U[0]!=0 || U[1]!=0 || U[2]!=0)
          {
            std::cout << "wall_velocity needs two flat walls (geometries channel and slab).\nAborting the program.\n";
            exit(1);
          }
        if(ConfinedStokes::needed(L))
          {
            confined=new ConfinedStokes(L,P.eta,sim_param->flow_solver);
            force.resize(3*n);
            g_e2.resize(3*confined->ne); lambda.resize(3*confined->ne);
          }
      }
    sigma.resize(9*n); u.resize(3*n);
    for(int a=0; a<3; a++) fk[a].resize(n);
    std::cout << "Backflow: Stokes flow " << (walls ? "between walls" : (confined ? "in a confined geometry" : "in a periodic box (FFT)"))
              << ", viscosity=" << P.eta << " Pa s, xi=" << P.xi << ".\n";
    flow_alignment_report(P.xi,L.bulk.S_eq);
  }

  ~FlowCpu() { delete walls; delete confined; }

  //Adds the flow terms to k_i (which holds Lambda*H on entry). Called outside parallel regions.
  void add(double * k_i, const double * Qij)
  {
    const LatticeView L=this->L;
    const FlowParams P=this->P;
    double * sg=sigma.data();

    #pragma omp parallel for schedule(static)
    for(long s=0; s<n; s++) if(L.point_type[s]==1) site_stress(L,P,Qij,k_i,s,&sg[9*s]);
    if(walls || confined)
      {
        #pragma omp parallel for schedule(static)
        for(long s=0; s<n; s++) if(L.point_type[s]!=1) boundary_stress(L,s,sg);
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

  //Stokes flow for a given force density (3 per site), e.g. for tests. With surface sites the force
  //off the bulk is ignored.
  void solve_force(const double * f)
  {
    for(long s=0; s<n; s++)
      for(int a=0; a<3; a++) fk[a][s]= L.point_type[s]==1 ? f[3*s+a] : 0.;
    solve();
  }

  //Velocity of the last solution (nm/us), 3 per site.
  const double * velocity() const { return u.data(); }

  //Confined geometries: the Lagrange-force solve (NULL otherwise).
  const ConfinedStokes * confined_solver() const { return confined; }

 private:

  //Velocity from the force density held in fk.
  void solve()
  {
    const LatticeView L=this->L;
    const FlowParams P=this->P;
    if(walls) { solve_walls(); return; }
    if(confined) { solve_confined(); return; }
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

  //Confined geometry (flow_confined.h): u = G(f + lambda) + U0, with lambda on E2 making u vanish there.
  void solve_confined()
  {
    const long ne=confined->ne;
    const long * e2=confined->e2.data();
    #pragma omp parallel for schedule(static)
    for(long s=0; s<n; s++)
      for(int a=0; a<3; a++) force[3*s+a]=fk[a][s].real();

    //u_G = G f on E2, and the total force (uniform mode of the transform)
    double total[3];
    periodic_velocity(force.data(),total);
    for(long i=0; i<ne; i++)
      for(int a=0; a<3; a++) g_e2[3*i+a]=u[3*e2[i]+a];

    std::vector<double> field;
    confined->solve(g_e2.data(),total,lambda.data(),[&](const double * p, double * out){
        field.assign(3*n,0.);
        for(long i=0; i<ne; i++)
          for(int a=0; a<3; a++) field[3*e2[i]+a]=p[3*i+a];
        periodic_velocity(field.data());
        for(long i=0; i<ne; i++)
          for(int a=0; a<3; a++) out[3*i+a]=u[3*e2[i]+a];
      });

    for(long i=0; i<ne; i++)
      for(int a=0; a<3; a++) force[3*e2[i]+a]+=lambda[3*i+a];
    periodic_velocity(force.data());
    double U0[3]={0,0,0};
    for(long i=0; i<ne; i++)
      for(int a=0; a<3; a++) U0[a]-=u[3*e2[i]+a];
    for(int a=0; a<3; a++) U0[a]/=ne;
    double * us=u.data();
    #pragma omp parallel for schedule(static)
    for(long s=0; s<n; s++)
      for(int a=0; a<3; a++) us[3*s+a]= L.point_type[s]==1 ? us[3*s+a]+U0[a] : 0.;
  }

  //u = G f: periodic Stokes flow (checkerboards kept, uniform mode removed) of the force f (3 per
  //site); total (if given): the sum of f.
  void periodic_velocity(const double * f, double * total=NULL)
  {
    const LatticeView L=this->L;
    const double eta=P.eta;
    for(long s=0; s<n; s++)
      for(int a=0; a<3; a++) fk[a][s]=f[3*s+a];
    fft(true);
    if(total) for(int a=0; a<3; a++) total[a]=fk[a][0].real();
    #pragma omp parallel for schedule(static)
    for(long s=0; s<n; s++)
      {
        const int m0=s%L.Nx, m1=(s/L.Nx)%L.Ny, m2=s/((long)L.Nx*L.Ny);
        double fr[3], fi[3];
        for(int a=0; a<3; a++) { fr[a]=fk[a][s].real(); fi[a]=fk[a][s].imag(); }
        stokes_mode(L,eta,m0,m1,m2,fr,fi,true);
        for(int a=0; a<3; a++) fk[a][s]=std::complex<double>(fr[a],fi[a]);
      }
    fft(false);
    double * us=u.data();
    #pragma omp parallel for schedule(static)
    for(long s=0; s<n; s++)
      for(int a=0; a<3; a++) us[3*s+a]=fk[a][s].real();
  }

  //In-place FFT of the three components of fk (flow_cpu.cpp): over every axis in a periodic box,
  //over the two in-plane axes between walls.
  void fft(bool forward);

  const LatticeView L;
  FlowParams P;
  const long n;
  WallStokes * walls;
  ConfinedStokes * confined;
  long n_plane;
  std::vector<double> sigma, u, scratch, force, g_e2, lambda;
  std::vector<std::complex<double>> fk[3];
};

#endif
