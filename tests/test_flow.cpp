// Backflow (flow stokes):
//  - the generated df/d(grad Q) (physics_ldg_gradient.h) matches finite differences of the energy;
//  - the FFT Stokes solver gives the exact discrete solution for a sinusoidal force and a
//    divergence-free velocity for a random one;
//  - a uniform state drives no flow;
//  - with flow the energy never increases (checks the signs of the stress against the coupling
//    terms of the Q equation), in 2D (defect pair) and 3D (random chiral state);
//  - the flow changes the dynamics but not the final state.
#include "driver.h"
#include "geometry.h"
#include "integrator.h"
#include "exec_cpu.h"
#include "integrator_factory.h"
#include "flow_cpu.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <random>
#include <vector>

static int failed=0;

static void check(const char * what, bool ok, const char * fmt, double value)
{
  char v[100];
  sprintf(v,fmt,value);
  fprintf(stderr,"  %-58s %-12s %s\n", what, v, ok ? "ok" : "FAILED");
  failed+=!ok;
}

static driver * setup(std::function<void(Simulation_Parameters&)> change=NULL)
{
  driver * d=new driver;
  char name[]="tests/inputs/flow_pair.txt";
  d->parse_input_file(name);
  if(change) change(d->sim_param);
  d->setup_LC();
  d->setup_Simulation();
  return d;
}

static void gradient_derivative()
{
  BulkParams P={};
  P.L1=6; P.L2=3; P.L3=2; P.Ls=1.5; P.Lq=5; P.Lq_tilde=2*5*0.1;
  std::mt19937 rng(3);
  std::uniform_real_distribution<double> U(-0.3,0.3);
  double QN[5], dQ[15], g[15], err=0, scale=0;
  const double v[3]={0,0,0};
  for(int trial=0; trial<20; trial++)
    {
      for(double& x : QN) x=U(rng);
      for(double& x : dQ) x=U(rng);
      ldg_gradient_derivative(P,QN,dQ,g);
      for(int m=0; m<15; m++)
        {
          const double h=1e-4, x0=dQ[m];
          dQ[m]=x0+h; const double fp=ldg_energy_calculation<double>(P,QN,dQ,v);
          dQ[m]=x0-h; const double fm=ldg_energy_calculation<double>(P,QN,dQ,v);
          dQ[m]=x0;
          err=fmax(err,fabs(g[m]-(fp-fm)/(2*h))); scale=fmax(scale,fabs(g[m]));
        }
    }
  check("df/d(grad Q) against finite differences (rel.)", err/scale<1e-8, "%.1e", err/scale);
}

static void stokes_solver()
{
  driver * d=setup([](Simulation_Parameters& p){ p.Nx=16; p.Ny=12; });
  GEOMETRY * g=d->LcS_Geometry;
  FlowCpu flow(g,&d->sim_param);
  const int Nx=16, Ny=12, N=Nx*Ny;
  const double eta=d->sim_param.viscosity, dy=d->sim_param.dy;
  std::vector<double> f(3*N,0.);
  for(int j=0;j<Ny;j++) for(int i=0;i<Nx;i++) f[3*(j*Nx+i)]=sin(2*M_PI*j/Ny);
  flow.solve_force(f.data());
  const double lap=(2-2*cos(2*M_PI/Ny))/(dy*dy);
  double err=0;
  for(int j=0;j<Ny;j++) for(int i=0;i<Nx;i++)
    {
      const double * u=&flow.velocity()[3*(j*Nx+i)];
      err=fmax(err,fabs(u[0]-sin(2*M_PI*j/Ny)/(eta*lap))+fabs(u[1])+fabs(u[2]));
    }
  check("Stokes: shear force sin(2 pi y/Ly), max error (nm/us)", err<1e-10, "%.1e", err);

  std::mt19937 rng(7);
  std::normal_distribution<double> G;
  for(double& x : f) x=G(rng);
  flow.solve_force(f.data());
  const double * u=flow.velocity();
  double div=0, umax=0;
  for(int j=0;j<Ny;j++) for(int i=0;i<Nx;i++)
    {
      const int s=j*Nx+i, xp=j*Nx+(i+1)%Nx, xm=j*Nx+(i+Nx-1)%Nx, yp=((j+1)%Ny)*Nx+i, ym=((j+Ny-1)%Ny)*Nx+i;
      div=fmax(div,fabs((u[3*xp]-u[3*xm])/(2*d->sim_param.dx)+(u[3*yp+1]-u[3*ym+1])/(2*dy)));
      umax=fmax(umax,fabs(u[3*s])+fabs(u[3*s+1]));
    }
  check("Stokes: random force, max |div u| / max |u| per nm", div/umax<1e-12, "%.1e", div/umax);
  delete d;
}

static double max_speed(driver * d)
{
  GEOMETRY * g=d->LcS_Geometry;
  const long N=(long)g->Nx*g->Ny*g->Nz;
  std::vector<double> k(5*N);
  CpuExec plain(g);
  plain.fill_ki(k.data(),d->Qij);
  FlowCpu flow(g,&d->sim_param);
  flow.add(k.data(),d->Qij);
  double m=0;
  for(long s=0; s<N; s++) m=fmax(m,sqrt(pow(flow.velocity()[3*s],2)+pow(flow.velocity()[3*s+1],2)+pow(flow.velocity()[3*s+2],2)));
  return m;
}

static void uniform_state()
{
  driver * d=setup([](Simulation_Parameters& p){ strcpy(p.initial_conditions,"homogeneous"); p.theta_i=90; p.ic_flag[2]=parameter_status::set; p.ic_flag[3]=parameter_status::set; });
  const double m=max_speed(d);
  check("uniform state: max |u| (nm/us)", m<1e-12, "%.1e", m);
  delete d;
}

//Energy at n intervals up to tf; returns the largest relative increase and the final Q.
static double energy_history(driver * d, double tf, int n, std::vector<double>& Q, double& E_end)
{
  Integrator * it=make_integrator<CpuExec>(d->LcS_Geometry,&d->sim_param);
  CpuExec plain(d->LcS_Geometry);
  double t=0, E=plain.energy(d->Qij), rise=0;
  const double E0=E;
  for(int m=1; m<=n; m++)
    {
      it->evolve(d->Qij,&t,tf*m/n);
      const double En=plain.energy(d->Qij);
      rise=fmax(rise,(En-E)/fabs(E0));
      E=En;
    }
  Q.assign(d->Qij,d->Qij+5L*d->LcS_Geometry->Nx*d->LcS_Geometry->Ny*d->LcS_Geometry->Nz);
  E_end=E;
  return rise;
}

static void energy_decreases()
{
  std::vector<double> Qf, Qn;
  double Ef, En;
  driver * d=setup([](Simulation_Parameters& p){ p.Nx=p.Ny=30; p.defect_x={8.5,20.5}; p.defect_y={15.5,15.5}; });
  const double u=max_speed(d);
  check("defect pair: backflow speed max |u| (nm/us)", u>1e-3, "%.2e", u);
  const double rise=energy_history(d,100,20,Qf,Ef);
  check("defect pair with flow: largest energy increase (rel.)", rise<1e-9, "%.1e", rise);
  delete d;

  driver * e=setup([](Simulation_Parameters& p){ strcpy(p.flow,"none"); p.Nx=p.Ny=30; p.defect_x={8.5,20.5}; p.defect_y={15.5,15.5}; });
  energy_history(e,100,20,Qn,En);
  double dq=0;
  for(size_t i=0;i<Qf.size();i++) dq=fmax(dq,fabs(Qf[i]-Qn[i]));
  check("defect pair at t=100: flow changes the dynamics, max |dQ|", dq>1e-3, "%.2e", dq);
  delete e;

  driver * r=setup([](Simulation_Parameters& p){
      p.dimension=3; p.Nx=p.Ny=p.Nz=10; p.grid_flag[2]=parameter_status::set; p.dz=p.dx; p.grid_spacing_flag[2]=parameter_status::set;
      p.p0=60; p.chirality_flag=chirality_status::p0;
      p.L2=3; p.Ls=1.5; p.Lq=6; p.plane_confinement=0;
      strcpy(p.initial_conditions,"random"); p.rng_seed=11; p.ic_flag[4]=parameter_status::set;
    });
  const double rise3=energy_history(r,20,20,Qf,Ef);
  check("3D random chiral state with flow: largest energy increase", rise3<1e-9, "%.1e", rise3);
  delete r;
}

//Long relaxation with and without flow: same equilibrium.
static void same_equilibrium()
{
  std::vector<double> Q[2];
  double E[2];
  for(int c=0; c<2; c++)
    {
      driver * d=setup([c](Simulation_Parameters& p){ if(c) strcpy(p.flow,"none"); p.Nx=p.Ny=16; p.Atol=p.Rtol=1e-5;
                                                     p.defect_x={4.5,11.5}; p.defect_y={8.5,8.5}; });
      energy_history(d,1500,3,Q[c],E[c]);
      delete d;
    }
  check("same final energy with and without flow (rel.)", fabs(E[0]-E[1])/fabs(E[1])<1e-6, "%.1e", fabs(E[0]-E[1])/fabs(E[1]));
}

int main()
{
  fprintf(stderr,"test_flow: backflow\n");
  gradient_derivative();
  stokes_solver();
  uniform_state();
  energy_decreases();
  same_equilibrium();
  return failed ? 1 : 0;
}
