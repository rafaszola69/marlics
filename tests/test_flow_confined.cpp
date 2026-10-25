// Backflow in confined geometries (flow_confined.h: disk, box, droplet) and between flat walls
// (flow_walls.h), which solve the same discrete problem:
//  - for a random force the velocity vanishes off the bulk, is divergence-free (central differences)
//    at the bulk sites and at the non-bulk sites next to them (no flux into the walls), and is the
//    exact constrained Stokes solution: <v, eta A u - f> = 0 for any such v (A: 3-point Laplacian);
//  - the solution operator is symmetric and positive;
//  - the force of an isotropic stress (a gradient) is balanced by the pressure: no flow;
//  - the direct and iterative (conjugate gradients) solves agree;
//  - a body force f = (-y, x) turns the fluid as the continuum solution, u_theta = r (R^2 - r^2)/(8 eta)
//    in a disk, (R^2 - r^2)/(10 eta) r_perp in a sphere, up to the staircase boundary (first order in h);
//  - with backflow the energy does not increase (homeotropic anchoring: with fg anchoring in the disk
//    the surface relaxation itself, without flow, can raise the energy for a while).
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
  fprintf(stderr,"  %-66s %-12s %s\n", what, v, ok ? "ok" : "FAILED");
  failed+=!ok;
}

static driver * setup(const char * input, std::function<void(Simulation_Parameters&)> change=NULL)
{
  driver * d=new driver;
  char name[100];
  strcpy(name,input);
  d->parse_input_file(name);
  if(change) change(d->sim_param);
  d->setup_LC();
  d->setup_Simulation();
  return d;
}

static double dot(const std::vector<double>& a, const std::vector<double>& b)
{
  double r=0;
  for(size_t i=0;i<a.size();i++) r+=a[i]*b[i];
  return r;
}

static double max_abs(const std::vector<double>& a)
{
  double m=0;
  for(double v : a) m=fmax(m,fabs(v));
  return m;
}

struct Run
{
  driver * d;
  LatticeView L;
  long n;
  Run(const char * input, std::function<void(Simulation_Parameters&)> change=NULL) : d(setup(input,change)), L(d->LcS_Geometry->view())
  {
    n=(long)L.Nx*L.Ny*L.Nz;
  }
  ~Run() { delete d; }
  bool bulk(long s) const { return L.point_type[s]==1; }
  bool next_to_bulk(long s) const
  {
    for(int a=0;a<3;a++)
      for(int dir=-1;dir<=1;dir+=2) if(bulk(periodic_neighbor(L,(int)s,a,dir))) return true;
    return false;
  }
  //Position (nm) from the centre of the box.
  double x(long s, int a) const
  {
    const int c[3]={(int)(s%L.Nx),(int)((s/L.Nx)%L.Ny),(int)(s/((long)L.Nx*L.Ny))}, N[3]={L.Nx,L.Ny,L.Nz};
    return (c[a]-(N[a]-1)/2.)/L.d_1[a];
  }
  std::vector<double> solve(FlowCpu& flow, const std::vector<double>& f)
  {
    flow.solve_force(f.data());
    return std::vector<double>(flow.velocity(),flow.velocity()+3*n);
  }
};

//Random forces: walls, divergence, variational characterization, symmetry and positivity; direct
//against iterative solve.
static void random_force(const char * input, bool compare_solvers)
{
  Run R(input);
  const LatticeView& L=R.L;
  const double eta=R.d->sim_param.viscosity;
  FlowCpu flow(R.d->LcS_Geometry,&R.d->sim_param);
  const char * name=R.d->sim_param.geometry;
  char what[300];

  std::mt19937 rng(17);
  std::normal_distribution<double> G;
  std::vector<double> f1(3*R.n,0.), f2(3*R.n,0.);
  for(long s=0;s<R.n;s++)
    if(R.bulk(s)) for(int c=0;c<3;c++) { f1[3*s+c]=G(rng); f2[3*s+c]=G(rng); }
  const std::vector<double> u1=R.solve(flow,f1), u2=R.solve(flow,f2);

  double wall=0, div=0;
  const double umax=max_abs(u1);
  for(long s=0;s<R.n;s++)
    {
      if(!R.bulk(s)) for(int c=0;c<3;c++) wall=fmax(wall,fabs(u1[3*s+c]));
      if(!R.bulk(s) && !R.next_to_bulk(s)) continue;
      double dv=0;
      for(int a=0;a<3;a++)
        dv+=0.5*(u1[3*periodic_neighbor(L,(int)s,a,1)+a]-u1[3*periodic_neighbor(L,(int)s,a,-1)+a])*L.d_1[a];
      div=fmax(div,fabs(dv));
    }
  sprintf(what,"%s: random force, max |u| off the bulk",name);
  check(what, wall==0, "%.1e", wall);
  sprintf(what,"%s: random force, max |div u| (bulk and walls) / max |u| per nm",name);
  check(what, div/umax<1e-12, "%.1e", div/umax);

  std::vector<double> r(3*R.n,0.);
  for(long s=0;s<R.n;s++)
    {
      if(!R.bulk(s)) continue;
      for(int c=0;c<3;c++)
        {
          double lap=0;
          for(int a=0;a<3;a++)
            lap+=(u1[3*periodic_neighbor(L,(int)s,a,1)+c]-2*u1[3*s+c]+u1[3*periodic_neighbor(L,(int)s,a,-1)+c])*L.d_1[a]*L.d_1[a];
          r[3*s+c]=-eta*lap-f1[3*s+c];
        }
    }
  const double orth=fabs(dot(u2,r))/sqrt(dot(u2,u2)*dot(f1,f1));
  sprintf(what,"%s: <v, eta A u - f> for divergence-free v (rel.)",name);
  check(what, orth<1e-11, "%.1e", orth);

  const double a12=dot(u1,f2), a21=dot(u2,f1);
  sprintf(what,"%s: symmetry |<u1,f2>-<u2,f1>| / |<u1,f2>|",name);
  check(what, fabs(a12-a21)<1e-11*fabs(a12), "%.1e", fabs(a12-a21)/fabs(a12));
  sprintf(what,"%s: dissipation <u,f> > 0",name);
  check(what, dot(u1,f1)>0, "%.3g", dot(u1,f1));

  if(compare_solvers)
    {
      Simulation_Parameters p=R.d->sim_param;
      p.flow_solver=2;
      FlowCpu iterative(R.d->LcS_Geometry,&p);
      const std::vector<double> v1=R.solve(iterative,f1);
      std::vector<double> diff(3*R.n);
      for(long k=0;k<3*R.n;k++) diff[k]=v1[k]-u1[k];
      sprintf(what,"%s: iterative (%d CG iterations) - direct solve, max / max |u|",name,iterative.confined_solver()->iterations);
      check(what, max_abs(diff)<1e-8*umax, "%.1e", max_abs(diff)/umax);
    }
}

//Isotropic stress Pi I (a smooth bump off the centre): its force, with the stress of the non-bulk
//sites from boundary_stress, is a gradient and moves no fluid.
static void gradient_force(const char * input)
{
  Run R(input);
  const LatticeView& L=R.L;
  FlowCpu flow(R.d->LcS_Geometry,&R.d->sim_param);
  const double w=0.3*L.Nx/L.d_1[0];
  std::vector<double> sigma(9*R.n,0.), f(3*R.n,0.), rot(3*R.n,0.);
  for(long s=0;s<R.n;s++)
    if(R.bulk(s))
      {
        double r2=0;
        for(int a=0;a<3;a++) if(axis_size(L,a)>1) r2+=pow(R.x(s,a)-0.2*w,2);
        for(int a=0;a<3;a++) sigma[9*s+4*a]=exp(-r2/(w*w));
      }
  for(long s=0;s<R.n;s++) if(!R.bulk(s)) boundary_stress(L,s,sigma.data());
  for(long s=0;s<R.n;s++) if(R.bulk(s)) site_force(L,sigma.data(),(int)s,&f[3*s]);
  //A force of the same size that is not a gradient: (-y, x).
  double rmax=0;
  for(long s=0;s<R.n;s++) if(R.bulk(s)) rmax=fmax(rmax,hypot(R.x(s,0),R.x(s,1)));
  for(long s=0;s<R.n;s++)
    if(R.bulk(s)) { rot[3*s]=-R.x(s,1)/rmax*max_abs(f); rot[3*s+1]=R.x(s,0)/rmax*max_abs(f); }
  const double ug=max_abs(R.solve(flow,f)), ur=max_abs(R.solve(flow,rot));
  char what[300];
  sprintf(what,"%s: isotropic stress, max |u| / that of a rotational force",R.d->sim_param.geometry);
  check(what, ug<1e-10*ur, "%.1e", ug/ur);
}

//Body force (-y, x) f0 in a disk (2D) or sphere: rigid-like rotation profile. Returns the largest
//error over the bulk relative to the largest speed, with R = (N-1)/2 h (half-way into the surface shell).
static double rotation_error(const char * input, int N)
{
  Run R(input,[&](Simulation_Parameters& p){ p.Nx=p.Ny=N; if(p.dimension!=2) p.Nz=N; });
  const LatticeView& L=R.L;
  const double eta=R.d->sim_param.viscosity, f0=1e-3, Rad=(N-1)/2./L.d_1[0];
  const bool sphere= L.Nz>1;
  FlowCpu flow(R.d->LcS_Geometry,&R.d->sim_param);
  std::vector<double> f(3*R.n,0.);
  for(long s=0;s<R.n;s++)
    if(R.bulk(s)) { f[3*s]=-f0*R.x(s,1); f[3*s+1]=f0*R.x(s,0); }
  const std::vector<double> u=R.solve(flow,f);
  double err=0, umax=0;
  for(long s=0;s<R.n;s++)
    {
      if(!R.bulk(s)) continue;
      const double x=R.x(s,0), y=R.x(s,1), z= sphere ? R.x(s,2) : 0, r2=x*x+y*y+z*z;
      const double g= sphere ? f0*(Rad*Rad-r2)/(10*eta) : f0*(Rad*Rad-r2)/(8*eta);
      const double ex[3]={-g*y,g*x,0};
      for(int a=0;a<3;a++) err=fmax(err,fabs(u[3*s+a]-ex[a]));
      umax=fmax(umax,hypot(ex[0],ex[1]));
    }
  return err/umax;
}

static double energy_rise(driver * d, double tf, int n)
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
  return rise;
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

static void energy_decreases()
{
  driver * d=setup("tests/inputs/flow_disk.txt",[](Simulation_Parameters& p){
      strcpy(p.initial_conditions,"defects"); p.defect_x={10.5,20.5}; p.defect_y={15.5,15.5}; p.defect_k={0.5,0.5};
      p.anchoring_type[0]="homeotropic"; p.Wo1[0]=2.0; });
  check("disk (homeotropic), two +1/2 defects: backflow speed max |u| (nm/us)", max_speed(d)>1e-3, "%.2e", max_speed(d));
  const double rise=energy_rise(d,50,20);
  check("disk (homeotropic), two +1/2 defects with flow: largest energy increase (rel.)", rise<1e-9, "%.1e", rise);
  delete d;

  driver * e=setup("tests/inputs/flow_droplet.txt");
  check("droplet, random chiral state: backflow speed max |u| (nm/us)", max_speed(e)>1e-3, "%.2e", max_speed(e));
  const double rise3=energy_rise(e,5,10);
  check("droplet, random chiral state with flow: largest energy increase", rise3<1e-9, "%.1e", rise3);
  delete e;
}

int main()
{
  fprintf(stderr,"test_flow_confined: backflow in confined geometries\n");
  random_force("tests/inputs/flow_disk.txt",true);
  random_force("tests/inputs/flow_box.txt",true);
  random_force("tests/inputs/flow_droplet.txt",true);
  random_force("tests/inputs/flow_channel.txt",false);
  random_force("tests/inputs/flow_slab.txt",false);
  gradient_force("tests/inputs/flow_disk.txt");
  gradient_force("tests/inputs/flow_box.txt");
  gradient_force("tests/inputs/flow_droplet.txt");
  gradient_force("tests/inputs/flow_channel.txt");
  gradient_force("tests/inputs/flow_slab.txt");

  const double e1=rotation_error("tests/inputs/flow_disk.txt",32), e2=rotation_error("tests/inputs/flow_disk.txt",64);
  check("disk N=32: rotation by a body force, max error / max u", e1<0.2, "%.3f", e1);
  check("disk N=64: rotation by a body force, max error / max u", e2<0.1, "%.3f", e2);
  check("disk: error ratio N=32 / N=64 (first order: about 2)", e1/e2>1.4, "%.2f", e1/e2);
  const double e3=rotation_error("tests/inputs/flow_droplet.txt",16);
  check("sphere N=16: rotation by a body force, max error / max u", e3<0.3, "%.3f", e3);

  energy_decreases();
  return failed ? 1 : 0;
}
