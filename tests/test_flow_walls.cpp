// Backflow between flat walls (flow stokes in the channel and slab geometries, flow_walls.h):
//  - pressure-driven (Poiseuille) and wall-driven (Couette) flows are the exact discrete profiles;
//  - for a random force the velocity vanishes on the walls, is divergence-free (central differences)
//    and is the exact constrained Stokes solution: <v, eta A u - f> = 0 for any divergence-free v
//    vanishing on the walls (A: 3-point Laplacian), i.e. u minimizes the dissipation minus the work;
//  - the solution operator is symmetric and positive (<u1,f2>=<u2,f1>, <u,f> > 0);
//  - a uniform state drives no flow, with activity, so the stress at the walls adds no spurious force;
//  - under shear (moving upper wall) the director aligns at the flow-aligning angle (for a uniaxial
//    Q, cos 2 theta = 3S/((2+S) xi)), uniformly across the channel when the walls do not anchor it;
//  - with resting walls the energy does not increase.
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
  fprintf(stderr,"  %-62s %-12s %s\n", what, v, ok ? "ok" : "FAILED");
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

//Sites, wall axis and lattice of a run.
struct Box
{
  LatticeView L;
  long n;
  int axis, M;
  long st;
  Box(driver * d) : L(d->LcS_Geometry->view())
  {
    n=(long)L.Nx*L.Ny*L.Nz;
    WallStokes::detect(L,axis);
    M=axis_size(L,axis);
    st=axis_step(L,axis);
  }
  int row(long s) const { return (int)((s/st)%M); }
  bool bulk(long s) const { return L.point_type[s]==1; }
};

static double dot(const std::vector<double>& a, const std::vector<double>& b)
{
  double r=0;
  for(size_t i=0;i<a.size();i++) r+=a[i]*b[i];
  return r;
}

static std::vector<double> solve(FlowCpu& flow, const Box& B, const std::vector<double>& f)
{
  flow.solve_force(f.data());
  return std::vector<double>(flow.velocity(),flow.velocity()+3*B.n);
}

//Poiseuille (uniform force along an in-plane axis a) and Couette (upper wall moving along a).
static void profiles(const char * input, int a)
{
  const double U=3.0, f0=0.01;
  double Uw[3]={0,0,0};
  Uw[a]=U;
  driver * d=setup(input,[&](Simulation_Parameters& p){ for(int c=0;c<3;c++) p.wall_velocity[c]=Uw[c]; });
  Box B(d);
  const double eta=d->sim_param.viscosity, h=1/B.L.d_1[B.axis], H=(B.M-1)*h;
  FlowCpu flow(d->LcS_Geometry,&d->sim_param);

  std::vector<double> f(3*B.n,0.);
  std::vector<double> u=solve(flow,B,f);
  double err=0;
  for(long s=0;s<B.n;s++)
    for(int c=0;c<3;c++) err=fmax(err,fabs(u[3*s+c]-(c==a ? U*B.row(s)/(B.M-1.) : 0)));
  char what[300];
  sprintf(what,"%s: Couette along %c, max error (nm/us)",d->sim_param.geometry,"xyz"[a]);
  check(what, err<1e-12*U, "%.1e", err);

  delete d;
  d=setup(input);
  Box C(d);
  for(long s=0;s<C.n;s++) f[3*s+a]=f0;
  FlowCpu still(d->LcS_Geometry,&d->sim_param);
  u=solve(still,C,f);
  err=0;
  double umax=0;
  for(long s=0;s<C.n;s++)
    {
      const double y=C.row(s)*h, exact=f0/(2*eta)*y*(H-y);
      for(int c=0;c<3;c++) err=fmax(err,fabs(u[3*s+c]-(c==a ? exact : 0)));
      umax=fmax(umax,exact);
    }
  sprintf(what,"%s: Poiseuille along %c, max error / max u",d->sim_param.geometry,"xyz"[a]);
  check(what, err<1e-12*umax, "%.1e", err/umax);
  delete d;
}

//Random forces: walls, divergence, variational characterization, symmetry and positivity.
static void random_force(const char * input)
{
  driver * d=setup(input);
  Box B(d);
  const LatticeView& L=B.L;
  const double eta=d->sim_param.viscosity;
  FlowCpu flow(d->LcS_Geometry,&d->sim_param);
  const char * name=d->sim_param.geometry;
  char what[300];

  std::mt19937 rng(17);
  std::normal_distribution<double> G;
  std::vector<double> f1(3*B.n,0.), f2(3*B.n,0.);
  for(long s=0;s<B.n;s++)
    if(B.bulk(s)) for(int c=0;c<3;c++) { f1[3*s+c]=G(rng); f2[3*s+c]=G(rng); }
  const std::vector<double> u1=solve(flow,B,f1), u2=solve(flow,B,f2);

  double wall=0, div=0, umax=0;
  for(long s=0;s<B.n;s++)
    {
      for(int c=0;c<3;c++) umax=fmax(umax,fabs(u1[3*s+c]));
      if(!B.bulk(s)) { for(int c=0;c<3;c++) wall=fmax(wall,fabs(u1[3*s+c])); continue; }
      double dv=0;
      for(int a=0;a<3;a++)
        dv+=0.5*(u1[3*periodic_neighbor(L,s,a,1)+a]-u1[3*periodic_neighbor(L,s,a,-1)+a])*L.d_1[a];
      div=fmax(div,fabs(dv));
    }
  sprintf(what,"%s: random force, max |u| on the walls",name);
  check(what, wall==0, "%.1e", wall);
  sprintf(what,"%s: random force, max |div u| / max |u| per nm",name);
  check(what, div/umax<1e-12, "%.1e", div/umax);

  //r = eta A u1 - f1 at bulk sites (A = minus the 3-point Laplacian, u=0 on the walls).
  std::vector<double> r(3*B.n,0.);
  for(long s=0;s<B.n;s++)
    {
      if(!B.bulk(s)) continue;
      for(int c=0;c<3;c++)
        {
          double lap=0;
          for(int a=0;a<3;a++)
            lap+=(u1[3*periodic_neighbor(L,s,a,1)+c]-2*u1[3*s+c]+u1[3*periodic_neighbor(L,s,a,-1)+c])*L.d_1[a]*L.d_1[a];
          r[3*s+c]=-eta*lap-f1[3*s+c];
        }
    }
  const double orth=fabs(dot(u2,r))/sqrt(dot(u2,u2)*dot(f1,f1));
  sprintf(what,"%s: <v, eta A u - f> for divergence-free v (rel.)",name);
  check(what, orth<1e-12, "%.1e", orth);

  const double a12=dot(u1,f2), a21=dot(u2,f1);
  sprintf(what,"%s: symmetry |<u1,f2>-<u2,f1>| / |<u1,f2>|",name);
  check(what, fabs(a12-a21)<1e-12*fabs(a12), "%.1e", fabs(a12-a21)/fabs(a12));
  sprintf(what,"%s: dissipation <u,f> > 0",name);
  check(what, dot(u1,f1)>0, "%.3g", dot(u1,f1));
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
  driver * d=setup("tests/inputs/flow_channel.txt",[](Simulation_Parameters& p){ p.activity=0.02; });
  const double m=max_speed(d);
  check("channel: uniform active state, max |u| (nm/us)", m<1e-12, "%.1e", m);
  delete d;
}

//Director angle in the middle of the channel after shearing from the director along x. A steady
//uniform state aligns at cos 2 theta = (l_n - l_m)/(xi (1 + l_n + l_m)), l_n > l_m the in-plane
//eigenvalues of Q (projection of the corotation term on the director; for a uniaxial Q, l_m=-S/2,
//this is 3S/((2+S) xi)). The plane confinement makes Q biaxial, so the uniaxial form is not used.
static void shear_alignment()
{
  const double U=20;
  driver * d=setup("tests/inputs/flow_channel.txt",[&](Simulation_Parameters& p){ p.Ny=24; p.wall_velocity[0]=U; p.Wo1[0]=p.Wo1[1]=0; });
  const double xi=d->sim_param.xi;
  Integrator * it=make_integrator<CpuExec>(d->LcS_Geometry,&d->sim_param);
  double t=0;
  it->evolve(d->Qij,&t,150.);
  const int Nx=d->LcS_Geometry->Nx, Ny=24;
  double amax=-10, amin=10, theta=0, ln=0, lm=0;
  for(int j=1;j<Ny-1;j++)
    for(int i=0;i<Nx;i++)
      {
        const double * q=&d->Qij[5*(j*Nx+i)];
        const double ang=0.5*atan2(2*q[1],q[0]-q[3]);
        amax=fmax(amax,ang); amin=fmin(amin,ang);
        if(j==Ny/2 && i==0)
          {
            const double r=sqrt(0.25*(q[0]-q[3])*(q[0]-q[3])+q[1]*q[1]);
            theta=ang; ln=0.5*(q[0]+q[3])+r; lm=0.5*(q[0]+q[3])-r;
          }
      }
  const double leslie=0.5*acos((ln-lm)/(xi*(1+ln+lm)));
  check("channel shear: director angle - aligning angle (degrees)", fabs(theta-leslie)*180/M_PI<0.02, "%.4f", (theta-leslie)*180/M_PI);
  check("channel shear: aligning angle (degrees), towards the extension axis", theta>0, "%.2f", theta*180/M_PI);
  check("channel shear: angle uniform across the channel (spread, degrees)", (amax-amin)*180/M_PI<0.2, "%.3f", (amax-amin)*180/M_PI);
  delete d;
}

//Largest relative energy increase over n intervals up to tf.
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

static void energy_decreases()
{
  driver * d=setup("tests/inputs/flow_channel.txt",[](Simulation_Parameters& p){
      p.Nx=30; strcpy(p.initial_conditions,"defects"); p.defect_x={8.5,20.5}; p.defect_y={9.5,9.5}; p.defect_k={0.5,-0.5}; });
  const double u=max_speed(d);
  check("channel, defect pair: backflow speed max |u| (nm/us)", u>1e-3, "%.2e", u);
  const double rise=energy_rise(d,50,20);
  check("channel, defect pair with flow: largest energy increase (rel.)", rise<1e-9, "%.1e", rise);
  delete d;

  driver * e=setup("tests/inputs/flow_slab.txt");
  const double rise3=energy_rise(e,10,20);
  check("slab, random chiral state with flow: largest energy increase", rise3<1e-9, "%.1e", rise3);
  delete e;
}

int main()
{
  fprintf(stderr,"test_flow_walls: backflow between flat walls\n");
  profiles("tests/inputs/flow_channel.txt",0);
  profiles("tests/inputs/flow_slab.txt",0);
  profiles("tests/inputs/flow_slab.txt",1);
  random_force("tests/inputs/flow_channel.txt");
  random_force("tests/inputs/flow_slab.txt");
  uniform_state();
  shear_alignment();
  energy_decreases();
  return failed ? 1 : 0;
}
