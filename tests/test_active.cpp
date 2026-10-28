// Activity (active stress -zeta*q, with flow stokes):
//  - the active part of the flow is linear in zeta;
//  - zeta>0 is extensile: a bend perturbation of a uniform state grows, while it decays faster
//    than in the passive case for zeta<0 (linear theory: the bend response is ~ (lambda+1) > 0);
//  - for a tumbling material (xi < 3S/(2+S)) contractility destabilizes splay (~ (lambda-1) < 0);
//  - a +1/2 defect self-propels along its axis, towards its head for zeta>0, and a -1/2 does
//    not: velocity at the core from the active stress in a window around each defect;
//  - in 3D (no plane confinement) an out-of-plane bend perturbation, n=(cos,0,sin), grows for
//    zeta>0 and the flow has a z component.
//  - with elastic_stress 0 only the active stress drives the flow.
// Chiral activity (-zeta_c (eps q - q eps)/2):
//  - its flow is the active flow -zeta q of the director turned by -45 degrees about z (exact);
//  - a +1/2 defect self-propels at -atan(zeta_c/zeta) from its axis, a -1/2 does not.
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
#include <vector>

static int failed=0;

static void check(const char * what, bool ok, const char * fmt, double value)
{
  char v[100];
  sprintf(v,fmt,value);
  fprintf(stderr,"  %-60s %-12s %s\n", what, v, ok ? "ok" : "FAILED");
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

static std::vector<double> velocity(driver * d)
{
  GEOMETRY * g=d->LcS_Geometry;
  const long n=(long)g->Nx*g->Ny*g->Nz;
  std::vector<double> k(5*n);
  CpuExec plain(g);
  plain.fill_ki(k.data(),d->Qij);
  FlowCpu flow(g,&d->sim_param);
  flow.add(k.data(),d->Qij);
  return std::vector<double>(flow.velocity(),flow.velocity()+3*n);
}

static void linearity()
{
  std::vector<double> u[3];
  const double zeta[3]={0,0.01,0.02};
  for(int c=0;c<3;c++)
    {
      driver * d=setup([&](Simulation_Parameters& p){ p.activity=zeta[c]; });
      u[c]=velocity(d);
      delete d;
    }
  double err=0, m=0;
  for(size_t i=0;i<u[0].size();i++)
    {
      err=fmax(err,fabs((u[2][i]-u[0][i])-2*(u[1][i]-u[0][i])));
      m=fmax(m,fabs(u[2][i]-u[0][i]));
    }
  check("active flow linear in zeta: |u(2z)-u(0)-2(u(z)-u(0))|/max", err/m<1e-10, "%.1e", err/m);
  //elastic_stress 0: no flow without activity, and the active flow alone, u(z)-u(0), with activity
  driver * d=setup([&](Simulation_Parameters& p){ p.elastic_stress=0; });
  std::vector<double> u0=velocity(d);
  delete d;
  d=setup([&](Simulation_Parameters& p){ p.elastic_stress=0; p.activity=zeta[1]; });
  std::vector<double> u1=velocity(d);
  delete d;
  double e0=0, e1=0;
  for(size_t i=0;i<u0.size();i++)
    {
      e0=fmax(e0,fabs(u0[i]));
      e1=fmax(e1,fabs(u1[i]-(u[1][i]-u[0][i])));
    }
  check("elastic_stress 0, zeta=0: no flow, max|u| / max active", e0/m<1e-12, "%.1e", e0/m);
  check("elastic_stress 0: flow = active part, error / max", e1/m<1e-10, "%.1e", e1/m);
}

//Amplitude ratio of a director perturbation (bend: varying along the director x; splay: along y).
static double growth(bool bend, double zeta, double xi)
{
  driver * d=setup([&](Simulation_Parameters& p){
      p.Nx=p.Ny=32; p.activity=zeta; p.xi=xi;
      strcpy(p.initial_conditions,"homogeneous"); p.theta_i=90; p.phi_i=0;
      p.ic_flag[2]=p.ic_flag[3]=parameter_status::set; });
  const int N=32;
  const double S=d->sim_param.S_eq;
  auto amplitude=[&](){
    double a=0;
    for(int j=0;j<N;j++) for(int i=0;i<N;i++)
      {
        const double * q=&d->Qij[5*(j*N+i)];
        a+=0.5*atan2(2*q[1],q[0]-q[3])*sin(2*M_PI*(bend ? i : j)/N);
      }
    return a; };
  for(int j=0;j<N;j++) for(int i=0;i<N;i++)
    {
      const double phi=0.01*sin(2*M_PI*(bend ? i : j)/N), n[2]={cos(phi),sin(phi)};
      double * q=&d->Qij[5*(j*N+i)];
      q[0]=0.5*S*(3*n[0]*n[0]-1); q[1]=1.5*S*n[0]*n[1]; q[2]=0; q[3]=0.5*S*(3*n[1]*n[1]-1); q[4]=0;
    }
  const double a0=amplitude();
  Integrator * it=make_integrator<CpuExec>(d->LcS_Geometry,&d->sim_param);
  double t=0;
  it->evolve(d->Qij,&t,15.);
  const double r=amplitude()/a0;
  delete d;
  return r;
}

static void instabilities()
{
  const double b_ext=growth(true,0.05,0.6), b_pas=growth(true,0,0.6), b_con=growth(true,-0.05,0.6);
  check("bend, extensile (zeta=+0.05): amplitude ratio after t=15", b_ext>1, "%.3g", b_ext);
  check("bend, contractile (zeta=-0.05): decays faster than passive", b_con<b_pas, "%.3g", b_con/b_pas);
  const double s_con=growth(false,-0.05,0.3), s_pas=growth(false,0,0.3);
  check("splay, tumbling (xi=0.3), contractile: grows", s_con>1, "%.3g", s_con);
  check("splay, tumbling (xi=0.3), passive: decays", s_pas<1, "%.3g", s_pas);
}

//Velocity at a defect core due to the active stress in a Gaussian window (radius R sites) around it.
static void core_velocity(driver * d, double x0, double y0, double zeta, double R, double u[2], double zeta_c=0)
{
  GEOMETRY * g=d->LcS_Geometry;
  const int Nx=g->Nx, Ny=g->Ny;
  const long n=(long)Nx*Ny;
  const double dx=d->sim_param.dx, dy=d->sim_param.dy;
  std::vector<double> sxx(n), sxy(n), syy(n), f(3*n,0.);
  for(int j=0;j<Ny;j++) for(int i=0;i<Nx;i++)
    {
      const long s=(long)j*Nx+i;
      const double w=exp(-((i-x0)*(i-x0)+(j-y0)*(j-y0))/(2*R*R));
      const double * Q=&d->Qij[5*s];
      sxx[s]=-w*2./3.*(zeta*Q[0]+zeta_c*Q[1]); sxy[s]=-w*2./3.*(zeta*Q[1]+zeta_c*0.5*(Q[3]-Q[0]));
      syy[s]=-w*2./3.*(zeta*Q[3]-zeta_c*Q[1]);
    }
  for(int j=0;j<Ny;j++) for(int i=0;i<Nx;i++)
    {
      const long s=(long)j*Nx+i, xp=(long)j*Nx+(i+1)%Nx, xm=(long)j*Nx+(i+Nx-1)%Nx;
      const long yp=(long)((j+1)%Ny)*Nx+i, ym=(long)((j+Ny-1)%Ny)*Nx+i;
      f[3*s]  =(sxx[xp]-sxx[xm])/(2*dx)+(sxy[yp]-sxy[ym])/(2*dy);
      f[3*s+1]=(sxy[xp]-sxy[xm])/(2*dx)+(syy[yp]-syy[ym])/(2*dy);
    }
  FlowCpu flow(g,&d->sim_param);
  flow.solve_force(f.data());
  u[0]=u[1]=0;
  for(int dj=0;dj<2;dj++) for(int di=0;di<2;di++)
    {
      const long s=(long)((int)y0+dj)*Nx+(int)x0+di;
      u[0]+=flow.velocity()[3*s]/4; u[1]+=flow.velocity()[3*s+1]/4;
    }
}

static void self_propulsion()
{
  //+1/2 at (32.5,64.5), -1/2 at (96.5,64.5): near the +1/2 the director is theta/2-pi/2, so its
  //head (where the director wraps around the core) points to +x.
  driver * d=setup([](Simulation_Parameters& p){ p.Nx=p.Ny=128;
                                                 p.defect_x={32.5,96.5}; p.defect_y={64.5,64.5}; });
  double up[2], um[2], upc[2];
  core_velocity(d,32.5,64.5,0.01,8,up);
  core_velocity(d,96.5,64.5,0.01,8,um);
  core_velocity(d,32.5,64.5,-0.01,8,upc);
  check("+1/2, extensile: moves along its axis towards its head (+x)", up[0]>0 && fabs(up[1])<1e-3*fabs(up[0]), "%.3g", up[0]);
  check("+1/2, contractile: opposite direction", upc[0]<0, "%.3g", upc[0]);
  const double ratio=sqrt(um[0]*um[0]+um[1]*um[1])/fabs(up[0]);
  check("-1/2: no self-propulsion, |u(-1/2)|/|u(+1/2)|", ratio<0.05, "%.1e", ratio);
  delete d;
}

//Turns the director of every site by beta about z (Q -> R Q R^T).
static void rotate_q(driver * d, double beta)
{
  GEOMETRY * g=d->LcS_Geometry;
  const long n=(long)g->Nx*g->Ny*g->Nz;
  const double c=cos(beta), s=sin(beta), R[3][3]={{c,-s,0},{s,c,0},{0,0,1}};
  for(long i=0;i<n;i++)
    {
      double * v=&d->Qij[5*i];
      const double Q[3][3]={{v[0],v[1],v[2]},{v[1],v[3],v[4]},{v[2],v[4],-v[0]-v[3]}};
      double P[3][3];
      for(int a=0;a<3;a++) for(int b=0;b<3;b++)
        {
          P[a][b]=0;
          for(int k=0;k<3;k++) for(int l=0;l<3;l++) P[a][b]+=R[a][k]*Q[k][l]*R[b][l];
        }
      v[0]=P[0][0]; v[1]=P[0][1]; v[2]=P[0][2]; v[3]=P[1][1]; v[4]=P[1][2];
    }
}

static void chirality()
{
  //active parts of the flow: chiral (zeta_c=z) of the pair, achiral (zeta=z) of the turned pair
  std::vector<double> u[4];
  for(int c=0;c<4;c++)
    {
      driver * d=setup([&](Simulation_Parameters& p){
          if(c==1) p.chiral_activity=0.02;
          if(c==3) p.activity=0.02; });
      if(c>=2) rotate_q(d,-M_PI/4);
      u[c]=velocity(d);
      delete d;
    }
  double err=0, m=0;
  for(size_t i=0;i<u[0].size();i++)
    {
      err=fmax(err,fabs((u[1][i]-u[0][i])-(u[3][i]-u[2][i])));
      m=fmax(m,fabs(u[1][i]-u[0][i]));
    }
  check("chiral flow = active flow of the director turned by -45 deg", err/m<1e-10, "%.1e", err/m);

  driver * d=setup([](Simulation_Parameters& p){ p.Nx=p.Ny=128;
                                                 p.defect_x={32.5,96.5}; p.defect_y={64.5,64.5}; });
  double up[2], uc[2], um[2];
  core_velocity(d,32.5,64.5,0.01,8,up,0.01);
  core_velocity(d,32.5,64.5,0.,8,uc,0.01);
  core_velocity(d,96.5,64.5,0.01,8,um,0.01);
  const double tilt=atan2(up[1],up[0])*180/M_PI;
  check("+1/2, zeta=zeta_c=0.01: direction from its axis (deg, -45)", fabs(tilt+45)<1, "%.2f", tilt);
  check("+1/2, chiral only: moves along -y (perpendicular to its axis)", uc[1]<0 && fabs(uc[0])<1e-3*fabs(uc[1]), "%.3g", uc[1]);
  const double ratio=sqrt(um[0]*um[0]+um[1]*um[1])/sqrt(up[0]*up[0]+up[1]*up[1]);
  check("-1/2 with chiral activity: no self-propulsion", ratio<0.05, "%.1e", ratio);
  delete d;
}

static void bend_3d()
{
  double r[2];
  const double zeta[2]={0.05,0.};
  for(int c=0; c<2; c++)
    {
      driver * d=setup([&](Simulation_Parameters& p){
          p.dimension=3; p.Nx=32; p.Ny=p.Nz=8; p.grid_flag[2]=parameter_status::set; p.dz=p.dx;
          p.grid_spacing_flag[2]=parameter_status::set; p.plane_confinement=0; p.activity=zeta[c];
          strcpy(p.initial_conditions,"homogeneous"); p.theta_i=90; p.phi_i=0;
          p.ic_flag[2]=p.ic_flag[3]=parameter_status::set; });
      const int N=32, M=8;
      const double S=d->sim_param.S_eq;
      auto amplitude=[&](){
        double a=0;
        for(int k=0;k<M;k++) for(int j=0;j<M;j++) for(int i=0;i<N;i++)
          {
            const double * q=&d->Qij[5*((k*M+j)*N+i)];
            a+=0.5*atan2(2*q[2],q[0]-(-q[0]-q[3]))*sin(2*M_PI*i/N);   //angle in the xz plane
          }
        return a; };
      for(int k=0;k<M;k++) for(int j=0;j<M;j++) for(int i=0;i<N;i++)
        {
          const double phi=0.01*sin(2*M_PI*i/N), n[3]={cos(phi),0,sin(phi)};
          double * q=&d->Qij[5*((k*M+j)*N+i)];
          q[0]=0.5*S*(3*n[0]*n[0]-1); q[1]=0; q[2]=1.5*S*n[0]*n[2]; q[3]=-0.5*S; q[4]=0;
        }
      const double a0=amplitude();
      if(c==0)
        {
          const std::vector<double> u=velocity(d);
          double uz=0, um=0;
          for(size_t s=0; s<u.size()/3; s++) { uz=fmax(uz,fabs(u[3*s+2])); um=fmax(um,fabs(u[3*s])+fabs(u[3*s+1])); }
          check("3D bend: the flow is along z (max |u_z| / max |u_x,u_y|)", uz>10*um, "%.2g", uz/um);
        }
      Integrator * it=make_integrator<CpuExec>(d->LcS_Geometry,&d->sim_param);
      double t=0;
      it->evolve(d->Qij,&t,15.);
      r[c]=amplitude()/a0;
      delete d;
    }
  check("3D out-of-plane bend, extensile: amplitude ratio after t=15", r[0]>1, "%.3g", r[0]);
  check("3D out-of-plane bend, passive: decays", r[1]<1, "%.3g", r[1]);
}

int main()
{
  fprintf(stderr,"test_active: activity\n");
  linearity();
  instabilities();
  self_propulsion();
  chirality();
  bend_3d();
  return failed ? 1 : 0;
}
