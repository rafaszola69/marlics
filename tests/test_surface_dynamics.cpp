// Surface dynamics (surface_dynamics local|variational):
//  - variational: the rate of every surface site is the exact gradient of the discrete energy,
//    fill_ki = -Lambda_s * P[dE/dQ] / area (P the traceless projection of test_energy_gradient), for
//    curved (disk, ellipse, sphere) and flat (box, slab) surfaces, every elastic term, chirality, plane
//    confinement and the three anchoring kinds; strong anchoring stays fixed;
//  - variational: the energy of a planar degenerate (fg) disk with two +1/2 defects never rises
//    (with the local scheme it rises for a while, see lattice.h), and the defects settle at
//    r = 5^(-1/4) R = 0.67 R, the strong-anchoring one-constant result (Duclos et al., Nat. Phys. 13,
//    58 (2017)), already at this coarse resolution (radius 16 sites), where the local scheme lets them
//    reach the wall (the two schemes agree from a radius of about 32 sites).
#include "driver.h"
#include "geometry.h"
#include "integrator.h"
#include "exec_cpu.h"
#include "integrator_factory.h"
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
  fprintf(stderr,"  %-58s %-12s %s\n", what, v, ok ? "ok" : "FAILED");
  failed+=!ok;
}

static driver * setup(const char * input, std::function<void(Simulation_Parameters&)> change)
{
  driver * d=new driver;
  char name[300];
  strcpy(name,input);
  d->parse_input_file(name);
  change(d->sim_param);
  d->setup_LC();
  d->setup_Simulation();
  return d;
}

static double total_energy(const GEOMETRY * g, std::vector<double>& e, const double * Q)
{
  g->Energy_calc(e.data(),Q);
  double E=0;
  for(double x : e) E+=x;
  return E;
}

//Largest difference between the surface rate and -Lambda_s P[dE/dQ]/area, relative to the largest rate.
static double surface_gradient_error(const char * input, std::function<void(Simulation_Parameters&)> change)
{
  driver * d=setup(input,[&](Simulation_Parameters& p){ change(p); p.surface_variational=true; strcpy(p.flow,"none"); });
  const GEOMETRY * g=d->LcS_Geometry;
  const LatticeView L=g->view();
  const long n=(long)g->Nx*g->Ny*g->Nz;
  std::vector<double> k(5*n), e(2*n);
  double * Q=d->Qij;
  g->fill_ki(k.data(),Q);

  const double h=1e-6;
  double err=0, kmax=0;
  const int step= L.n_surface>80 ? L.n_surface/80 : 1;
  for(int m=0; m<L.n_surface; m+=step)
    {
      const SurfaceSite& S=L.surface[m];
      double gr[5];
      for(int c=0;c<5;c++)
        {
          const double q=Q[5*S.site+c];
          Q[5*S.site+c]=q+h; const double Ep=total_energy(g,e,Q);
          Q[5*S.site+c]=q-h; const double Em=total_energy(g,e,Q);
          Q[5*S.site+c]=q;
          gr[c]=(Ep-Em)/(2*h);
        }
      const double c=L.bc[S.bc].kind==BC_STRONG ? 0. : L.bulk.Lambda_s/S.area;
      const double expect[5]={-c*(2*gr[0]-gr[3])/3, -c*gr[1]/2, -c*gr[2]/2, -c*(2*gr[3]-gr[0])/3, -c*gr[4]/2};
      for(int l=0;l<5;l++)
        {
          err=fmax(err,fabs(k[5*S.site+l]-expect[l]));
          kmax=fmax(kmax,fabs(k[5*S.site+l]));
        }
    }
  return err/kmax;
}

//Planar degenerate disk with two +1/2 defects relaxed with DP5 up to tf: the largest rise of the
//energy between snapshots (relative to the energy released), and the distance of the defect cores
//(smallest tr Q^2 inside the disk) to the centre, over the radius.
static double disk_relaxation(bool variational, double tf, double& defect_r)
{
  driver * d=setup("tests/inputs/flow_disk.txt",[&](Simulation_Parameters& p)
    {
      strcpy(p.flow,"none");
      p.surface_variational=variational;
      strcpy(p.initial_conditions,"defects");
      p.defect_x={10.5,20.5}; p.defect_y={15.5,15.5}; p.defect_k={0.5,0.5};
    });
  Integrator * it=make_integrator<CpuExec>(d->LcS_Geometry,&d->sim_param);
  CpuExec plain(d->LcS_Geometry);
  double t=0, E=0, E0=0, rise=0;
  for(int m=1; m<=40; m++)
    {
      it->evolve(d->Qij,&t,tf*m/40);
      const double En=plain.energy(d->Qij);
      if(m==1) E0=En;
      else rise=fmax(rise,En-E);
      E=En;
    }
  const int N=d->sim_param.Nx;
  const double c=(N-1)/2., R=N/2.;
  double qmin=1e30;
  for(int j=0;j<N;j++) for(int i=0;i<N;i++)
    {
      const double r=hypot(i-c,j-c);
      if(r>=R-1) continue;
      const double * q=&d->Qij[5*(j*N+i)];
      const double q2=q[0]*q[0]+q[3]*q[3]+(q[0]+q[3])*(q[0]+q[3])+2*(q[1]*q[1]+q[2]*q[2]+q[4]*q[4]);
      if(q2<qmin){ qmin=q2; defect_r=r/R; }
    }
  return rise/fabs(E-E0);
}

int main()
{
  fprintf(stderr,"test_surface_dynamics\n");
  auto all_terms=[](Simulation_Parameters& p){ p.L2=3.0; p.L3=2.0; p.Ls=1.0; };
  check("disk, homeotropic, all elastic terms, chiral, plane W",
        surface_gradient_error("tests/inputs/2d_disk.txt",all_terms)<1e-6, "%.2e",
        surface_gradient_error("tests/inputs/2d_disk.txt",all_terms));
  const double ell=surface_gradient_error("tests/inputs/2d_ellipse.txt",all_terms);
  check("ellipse 70 x 45 nm, homeotropic, all terms", ell<1e-6, "%.2e", ell);
  const double ellfg=surface_gradient_error("tests/inputs/2d_ellipse.txt",[](Simulation_Parameters& p){ p.anchoring_type[0]="fg"; });
  check("ellipse, fg", ellfg<1e-6, "%.2e", ellfg);
  const double hex=surface_gradient_error("tests/inputs/2d_polygon.txt",all_terms);
  check("hexagon, homeotropic, all terms", hex<1e-6, "%.2e", hex);
  const double sq=surface_gradient_error("tests/inputs/2d_polygon.txt",[](Simulation_Parameters& p){ p.anchoring_type[0]="fg"; p.polygon_sides=4; });
  check("square (polygon), fg", sq<1e-6, "%.2e", sq);
  const double tri=surface_gradient_error("tests/inputs/2d_polygon.txt",[](Simulation_Parameters& p){ p.anchoring_type[0]="fg"; p.polygon_sides=3; });
  check("triangle (polygon), fg", tri<1e-6, "%.2e", tri);
  const double box=surface_gradient_error("tests/inputs/2d_box.txt",all_terms);
  check("box, fg and rp walls with corners, all terms", box<1e-6, "%.2e", box);
  const double fgdisk=surface_gradient_error("tests/inputs/flow_disk.txt",[](Simulation_Parameters& p){ strcpy(p.initial_conditions,"random"); });
  check("disk, fg", fgdisk<1e-6, "%.2e", fgdisk);
  const double sphere=surface_gradient_error("tests/inputs/droplet.txt",all_terms);
  check("sphere, homeotropic, all terms, chiral", sphere<1e-6, "%.2e", sphere);
  const double slab=surface_gradient_error("tests/inputs/slab.txt",all_terms);
  check("slab, homeotropic and rp walls, all terms", slab<1e-6, "%.2e", slab);
  const double strong=surface_gradient_error("tests/inputs/slab.txt",[](Simulation_Parameters& p)
    { p.anchoring_type[1]="strong"; });
  check("slab with a strong wall (stays fixed)", strong<1e-6, "%.2e", strong);

  double r_local, r_var;
  const double rise_local=disk_relaxation(false,40,r_local);
  fprintf(stderr,"  fg disk, local scheme: largest energy rise %.2e of the energy released\n", rise_local);
  const double rise_var=disk_relaxation(true,300,r_var);
  check("fg disk, variational: energy never rises", rise_var<=0, "%.2e", rise_var);
  check("fg disk, variational: defects at r/R = 0.67 +- 0.05", fabs(r_var-0.67)<0.05, "%.3f", r_var);

  return failed ? 1 : 0;
}
