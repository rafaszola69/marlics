// GPU against CPU (make test-gpu):
//  - fill_ki and Energy_calc per site, on the initial condition of several inputs covering all
//    geometries the GPU supports, the four anchoring types, every elastic term and the field;
//  - FIRE and DP5 relaxations: same final energy and state on both devices.
#include "gpu.cu"
#include "exec_cpu.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

static int failed=0;

static void check(const char * what, double diff, double tol)
{
  const bool ok= diff<=tol;
  if(!ok) failed++;
  fprintf(stderr,"  %-58s %.2e  %s\n", what, diff, ok ? "ok" : "FAILED");
}

static double rel_diff(const std::vector<double>& a, const std::vector<double>& b)
{
  double d=0, m=0;
  for(size_t i=0; i<a.size(); i++) { d=fmax(d,fabs(a[i]-b[i])); m=fmax(m,fabs(a[i])); }
  return m>0 ? d/m : d;
}

static driver * setup(const char * input)
{
  driver * d=new driver;
  char name[300];
  strcpy(name,input);
  d->parse_input_file(name);
  d->setup_LC();
  d->setup_Simulation();
  return d;
}

static void compare_forces(const char * input)
{
  driver * d=setup(input);
  GEOMETRY * g=d->LcS_Geometry;
  const long N=(long)g->Nx*g->Ny*g->Nz;
  std::vector<double> kc(5*N), kg(5*N), ec(2*N), eg(2*N);

  #pragma omp parallel
  {
    g->fill_ki(kc.data(),d->Qij);
    #pragma omp barrier
    g->Energy_calc(ec.data(),d->Qij);
  }

  CudaExec exec(g);
  double * Q=exec.attach(d->Qij);
  double * k=exec.alloc(5*N);
  exec.fill_ki(k,Q);
  CUDA_CHECK(cudaMemcpy(kg.data(),k,5*N*sizeof(double),cudaMemcpyDeviceToHost));
  double * e=exec.alloc(2*N);
  LatticeView L=exec.view();
  gpu::bulk_energies<<<gpu::blocks(N),gpu::threads>>>(L,Q,e);
  if(L.n_surface>0) gpu::surface_energies<<<gpu::blocks(L.n_surface),gpu::threads>>>(L,Q,e);
  CUDA_CHECK(cudaMemcpy(eg.data(),e,2*N*sizeof(double),cudaMemcpyDeviceToHost));
  exec.release(k);
  exec.release(e);

  char what[300];
  sprintf(what,"%s fill_ki (max rel. diff)",input);
  check(what,rel_diff(kc,kg),1e-12);
  sprintf(what,"%s Energy_calc (max rel. diff)",input);
  check(what,rel_diff(ec,eg),1e-12);
  delete d;
}

//Runs the input's integrator to its tf on the CPU and on the GPU from the same initial condition.
static void compare_relaxation(const char * input, const char * integrator)
{
  double E[2];
  std::vector<double> Q[2];
  for(int dev=0; dev<2; dev++)
    {
      driver * d=setup(input);
      strcpy(d->sim_param.integrator_type,integrator);
      GEOMETRY * g=d->LcS_Geometry;
      Integrator * it= dev==0 ? make_integrator<CpuExec>(g,&d->sim_param) : make_integrator<CudaExec>(g,&d->sim_param);
      double t=d->sim_param.ti;
      while(t<d->sim_param.tf && it->evolve(d->Qij,&t,d->sim_param.tf)) {}
      CpuExec cpu(g);
      E[dev]=cpu.energy(d->Qij);
      Q[dev].assign(d->Qij,d->Qij+5L*g->Nx*g->Ny*g->Nz);
      delete d;
    }
  char what[300];
  sprintf(what,"%s %s final energy (rel. diff)",input,integrator);
  check(what,fabs(E[0]-E[1])/fabs(E[0]),1e-6);
  sprintf(what,"%s %s final Q (max abs. diff)",input,integrator);
  double d=0;
  for(size_t i=0; i<Q[0].size(); i++) d=fmax(d,fabs(Q[0][i]-Q[1][i]));
  check(what,d,1e-4);
}

int main()
{
  fprintf(stderr,"test_gpu: GPU against CPU\n");
  const char * inputs[]={"tests/gpu/sphere_homeotropic.txt","tests/gpu/sphere_rp.txt","tests/gpu/sphere_fg.txt",
                         "tests/gpu/sphere_strong.txt","tests/gpu/hybrid_sphere.txt","tests/gpu/hemisphere.txt",
                         "tests/gpu/slab.txt","tests/gpu/bulk.txt","tests/gpu/custom.txt",
                         "tests/inputs/2d_disk.txt","tests/inputs/2d_channel.txt","tests/inputs/2d_box.txt"};
  for(const char * in : inputs) compare_forces(in);

  compare_relaxation("tests/inputs/droplet.txt","FIRE");
  compare_relaxation("tests/inputs/slab.txt","FIRE");
  compare_relaxation("tests/inputs/slab.txt","DP5");

  fprintf(stderr, failed ? "SOME GPU CHECKS FAILED\n" : "All GPU checks passed\n");
  return failed ? 1 : 0;
}
