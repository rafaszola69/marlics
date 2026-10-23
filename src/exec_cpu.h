#ifndef EXEC_CPU_H
#define EXEC_CPU_H
// CPU (OpenMP) executor for the integrators. The integrators write their element-wise work once,
// as lambdas over the array index (MARLICS_LAMBDA), and an executor runs it: this one with OpenMP
// loops, a CUDA one with kernels. Arrays allocated by the executor live where it computes.
//
// Reductions are deterministic: each thread sums a fixed contiguous block and the partial results
// are combined in thread order, so a run is reproducible for a given number of threads.
#include <omp.h>
#include <cstdlib>
#include <vector>
#include "driver.h"
#include "physics.h"
#include "geometry.h"

class CpuExec
{
 public:

  CpuExec(GEOMETRY * g, const Simulation_Parameters * sim_param=NULL) : geometry(g), n_sites(g->Nx*g->Ny*g->Nz), energy_buffer(NULL) {}
  ~CpuExec() { free(energy_buffer); }

  static const char * name() { return "CPU"; }

  double * alloc(long n)
  {
    double * p=(double *)calloc(n, sizeof(double));
    if(p==NULL){ ERROr }
    return p;
  }
  void release(double * p) { free(p); }

  //The integrators work on the executor's copy of Q: on the CPU it is the driver's array itself.
  double * attach(double * Qij) { return Qij; }
  void detach(double * Qij, const double * work) {}

  template<class F> void for_each(long n, F f) const
  {
    #pragma omp parallel for schedule(static)
    for(long i=0; i<n; i++) f(i);
  }

  //Sum of f(i). f may also update the arrays at index i: it is called exactly once per index.
  template<class F> double sum(long n, F f) const
  {
    std::vector<double> part(omp_get_max_threads(), 0.);
    #pragma omp parallel
    {
      const int t=omp_get_thread_num(), T=omp_get_num_threads();
      double s=0;
      for(long i=n*t/T; i<n*(t+1)/T; i++) s+=f(i);
      part[t]=s;
    }
    double s=0;
    for(double p : part) s+=p;
    return s;
  }

  template<class F> double max(long n, F f) const
  {
    std::vector<double> part(omp_get_max_threads(), 0.);
    #pragma omp parallel
    {
      const int t=omp_get_thread_num(), T=omp_get_num_threads();
      double m=0;
      for(long i=n*t/T; i<n*(t+1)/T; i++) { const double x=f(i); m= x>m ? x : m; }
      part[t]=m;
    }
    double m=0;
    for(double p : part) m= p>m ? p : m;
    return m;
  }

  //Two sums and a maximum in one pass: f(i, a, b, c) sets the three terms of index i.
  template<class F> SumSumMax sum_sum_max(long n, F f) const
  {
    const int nt=omp_get_max_threads();
    std::vector<SumSumMax> part(nt, SumSumMax{0.,0.,0.});
    #pragma omp parallel
    {
      const int t=omp_get_thread_num(), T=omp_get_num_threads();
      SumSumMax r={0.,0.,0.};
      for(long i=n*t/T; i<n*(t+1)/T; i++)
        {
          double a, b, c;
          f(i,a,b,c);
          r.sum0+=a; r.sum1+=b; r.max= c>r.max ? c : r.max;
        }
      part[t]=r;
    }
    SumSumMax r={0.,0.,0.};
    for(const SumSumMax & p : part) { r.sum0+=p.sum0; r.sum1+=p.sum1; r.max= p.max>r.max ? p.max : r.max; }
    return r;
  }

  void fill_ki(double * k_i, const double * Qij) const
  {
    #pragma omp parallel
    geometry->fill_ki(k_i,Qij);
  }

  //Total energy (10^-21 J).
  double energy(const double * Qij)
  {
    if(energy_buffer==NULL) energy_buffer=alloc(2*n_sites);
    double * e=energy_buffer;
    #pragma omp parallel
    geometry->Energy_calc(e,Qij);
    return sum(2*n_sites, MARLICS_LAMBDA (long i) { return e[i]; });
  }

 private:

  GEOMETRY * geometry;
  const long n_sites;
  double * energy_buffer;
};

#endif
