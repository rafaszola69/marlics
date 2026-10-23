#ifndef EXEC_CUDA_H
#define EXEC_CUDA_H
// CUDA executor for the integrators: runs their element-wise lambdas (MARLICS_LAMBDA) as kernels
// on arrays that live on the GPU. Same interface as CpuExec (exec_cpu.h).
//
// Reductions are deterministic: a fixed grid of blocks, a fixed tree inside each block, and the
// block results added on the host in block order.
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cuda_runtime.h>
#include "driver.h"
#include "physics.h"
#include "lattice.h"
#include "geometry.h"

#define MARLICS_LAMBDA [=] MARLICS_HD

#define CUDA_CHECK(call) do { cudaError_t err_=(call); if(err_!=cudaSuccess) { \
  fprintf(stderr,"CUDA error %s at %s:%d\n",cudaGetErrorString(err_),__FILE__,__LINE__); exit(20);} } while(0)

struct SumSumMax { double sum0, sum1, max; };

namespace gpu
{
  const int threads=256;
  const int reduce_blocks=512;

  template<class F> __global__ void for_each(long n, F f)
  {
    for(long i=blockIdx.x*(long)blockDim.x+threadIdx.x; i<n; i+=(long)blockDim.x*gridDim.x) f(i);
  }

  //Per-block sums (MAX=false) or maxima (MAX=true) of f(i), written to out[blockIdx.x].
  template<bool MAX, class F> __global__ void reduce(long n, F f, double * out)
  {
    __shared__ double s[threads];
    double acc=0;
    for(long i=blockIdx.x*(long)blockDim.x+threadIdx.x; i<n; i+=(long)blockDim.x*gridDim.x)
      {
        const double x=f(i);
        acc= MAX ? (x>acc ? x : acc) : acc+x;
      }
    s[threadIdx.x]=acc;
    __syncthreads();
    for(int w=threads/2; w>0; w/=2)
      {
        if(threadIdx.x<w) s[threadIdx.x]= MAX ? fmax(s[threadIdx.x],s[threadIdx.x+w]) : s[threadIdx.x]+s[threadIdx.x+w];
        __syncthreads();
      }
    if(threadIdx.x==0) out[blockIdx.x]=s[0];
  }

  template<class F> __global__ void reduce_sum_sum_max(long n, F f, SumSumMax * out)
  {
    __shared__ double s0[threads], s1[threads], m[threads];
    double a0=0, a1=0, am=0;
    for(long i=blockIdx.x*(long)blockDim.x+threadIdx.x; i<n; i+=(long)blockDim.x*gridDim.x)
      {
        double a, b, c;
        f(i,a,b,c);
        a0+=a; a1+=b; am= c>am ? c : am;
      }
    s0[threadIdx.x]=a0; s1[threadIdx.x]=a1; m[threadIdx.x]=am;
    __syncthreads();
    for(int w=threads/2; w>0; w/=2)
      {
        if(threadIdx.x<w)
          {
            s0[threadIdx.x]+=s0[threadIdx.x+w];
            s1[threadIdx.x]+=s1[threadIdx.x+w];
            m[threadIdx.x]=fmax(m[threadIdx.x],m[threadIdx.x+w]);
          }
        __syncthreads();
      }
    if(threadIdx.x==0) { out[blockIdx.x].sum0=s0[0]; out[blockIdx.x].sum1=s1[0]; out[blockIdx.x].max=m[0]; }
  }

  __global__ void bulk_forces(LatticeView L, const double * Qij, double * k_i)
  {
    const long n=(long)L.Nx*L.Ny*L.Nz;
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<n; s+=(long)blockDim.x*gridDim.x)
      {
        if(L.point_type[s]==1) bulk_site_force(L,Qij,s,&k_i[5*s]);
        else for(int ll=0; ll<=4; ll++) k_i[5*s+ll]=0.;
      }
  }

  __global__ void surface_forces(LatticeView L, const double * Qij, double * k_i)
  {
    for(int n=blockIdx.x*blockDim.x+threadIdx.x; n<L.n_surface; n+=blockDim.x*gridDim.x)
      surface_site_force(L,Qij,L.surface[n],&k_i[5*L.surface[n].site]);
  }

  __global__ void bulk_energies(LatticeView L, const double * Qij, double * e)
  {
    const long n=(long)L.Nx*L.Ny*L.Nz;
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<n; s+=(long)blockDim.x*gridDim.x)
      {
        e[2*s+0]= L.point_type[s]==0 ? 0. : L.dV*bulk_site_energy(L,Qij,s);
        e[2*s+1]=0.;
      }
  }

  __global__ void surface_energies(LatticeView L, const double * Qij, double * e)
  {
    for(int n=blockIdx.x*blockDim.x+threadIdx.x; n<L.n_surface; n+=blockDim.x*gridDim.x)
      e[2*L.surface[n].site+1]=surface_site_energy(L,Qij,L.surface[n]);
  }

  inline int blocks(long n) { const long b=(n+threads-1)/threads; return b<65535 ? (int)b : 65535; }
}

class CudaExec
{
 public:

  //Copies the lattice description (point types, surface sites, anchoring parameters) to the GPU.
  CudaExec(GEOMETRY * g) : n_sites((long)g->Nx*g->Ny*g->Nz), Q_dev(NULL), energy_buffer(NULL)
  {
    const LatticeView host=g->view();
    L=host;

    int * pt;
    CUDA_CHECK(cudaMalloc(&pt, n_sites*sizeof(int)));
    CUDA_CHECK(cudaMemcpy(pt, host.point_type, n_sites*sizeof(int), cudaMemcpyHostToDevice));
    L.point_type=pt;

    SurfaceSite * surf=NULL;
    if(host.n_surface>0)
      {
        CUDA_CHECK(cudaMalloc(&surf, host.n_surface*sizeof(SurfaceSite)));
        CUDA_CHECK(cudaMemcpy(surf, host.surface, host.n_surface*sizeof(SurfaceSite), cudaMemcpyHostToDevice));
      }
    L.surface=surf;

    BoundaryParams * bc=NULL;
    if(host.n_bc>0)
      {
        CUDA_CHECK(cudaMalloc(&bc, host.n_bc*sizeof(BoundaryParams)));
        CUDA_CHECK(cudaMemcpy(bc, host.bc, host.n_bc*sizeof(BoundaryParams), cudaMemcpyHostToDevice));
      }
    L.bc=bc;

    CUDA_CHECK(cudaMalloc(&partial, gpu::reduce_blocks*sizeof(SumSumMax)));

    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop,0));
    std::cout << "GPU: " << prop.name << ", " << prop.totalGlobalMem/1048576 << " MiB.\n";
  }

  ~CudaExec()
  {
    cudaFree((void *)L.point_type);
    cudaFree((void *)L.surface);
    cudaFree((void *)L.bc);
    cudaFree(partial);
    cudaFree(Q_dev);
    cudaFree(energy_buffer);
  }

  static const char * name() { return "GPU"; }

  double * alloc(long n)
  {
    double * p;
    if(cudaMalloc(&p, n*sizeof(double))!=cudaSuccess)
      {
        std::cout << "Not enough GPU memory for " << n*sizeof(double)/1048576 << " MiB more.\nAborting the program.\n";
        exit(4);
      }
    CUDA_CHECK(cudaMemset(p, 0, n*sizeof(double)));
    return p;
  }
  void release(double * p) { cudaFree(p); }

  //The integrators work on a copy of Q on the GPU, sent at the start of evolve and brought back
  //at its end (for the snapshot).
  double * attach(double * Qij)
  {
    if(Q_dev==NULL) Q_dev=alloc(5*n_sites);
    CUDA_CHECK(cudaMemcpy(Q_dev, Qij, 5*n_sites*sizeof(double), cudaMemcpyHostToDevice));
    return Q_dev;
  }
  void detach(double * Qij, const double * work)
  {
    CUDA_CHECK(cudaMemcpy(Qij, work, 5*n_sites*sizeof(double), cudaMemcpyDeviceToHost));
  }

  template<class F> void for_each(long n, F f) const
  {
    gpu::for_each<<<gpu::blocks(n),gpu::threads>>>(n,f);
    CUDA_CHECK(cudaGetLastError());
  }

  template<class F> double sum(long n, F f) const { return reduce<false>(n,f); }
  template<class F> double max(long n, F f) const { return reduce<true>(n,f); }

  template<class F> SumSumMax sum_sum_max(long n, F f) const
  {
    SumSumMax part[gpu::reduce_blocks];
    gpu::reduce_sum_sum_max<<<gpu::reduce_blocks,gpu::threads>>>(n,f,(SumSumMax *)partial);
    CUDA_CHECK(cudaMemcpy(part, partial, sizeof(part), cudaMemcpyDeviceToHost));
    SumSumMax r={0.,0.,0.};
    for(const SumSumMax & p : part) { r.sum0+=p.sum0; r.sum1+=p.sum1; r.max= p.max>r.max ? p.max : r.max; }
    return r;
  }

  void fill_ki(double * k_i, const double * Qij) const
  {
    gpu::bulk_forces<<<gpu::blocks(n_sites),gpu::threads>>>(L,Qij,k_i);
    if(L.n_surface>0) gpu::surface_forces<<<gpu::blocks(L.n_surface),gpu::threads>>>(L,Qij,k_i);
    CUDA_CHECK(cudaGetLastError());
  }

  //Total energy (10^-21 J).
  double energy(const double * Qij)
  {
    if(energy_buffer==NULL) energy_buffer=alloc(2*n_sites);
    double * e=energy_buffer;
    gpu::bulk_energies<<<gpu::blocks(n_sites),gpu::threads>>>(L,Qij,e);
    if(L.n_surface>0) gpu::surface_energies<<<gpu::blocks(L.n_surface),gpu::threads>>>(L,Qij,e);
    CUDA_CHECK(cudaGetLastError());
    return sum(2*n_sites, MARLICS_LAMBDA (long i) { return e[i]; });
  }

  //Device copy of the lattice description (tests).
  const LatticeView & view() const { return L; }

 private:

  template<bool MAX, class F> double reduce(long n, F f) const
  {
    double part[gpu::reduce_blocks];
    gpu::reduce<MAX><<<gpu::reduce_blocks,gpu::threads>>>(n,f,(double *)partial);
    CUDA_CHECK(cudaMemcpy(part, partial, sizeof(part), cudaMemcpyDeviceToHost));
    double r=0;
    for(double p : part) r= MAX ? (p>r ? p : r) : r+p;
    return r;
  }

  const long n_sites;
  LatticeView L;
  void * partial;
  double * Q_dev;
  double * energy_buffer;
};

#endif
