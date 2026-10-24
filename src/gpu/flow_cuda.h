#ifndef FLOW_CUDA_H
#define FLOW_CUDA_H
// Backflow on the GPU: the per-site functions of flow_physics.h as kernels, and cuFFT for the Stokes
// solve. Same steps as FlowCpu (flow_cpu.h).
#include <cufft.h>
#include <iostream>
#include <cstdlib>
#include "driver.h"
#include "flow_physics.h"

#define CUFFT_CHECK(call) do { cufftResult r_=(call); if(r_!=CUFFT_SUCCESS) { \
  fprintf(stderr,"cuFFT error %d at %s:%d\n",(int)r_,__FILE__,__LINE__); exit(21);} } while(0)

namespace gpu
{
  __global__ void flow_stress(LatticeView L, FlowParams P, const double * Qij, const double * k_i, double * sigma, long n)
  {
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<n; s+=(long)blockDim.x*gridDim.x)
      site_stress(L,P,Qij,k_i,s,&sigma[9*s]);
  }

  __global__ void flow_force(LatticeView L, const double * sigma, cufftDoubleComplex * fk, long n)
  {
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<n; s+=(long)blockDim.x*gridDim.x)
      {
        double f[3];
        site_force(L,sigma,s,f);
        for(int a=0; a<3; a++) { fk[a*n+s].x=f[a]; fk[a*n+s].y=0; }
      }
  }

  __global__ void flow_stokes(LatticeView L, double eta, cufftDoubleComplex * fk, long n)
  {
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<n; s+=(long)blockDim.x*gridDim.x)
      {
        const int m0=s%L.Nx, m1=(s/L.Nx)%L.Ny, m2=s/((long)L.Nx*L.Ny);
        double fr[3], fi[3];
        for(int a=0; a<3; a++) { fr[a]=fk[a*n+s].x; fi[a]=fk[a*n+s].y; }
        stokes_mode(L,eta,m0,m1,m2,fr,fi);
        for(int a=0; a<3; a++) { fk[a*n+s].x=fr[a]; fk[a*n+s].y=fi[a]; }
      }
  }

  __global__ void flow_velocity(const cufftDoubleComplex * fk, double * u, long n)
  {
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<n; s+=(long)blockDim.x*gridDim.x)
      for(int a=0; a<3; a++) u[3*s+a]=fk[a*n+s].x/n;
  }

  __global__ void flow_rhs(LatticeView L, FlowParams P, const double * Qij, const double * u, double * k_i, long n)
  {
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<n; s+=(long)blockDim.x*gridDim.x)
      {
        double r[5];
        site_flow_rhs(L,P,Qij,u,s,r);
        for(int c=0; c<5; c++) k_i[5*s+c]+=r[c];
      }
  }
}

class FlowCuda
{
 public:

  //L: the lattice description with device pointers (CudaExec::view()).
  FlowCuda(const LatticeView & lattice, const Simulation_Parameters * sim_param) : L(lattice), n((long)L.Nx*L.Ny*L.Nz)
  {
    P.eta=sim_param->viscosity;
    P.xi=sim_param->xi;
    if(!L.periodic[0] || !L.periodic[1] || !L.periodic[2] || L.n_surface>0)
      {
        std::cout << "The flow (backflow) is only available in periodic domains (geometry bulk).\nAborting the program.\n";
        exit(1);
      }
    CUDA_CHECK(cudaMalloc(&sigma, 9*n*sizeof(double)));
    CUDA_CHECK(cudaMalloc(&u, 3*n*sizeof(double)));
    CUDA_CHECK(cudaMalloc(&fk, 3*n*sizeof(cufftDoubleComplex)));

    //Three transforms (one per component) of the dimensions larger than 1, row-major (z,y,x).
    int dims[3], rank=0;
    if(L.Nz>1) dims[rank++]=L.Nz;
    if(L.Ny>1) dims[rank++]=L.Ny;
    dims[rank++]=L.Nx;
    CUFFT_CHECK(cufftPlanMany(&plan, rank, dims, NULL, 1, (int)n, NULL, 1, (int)n, CUFFT_Z2Z, 3));
    std::cout << "Backflow on the GPU: Stokes flow in a periodic box (cuFFT), viscosity=" << P.eta << " Pa s, xi=" << P.xi << ".\n";
    flow_alignment_report(P.xi,L.bulk.S_eq);
  }

  ~FlowCuda()
  {
    cufftDestroy(plan);
    cudaFree(sigma); cudaFree(u); cudaFree(fk);
  }

  void add(double * k_i, const double * Qij)
  {
    const int b=gpu::blocks(n), t=gpu::threads;
    gpu::flow_stress<<<b,t>>>(L,P,Qij,k_i,sigma,n);
    gpu::flow_force<<<b,t>>>(L,sigma,fk,n);
    CUFFT_CHECK(cufftExecZ2Z(plan,fk,fk,CUFFT_FORWARD));
    gpu::flow_stokes<<<b,t>>>(L,P.eta,fk,n);
    CUFFT_CHECK(cufftExecZ2Z(plan,fk,fk,CUFFT_INVERSE));
    gpu::flow_velocity<<<b,t>>>(fk,u,n);
    gpu::flow_rhs<<<b,t>>>(L,P,Qij,u,k_i,n);
    CUDA_CHECK(cudaGetLastError());
  }

  //Velocity of the last call (device pointer, nm/us, 3 per site).
  const double * velocity() const { return u; }

 private:

  const LatticeView L;
  FlowParams P;
  const long n;
  double * sigma;
  double * u;
  cufftDoubleComplex * fk;
  cufftHandle plan;
};

#endif
