#ifndef FLOW_CUDA_H
#define FLOW_CUDA_H
// Backflow on the GPU: the per-site functions of flow_physics.h and flow_walls.h as kernels, and
// cuFFT for the Stokes solve. Same steps as FlowCpu (flow_cpu.h).
#include <cufft.h>
#include <iostream>
#include <cstdlib>
#include "driver.h"
#include "flow_physics.h"
#include "flow_walls.h"

#define CUFFT_CHECK(call) do { cufftResult r_=(call); if(r_!=CUFFT_SUCCESS) { \
  fprintf(stderr,"cuFFT error %d at %s:%d\n",(int)r_,__FILE__,__LINE__); exit(21);} } while(0)

namespace gpu
{
  __global__ void flow_stress(LatticeView L, FlowParams P, const double * Qij, const double * k_i, double * sigma, long n)
  {
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<n; s+=(long)blockDim.x*gridDim.x)
      if(L.point_type[s]==1) site_stress(L,P,Qij,k_i,s,&sigma[9*s]);
  }

  __global__ void flow_wall_stress(LatticeView L, WallFlowView W, double * sigma, long n)
  {
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<n; s+=(long)blockDim.x*gridDim.x)
      if(L.point_type[s]!=1) wall_stress(L,W,s,sigma);
  }

  __global__ void flow_force(LatticeView L, const double * sigma, cufftDoubleComplex * fk, long n)
  {
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<n; s+=(long)blockDim.x*gridDim.x)
      {
        double f[3]={0,0,0};
        if(L.point_type[s]==1) site_force(L,sigma,s,f);
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

  __global__ void flow_wall_modes(LatticeView L, WallFlowView W, double eta, cufftDoubleComplex * fk, double * scratch, long n, long n_plane)
  {
    double * const F[3]={(double *)fk,(double *)(fk+n),(double *)(fk+2*n)};
    for(long t=blockIdx.x*(long)blockDim.x+threadIdx.x; t<n_plane; t+=(long)blockDim.x*gridDim.x)
      wall_stokes_mode(L,W,eta,t,F,&scratch[5*W.M*t],n_plane);
  }

  __global__ void flow_wall_velocity(LatticeView L, WallFlowView W, const cufftDoubleComplex * fk, double * u, long n, long n_plane)
  {
    const long st=axis_step(L,W.axis);
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<n; s+=(long)blockDim.x*gridDim.x)
      {
        const int k=(int)((s/st)%W.M);
        for(int a=0; a<3; a++) u[3*s+a]= k==W.M-1 ? W.U[a] : (k==0 ? 0. : fk[a*n+s].x/n_plane);
      }
  }

  __global__ void flow_rhs(LatticeView L, FlowParams P, const double * Qij, const double * u, double * k_i, long n)
  {
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<n; s+=(long)blockDim.x*gridDim.x)
      {
        if(L.point_type[s]!=1) continue;
        double r[5];
        site_flow_rhs(L,P,Qij,u,s,r);
        for(int c=0; c<5; c++) k_i[5*s+c]+=r[c];
      }
  }
}

class FlowCuda
{
 public:

  //L: the lattice description with device pointers (CudaExec::view()); host: the same with host
  //pointers (GEOMETRY::view()).
  FlowCuda(const LatticeView & lattice, const LatticeView & host, const Simulation_Parameters * sim_param) : L(lattice), n((long)L.Nx*L.Ny*L.Nz),
    walls(NULL), Sinv(NULL), scratch(NULL), have_plan(true)
  {
    P.eta=sim_param->viscosity;
    P.xi=sim_param->xi;
    P.zeta=sim_param->activity;
    CUDA_CHECK(cudaMalloc(&sigma, 9*n*sizeof(double)));
    CUDA_CHECK(cudaMalloc(&u, 3*n*sizeof(double)));
    CUDA_CHECK(cudaMalloc(&fk, 3*n*sizeof(cufftDoubleComplex)));

    int axis;
    if(WallStokes::detect(host,axis))
      {
        walls=new WallStokes(host,axis,P.eta,sim_param->wall_velocity);
        W=walls->W;
        n_plane=n/W.M;
        CUDA_CHECK(cudaMalloc(&Sinv, walls->Sinv.size()*sizeof(double)));
        CUDA_CHECK(cudaMemcpy(Sinv, walls->Sinv.data(), walls->Sinv.size()*sizeof(double), cudaMemcpyHostToDevice));
        W.Sinv=Sinv;
        CUDA_CHECK(cudaMalloc(&scratch, 5*W.M*n_plane*sizeof(double)));

        //Transforms over the in-plane axes larger than 1 (fast: smaller stride), one batch entry per
        //row along the wall normal; one execution per component.
        int in[2], na=0;
        for(int a=0; a<3; a++) if(a!=axis && axis_size(L,a)>1) in[na++]=a;
        if(na==0) have_plan=false;
        else
          {
            const int f=in[0];
            int dims[2], embed[2], rank=0;
            if(na==2) { dims[rank]=axis_size(L,in[1]); embed[rank++]=axis_size(L,in[1]); }
            dims[rank]=axis_size(L,f);
            embed[rank++]= na==2 ? (int)(axis_step(L,in[1])/axis_step(L,f)) : axis_size(L,f);
            CUFFT_CHECK(cufftPlanMany(&plan, rank, dims, embed, (int)axis_step(L,f), (int)axis_step(L,axis),
                                      embed, (int)axis_step(L,f), (int)axis_step(L,axis), CUFFT_Z2Z, W.M));
          }
      }
    else
      {
        //Three transforms (one per component) of the dimensions larger than 1, row-major (z,y,x).
        int dims[3], rank=0;
        if(L.Nz>1) dims[rank++]=L.Nz;
        if(L.Ny>1) dims[rank++]=L.Ny;
        dims[rank++]=L.Nx;
        CUFFT_CHECK(cufftPlanMany(&plan, rank, dims, NULL, 1, (int)n, NULL, 1, (int)n, CUFFT_Z2Z, 3));
      }
    std::cout << "Backflow on the GPU: Stokes flow " << (walls ? "between walls" : "in a periodic box (cuFFT)") << ", viscosity=" << P.eta << " Pa s, xi=" << P.xi << ".\n";
    flow_alignment_report(P.xi,L.bulk.S_eq);
  }

  ~FlowCuda()
  {
    if(have_plan) cufftDestroy(plan);
    cudaFree(sigma); cudaFree(u); cudaFree(fk);
    cudaFree(Sinv); cudaFree(scratch);
    delete walls;
  }

  void add(double * k_i, const double * Qij)
  {
    const int b=gpu::blocks(n), t=gpu::threads;
    gpu::flow_stress<<<b,t>>>(L,P,Qij,k_i,sigma,n);
    if(walls) gpu::flow_wall_stress<<<b,t>>>(L,W,sigma,n);
    gpu::flow_force<<<b,t>>>(L,sigma,fk,n);
    if(walls)
      {
        for(int a=0; a<3 && have_plan; a++) CUFFT_CHECK(cufftExecZ2Z(plan,fk+a*n,fk+a*n,CUFFT_FORWARD));
        gpu::flow_wall_modes<<<gpu::blocks(n_plane),t>>>(L,W,P.eta,fk,scratch,n,n_plane);
        for(int a=0; a<3 && have_plan; a++) CUFFT_CHECK(cufftExecZ2Z(plan,fk+a*n,fk+a*n,CUFFT_INVERSE));
        gpu::flow_wall_velocity<<<b,t>>>(L,W,fk,u,n,n_plane);
      }
    else
      {
        CUFFT_CHECK(cufftExecZ2Z(plan,fk,fk,CUFFT_FORWARD));
        gpu::flow_stokes<<<b,t>>>(L,P.eta,fk,n);
        CUFFT_CHECK(cufftExecZ2Z(plan,fk,fk,CUFFT_INVERSE));
        gpu::flow_velocity<<<b,t>>>(fk,u,n);
      }
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
  bool have_plan;
  WallStokes * walls;    //host setup between walls (NULL in a periodic box)
  WallFlowView W;        //its device view
  long n_plane;
  double * Sinv;
  double * scratch;
};

#endif
