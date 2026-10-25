#ifndef FLOW_CUDA_H
#define FLOW_CUDA_H
// Backflow on the GPU: the per-site functions of flow_physics.h and flow_walls.h as kernels, and
// cuFFT for the Stokes solve. Same steps as FlowCpu (flow_cpu.h); in confined geometries the
// Lagrange-force solve of flow_confined.h runs on the host on vectors of the E2 sites (direct: the
// matrix-vector product on the device).
#include <cufft.h>
#include <iostream>
#include <cstdlib>
#include "driver.h"
#include "flow_physics.h"
#include "flow_walls.h"
#include "flow_confined.h"
#include <vector>

#define CUFFT_CHECK(call) do { cufftResult r_=(call); if(r_!=CUFFT_SUCCESS) { \
  fprintf(stderr,"cuFFT error %d at %s:%d\n",(int)r_,__FILE__,__LINE__); exit(21);} } while(0)

namespace gpu
{
  __global__ void flow_stress(LatticeView L, FlowParams P, const double * Qij, const double * k_i, double * sigma, long n)
  {
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<n; s+=(long)blockDim.x*gridDim.x)
      if(L.point_type[s]==1) site_stress(L,P,Qij,k_i,s,&sigma[9*s]);
  }

  __global__ void flow_boundary_stress(LatticeView L, double * sigma, long n)
  {
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<n; s+=(long)blockDim.x*gridDim.x)
      if(L.point_type[s]!=1) boundary_stress(L,s,sigma);
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

  __global__ void flow_stokes(LatticeView L, double eta, cufftDoubleComplex * fk, long n, bool keep_checkerboards=false)
  {
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<n; s+=(long)blockDim.x*gridDim.x)
      {
        const int m0=s%L.Nx, m1=(s/L.Nx)%L.Ny, m2=s/((long)L.Nx*L.Ny);
        double fr[3], fi[3];
        for(int a=0; a<3; a++) { fr[a]=fk[a*n+s].x; fi[a]=fk[a*n+s].y; }
        stokes_mode(L,eta,m0,m1,m2,fr,fi,keep_checkerboards);
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

  //Confined geometries. Real field (3 per site) into fk; force (if not NULL) keeps a copy.
  __global__ void flow_load(const double * f, cufftDoubleComplex * fk, long n)
  {
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<n; s+=(long)blockDim.x*gridDim.x)
      for(int a=0; a<3; a++) { fk[a*n+s].x=f[3*s+a]; fk[a*n+s].y=0; }
  }

  __global__ void flow_save_force(const cufftDoubleComplex * fk, double * force, long n)
  {
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<n; s+=(long)blockDim.x*gridDim.x)
      for(int a=0; a<3; a++) force[3*s+a]=fk[a*n+s].x;
  }

  //fk = 0, then p on the E2 sites (for E G E^T p).
  __global__ void flow_clear(cufftDoubleComplex * fk, long n)
  {
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<3*n; s+=(long)blockDim.x*gridDim.x) { fk[s].x=0; fk[s].y=0; }
  }

  __global__ void flow_add_e2(const long * e2, long ne, const double * p, cufftDoubleComplex * fk, long n)
  {
    for(long i=blockIdx.x*(long)blockDim.x+threadIdx.x; i<ne; i+=(long)blockDim.x*gridDim.x)
      for(int a=0; a<3; a++) fk[a*n+e2[i]].x+=p[3*i+a];
  }

  //Velocity on E2 after the inverse transform.
  __global__ void flow_gather_e2(const long * e2, long ne, const cufftDoubleComplex * fk, double * g, long n)
  {
    for(long i=blockIdx.x*(long)blockDim.x+threadIdx.x; i<ne; i+=(long)blockDim.x*gridDim.x)
      for(int a=0; a<3; a++) g[3*i+a]=fk[a*n+e2[i]].x/n;
  }

  __global__ void flow_confined_velocity(LatticeView L, const cufftDoubleComplex * fk, double U0x, double U0y, double U0z, double * u, long n)
  {
    const double U0[3]={U0x,U0y,U0z};
    for(long s=blockIdx.x*(long)blockDim.x+threadIdx.x; s<n; s+=(long)blockDim.x*gridDim.x)
      for(int a=0; a<3; a++) u[3*s+a]= L.point_type[s]==1 ? fk[a*n+s].x/n+U0[a] : 0.;
  }

  //out = C w, C dense m x m: one warp per row.
  __global__ void flow_matvec(const double * C, const double * w, double * out, long m)
  {
    const long warp=(blockIdx.x*(long)blockDim.x+threadIdx.x)/32, lane=threadIdx.x%32;
    const long warps=(long)gridDim.x*blockDim.x/32;
    for(long r=warp; r<m; r+=warps)
      {
        double v=0;
        for(long c=lane; c<m; c+=32) v+=C[r*m+c]*w[c];
        for(int o=16; o>0; o/=2) v+=__shfl_down_sync(0xffffffff,v,o);
        if(lane==0) out[r]=v;
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
    walls(NULL), Sinv(NULL), scratch(NULL), have_plan(true), confined(NULL), force(NULL), e2(NULL), vec(NULL), vec2(NULL), Cinv(NULL)
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
        const double * U=sim_param->wall_velocity;
        if(U[0]!=0 || U[1]!=0 || U[2]!=0)
          {
            std::cout << "wall_velocity needs two flat walls (geometries channel and slab).\nAborting the program.\n";
            exit(1);
          }
        if(ConfinedStokes::needed(host))
          {
            confined=new ConfinedStokes(host,P.eta,sim_param->flow_solver);
            const long ne=confined->ne, nc=3*ne;
            CUDA_CHECK(cudaMalloc(&force, 3*n*sizeof(double)));
            CUDA_CHECK(cudaMalloc(&e2, ne*sizeof(long)));
            CUDA_CHECK(cudaMemcpy(e2, confined->e2.data(), ne*sizeof(long), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMalloc(&vec, nc*sizeof(double)));
            CUDA_CHECK(cudaMalloc(&vec2, nc*sizeof(double)));
            if(confined->direct)
              {
                CUDA_CHECK(cudaMalloc(&Cinv, nc*nc*sizeof(double)));
                CUDA_CHECK(cudaMemcpy(Cinv, confined->Cinv.data(), nc*nc*sizeof(double), cudaMemcpyHostToDevice));
                std::vector<double>().swap(confined->Cinv);
              }
            h_g.resize(nc); h_lambda.resize(nc);
          }
        //Three transforms (one per component) of the dimensions larger than 1, row-major (z,y,x).
        int dims[3], rank=0;
        if(L.Nz>1) dims[rank++]=L.Nz;
        if(L.Ny>1) dims[rank++]=L.Ny;
        dims[rank++]=L.Nx;
        CUFFT_CHECK(cufftPlanMany(&plan, rank, dims, NULL, 1, (int)n, NULL, 1, (int)n, CUFFT_Z2Z, 3));
      }
    std::cout << "Backflow on the GPU: Stokes flow " << (walls ? "between walls" : (confined ? "in a confined geometry" : "in a periodic box (cuFFT)")) << ", viscosity=" << P.eta << " Pa s, xi=" << P.xi << ".\n";
    flow_alignment_report(P.xi,L.bulk.S_eq);
  }

  ~FlowCuda()
  {
    if(have_plan) cufftDestroy(plan);
    cudaFree(sigma); cudaFree(u); cudaFree(fk);
    cudaFree(Sinv); cudaFree(scratch);
    cudaFree(force); cudaFree(e2); cudaFree(vec); cudaFree(vec2); cudaFree(Cinv);
    delete walls;
    delete confined;
  }

  void add(double * k_i, const double * Qij)
  {
    const int b=gpu::blocks(n), t=gpu::threads;
    gpu::flow_stress<<<b,t>>>(L,P,Qij,k_i,sigma,n);
    if(walls || confined) gpu::flow_boundary_stress<<<b,t>>>(L,sigma,n);
    gpu::flow_force<<<b,t>>>(L,sigma,fk,n);
    solve();
    gpu::flow_rhs<<<b,t>>>(L,P,Qij,u,k_i,n);
    CUDA_CHECK(cudaGetLastError());
  }

  //Stokes flow for a force density on the device (3 per site), e.g. for tests; with surface sites
  //the force off the bulk must be zero.
  void solve_force(const double * f)
  {
    gpu::flow_load<<<gpu::blocks(n),gpu::threads>>>(f,fk,n);
    solve();
    CUDA_CHECK(cudaGetLastError());
  }

  //Velocity of the last call (device pointer, nm/us, 3 per site).
  const double * velocity() const { return u; }

 private:

  //Velocity u from the force held in fk.
  void solve()
  {
    const int b=gpu::blocks(n), t=gpu::threads;
    if(confined) { solve_confined(); return; }
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
  }

  //Periodic Stokes solve of fk in place (checkerboards kept); total (if given): the sum of the force.
  void periodic_solve(double * total=NULL)
  {
    CUFFT_CHECK(cufftExecZ2Z(plan,fk,fk,CUFFT_FORWARD));
    if(total)
      for(int a=0; a<3; a++)
        {
          cufftDoubleComplex z;
          CUDA_CHECK(cudaMemcpy(&z, fk+a*n, sizeof(z), cudaMemcpyDeviceToHost));
          total[a]=z.x;
        }
    gpu::flow_stokes<<<gpu::blocks(n),gpu::threads>>>(L,P.eta,fk,n,true);
    CUFFT_CHECK(cufftExecZ2Z(plan,fk,fk,CUFFT_INVERSE));
  }

  //E2 velocities of fk into the host vector out.
  void gather(double * out)
  {
    const long ne=confined->ne;
    gpu::flow_gather_e2<<<gpu::blocks(ne),gpu::threads>>>(e2,ne,fk,vec,n);
    CUDA_CHECK(cudaMemcpy(out, vec, 3*ne*sizeof(double), cudaMemcpyDeviceToHost));
  }

  //Confined geometry: the steps of FlowCpu::solve_confined.
  void solve_confined()
  {
    const int b=gpu::blocks(n), t=gpu::threads;
    const long ne=confined->ne, nc=3*ne;
    gpu::flow_save_force<<<b,t>>>(fk,force,n);
    double total[3];
    periodic_solve(total);
    gather(h_g.data());

    const ConfinedStokes::Apply Cmul=[&](const double * p, double * out){
        CUDA_CHECK(cudaMemcpy(vec2, p, nc*sizeof(double), cudaMemcpyHostToDevice));
        gpu::flow_clear<<<gpu::blocks(3*n),t>>>(fk,n);
        gpu::flow_add_e2<<<gpu::blocks(ne),t>>>(e2,ne,vec2,fk,n);
        periodic_solve();
        gather(out);
      };
    const ConfinedStokes::Apply Dmul=[&](const double * w, double * out){
        CUDA_CHECK(cudaMemcpy(vec2, w, nc*sizeof(double), cudaMemcpyHostToDevice));
        gpu::flow_matvec<<<gpu::blocks(32*nc),t>>>(Cinv,vec2,vec,nc);
        CUDA_CHECK(cudaMemcpy(out, vec, nc*sizeof(double), cudaMemcpyDeviceToHost));
      };
    confined->solve(h_g.data(),total,h_lambda.data(),Cmul,confined->direct ? Dmul : ConfinedStokes::Apply());

    CUDA_CHECK(cudaMemcpy(vec2, h_lambda.data(), nc*sizeof(double), cudaMemcpyHostToDevice));
    gpu::flow_load<<<b,t>>>(force,fk,n);
    gpu::flow_add_e2<<<gpu::blocks(ne),t>>>(e2,ne,vec2,fk,n);
    periodic_solve();
    gather(h_g.data());
    double U0[3]={0,0,0};
    for(long i=0; i<ne; i++)
      for(int a=0; a<3; a++) U0[a]-=h_g[3*i+a];
    for(int a=0; a<3; a++) U0[a]/=ne;
    gpu::flow_confined_velocity<<<b,t>>>(L,fk,U0[0],U0[1],U0[2],u,n);
  }

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
  ConfinedStokes * confined;   //host setup in confined geometries (NULL otherwise)
  double * force;              //device: force (3 per site)
  long * e2;                   //device: E2 sites
  double * vec, * vec2;        //device: vectors of the E2 sites
  double * Cinv;               //device: direct inverse
  std::vector<double> h_g, h_lambda;
};

#endif
