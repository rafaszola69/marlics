#ifndef LATTICE_H
#define LATTICE_H
// Per-site right-hand side (fill_ki) and energy, written once for the CPU (OpenMP loops in
// GEOMETRY) and the GPU (CUDA kernels). A geometry only describes its lattice: the point types,
// which axes are periodic and, for each surface site, its normal, area, anchoring condition and
// the neighbours used for its first derivatives.
#include "physics.h"
#include "physics_ldg.h"
#include "physics_rp.h"
#include "physics_fg.h"
#include "physics_homeotropic.h"
#include "physics_ldg_gradient.h"

//A surface site and its stencil: dQ along axis a is w[a]*(Q[np[a]]-Q[nm[a]])/d_a, with w=0.5 for
//a central difference and 1 for a one-sided one.
struct SurfaceSite
{
  int site;
  int bc;
  double v[3];
  double area;
  int nm[3], np[3];
  double w[3];
};

//Everything the per-site functions need. On the GPU the pointers are device pointers.
struct LatticeView
{
  int Nx, Ny, Nz;
  int periodic[3];
  double d_1[3];      //1/dx, 1/dy, 1/dz
  double dV;
  const int * point_type;
  int n_surface;
  const SurfaceSite * surface;
  BulkParams bulk;
  int n_bc;
  const BoundaryParams * bc;
};

template<class T> MARLICS_HD inline void boundary_derivative(const BulkParams& P, const BoundaryParams& B, const T QN[5], const T dQ[15], const T v[3], double out[5])
{
  switch(B.kind)
    {
    case BC_RP:
      out[0]=rp_functional_derivative_00(P,B,QN,dQ,v); out[1]=rp_functional_derivative_01(P,B,QN,dQ,v);
      out[2]=rp_functional_derivative_02(P,B,QN,dQ,v); out[3]=rp_functional_derivative_11(P,B,QN,dQ,v);
      out[4]=rp_functional_derivative_12(P,B,QN,dQ,v);
      break;
    case BC_FG:
      out[0]=fg_functional_derivative_00(P,B,QN,dQ,v); out[1]=fg_functional_derivative_01(P,B,QN,dQ,v);
      out[2]=fg_functional_derivative_02(P,B,QN,dQ,v); out[3]=fg_functional_derivative_11(P,B,QN,dQ,v);
      out[4]=fg_functional_derivative_12(P,B,QN,dQ,v);
      break;
    case BC_HOMEOTROPIC:
      out[0]=homeotropic_functional_derivative_00(P,B,QN,dQ,v); out[1]=homeotropic_functional_derivative_01(P,B,QN,dQ,v);
      out[2]=homeotropic_functional_derivative_02(P,B,QN,dQ,v); out[3]=homeotropic_functional_derivative_11(P,B,QN,dQ,v);
      out[4]=homeotropic_functional_derivative_12(P,B,QN,dQ,v);
      break;
    default:  //strong anchoring: surface sites stay fixed
      for(int ll=0; ll<=4; ll++) out[ll]=0.;
    }
}

MARLICS_HD inline double boundary_energy(const BulkParams& P, const BoundaryParams& B, const double QN[5], const double dQ[15], const double v[3])
{
  switch(B.kind)
    {
    case BC_RP:          return rp_energy_calculation(P,B,QN,dQ,v);
    case BC_FG:          return fg_energy_calculation(P,B,QN,dQ,v);
    case BC_HOMEOTROPIC: return homeotropic_energy_calculation(P,B,QN,dQ,v);
    default:             return 0.;
    }
}

MARLICS_HD inline int lattice_index(const LatticeView& L, int i, int j, int k) { return (k*L.Ny+j)*L.Nx+i; }

//Copy of a double array in the evaluation type Real (no-op for double).
template<class Real, int n> struct RealArray
{
  Real x[n];
  MARLICS_HD RealArray(const double * a) { for(int i=0; i<n; i++) x[i]=Real(a[i]); }
};

//Bulk site (point_type 1): central differences, wrapping around periodic axes. The differences are
//taken in double; Real is the type the functional derivative is evaluated in (float for the GPU
//mixed precision).
template<class Real=double>
MARLICS_HD inline void bulk_site_force(const LatticeView& L, const double * Qij, int s, double out[5])
{
  const int Nx=L.Nx, Ny=L.Ny, Nz=L.Nz;
  const int i=s%Nx, j=(s/Nx)%Ny, k=s/(Nx*Ny);
  const int ip1= L.periodic[0] ? (i+1)%Nx : i+1, im1= L.periodic[0] ? (i+Nx-1)%Nx : i-1;
  const int jp1= L.periodic[1] ? (j+1)%Ny : j+1, jm1= L.periodic[1] ? (j+Ny-1)%Ny : j-1;
  const int kp1= L.periodic[2] ? (k+1)%Nz : k+1, km1= L.periodic[2] ? (k+Nz-1)%Nz : k-1;
  const double dx_1=L.d_1[0], dy_1=L.d_1[1], dz_1=L.d_1[2];
  double QN[5], dQ[15], ddQ[30];
  const double v[3]={0.,0.,0.};

  #define Q_AT(ii,jj,kk,ll) Qij[5*(((kk)*Ny+(jj))*Nx+(ii))+(ll)]
  for(int ll=0; ll<=4;ll++)
    {
      QN[ll]=Qij[5*s+ll];

      dQ[ll]   = 0.5*(Q_AT(ip1,j,k,ll)-Q_AT(im1,j,k,ll))*dx_1;
      dQ[5+ll] = 0.5*(Q_AT(i,jp1,k,ll)-Q_AT(i,jm1,k,ll))*dy_1;
      dQ[10+ll]= 0.5*(Q_AT(i,j,kp1,ll)-Q_AT(i,j,km1,ll))*dz_1;

      ddQ[ll]   = (Q_AT(ip1,j,k,ll)-2.0*QN[ll]+Q_AT(im1,j,k,ll))*dx_1*dx_1;
      ddQ[5+ll] = 0.25*(Q_AT(ip1,jp1,k,ll)+Q_AT(im1,jm1,k,ll)-Q_AT(im1,jp1,k,ll)-Q_AT(ip1,jm1,k,ll))*dx_1*dy_1;
      ddQ[10+ll]= 0.25*(Q_AT(ip1,j,kp1,ll)+Q_AT(im1,j,km1,ll)-Q_AT(im1,j,kp1,ll)-Q_AT(ip1,j,km1,ll))*dx_1*dz_1;
      ddQ[15+ll]= (Q_AT(i,jp1,k,ll)-2.0*QN[ll]+Q_AT(i,jm1,k,ll))*dy_1*dy_1;
      ddQ[20+ll]= 0.25*(Q_AT(i,jp1,kp1,ll)+Q_AT(i,jm1,km1,ll)-Q_AT(i,jm1,kp1,ll)-Q_AT(i,jp1,km1,ll))*dz_1*dy_1;
      ddQ[25+ll]= (Q_AT(i,j,kp1,ll)-2.0*QN[ll]+Q_AT(i,j,km1,ll))*dz_1*dz_1;
    }
  #undef Q_AT

  const RealArray<Real,5> q(QN);
  const RealArray<Real,15> d(dQ);
  const RealArray<Real,30> dd(ddQ);
  const RealArray<Real,3> n(v);
  out[0]=ldg_functional_derivative_00<Real>(L.bulk,q.x,d.x,dd.x,n.x);
  out[1]=ldg_functional_derivative_01<Real>(L.bulk,q.x,d.x,dd.x,n.x);
  out[2]=ldg_functional_derivative_02<Real>(L.bulk,q.x,d.x,dd.x,n.x);
  out[3]=ldg_functional_derivative_11<Real>(L.bulk,q.x,d.x,dd.x,n.x);
  out[4]=ldg_functional_derivative_12<Real>(L.bulk,q.x,d.x,dd.x,n.x);

  //Plane confinement f=W*Q_zz=-W*(Q_00+Q_11).
  if(L.bulk.plane_W!=0)
    {
      out[0]+=L.bulk.Lambda*L.bulk.plane_W/3.;
      out[3]+=L.bulk.Lambda*L.bulk.plane_W/3.;
    }
}

//Neighbour of site s along axis in direction dir (+1/-1). Outside points and non-periodic domain
//edges return the site itself (zero one-sided difference).
MARLICS_HD inline int lattice_neighbor(const LatticeView& L, int s, int axis, int dir)
{
  int c[3]={s%L.Nx, (s/L.Nx)%L.Ny, s/(L.Nx*L.Ny)};
  const int n[3]={L.Nx,L.Ny,L.Nz};

  c[axis]+=dir;
  if(c[axis]<0 || c[axis]>=n[axis])
    {
      if(!L.periodic[axis]) return s;
      c[axis]=(c[axis]+n[axis])%n[axis];
    }
  const int t=lattice_index(L,c[0],c[1],c[2]);
  return L.point_type[t]==0 ? s : t;
}

//Bulk energy density at site s, averaged over the 8 octants of one-sided differences.
//Its exact gradient is the compact stencil of bulk_site_force (3-point second derivatives and
//4-corner mixed derivatives), so the time integrators decrease it, unlike a central-difference energy.
MARLICS_HD inline double bulk_site_energy(const LatticeView& L, const double * Qij, int s)
{
  const int nx[2]={lattice_neighbor(L,s,0,-1),lattice_neighbor(L,s,0,1)};
  const int ny[2]={lattice_neighbor(L,s,1,-1),lattice_neighbor(L,s,1,1)};
  const int nz[2]={lattice_neighbor(L,s,2,-1),lattice_neighbor(L,s,2,1)};
  const double sg[2]={-1.,1.};
  const double v[3]={0.,0.,0.};
  double QN[5], dQ[15], f=0;

  for(int ll=0; ll<=4;ll++) QN[ll]=Qij[5*s+ll];
  for(int a=0; a<2; a++)
    for(int b=0; b<2; b++)
      for(int c=0; c<2; c++)
        {
          for(int ll=0; ll<=4;ll++)
            {
              dQ[ll]   =sg[a]*(Qij[5*nx[a]+ll]-QN[ll])*L.d_1[0];
              dQ[5+ll] =sg[b]*(Qij[5*ny[b]+ll]-QN[ll])*L.d_1[1];
              dQ[10+ll]=sg[c]*(Qij[5*nz[c]+ll]-QN[ll])*L.d_1[2];
            }
          f+=ldg_energy_calculation(L.bulk,QN,dQ,v);
        }
  return 0.125*f-L.bulk.plane_W*(QN[0]+QN[3]);
}

//Gradient of the bulk part of the discrete energy, sum over sites of dV*bulk_site_energy, with
//respect to the independent components of Q at site t, divided by dV: the octants of t and those of
//its neighbours that reach t.
MARLICS_HD inline void bulk_energy_gradient(const LatticeView& L, const double * Qij, int t, double g[5])
{
  const double sg[2]={-1.,1.};
  const double v[3]={0.,0.,0.};
  double QN[5], dQ[15], fQ[5], fD[15];

  for(int ll=0; ll<=4; ll++) g[ll]=0.;
  for(int s_axis=-1; s_axis<3; s_axis++)
    for(int s_dir=0; s_dir<2; s_dir++)
      {
        //s_axis=-1: site t itself (once); otherwise its neighbour along s_axis, whose octants pointing
        //back to t depend on Q_t through d_{s_axis}Q.
        if(s_axis<0 && s_dir==1) continue;
        const int s= s_axis<0 ? t : lattice_neighbor(L,t,s_axis,2*s_dir-1);
        if(s_axis>=0 && (s==t || lattice_neighbor(L,s,s_axis,1-2*s_dir)!=t)) continue;
        int nb[3][2];
        for(int a=0; a<3; a++) for(int d=0; d<2; d++) nb[a][d]=lattice_neighbor(L,s,a,2*d-1);
        for(int ll=0; ll<=4; ll++) QN[ll]=Qij[5*s+ll];
        for(int o=0; o<8; o++)
          {
            const int oc[3]={o&1, (o>>1)&1, (o>>2)&1};
            if(s_axis>=0 && oc[s_axis]!=1-s_dir) continue;
            for(int a=0; a<3; a++)
              for(int ll=0; ll<=4; ll++) dQ[5*a+ll]=sg[oc[a]]*(Qij[5*nb[a][oc[a]]+ll]-QN[ll])*L.d_1[a];
            ldg_gradient_derivative(L.bulk,QN,dQ,fD);
            if(s_axis<0)
              {
                ldg_value_derivative(L.bulk,QN,dQ,fQ);
                for(int ll=0; ll<=4; ll++) g[ll]+=0.125*fQ[ll];
                for(int a=0; a<3; a++)
                  if(nb[a][oc[a]]!=t)
                    for(int ll=0; ll<=4; ll++) g[ll]-=0.125*sg[oc[a]]*L.d_1[a]*fD[5*a+ll];
              }
            else
              for(int ll=0; ll<=4; ll++) g[ll]+=0.125*sg[oc[s_axis]]*L.d_1[s_axis]*fD[5*s_axis+ll];
          }
      }
  g[0]-=L.bulk.plane_W;
  g[3]-=L.bulk.plane_W;
}

//Surface site. Local dynamics: first derivatives from its stencil, then the natural boundary
//condition (the published Marlics scheme; its elastic part, a one-sided normal derivative, is not the
//gradient of the discrete energy, which can then rise for a while, mostly with degenerate anchoring).
//Variational dynamics: the exact gradient of the discrete energy, the bulk energy the site carries
//and its anchoring energy, per unit area: L nu.grad Q + df_s/dQ up to O(h) on a flat wall.
template<class Real=double>
MARLICS_HD inline void surface_site_force(const LatticeView& L, const double * Qij, const SurfaceSite& S, double out[5])
{
  double QN[5], dQ[15];

  if(L.bulk.surface_variational && L.bc[S.bc].kind!=BC_STRONG)
    {
      BulkParams P=L.bulk;
      P.L1=P.L2=P.L3=P.Ls=P.Lq_tilde=0.;
      for(int ll=0; ll<=4; ll++) QN[ll]=Qij[5*S.site+ll];
      for(int m=0; m<15; m++) dQ[m]=0.;
      boundary_derivative<double>(P,L.bc[S.bc],QN,dQ,S.v,out);
      double g[5];
      bulk_energy_gradient(L,Qij,S.site,g);
      const double c=L.bulk.Lambda_s*L.dV/S.area;
      out[0]-=c*(2.*g[0]-g[3])/3.;
      out[3]-=c*(2.*g[3]-g[0])/3.;
      out[1]-=0.5*c*g[1];
      out[2]-=0.5*c*g[2];
      out[4]-=0.5*c*g[4];
      return;
    }

  for(int ll=0; ll<=4;ll++)
    {
      QN[ll]=Qij[5*S.site+ll];
      for(int a=0; a<3; a++) dQ[5*a+ll]=S.w[a]*(Qij[5*S.np[a]+ll]-Qij[5*S.nm[a]+ll])*L.d_1[a];
    }
  const RealArray<Real,5> q(QN);
  const RealArray<Real,15> d(dQ);
  const RealArray<Real,3> n(S.v);
  boundary_derivative<Real>(L.bulk,L.bc[S.bc],q.x,d.x,n.x,out);
}

//Anchoring energy of a surface site (10^-21 J): density times the area the site represents.
MARLICS_HD inline double surface_site_energy(const LatticeView& L, const double * Qij, const SurfaceSite& S)
{
  const double dQ[15]={0.};
  return S.area*boundary_energy(L.bulk,L.bc[S.bc],&Qij[5*S.site],dQ,S.v);
}

#endif
