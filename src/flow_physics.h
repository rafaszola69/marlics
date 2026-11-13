#ifndef FLOW_PHYSICS_H
#define FLOW_PHYSICS_H
// Coupling between Q and the flow (backflow), per site, shared by the CPU and GPU code.
// Beris-Edwards nematodynamics in the form of Denniston, Orlandini and Yeomans, PRE 63, 056702 (2001):
//   dQ/dt + u.grad Q - S(W,Q) = Gamma H,     eta lap u - grad p + div sigma = 0,  div u = 0  (Stokes limit)
// written there for q=2Q/3 (Marlics uses Q=(S/2)(3nn-I)): h=(3/2)H, with H the molecular field that
// fill_ki already gives (k=Lambda*H at bulk sites), and dQ/dt = (3/2) dq/dt.
// Units: u in nm/us, sigma in 10^6 Pa, eta in Pa s, time in us.
// Periodic domains (bulk geometry), two flat walls (flow_walls.h) or any other geometry
// (flow_confined.h): every derivative is a central difference, wrapping around periodic axes; with
// surface sites the flow terms act at bulk sites only, the velocity vanishes elsewhere, and the
// stress at the non-bulk sites next to the bulk is extrapolated from it (boundary_stress).
#include <cstdio>
#include "physics.h"
#include "lattice.h"
#include "physics_ldg_gradient.h"

struct FlowParams
{
  double eta;   //isotropic viscosity (Pa s)
  double xi;    //flow-aligning parameter
  double zeta;  //activity: active stress -zeta*q (10^6 Pa)
  double zeta_c;//chiral activity: active stress -zeta_c*(eps.q-q.eps)/2, eps_ij=eps_ijz (10^6 Pa)
  double el;    //1, or 0 to drop the passive (elastic and flow-alignment) stress
  const double * xis; //flow-aligning parameter per site (xi_wall near the walls), or NULL: xi everywhere
  double xi_core;     //flow-aligning parameter at the sites with tr Q^2 < q2_core (defect cores)
  double q2_core;     //tr Q^2 = (3/2) S^2 for uniaxial Q; 0: no core mask
};

//Flow-aligning parameter at site s: xi_core in the defect cores (low order), else xis or xi.
MARLICS_HD inline double site_xi(const FlowParams& F, long s, const double * Qij)
{
  if(F.q2_core>0)
    {
      const double * v=&Qij[5*s];
      const double zz=-v[0]-v[3];
      if(v[0]*v[0]+v[3]*v[3]+zz*zz+2*(v[1]*v[1]+v[2]*v[2]+v[4]*v[4])<F.q2_core) return F.xi_core;
    }
  return F.xis ? F.xis[s] : F.xi;
}

//Host only: whether the liquid crystal aligns with a shear flow, and at which angle from the flow
//direction (cos 2 theta = 3S/((2+S) xi), Denniston, Orlandini and Yeomans 2001), or tumbles.
inline void flow_alignment_report(double xi, double S)
{
  const double xi_c=3*S/(2+S);
  if(xi>xi_c)
    printf("xi=%g > 3S/(2+S)=%.3f (S=%.3f): flow aligning, at %.1f degrees from the flow direction.\n",
           xi, xi_c, S, 0.5*acos(xi_c/xi)*180/3.141592653589793);
  else
    printf("xi=%g < 3S/(2+S)=%.3f (S=%.3f): flow tumbling (no aligning angle).\n", xi, xi_c, S);
}

//Neighbour of site s along axis at distance dir (+1/-1), periodic.
MARLICS_HD inline int periodic_neighbor(const LatticeView& L, int s, int axis, int dir)
{
  int c[3]={s%L.Nx, (s/L.Nx)%L.Ny, s/(L.Nx*L.Ny)};
  const int n[3]={L.Nx,L.Ny,L.Nz};
  c[axis]=(c[axis]+dir+n[axis])%n[axis];
  return lattice_index(L,c[0],c[1],c[2]);
}

//Central first derivatives of an array with m components per site: d[m*axis+c].
MARLICS_HD inline void central_gradient(const LatticeView& L, const double * f, int m, int s, double * d)
{
  for(int a=0; a<3; a++)
    {
      const int p=periodic_neighbor(L,s,a,1), q=periodic_neighbor(L,s,a,-1);
      for(int c=0; c<m; c++) d[m*a+c]=0.5*(f[m*p+c]-f[m*q+c])*L.d_1[a];
    }
}

MARLICS_HD inline void full_tensor(const double * v, double M[3][3])
{
  M[0][0]=v[0]; M[0][1]=v[1]; M[0][2]=v[2];
  M[1][0]=v[1]; M[1][1]=v[3]; M[1][2]=v[4];
  M[2][0]=v[2]; M[2][1]=v[4]; M[2][2]=-v[0]-v[3];
}

//Stress at site s (sigma[3*i+j], force density f_i = d_j sigma_ij), from Q and k=Lambda*H at s:
//  sigma = -xi (h A + A h) + 2 xi A (q:h) - T + q h - h q - zeta q - zeta_c (eps q - q eps)/2,
//  A=q+I/3, T_ij = d_i Q_c df/d(d_j Q_c) (distortion stress; isotropic parts go into the pressure).
//The passive terms (all but the last two) are multiplied by el (elastic_stress 0 drops them).
//The chiral active stress (Hoffmann, Schakenraad, Merks and Giomi, Soft Matter 16, 764 (2020):
//-2 tau eps.Q) is symmetric, with eps_ij=eps_ijz (eps_xy=1) the rotation about z: for q in the xy
//plane it is -zeta_c eps.q, the stress -zeta_c q of the director turned by -45 degrees, so a +1/2
//defect self-propels at -atan(zeta_c/zeta) from its axis.
MARLICS_HD inline void site_stress(const LatticeView& L, const FlowParams& F, const double * Qij, const double * k_i, int s, double sigma[9])
{
  double dQ[15], g[15], q[3][3], h[3][3], A[3][3], qh=0;
  const double xi=site_xi(F,s,Qij);
  central_gradient(L,Qij,5,s,dQ);
  ldg_gradient_derivative(L.bulk,&Qij[5*s],dQ,g);

  double hv[5];
  for(int c=0; c<5; c++) hv[c]=1.5*k_i[5*s+c]/L.bulk.Lambda;
  full_tensor(&Qij[5*s],q);
  full_tensor(hv,h);
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      {
        q[i][j]*=2./3.;
        A[i][j]=q[i][j]+(i==j ? 1./3. : 0.);
        qh+=q[i][j]*h[i][j];
      }

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      {
        double hA=0, Ah=0, qhm=0, hq=0, T=0;
        for(int l=0; l<3; l++)
          {
            hA+=h[i][l]*A[l][j]; Ah+=A[i][l]*h[l][j];
            qhm+=q[i][l]*h[l][j]; hq+=h[i][l]*q[l][j];
          }
        for(int c=0; c<5; c++) T+=dQ[5*i+c]*g[5*j+c];
        //(eps q)_ij - (q eps)_ij with eps_01=1, eps_10=-1
        const double eq=(i==0 ? q[1][j] : i==1 ? -q[0][j] : 0.)-(j==0 ? -q[i][1] : j==1 ? q[i][0] : 0.);
        sigma[3*i+j]=F.el*(-xi*(hA+Ah)+2*xi*A[i][j]*qh-T+qhm-hq)-F.zeta*q[i][j]-0.5*F.zeta_c*eq;
      }
}

//Stress at a non-bulk site s (the force at bulk sites is a central difference of the stress, which
//reaches the non-bulk sites next to the bulk): the mean, over its bulk face neighbours x, of the
//linear extrapolation 2 sigma(x) - sigma(2x-s) (sigma(x) when 2x-s is not bulk). A uniform stress
//then drives no flow, and an isotropic one is balanced by the pressure. Zero without bulk neighbours.
MARLICS_HD inline void boundary_stress(const LatticeView& L, long s, double * sigma)
{
  double acc[9]={0,0,0,0,0,0,0,0,0};
  int count=0;
  for(int a=0; a<3; a++)
    for(int dir=-1; dir<=1; dir+=2)
      {
        const int x=periodic_neighbor(L,(int)s,a,dir);
        if(L.point_type[x]!=1) continue;
        const int x2=periodic_neighbor(L,x,a,dir);
        const bool two= L.point_type[x2]==1;
        for(int c=0; c<9; c++) acc[c]+= two ? 2*sigma[9*(long)x+c]-sigma[9*(long)x2+c] : sigma[9*(long)x+c];
        count++;
      }
  for(int c=0; c<9; c++) sigma[9*s+c]= count ? acc[c]/count : 0.;
}

//Force density f_i = d_j sigma_ij at site s.
MARLICS_HD inline void site_force(const LatticeView& L, const double * sigma, int s, double f[3])
{
  for(int i=0; i<3; i++) f[i]=0;
  for(int j=0; j<3; j++)
    {
      const int p=periodic_neighbor(L,s,j,1), q=periodic_neighbor(L,s,j,-1);
      for(int i=0; i<3; i++) f[i]+=0.5*(sigma[9*p+3*i+j]-sigma[9*q+3*i+j])*L.d_1[j];
    }
}

//Stokes problem for one Fourier mode (m0,m1,m2), with the symbols of the central first difference
//(i s_a) and of the 3-point Laplacian (lap), so that the velocity is divergence-free in the same
//discrete sense as the rest of the code:  u = (f - s (s.f)/|s|^2) / (eta |lap|).
//The uniform mode carries no velocity; modes with s=0 and lap>0 (checkerboards) carry none either in a
//periodic box, while keep_checkerboards gives them u = f/(eta lap) (they feel no pressure), which makes
//this the exact inverse of the periodic Stokes operator apart from the uniform mode (flow_confined.h).
//fr/fi: real and imaginary parts.
MARLICS_HD inline void stokes_mode(const LatticeView& L, double eta, int m0, int m1, int m2, double fr[3], double fi[3], bool keep_checkerboards=false)
{
  const int m[3]={m0,m1,m2}, n[3]={L.Nx,L.Ny,L.Nz};
  const double two_pi=6.283185307179586;
  double s[3], s2=0, lap=0;
  for(int a=0; a<3; a++)
    {
      const double ka=two_pi*m[a]/n[a];
      s[a]=sin(ka)*L.d_1[a];
      s2+=s[a]*s[a];
      lap+=(2-2*cos(ka))*L.d_1[a]*L.d_1[a];
    }
  if(lap==0 || (s2<1e-12*lap && !keep_checkerboards))
    {
      for(int a=0; a<3; a++) fr[a]=fi[a]=0;
      return;
    }
  if(s2<1e-12*lap)
    {
      for(int a=0; a<3; a++) { fr[a]/=eta*lap; fi[a]/=eta*lap; }
      return;
    }
  const double sfr=s[0]*fr[0]+s[1]*fr[1]+s[2]*fr[2], sfi=s[0]*fi[0]+s[1]*fi[1]+s[2]*fi[2];
  for(int a=0; a<3; a++)
    {
      fr[a]=(fr[a]-s[a]*sfr/s2)/(eta*lap);
      fi[a]=(fi[a]-s[a]*sfi/s2)/(eta*lap);
    }
}

//Flow terms of dQ/dt at site s: -u.grad Q + (3/2) S(W,q), with W_ab = d_b u_a and
//  S = (xi D + Omega) A + A (xi D - Omega) - 2 xi A tr(q W).
MARLICS_HD inline void site_flow_rhs(const LatticeView& L, const FlowParams& F, const double * Qij, const double * u, int s, double out[5])
{
  double dQ[15], du[9], W[3][3], D[3][3], O[3][3], q[3][3], A[3][3], qW=0;
  const double xi=site_xi(F,s,Qij);
  central_gradient(L,Qij,5,s,dQ);
  central_gradient(L,u,3,s,du);                       //du[3*b+a] = d_b u_a
  full_tensor(&Qij[5*s],q);
  for(int a=0; a<3; a++)
    for(int b=0; b<3; b++) W[a][b]=du[3*b+a];
  for(int a=0; a<3; a++)
    for(int b=0; b<3; b++)
      {
        D[a][b]=0.5*(W[a][b]+W[b][a]);
        O[a][b]=0.5*(W[a][b]-W[b][a]);
        q[a][b]*=2./3.;
        A[a][b]=q[a][b]+(a==b ? 1./3. : 0.);
      }
  for(int a=0; a<3; a++)
    for(int b=0; b<3; b++) qW+=q[a][b]*W[b][a];

  const int ia[5]={0,0,0,1,1}, ib[5]={0,1,2,1,2};
  const double * us=&u[3*s];
  for(int c=0; c<5; c++)
    {
      const int a=ia[c], b=ib[c];
      double S=-2*xi*A[a][b]*qW;
      for(int l=0; l<3; l++)
        S+=(xi*D[a][l]+O[a][l])*A[l][b]+A[a][l]*(xi*D[l][b]-O[l][b]);
      const double adv=us[0]*dQ[c]+us[1]*dQ[5+c]+us[2]*dQ[10+c];
      out[c]=-adv+1.5*S;
    }
}

#endif
