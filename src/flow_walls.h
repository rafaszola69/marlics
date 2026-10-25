#ifndef FLOW_WALLS_H
#define FLOW_WALLS_H
// Stokes flow between two flat walls (geometries channel and slab): Fourier modes along the two
// periodic axes, and for each mode a one-dimensional problem along the wall normal, with no-slip
// at the two rows of surface sites (k=0 and k=M-1; the upper wall may move, wall_velocity).
// This is the problem of flow_confined.h for two flat walls: for each mode the unknowns are u at the
// bulk rows k=1..M-2 (m of them) and p at every row, with central first differences for div u and
// grad p and the 3-point Laplacian; the continuity equation holds at every row, the walls included:
//   T u_a - (B^H p)_a = f_a,   B u = 0,   T = eta (K + lap),  B u = i s.u_inplane + Dz u_normal,
// K the 3-point -d^2/dz^2 with u=0 at the walls, Dz the central difference (rows 1..M-2 -> 0..M-1).
// The pressure is eliminated: S p = B T^-1 f, S = B T^-1 B^H = |s|^2 T^-1 + Dz T^-1 Dz^T (M x M), then
// u = T^-1 f + T^-1 B^H p. S depends only on the mode and the geometry: its inverse is computed once
// per |mode| (key). A force that is a discrete gradient (with the stress of boundary_stress at the
// walls) is balanced by the pressure and drives no flow.
// Modes with s=0 (the in-plane mean and checkerboards): no normal flow, u_inplane = T^-1 f (+ the
// moving wall for the mean).
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "physics.h"
#include "lattice.h"

struct WallFlowView
{
  int axis;            //wall-normal axis
  int pa, pb;          //in-plane (periodic) axes
  int M;               //sites along the wall normal: walls at 0 and M-1
  int Ka;              //keys along pa (N_pa/2+1)
  const double * Sinv; //M^2 per key, row-major
  double U[3];         //velocity of the wall at M-1 (nm/us)
};

MARLICS_HD inline long axis_step(const LatticeView& L, int a) { return a==0 ? 1 : (a==1 ? (long)L.Nx : (long)L.Nx*L.Ny); }
MARLICS_HD inline int axis_size(const LatticeView& L, int a) { return a==0 ? L.Nx : (a==1 ? L.Ny : L.Nz); }

//A complex column along the wall normal: element k at p[2*st*k] (real), p[2*st*k+1] (imaginary).
struct WallColumn
{
  double * p;
  long st;
  MARLICS_HD double& re(int k) const { return p[2*st*k]; }
  MARLICS_HD double& im(int k) const { return p[2*st*k+1]; }
};

//Solves T x = r on the rows 1..m of x (in place), cp: scratch of M doubles.
MARLICS_HD inline void wall_thomas(WallColumn x, double * cp, double diag, double off, int m)
{
  cp[1]=off/diag;
  x.re(1)/=diag; x.im(1)/=diag;
  for(int k=2; k<=m; k++)
    {
      const double den=diag-off*cp[k-1];
      cp[k]=off/den;
      x.re(k)=(x.re(k)-off*x.re(k-1))/den;
      x.im(k)=(x.im(k)-off*x.im(k-1))/den;
    }
  for(int k=m-1; k>=1; k--)
    {
      x.re(k)-=cp[k]*x.re(k+1);
      x.im(k)-=cp[k]*x.im(k+1);
    }
}

//Symbols of the in-plane mode (ia,ib): s[3] (central difference, zero along the normal), lap.
MARLICS_HD inline void wall_symbols(const LatticeView& L, const WallFlowView& W, int ia, int ib, double s[3], double& lap)
{
  const double two_pi=6.283185307179586;
  const int Na=axis_size(L,W.pa), Nb=axis_size(L,W.pb);
  const double ka=two_pi*ia/Na, kb=two_pi*ib/Nb;
  s[0]=s[1]=s[2]=0;
  s[W.pa]=sin(ka)*L.d_1[W.pa];
  s[W.pb]=sin(kb)*L.d_1[W.pb];
  lap=(2-2*cos(ka))*L.d_1[W.pa]*L.d_1[W.pa]+(2-2*cos(kb))*L.d_1[W.pb]*L.d_1[W.pb];
}

//Velocity of in-plane mode t (0..N_pa*N_pb-1) from its force, in place in F[3] (interleaved complex,
//the in-plane Fourier transform of the force, zero on the wall rows). scratch: 5*M doubles.
//n_plane: number of sites of a plane (the mean mode holds n_plane times the mean).
MARLICS_HD inline void wall_stokes_mode(const LatticeView& L, const WallFlowView& W, double eta, long t, double * const F[3], double * scratch, long n_plane)
{
  const int Na=axis_size(L,W.pa), Nb=axis_size(L,W.pb), M=W.M, m=M-2;
  const int ia=(int)(t%Na), ib=(int)(t/Na);
  const long base=ia*axis_step(L,W.pa)+ib*axis_step(L,W.pb), st=axis_step(L,W.axis);
  WallColumn u[3];
  for(int a=0; a<3; a++) u[a]=WallColumn{F[a]+2*base,st};

  double s[3], lap;
  wall_symbols(L,W,ia,ib,s,lap);
  const double h_1=L.d_1[W.axis], off=-eta*h_1*h_1, diag=eta*(2*h_1*h_1+lap);
  double * cp=scratch;

  if((ia==0 || 2*ia==Na) && (ib==0 || 2*ib==Nb))   //s=0: in-plane mean (shear and pressure-driven flows), checkerboards
    {
      for(int a=0; a<3; a++)
        {
          if(a==W.axis) { for(int k=0; k<M; k++) u[a].re(k)=u[a].im(k)=0; continue; }
          if(ia==0 && ib==0) u[a].re(m)-=off*W.U[a]*n_plane;
          wall_thomas(u[a],cp,diag,off,m);
        }
      return;
    }

  //g = T^-1 f
  for(int a=0; a<3; a++) wall_thomas(u[a],cp,diag,off,m);

  //r = B g = i s.g + Dz g_normal at every row (g=0 on the walls; scratch column R), then p = S^-1 r (column P).
  WallColumn R{scratch+M,1}, P{scratch+3*M,1};
  const WallColumn& gz=u[W.axis];
  const int key=(ia<Na-ia ? ia : Na-ia)+W.Ka*(ib<Nb-ib ? ib : Nb-ib);
  const double * Si=W.Sinv+(long)key*M*M;
  for(int k=0; k<M; k++)
    {
      double rr=0, ri=0;
      if(k>0 && k<M-1)
        for(int a=0; a<3; a++) { rr-=s[a]*u[a].im(k); ri+=s[a]*u[a].re(k); }
      const double zp_r= k<M-1 ? gz.re(k+1) : 0, zp_i= k<M-1 ? gz.im(k+1) : 0;
      const double zm_r= k>0 ? gz.re(k-1) : 0, zm_i= k>0 ? gz.im(k-1) : 0;
      R.re(k)=rr+0.5*h_1*(zp_r-zm_r);
      R.im(k)=ri+0.5*h_1*(zp_i-zm_i);
    }
  for(int k=0; k<M; k++)
    {
      double pr=0, pi=0;
      const double * row=Si+(long)k*M;
      for(int l=0; l<M; l++) { pr+=row[l]*R.re(l); pi+=row[l]*R.im(l); }
      P.re(k)=pr; P.im(k)=pi;
    }

  //u_normal = g_normal + T^-1 Dz p
  for(int k=1; k<=m; k++)
    {
      R.re(k)=0.5*h_1*(P.re(k+1)-P.re(k-1));
      R.im(k)=0.5*h_1*(P.im(k+1)-P.im(k-1));
    }
  wall_thomas(R,cp,diag,off,m);
  for(int k=1; k<=m; k++) { gz.re(k)+=R.re(k); gz.im(k)+=R.im(k); }

  //u_inplane = g_inplane + i s T^-1 p
  wall_thomas(P,cp,diag,off,m);
  for(int a=0; a<3; a++)
    {
      if(a==W.axis) continue;
      for(int k=1; k<=m; k++)
        {
          u[a].re(k)-=s[a]*P.im(k);
          u[a].im(k)+=s[a]*P.re(k);
        }
    }
}

//Host side: checks the geometry (one non-periodic axis, surface sites exactly on its first and last
//rows) and computes S^-1 for every key.
class WallStokes
{
 public:

  //Whether the domain is two flat walls: one non-periodic axis (at least 4 sites across), surface
  //sites exactly on its first and last rows, bulk sites elsewhere.
  static bool detect(const LatticeView& L, int& axis)
  {
    int nonperiodic=0;
    axis=-1;
    for(int a=0; a<3; a++) if(!L.periodic[a]) { nonperiodic++; axis=a; }
    if(nonperiodic!=1 || axis_size(L,axis)<4) return false;
    const long n=(long)L.Nx*L.Ny*L.Nz, st=axis_step(L,axis);
    const int M=axis_size(L,axis);
    for(long s=0; s<n; s++)
      {
        const int k=(int)((s/st)%M);
        if( (k==0 || k==M-1) ? L.point_type[s]<2 : L.point_type[s]!=1 ) return false;
      }
    return true;
  }

  WallStokes(const LatticeView& L, int axis, double eta, const double U[3])
  {
    W.axis=axis;
    W.pa= axis==0 ? 1 : 0;
    W.pb= axis==2 ? 1 : 2;
    W.M=axis_size(L,axis);
    W.Ka=axis_size(L,W.pa)/2+1;
    for(int a=0; a<3; a++) W.U[a]=U[a];
    if(U[axis]!=0)
      {
        std::cout << "wall_velocity must be parallel to the walls.\nAborting the program.\n";
        exit(1);
      }
    const int Kb=axis_size(L,W.pb)/2+1, M=W.M;
    const long keys=(long)W.Ka*Kb;
    Sinv.assign(keys*M*M,0.);
    std::cout << "Stokes flow between walls (normal along " << "xyz"[axis] << "): " << keys << " pressure matrices, "
              << keys*M*M*8/1048576. << " MiB.\n";

    #pragma omp parallel for schedule(dynamic)
    for(long key=0; key<keys; key++)
      {
        double s[3], lap;
        wall_symbols(L,W,(int)(key%W.Ka),(int)(key/W.Ka),s,lap);
        const double s2=s[0]*s[0]+s[1]*s[1]+s[2]*s[2];
        if(s2<1e-12*lap || lap==0) continue;   //s=0: no pressure
        schur_inverse(L.d_1[axis],eta,lap,s2,M,&Sinv[key*M*M]);
      }
    W.Sinv=Sinv.data();
  }

  WallFlowView W;
  std::vector<double> Sinv;

 private:

  //S^-1 for one key: S = s2 T^-1 + Dz T^-1 Dz^T (M x M; T^-1 acts on the m=M-2 bulk rows, the
  //in-plane term vanishes on the wall rows), T = eta (K + lap).
  static void schur_inverse(double h_1, double eta, double lap, double s2, int M, double * out)
  {
    const int m=M-2;
    std::vector<double> X(m*m), Y((long)m*M), cp(m+2);
    const double off=-eta*h_1*h_1, diag=eta*(2*h_1*h_1+lap);
    #define AT(A,i,j) A[(long)(i)*ld+(j)]
    long ld=m;
    //X = T^-1, column by column (T is symmetric, so X is too). X(i,j): bulk rows i+1, j+1.
    std::vector<double> col(2*(m+2));
    for(int j=0; j<m; j++)
      {
        for(double& c : col) c=0;
        col[2*(j+1)]=1;
        wall_thomas(WallColumn{col.data(),1},cp.data(),diag,off,m);
        for(int i=0; i<m; i++) AT(X,i,j)=col[2*(i+1)];
      }
    //Dz(i,k) = (h_1/2)(delta_{k,i+1} - delta_{k,i-1}), row i=0..M-1, bulk row k=1..m (column k-1).
    //Y = X Dz^T (m x M): Y(r,i) = (h_1/2)(X(r,i) - X(r,i-2)).
    ld=M;
    for(int r=0; r<m; r++)
      for(int i=0; i<M; i++)
        {
          const double a= i<m ? X[(long)r*m+i] : 0, b= i>=2 ? X[(long)r*m+i-2] : 0;
          AT(Y,r,i)=0.5*h_1*(a-b);
        }
    //S(i,j) = s2 X(i-1,j-1) + (h_1/2)(Y(i,j) - Y(i-2,j)).
    for(int i=0; i<M; i++)
      for(int j=0; j<M; j++)
        {
          double v= (i>=1 && i<=m && j>=1 && j<=m) ? s2*X[(long)(i-1)*m+j-1] : 0;
          v+=0.5*h_1*((i<m ? AT(Y,i,j) : 0)-(i>=2 ? AT(Y,i-2,j) : 0));
          AT(out,i,j)=v;
        }
    const int n=M;
    //Inverse by Cholesky: S = C C^T, then S^-1 = C^-T C^-1.
    std::vector<double> C((long)n*n,0.), Ci((long)n*n,0.);
    for(int j=0; j<n; j++)
      {
        double d=AT(out,j,j);
        for(int l=0; l<j; l++) d-=AT(C,j,l)*AT(C,j,l);
        AT(C,j,j)=sqrt(d);
        for(int i=j+1; i<n; i++)
          {
            double v=AT(out,i,j);
            for(int l=0; l<j; l++) v-=AT(C,i,l)*AT(C,j,l);
            AT(C,i,j)=v/AT(C,j,j);
          }
      }
    for(int j=0; j<n; j++)
      for(int i=j; i<n; i++)
        {
          double v= i==j ? 1. : 0.;
          for(int l=j; l<i; l++) v-=AT(C,i,l)*AT(Ci,l,j);
          AT(Ci,i,j)=v/AT(C,i,i);
        }
    for(int i=0; i<n; i++)
      for(int j=0; j<=i; j++)
        {
          double v=0;
          for(int l=i; l<n; l++) v+=AT(Ci,l,i)*AT(Ci,l,j);
          AT(out,i,j)=AT(out,j,i)=v;
        }
    #undef AT
  }
};

#endif
