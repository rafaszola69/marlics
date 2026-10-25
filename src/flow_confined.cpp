// Setup and Lagrange-force solve of the Stokes flow in general geometries (flow_confined.h).
#include "flow_confined.h"
#include "flow_physics.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>

static inline double dot(const double * a, const double * b, long m)
{
  double v=0;
  for(long k=0; k<m; k++) v+=a[k]*b[k];
  return v;
}

bool ConfinedStokes::needed(const LatticeView& L)
{
  const long n=(long)L.Nx*L.Ny*L.Nz;
  for(long s=0; s<n; s++) if(L.point_type[s]!=1) return true;
  return false;
}

ConfinedStokes::ConfinedStokes(const LatticeView& L, double eta, int solver) : ne(0), direct(false), iterations(0)
{
  const long n=(long)L.Nx*L.Ny*L.Nz;
  std::vector<char> e1(n,0), in2(n,0);
  auto bulk=[&](long s){ return L.point_type[s]==1; };
  for(long s=0; s<n; s++)
    if(!bulk(s))
      for(int a=0; a<3; a++)
        for(int dir=-1; dir<=1; dir+=2) if(bulk(periodic_neighbor(L,(int)s,a,dir))) e1[s]=1;
  for(long s=0; s<n; s++)
    if(!bulk(s))
      {
        bool near=e1[s];
        for(int a=0; a<3; a++)
          for(int dir=-1; dir<=1; dir+=2) if(e1[periodic_neighbor(L,(int)s,a,dir)]) near=true;
        in2[s]=near;
      }
  std::vector<long> index(n,-1);
  for(long s=0; s<n; s++) if(in2[s]) { index[s]=(long)e2.size(); e2.push_back(s); }
  ne=(long)e2.size();
  e2_nb.resize(6*ne);
  for(long i=0; i<ne; i++)
    for(int a=0; a<3; a++)
      for(int q=0; q<2; q++) e2_nb[6*i+2*a+q]=index[periodic_neighbor(L,(int)e2[i],a,q ? 1 : -1)];
  for(int a=0; a<3; a++) d2.push_back(L.d_1[a]*L.d_1[a]);

  const long nc=3*ne;
  direct= solver==1 || (solver==0 && nc<=max_direct);
  std::cout << "Stokes flow in a confined geometry: no-slip on " << ne << " sites (two layers around the bulk), "
            << nc << " Lagrange forces, " << (direct ? "direct" : "iterative (conjugate gradients)") << " solve.\n";
  //Response on E2 to a unit force along a on every E2 site.
  {
    std::vector<std::complex<double>> F(3*n);
    g_uniform.assign(3*nc,0.);
    for(int b=0; b<3; b++)
      {
        for(long s=0; s<3*n; s++) F[s]=0;
        for(long i=0; i<ne; i++) F[b*n+e2[i]]=1;
        fft_periodic(F.data(),L.Nx,L.Ny,L.Nz,3,true);
        for(long s=0; s<n; s++)
          {
            double fr[3], fi[3];
            for(int a=0; a<3; a++) { fr[a]=F[a*n+s].real(); fi[a]=F[a*n+s].imag(); }
            stokes_mode(L,eta,(int)(s%L.Nx),(int)((s/L.Nx)%L.Ny),(int)(s/((long)L.Nx*L.Ny)),fr,fi,true);
            for(int a=0; a<3; a++) F[a*n+s]=std::complex<double>(fr[a],fi[a]);
          }
        fft_periodic(F.data(),L.Nx,L.Ny,L.Nz,3,false);
        for(long i=0; i<ne; i++)
          for(int a=0; a<3; a++) g_uniform[b*nc+3*i+a]=F[a*n+e2[i]].real();
      }
  }
  if(direct) build_direct(L,eta);
  if(!direct) last.assign(nc,0.);
}

void ConfinedStokes::project(double * v) const
{
  for(int a=0; a<3; a++)
    {
      double m=0;
      for(long i=0; i<ne; i++) m+=v[3*i+a];
      m/=ne;
      for(long i=0; i<ne; i++) v[3*i+a]-=m;
    }
}

void ConfinedStokes::precondition(const double * r, double * z) const
{
  #pragma omp parallel for schedule(static)
  for(long i=0; i<ne; i++)
    for(int c=0; c<3; c++)
      {
        double v=0;
        for(int a=0; a<3; a++)
          {
            v+=2*d2[a]*r[3*i+c];
            for(int q=0; q<2; q++)
              {
                const long j=e2_nb[6*i+2*a+q];
                if(j>=0) v-=d2[a]*r[3*j+c];
              }
          }
        z[3*i+c]=v;
      }
}

void ConfinedStokes::solve(double * g, const double total[3], double * lambda, const Apply& Cmul, const Apply& Dmul)
{
  const long nc=3*ne;
  double l0[3];
  for(int b=0; b<3; b++) l0[b]=-total[b]/ne;
  std::vector<double> w(nc);
  for(long r=0; r<nc; r++) w[r]=-(g[r]+l0[0]*g_uniform[r]+l0[1]*g_uniform[nc+r]+l0[2]*g_uniform[2*nc+r]);
  project(w.data());
  solve_projected(w.data(),lambda,Cmul,Dmul);
  for(long i=0; i<ne; i++)
    for(int a=0; a<3; a++) lambda[3*i+a]+=l0[a];
}

void ConfinedStokes::solve_projected(const double * w_in, double * lambda, const Apply& Cmul, const Apply& Dmul)
{
  const long nc=3*ne;
  std::vector<double> w(w_in,w_in+nc);
  if(direct)
    {
      if(Dmul) Dmul(w.data(),lambda);
      else
        {
          const double * C=Cinv.data();
          #pragma omp parallel for schedule(static)
          for(long r=0; r<nc; r++) lambda[r]=dot(&C[r*nc],w.data(),nc);
        }
      project(lambda);
      return;
    }

  //Projected, preconditioned conjugate gradients from the previous solution.
  const double tol=1e-10;
  const int max_iterations=5000;
  std::vector<double> x=last, r(nc), z(nc), p(nc), Ap(nc);
  const double wn=sqrt(dot(w.data(),w.data(),nc));
  if(wn==0)
    {
      for(long k=0; k<nc; k++) lambda[k]=0;
      return;
    }
  Cmul(x.data(),Ap.data());
  project(Ap.data());
  for(long k=0; k<nc; k++) r[k]=w[k]-Ap[k];
  precondition(r.data(),z.data());
  project(z.data());
  p=z;
  double rz=dot(r.data(),z.data(),nc);
  int it=0;
  for(; it<max_iterations && sqrt(dot(r.data(),r.data(),nc))>tol*wn; it++)
    {
      Cmul(p.data(),Ap.data());
      project(Ap.data());
      const double alpha=rz/dot(p.data(),Ap.data(),nc);
      for(long k=0; k<nc; k++) { x[k]+=alpha*p[k]; r[k]-=alpha*Ap[k]; }
      precondition(r.data(),z.data());
      project(z.data());
      const double rzn=dot(r.data(),z.data(),nc);
      for(long k=0; k<nc; k++) p[k]=z[k]+rzn/rz*p[k];
      rz=rzn;
    }
  if(it==max_iterations) std::cout << "Warning: the confined Stokes solve did not converge in " << it << " iterations.\n";
  iterations=it;
  last=x;
  for(long k=0; k<nc; k++) lambda[k]=x[k];
}

void ConfinedStokes::build_direct(const LatticeView& L, double eta)
{
  const auto t0=std::chrono::steady_clock::now();
  const int N[3]={L.Nx,L.Ny,L.Nz};
  const long n=(long)L.Nx*L.Ny*L.Nz, nc=3*ne;

  //Periodic Green's function: Gr[(3b+a)n+x] = velocity a at x for a unit force along b at the origin.
  std::vector<double> Gr(9*n);
  {
    std::vector<std::complex<double>> F(3*n);
    for(int b=0; b<3; b++)
      {
        for(long s=0; s<n; s++)
          {
            double fr[3]={0,0,0}, fi[3]={0,0,0};
            fr[b]=1;
            stokes_mode(L,eta,(int)(s%L.Nx),(int)((s/L.Nx)%L.Ny),(int)(s/((long)L.Nx*L.Ny)),fr,fi,true);
            for(int a=0; a<3; a++) F[a*n+s]=std::complex<double>(fr[a],fi[a]);
          }
        fft_periodic(F.data(),L.Nx,L.Ny,L.Nz,3,false);
        for(int a=0; a<3; a++)
          for(long s=0; s<n; s++) Gr[(3*b+a)*n+s]=F[a*n+s].real();
      }
  }

  //C = E G E^T, then P C P.
  Cinv.assign(nc*nc,0.);
  double * C=Cinv.data();
  #pragma omp parallel for schedule(static)
  for(long i=0; i<ne; i++)
    {
      const long si=e2[i];
      const int ci[3]={(int)(si%L.Nx),(int)((si/L.Nx)%L.Ny),(int)(si/((long)L.Nx*L.Ny))};
      for(long j=0; j<ne; j++)
        {
          const long sj=e2[j];
          const int cj[3]={(int)(sj%L.Nx),(int)((sj/L.Nx)%L.Ny),(int)(sj/((long)L.Nx*L.Ny))};
          int d[3];
          for(int a=0; a<3; a++) d[a]=((ci[a]-cj[a])%N[a]+N[a])%N[a];
          const long t=lattice_index(L,d[0],d[1],d[2]);
          for(int a=0; a<3; a++)
            for(int b=0; b<3; b++) C[(3*i+a)*nc+3*j+b]=Gr[(3*b+a)*n+t];
        }
    }
  #pragma omp parallel for schedule(static)
  for(long r=0; r<nc; r++) project(&C[r*nc]);
  {
    //Columns: mean over the rows of each component.
    std::vector<double> mean(3*nc,0.);
    for(long i=0; i<ne; i++)
      for(int a=0; a<3; a++)
        {
          const double * row=&C[(3*i+a)*nc];
          double * m=&mean[a*nc];
          for(long c=0; c<nc; c++) m[c]+=row[c];
        }
    #pragma omp parallel for schedule(static)
    for(long r=0; r<nc; r++)
      for(long c=0; c<nc; c++) C[r*nc+c]-=mean[(r%3)*nc+c]/ne;
    //Symmetrize, removing rounding differences.
    #pragma omp parallel for schedule(static)
    for(long r=0; r<nc; r++)
      for(long c=0; c<r; c++)
        {
          const double v=0.5*(C[r*nc+c]+C[c*nc+r]);
          C[r*nc+c]=C[c*nc+r]=v;
        }
  }
  double scale=0;
  for(long r=0; r<nc; r++) scale+=C[r*nc+r];
  scale/=nc;

  //Complete C on the directions it does not see: the mean of each component (taken by U0)...
  for(long i=0; i<ne; i++)
    for(long j=0; j<ne; j++)
      for(int a=0; a<3; a++) C[(3*i+a)*nc+3*j+a]+=scale/ne;

  //...and the pressure fields that move no fluid: lambda = grad psi (central differences) with
  //grad psi = 0 off E2. Those psi are constant on the classes of the relation psi(x+e) = psi(x-e) for
  //x off E2 (union-find); each class K gives grad 1_K, nonzero only on E2 (pressure jumps across the
  //walls, pressure at sites surrounded by E2).
  std::vector<long> index(n,-1), parent(n);
  for(long i=0; i<ne; i++) index[e2[i]]=i;
  for(long s=0; s<n; s++) parent[s]=s;
  std::function<long(long)> root=[&](long s){ while(parent[s]!=s) { parent[s]=parent[parent[s]]; s=parent[s]; } return s; };
  for(long s=0; s<n; s++)
    if(index[s]<0)
      for(int a=0; a<3; a++)
        {
          const long p=root(periodic_neighbor(L,(int)s,a,1)), q=root(periodic_neighbor(L,(int)s,a,-1));
          if(p!=q) parent[p]=q;
        }
  std::vector<long> cls(n,-1);
  long classes=0;
  std::vector<long> class_of_root(n,-1);
  for(long s=0; s<n; s++)
    {
      const long r=root(s);
      if(class_of_root[r]<0) class_of_root[r]=classes++;
      cls[s]=class_of_root[r];
    }
  //Sparse gradients of the class indicators, on E2.
  std::vector<std::vector<std::pair<long,double>>> grad(classes);
  for(long i=0; i<ne; i++)
    for(int a=0; a<3; a++)
      {
        const long p=cls[periodic_neighbor(L,(int)e2[i],a,1)], q=cls[periodic_neighbor(L,(int)e2[i],a,-1)];
        if(p==q) continue;
        grad[p].push_back({3*i+a,0.5*L.d_1[a]});
        grad[q].push_back({3*i+a,-0.5*L.d_1[a]});
      }
  std::vector<std::vector<double>> nulls;
  for(long k=0; k<classes; k++)
    {
      if(grad[k].empty()) continue;
      std::vector<double> nu(nc,0.);
      for(const auto& e : grad[k]) nu[e.first]+=e.second;
      const double n0=sqrt(dot(nu.data(),nu.data(),nc));
      if(n0==0) continue;
      for(int pass=0; pass<2; pass++)
        for(const auto& v : nulls)
          {
            const double c=dot(nu.data(),v.data(),nc);
            for(long r=0; r<nc; r++) nu[r]-=c*v[r];
          }
      const double n1=sqrt(dot(nu.data(),nu.data(),nc));
      if(n1<1e-8*n0) continue;
      for(long r=0; r<nc; r++) nu[r]/=n1;
      nulls.push_back(nu);
    }
  for(const auto& v : nulls)
    {
      #pragma omp parallel for schedule(static)
      for(long r=0; r<nc; r++)
        if(v[r]!=0)
          for(long c=0; c<nc; c++) C[r*nc+c]+=scale*v[r]*v[c];
    }

  if(!spd_inverse(C,nc))
    {
      std::cout << "The Lagrange-force matrix is not positive definite: using the iterative solve.\n";
      direct=false;
      std::vector<double>().swap(Cinv);
      return;
    }
  const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
  std::cout << "Lagrange-force matrix (" << nulls.size() << " pressure modes that move no fluid) inverted in "
            << seconds << " s, " << nc*nc*8/1048576. << " MiB.\n";
}

//s[r][c] = sum over k<K of a[r][k] b[c][k], for 4 rows a and 4 rows b (register tile).
static inline void tile_dot(const double * const a[4], const double * const b[4], long K, double s[4][4])
{
  double t[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
  const double * a0=a[0], * a1=a[1], * a2=a[2], * a3=a[3], * b0=b[0], * b1=b[1], * b2=b[2], * b3=b[3];
  for(long k=0; k<K; k++)
    {
      const double x0=a0[k], x1=a1[k], x2=a2[k], x3=a3[k], y0=b0[k], y1=b1[k], y2=b2[k], y3=b3[k];
      t[0][0]+=x0*y0; t[0][1]+=x0*y1; t[0][2]+=x0*y2; t[0][3]+=x0*y3;
      t[1][0]+=x1*y0; t[1][1]+=x1*y1; t[1][2]+=x1*y2; t[1][3]+=x1*y3;
      t[2][0]+=x2*y0; t[2][1]+=x2*y1; t[2][2]+=x2*y2; t[2][3]+=x2*y3;
      t[3][0]+=x3*y0; t[3][1]+=x3*y1; t[3][2]+=x3*y2; t[3][3]+=x3*y3;
    }
  for(int r=0; r<4; r++) for(int c=0; c<4; c++) s[r][c]=t[r][c];
}

//Rows r0..r0+3 of A from column k0 (rows past n repeat the last one; their results are ignored).
static inline void tile_rows(const double * A, long n, long r0, long k0, const double * p[4])
{
  for(int q=0; q<4; q++) p[q]=&A[std::min(r0+q,n-1)*n+k0];
}

bool spd_inverse(double * A, long n)
{
  const long nb=64;
  std::vector<double> diag(n);
  double dmax=0;
  for(long i=0; i<n; i++) dmax=std::max(dmax,A[i*n+i]);

  //1) Cholesky A = L L^T, left-looking by blocks of nb columns: strictly lower part in A, diagonal
  //   in diag. For each block J: subtract L[i][0:j0] L[J][0:j0]^T (tiles), then factor the block.
  for(long j0=0; j0<n; j0+=nb)
    {
      const long j1=std::min(n,j0+nb);
      if(j0>0)
        {
          #pragma omp parallel for schedule(dynamic)
          for(long i0=j0; i0<n; i0+=4)
            for(long jt=j0; jt<j1 && jt<=i0+3; jt+=4)
              {
                const double * a[4], * b[4];
                double t[4][4];
                tile_rows(A,n,i0,0,a); tile_rows(A,n,jt,0,b);
                tile_dot(a,b,j0,t);
                for(int r=0; r<4 && i0+r<n; r++)
                  for(int c=0; c<4 && jt+c<j1; c++) A[(i0+r)*n+jt+c]-=t[r][c];
              }
        }
      for(long i=j0; i<j1; i++)
        for(long j=j0; j<=i; j++)
          {
            const double v=A[i*n+j]-dot(&A[i*n+j0],&A[j*n+j0],j-j0);
            if(i==j)
              {
                if(!(v>1e-12*dmax)) return false;
                diag[i]=sqrt(v);
              }
            else A[i*n+j]=v/diag[j];
          }
      #pragma omp parallel for schedule(static)
      for(long i=j1; i<n; i++)
        for(long j=j0; j<j1; j++) A[i*n+j]=(A[i*n+j]-dot(&A[i*n+j0],&A[j*n+j0],j-j0))/diag[j];
    }

  //2) U = L^-T in the upper part, diagonal included: row j of U is column j of L^-1,
  //   U_ji = -(sum over k in [j,i) of L_ik U_jk)/L_ii. Blocks of nb rows of U per thread; for tiles of
  //   4 rows i of L and 4 rows j of U the common range k in [jt+4, i0) goes through tile_dot.
  #pragma omp parallel for schedule(dynamic)
  for(long j0=0; j0<n; j0+=nb)
    {
      const long j1=std::min(n,j0+nb);
      for(long j=j0; j<j1; j++) A[j*n+j]=1/diag[j];
      for(long i0=j0+1; i0<n; i0+=4)
        for(long jt=j0; jt<j1 && jt<i0+4; jt+=4)
          {
            const long kc=jt+4;
            double t[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
            if(kc<i0)
              {
                const double * a[4], * b[4];
                tile_rows(A,n,i0,kc,a); tile_rows(A,n,jt,kc,b);
                tile_dot(a,b,i0-kc,t);
              }
            for(int r=0; r<4 && i0+r<n; r++)
              {
                const long i=i0+r;
                for(int c=0; c<4 && jt+c<j1; c++)
                  {
                    const long j=jt+c;
                    if(j>=i) continue;
                    double v;
                    if(kc<i0) v=t[r][c]+dot(&A[i*n+j],&A[j*n+j],kc-j)+dot(&A[i*n+i0],&A[j*n+i0],i-i0);
                    else v=dot(&A[i*n+j],&A[j*n+j],i-j);
                    A[j*n+i]=-v/diag[i];
                  }
              }
          }
    }

  //3) A^-1 = U U^T: (i,j) for j <= i is the sum over k >= i of U_ik U_jk; tiles over the common range
  //   k >= it+4, then the few remaining terms. Lower part in A, diagonal in diag, then mirrored.
  #pragma omp parallel for schedule(dynamic)
  for(long b0=0; b0<n; b0+=nb)
    {
      const long b1=std::min(n,b0+nb);
      for(long jt=0; jt<b1; jt+=4)
        for(long it=b0; it<b1; it+=4)
          {
            if(jt>it+3) continue;
            const long kc=std::min(n,it+4);
            double t[4][4];
            const double * a[4], * b[4];
            tile_rows(A,n,it,kc,a); tile_rows(A,n,jt,kc,b);
            tile_dot(a,b,n-kc,t);
            for(int r=0; r<4 && it+r<b1; r++)
              {
                const long i=it+r;
                for(int c=0; c<4 && jt+c<=i; c++)
                  {
                    const long j=jt+c;
                    const double v=t[r][c]+dot(&A[i*n+i],&A[j*n+i],kc-i);
                    if(i==j) diag[i]=v;
                    else A[i*n+j]=v;
                  }
              }
          }
    }
  #pragma omp parallel for schedule(static)
  for(long i=0; i<n; i++)
    {
      A[i*n+i]=diag[i];
      for(long j=0; j<i; j++) A[j*n+i]=A[i*n+j];
    }
  return true;
}
