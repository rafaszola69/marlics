#ifndef FLOW_CONFINED_H
#define FLOW_CONFINED_H
// Stokes flow in any geometry with surface sites other than two flat walls (disk, box, droplets,
// custom geometries, ...). The discrete problem is the one of flow_walls.h:
//   velocity u at the bulk sites, zero at the others; pressure p at the bulk sites and at the
//   non-bulk sites next to them (E1);  eta A u + grad p = f at the bulk sites,  div u = 0 at the bulk
//   and E1 sites;  central differences and the 3-point Laplacian along each axis.
// The continuity equation at E1 keeps the fluid in (no flux into the walls), and since the pressure
// lives where the stress of boundary_stress does, a force that is a discrete gradient is balanced by
// the pressure and moves no fluid.
// Solution: the whole box is a periodic fluid, with Lagrange forces lambda on E2 (E1 and its non-bulk
// neighbours, two layers) that make u vanish there:
//   u = G (f + lambda) + U0,
// G the periodic Stokes solution (FFT, uniform mode removed: stokes_mode with checkerboards kept),
// U0 a uniform velocity. With u=0 on two layers the bulk equations hold exactly and the fluid outside
// stays at rest. lambda solves C lambda = -P g, g = (G f) on E2, C = P E G E^T P, P removing the mean
// of each component over E2 (absorbed by U0). C is symmetric positive semidefinite; its null vectors
// are pressure jumps across E2, which move no fluid.
//  - direct: C, completed on its null space, is inverted once (Cholesky); then each flow costs two
//    FFT solves and a matrix-vector product;
//  - iterative (large 3D droplets): conjugate gradients, with the Laplacian restricted to E2 as
//    preconditioner, from the previous lambda; each iteration is one FFT solve.
#include <complex>
#include <functional>
#include <vector>
#include "lattice.h"

//In-place FFT of howmany consecutive fields of a periodic Nx x Ny x Nz box (flow_cpu.cpp); the inverse
//is scaled by 1/(Nx Ny Nz).
void fft_periodic(std::complex<double> * a, int Nx, int Ny, int Nz, int howmany, bool forward);

//Inverse of a dense symmetric positive definite n x n matrix (row-major), in place, through its
//Cholesky factor (flow_confined.cpp, OpenMP). False if the matrix is not positive definite.
bool spd_inverse(double * A, long n);

class ConfinedStokes
{
 public:

  //solver: 0 automatic (direct when the matrix has at most max_direct rows), 1 direct, 2 iterative.
  ConfinedStokes(const LatticeView& L, double eta, int solver);

  //Whether the geometry has non-bulk sites (otherwise the periodic solver applies).
  static bool needed(const LatticeView& L);

  //Lagrange forces lambda (3 per E2 site, site-major) from g = (G f) on E2 and the total force f.
  //lambda = lambda0 + x: lambda0 uniform on E2 balances the total force (G drops the uniform mode of
  //the force, which would otherwise act as a uniform body force in the bulk), x the solve.
  //Cmul(p, out): out = E G E^T p (one FFT solve), for the iterative solve. Dmul(w, out): out = Cinv w,
  //for the direct solve (on the host when not given).
  typedef std::function<void(const double *, double *)> Apply;
  void solve(double * g, const double total[3], double * lambda, const Apply& Cmul, const Apply& Dmul=Apply());

  //Removes the mean of each component over E2.
  void project(double * v) const;

  std::vector<long> e2;        //E2 sites
  long ne;                     //number of E2 sites
  bool direct;
  std::vector<double> Cinv;    //direct: (3 ne)^2, row-major
  std::vector<double> g_uniform; //(G 1_E2 e_a) on E2, a=0,1,2: 3 x 3 ne
  int iterations;              //iterative: conjugate-gradient iterations of the last solve

  static const long max_direct=10000;

 private:

  void precondition(const double * r, double * z) const;
  void solve_projected(const double * w, double * lambda, const Apply& Cmul, const Apply& Dmul);
  void build_direct(const LatticeView& L, double eta);

  std::vector<long> e2_nb;     //E2 neighbours of each E2 site (6 per site, -1 outside E2)
  std::vector<double> d2;      //1/h^2 along each axis
  std::vector<double> last;    //iterative: lambda of the previous solve
};

#endif
