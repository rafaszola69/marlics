#ifndef PHYSICS_H
#define PHYSICS_H
// Plain parameter structs for the physics shared by the CPU and GPU code. They hold no pointers,
// so they can be copied to the GPU as they are.

#include <cmath>

#ifdef __CUDACC__
#define MARLICS_HD __host__ __device__
#else
#define MARLICS_HD
#endif

//Element-wise work of the integrators, run by an executor (exec_cpu.h, gpu/exec_cuda.h).
#define MARLICS_LAMBDA [=] MARLICS_HD

//Result of the executors' combined reduction (two sums and a maximum).
struct SumSumMax { double sum0, sum1, max; };

//Landau-de Gennes bulk parameters (see ENERGY for units).
struct BulkParams
{
  double sigma, a, bb, cc;
  double L1, L2, L3, Lq, Ls, Lq_tilde;
  double Lambda, Lambda_s;
  double S_eq;
  double deltaepslon, elecfieldx, elecfieldy, elecfieldz;
};

enum BoundaryKind { BC_STRONG=0, BC_RP=1, BC_FG=2, BC_HOMEOTROPIC=3 };

//One anchoring condition. Q0 is the preferred Q tensor (Rapini-Papoular only).
struct BoundaryParams
{
  int kind;
  double Wo1;
  double Q0[5];
};

#endif
