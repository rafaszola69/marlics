#ifndef PHYSICS_FG_H
#define PHYSICS_FG_H
// Boundary_Fg functional derivatives and energy density,
// shared by the CPU and GPU code. T is the floating
// point type of the evaluation (double, or float for the GPU mixed precision). Formulas moved unchanged from boundary_fournier_galatola.cpp.
#include "physics.h"

#define QN00 QN[0] 
#define QN01 QN[1]
#define QN02 QN[2]
#define QN11 QN[3]
#define QN12 QN[4]
#define  Q_00_0 dQ[0]
#define  Q_01_0 dQ[1]
#define  Q_02_0 dQ[2]
#define  Q_11_0 dQ[3]
#define  Q_12_0 dQ[4]
#define  Q_00_1 dQ[5]
#define  Q_01_1 dQ[6]
#define  Q_02_1 dQ[7]
#define  Q_11_1 dQ[8]
#define  Q_12_1 dQ[9]
#define  Q_00_2 dQ[10]
#define  Q_01_2 dQ[11]
#define  Q_02_2 dQ[12]
#define  Q_11_2 dQ[13]
#define  Q_12_2 dQ[14]
#define  Q_00_00 ddQ[0]
#define  Q_01_00 ddQ[1]
#define  Q_02_00 ddQ[2]
#define  Q_11_00 ddQ[3]
#define  Q_12_00 ddQ[4]
#define  Q_00_01 ddQ[5]
#define  Q_01_01 ddQ[6]
#define  Q_02_01 ddQ[7]
#define  Q_11_01 ddQ[8]
#define  Q_12_01 ddQ[9]
#define  Q_00_02 ddQ[10]
#define  Q_01_02 ddQ[11]
#define  Q_02_02 ddQ[12]
#define  Q_11_02 ddQ[13]
#define  Q_12_02 ddQ[14]
#define  Q_00_11 ddQ[15]
#define  Q_01_11 ddQ[16]
#define  Q_02_11 ddQ[17]
#define  Q_11_11 ddQ[18]
#define  Q_12_11 ddQ[19]
#define  Q_00_12 ddQ[20]
#define  Q_01_12 ddQ[21]
#define  Q_02_12 ddQ[22]
#define  Q_11_12 ddQ[23]
#define  Q_12_12 ddQ[24]
#define  Q_00_22 ddQ[25]
#define  Q_01_22 ddQ[26]
#define  Q_02_22 ddQ[27]
#define  Q_11_22 ddQ[28]
#define  Q_12_22 ddQ[29]

template<class T> MARLICS_HD inline T fg_functional_derivative_00(const BulkParams& P, const BoundaryParams& B, const T QN[5], const T dQ[15], const T v[3])
{
  const T sigma=T(P.sigma), a=T(P.a), bb=T(P.bb), cc=T(P.cc), L1=T(P.L1), L2=T(P.L2), L3=T(P.L3), Lq=T(P.Lq), Ls=T(P.Ls), Lq_tilde=T(P.Lq_tilde), Lambda=T(P.Lambda), Lambda_s=T(P.Lambda_s), S_eq=T(P.S_eq), deltaepslon=T(P.deltaepslon), elecfieldx=T(P.elecfieldx), elecfieldy=T(P.elecfieldy), elecfieldz=T(P.elecfieldz);
  const T Wo1=T(B.Wo1);
 
 return Lambda_s*((-3*L1*Q_00_0*v[0] - 2*L2*Q_00_0*v[0] - 2*Ls*Q_00_0*v[0] - 3*L3*QN00*Q_00_0*v[0] - 3*L3*QN01*Q_00_1*v[0] - 3*L3*QN02*Q_00_2*v[0] - 2*L2*Q_01_1*v[0] + Ls*Q_01_1*v[0] - 2*L2*Q_02_2*v[0] + Ls*Q_02_2*v[0] - 8*S_eq*Wo1*(v[0]*v[0]*v[0]*v[0]) - 16*Wo1*QN00*(v[0]*v[0]*v[0]*v[0]) + 8*S_eq*Wo1*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0]) + 16*Wo1*QN00*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0]) - 2*S_eq*Wo1*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0]*v[0]*v[0]) - 4*Wo1*QN00*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0]*v[0]*v[0]) + 3*Lq_tilde*QN02*v[1] - 3*L3*QN01*Q_00_0*v[1] - 3*L1*Q_00_1*v[1] - 3*L3*QN11*Q_00_1*v[1] - 3*L3*QN12*Q_00_2*v[1] + L2*Q_01_0*v[1] - 2*Ls*Q_01_0*v[1] + L2*Q_11_1*v[1] + Ls*Q_11_1*v[1] + L2*Q_12_2*v[1] + Ls*Q_12_2*v[1] - 20*Wo1*QN01*(v[0]*v[0]*v[0])*v[1] + 32*Wo1*QN01*(v[0]*v[0]*v[0]*v[0]*v[0])*v[1] - 8*Wo1*QN01*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0]*v[0])*v[1] - 4*S_eq*Wo1*v[0]*v[0]*v[1]*v[1] - 4*Wo1*QN00*v[0]*v[0]*v[1]*v[1] - 4*Wo1*QN11*v[0]*v[0]*v[1]*v[1] + 12*S_eq*Wo1*(v[0]*v[0]*v[0]*v[0])*v[1]*v[1] + 8*Wo1*QN00*(v[0]*v[0]*v[0]*v[0])*v[1]*v[1] + 16*Wo1*QN11*(v[0]*v[0]*v[0]*v[0])*v[1]*v[1] - 5*S_eq*Wo1*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0])*v[1]*v[1] - 6*Wo1*QN00*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0])*v[1]*v[1] - 4*Wo1*QN11*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0])*v[1]*v[1] + 4*Wo1*QN01*v[0]*(v[1]*v[1]*v[1]) + 16*Wo1*QN01*(v[0]*v[0]*v[0])*(v[1]*v[1]*v[1]) - 12*Wo1*QN01*(v[0]*v[0]*v[0]*v[0]*v[0])*(v[1]*v[1]*v[1]) + 4*S_eq*Wo1*(v[1]*v[1]*v[1]*v[1]) + 8*Wo1*QN11*(v[1]*v[1]*v[1]*v[1]) - 8*Wo1*QN00*v[0]*v[0]*(v[1]*v[1]*v[1]*v[1]) + 8*Wo1*QN11*v[0]*v[0]*(v[1]*v[1]*v[1]*v[1]) - 3*S_eq*Wo1*(v[0]*v[0]*v[0]*v[0])*(v[1]*v[1]*v[1]*v[1]) - 6*Wo1*QN11*(v[0]*v[0]*v[0]*v[0])*(v[1]*v[1]*v[1]*v[1]) - 16*Wo1*QN01*v[0]*(v[1]*v[1]*v[1]*v[1]*v[1]) - 4*S_eq*Wo1*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) - 8*Wo1*QN11*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) + S_eq*Wo1*v[0]*v[0]*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) + 2*Wo1*QN00*v[0]*v[0]*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) + 4*Wo1*QN01*v[0]*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) + S_eq*Wo1*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) + 2*Wo1*QN11*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) - (3*Lq_tilde*QN01 + 3*L3*QN02*Q_00_0 + 3*L3*QN12*Q_00_1 + (3*L1 + L2 + Ls)*Q_00_2 - 3*L3*(QN00 + QN11)*Q_00_2 - (L2 - 2*Ls)*Q_02_0 + (L2 + Ls)*(Q_11_2 - Q_12_1) + 4*Wo1*(-(QN12*v[1]*(-2*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0]) + (v[0]*v[0]*v[0]*v[0])*(8 - 3*v[1]*v[1]) + v[1]*v[1]*pow(-2 + v[1]*v[1],2) + v[0]*v[0]*(-2 + 4*v[1]*v[1]))) + QN02*v[0]*(2*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0]) + v[0]*v[0]*(5 - 4*v[1]*v[1]) + (v[0]*v[0]*v[0]*v[0])*(-8 + 3*v[1]*v[1]) - v[1]*v[1]*(1 - 4*v[1]*v[1] + (v[1]*v[1]*v[1]*v[1])))))*v[2] + Wo1*(-(S_eq*(-2 + v[0]*v[0] + v[1]*v[1])*(5*(v[0]*v[0]*v[0]*v[0]) + 4*v[1]*v[1] - 4*(v[1]*v[1]*v[1]*v[1]) + v[0]*v[0]*(-2 + v[1]*v[1]))) + 2*(2*QN01*v[0]*v[1]*(1 + 4*v[0]*v[0] - 3*(v[0]*v[0]*v[0]*v[0]) - 8*v[1]*v[1] + 3*(v[1]*v[1]*v[1]*v[1])) + 2*QN11*(v[0]*v[0] - 4*(v[0]*v[0]*v[0]*v[0]) + (v[0]*v[0]*v[0]*v[0]*v[0]*v[0]) - 2*(v[1]*v[1]*v[1]*v[1]) + (v[1]*v[1]*v[1]*v[1]*v[1]*v[1])) - QN00*((v[0]*v[0]*v[0]*v[0]*v[0]*v[0]) + (v[0]*v[0]*v[0]*v[0])*(4 - 3*v[1]*v[1]) - 3*v[0]*v[0]*v[1]*v[1]*(-4 + v[1]*v[1]) + v[1]*v[1]*pow(-2 + v[1]*v[1],2))))*v[2]*v[2] - 4*Wo1*(QN12*v[1]*(2 + 3*v[0]*v[0] - 3*v[1]*v[1])*(-2 + v[0]*v[0] + v[1]*v[1]) + QN02*v[0]*(-1 - 4*v[0]*v[0] + 3*(v[0]*v[0]*v[0]*v[0]) + 8*v[1]*v[1] - 3*(v[1]*v[1]*v[1]*v[1])))*(v[2]*v[2]*v[2]) + Wo1*(S_eq*(4 - 3*(v[0]*v[0]*v[0]*v[0]) + 3*(-4 + v[0]*v[0])*v[1]*v[1] + 6*(v[1]*v[1]*v[1]*v[1])) + 2*(2*QN01*v[0]*v[1]*(-4 + 3*v[1]*v[1]) + QN11*(-4 - 4*v[0]*v[0] + 3*(v[0]*v[0]*v[0]*v[0]) + 4*v[1]*v[1]) + QN00*(-4 + 3*(v[0]*v[0]*v[0]*v[0]) + 8*v[1]*v[1] - 3*(v[1]*v[1]*v[1]*v[1]) + v[0]*v[0]*(-8 + 3*v[1]*v[1]))))*(v[2]*v[2]*v[2]*v[2]) + 4*Wo1*(QN02*v[0] + QN12*v[1])*(-4 + 3*v[1]*v[1])*(v[2]*v[2]*v[2]*v[2]*v[2]) + Wo1*(4*QN01*v[0]*v[1] + 2*QN00*(4 + v[0]*v[0] - 3*v[1]*v[1]) - 4*QN11*(-2 + v[1]*v[1]) + S_eq*(-4 + v[0]*v[0] + 4*v[1]*v[1]))*(v[2]*v[2]*v[2]*v[2]*v[2]*v[2]) + 4*Wo1*(QN02*v[0] + QN12*v[1])*(v[2]*v[2]*v[2]*v[2]*v[2]*v[2]*v[2]) + Wo1*(S_eq - 2*(QN00 + QN11))*(v[2]*v[2]*v[2]*v[2]*v[2]*v[2]*v[2]*v[2]))/T(3.));
}

template<class T> MARLICS_HD inline T fg_functional_derivative_01(const BulkParams& P, const BoundaryParams& B, const T QN[5], const T dQ[15], const T v[3])
{
  const T sigma=T(P.sigma), a=T(P.a), bb=T(P.bb), cc=T(P.cc), L1=T(P.L1), L2=T(P.L2), L3=T(P.L3), Lq=T(P.Lq), Ls=T(P.Ls), Lq_tilde=T(P.Lq_tilde), Lambda=T(P.Lambda), Lambda_s=T(P.Lambda_s), S_eq=T(P.S_eq), deltaepslon=T(P.deltaepslon), elecfieldx=T(P.elecfieldx), elecfieldy=T(P.elecfieldy), elecfieldz=T(P.elecfieldz);
  const T Wo1=T(B.Wo1);
 
 return Lambda_s*((-(Lq_tilde*QN02*v[0]) - Ls*Q_00_1*v[0] - 2*L1*Q_01_0*v[0] - L2*Q_01_0*v[0] - Ls*Q_01_0*v[0] - 2*L3*QN00*Q_01_0*v[0] - 2*L3*QN01*Q_01_1*v[0] - 2*L3*QN02*Q_01_2*v[0] - L2*Q_11_1*v[0] - L2*Q_12_2*v[0] - 4*Wo1*QN01*(v[0]*v[0]*v[0]*v[0]) + Lq_tilde*QN12*v[1] - L2*Q_00_0*v[1] - 2*L3*QN01*Q_01_0*v[1] - 2*L1*Q_01_1*v[1] - L2*Q_01_1*v[1] - Ls*Q_01_1*v[1] - 2*L3*QN11*Q_01_1*v[1] - 2*L3*QN12*Q_01_2*v[1] - L2*Q_02_2*v[1] - Ls*Q_11_0*v[1] - 8*S_eq*Wo1*(v[0]*v[0]*v[0])*v[1] - 12*Wo1*QN00*(v[0]*v[0]*v[0])*v[1] - 4*Wo1*QN11*(v[0]*v[0]*v[0])*v[1] + 8*S_eq*Wo1*(v[0]*v[0]*v[0]*v[0]*v[0])*v[1] + 16*Wo1*QN00*(v[0]*v[0]*v[0]*v[0]*v[0])*v[1] - 2*S_eq*Wo1*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0]*v[0])*v[1] - 4*Wo1*QN00*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0]*v[0])*v[1] - 24*Wo1*QN01*v[0]*v[0]*v[1]*v[1] + 32*Wo1*QN01*(v[0]*v[0]*v[0]*v[0])*v[1]*v[1] - 8*Wo1*QN01*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0])*v[1]*v[1] - 8*S_eq*Wo1*v[0]*(v[1]*v[1]*v[1]) - 4*Wo1*QN00*v[0]*(v[1]*v[1]*v[1]) - 12*Wo1*QN11*v[0]*(v[1]*v[1]*v[1]) + 16*S_eq*Wo1*(v[0]*v[0]*v[0])*(v[1]*v[1]*v[1]) + 16*Wo1*QN00*(v[0]*v[0]*v[0])*(v[1]*v[1]*v[1]) + 16*Wo1*QN11*(v[0]*v[0]*v[0])*(v[1]*v[1]*v[1]) - 6*S_eq*Wo1*(v[0]*v[0]*v[0]*v[0]*v[0])*(v[1]*v[1]*v[1]) - 8*Wo1*QN00*(v[0]*v[0]*v[0]*v[0]*v[0])*(v[1]*v[1]*v[1]) - 4*Wo1*QN11*(v[0]*v[0]*v[0]*v[0]*v[0])*(v[1]*v[1]*v[1]) - 4*Wo1*QN01*(v[1]*v[1]*v[1]*v[1]) + 32*Wo1*QN01*v[0]*v[0]*(v[1]*v[1]*v[1]*v[1]) - 16*Wo1*QN01*(v[0]*v[0]*v[0]*v[0])*(v[1]*v[1]*v[1]*v[1]) + 8*S_eq*Wo1*v[0]*(v[1]*v[1]*v[1]*v[1]*v[1]) + 16*Wo1*QN11*v[0]*(v[1]*v[1]*v[1]*v[1]*v[1]) - 6*S_eq*Wo1*(v[0]*v[0]*v[0])*(v[1]*v[1]*v[1]*v[1]*v[1]) - 4*Wo1*QN00*(v[0]*v[0]*v[0])*(v[1]*v[1]*v[1]*v[1]*v[1]) - 8*Wo1*QN11*(v[0]*v[0]*v[0])*(v[1]*v[1]*v[1]*v[1]*v[1]) - 8*Wo1*QN01*v[0]*v[0]*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) - 2*S_eq*Wo1*v[0]*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) - 4*Wo1*QN11*v[0]*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) + (Lq_tilde*(QN00 - QN11) - 2*L3*QN02*Q_01_0 - 2*L3*QN12*Q_01_1 + 2*(-L1 + L3*(QN00 + QN11))*Q_01_2 - Ls*(Q_02_1 + Q_12_0) + 4*Wo1*(-(QN12*v[0]*(v[0]*v[0] + (5 - 8*v[0]*v[0] + 2*(v[0]*v[0]*v[0]*v[0]))*v[1]*v[1] + 4*(-2 + v[0]*v[0])*(v[1]*v[1]*v[1]*v[1]) + 2*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]))) - QN02*v[1]*(v[1]*v[1] + v[0]*v[0]*(5 - 8*v[1]*v[1] + 2*((v[0]*v[0]*v[0]*v[0]) + (v[1]*v[1]*v[1]*v[1]) + 2*v[0]*v[0]*(-2 + v[1]*v[1]))))))*v[2] - 2*Wo1*(2*QN01*(v[1]*v[1] + v[0]*v[0]*(1 + 4*v[1]*v[1]*(-2 + v[0]*v[0] + v[1]*v[1]))) + v[0]*v[1]*(S_eq*(-2 + v[0]*v[0] + v[1]*v[1])*(-2 + 3*v[0]*v[0] + 3*v[1]*v[1]) + 2*(QN00*(-1 + (v[0]*v[0]*v[0]*v[0]) + 4*v[1]*v[1] - (v[1]*v[1]*v[1]*v[1])) + QN11*(-1 + 4*v[0]*v[0] - (v[0]*v[0]*v[0]*v[0]) + (v[1]*v[1]*v[1]*v[1])))))*v[2]*v[2] - 4*Wo1*(QN02*v[1]*(1 + 4*v[0]*v[0]*(-2 + v[0]*v[0] + v[1]*v[1])) + QN12*v[0]*(1 + 4*v[1]*v[1]*(-2 + v[0]*v[0] + v[1]*v[1])))*(v[2]*v[2]*v[2]) + 2*Wo1*v[0]*v[1]*(S_eq*(4 - 3*v[0]*v[0] - 3*v[1]*v[1]) + 2*(-2*QN01*v[0]*v[1] + QN11*(-4 + 2*v[0]*v[0] + v[1]*v[1]) + QN00*(-4 + v[0]*v[0] + 2*v[1]*v[1])))*(v[2]*v[2]*v[2]*v[2]) - 8*Wo1*v[0]*v[1]*(QN02*v[0] + QN12*v[1])*(v[2]*v[2]*v[2]*v[2]*v[2]) - 2*Wo1*(S_eq - 2*(QN00 + QN11))*v[0]*v[1]*(v[2]*v[2]*v[2]*v[2]*v[2]*v[2]))/T(2.));
}

template<class T> MARLICS_HD inline T fg_functional_derivative_02(const BulkParams& P, const BoundaryParams& B, const T QN[5], const T dQ[15], const T v[3])
{
  const T sigma=T(P.sigma), a=T(P.a), bb=T(P.bb), cc=T(P.cc), L1=T(P.L1), L2=T(P.L2), L3=T(P.L3), Lq=T(P.Lq), Ls=T(P.Ls), Lq_tilde=T(P.Lq_tilde), Lambda=T(P.Lambda), Lambda_s=T(P.Lambda_s), S_eq=T(P.S_eq), deltaepslon=T(P.deltaepslon), elecfieldx=T(P.elecfieldx), elecfieldy=T(P.elecfieldy), elecfieldz=T(P.elecfieldz);
  const T Wo1=T(B.Wo1);
 
 return Lambda_s*((Lq_tilde*QN01*v[0] + L2*Q_00_2*v[0] - Ls*Q_00_2*v[0] - 2*L1*Q_02_0*v[0] - L2*Q_02_0*v[0] - Ls*Q_02_0*v[0] - 2*L3*QN00*Q_02_0*v[0] - 2*L3*QN01*Q_02_1*v[0] - 2*L3*QN02*Q_02_2*v[0] + L2*Q_11_2*v[0] - L2*Q_12_1*v[0] - 4*Wo1*QN02*(v[0]*v[0]*v[0]*v[0]) - 2*Lq_tilde*QN00*v[1] - Lq_tilde*QN11*v[1] - Ls*Q_01_2*v[1] - 2*L3*QN01*Q_02_0*v[1] - 2*L1*Q_02_1*v[1] - 2*L3*QN11*Q_02_1*v[1] - 2*L3*QN12*Q_02_2*v[1] - Ls*Q_12_0*v[1] - 4*Wo1*QN12*(v[0]*v[0]*v[0])*v[1] - 4*Wo1*QN02*v[0]*v[0]*v[1]*v[1] - 4*Wo1*QN12*v[0]*(v[1]*v[1]*v[1]) - (Lq_tilde*QN12 + (L2 - Ls)*Q_00_0 + L2*Q_01_1 + (2*L1 + L2 + Ls)*Q_02_2 - Ls*Q_11_0 + 8*S_eq*Wo1*(v[0]*v[0]*v[0]) + 2*(L3*QN02*Q_02_0 + L3*QN12*Q_02_1 - L3*QN00*Q_02_2 - L3*QN11*Q_02_2 + 4*Wo1*QN00*(v[0]*v[0]*v[0]) - 2*Wo1*QN11*(v[0]*v[0]*v[0]) - 4*S_eq*Wo1*(v[0]*v[0]*v[0]*v[0]*v[0]) - 8*Wo1*QN00*(v[0]*v[0]*v[0]*v[0]*v[0]) + S_eq*Wo1*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0]*v[0]) + 2*Wo1*QN00*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0]*v[0]) + 2*Wo1*QN01*v[0]*v[0]*(5 - 8*v[0]*v[0] + 2*(v[0]*v[0]*v[0]*v[0]))*v[1] + Wo1*v[0]*(4*S_eq + 2*QN11 - 8*(S_eq + QN00 + QN11)*v[0]*v[0] + (3*S_eq + 4*QN00 + 2*QN11)*(v[0]*v[0]*v[0]*v[0]))*v[1]*v[1] + 2*Wo1*QN01*(1 - 8*v[0]*v[0] + 4*(v[0]*v[0]*v[0]*v[0]))*(v[1]*v[1]*v[1]) + Wo1*v[0]*(-4*(S_eq + 2*QN11) + (3*S_eq + 2*QN00 + 4*QN11)*v[0]*v[0])*(v[1]*v[1]*v[1]*v[1]) + 4*Wo1*QN01*v[0]*v[0]*(v[1]*v[1]*v[1]*v[1]*v[1]) + Wo1*(S_eq + 2*QN11)*v[0]*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1])))*v[2] - 4*Wo1*(QN02*(v[1]*v[1] + 2*v[0]*v[0]*(-3 + v[0]*v[0] + v[1]*v[1])*(-1 + v[0]*v[0] + v[1]*v[1])) + QN12*v[0]*v[1]*(5 - 8*v[1]*v[1] + 2*((v[0]*v[0]*v[0]*v[0]) + (v[1]*v[1]*v[1]*v[1]) + 2*v[0]*v[0]*(-2 + v[1]*v[1]))))*v[2]*v[2] - 2*Wo1*(S_eq*v[0]*(-2 + v[0]*v[0] + v[1]*v[1])*(-2 + 3*v[0]*v[0] + 3*v[1]*v[1]) + 2*(QN00*v[0]*(-2 + (v[0]*v[0]*v[0]*v[0]) + 4*v[1]*v[1] - (v[1]*v[1]*v[1]*v[1])) + QN11*v[0]*(-3 + 4*v[0]*v[0] - (v[0]*v[0]*v[0]*v[0]) + (v[1]*v[1]*v[1]*v[1])) + QN01*v[1]*(1 + 4*v[0]*v[0]*(-2 + v[0]*v[0] + v[1]*v[1]))))*(v[2]*v[2]*v[2]) - 4*Wo1*(4*QN12*v[0]*v[1]*(-2 + v[0]*v[0] + v[1]*v[1]) + QN02*(1 + 4*v[0]*v[0]*(-2 + v[0]*v[0] + v[1]*v[1])))*(v[2]*v[2]*v[2]*v[2]) + 2*Wo1*v[0]*(S_eq*(4 - 3*v[0]*v[0] - 3*v[1]*v[1]) + 2*(-2*QN01*v[0]*v[1] + QN11*(-4 + 2*v[0]*v[0] + v[1]*v[1]) + QN00*(-4 + v[0]*v[0] + 2*v[1]*v[1])))*(v[2]*v[2]*v[2]*v[2]*v[2]) - 8*Wo1*v[0]*(QN02*v[0] + QN12*v[1])*(v[2]*v[2]*v[2]*v[2]*v[2]*v[2]) - 2*Wo1*(S_eq - 2*(QN00 + QN11))*v[0]*(v[2]*v[2]*v[2]*v[2]*v[2]*v[2]*v[2]))/T(2.));
}

template<class T> MARLICS_HD inline T fg_functional_derivative_11(const BulkParams& P, const BoundaryParams& B, const T QN[5], const T dQ[15], const T v[3])
{
  const T sigma=T(P.sigma), a=T(P.a), bb=T(P.bb), cc=T(P.cc), L1=T(P.L1), L2=T(P.L2), L3=T(P.L3), Lq=T(P.Lq), Ls=T(P.Ls), Lq_tilde=T(P.Lq_tilde), Lambda=T(P.Lambda), Lambda_s=T(P.Lambda_s), S_eq=T(P.S_eq), deltaepslon=T(P.deltaepslon), elecfieldx=T(P.elecfieldx), elecfieldy=T(P.elecfieldy), elecfieldz=T(P.elecfieldz);
  const T Wo1=T(B.Wo1);
 
 return Lambda_s*((-3*Lq_tilde*QN12*v[0] + L2*Q_00_0*v[0] + Ls*Q_00_0*v[0] + L2*Q_01_1*v[0] - 2*Ls*Q_01_1*v[0] + L2*Q_02_2*v[0] + Ls*Q_02_2*v[0] - 3*L1*Q_11_0*v[0] - 3*L3*QN00*Q_11_0*v[0] - 3*L3*QN01*Q_11_1*v[0] - 3*L3*QN02*Q_11_2*v[0] + 4*S_eq*Wo1*(v[0]*v[0]*v[0]*v[0]) + 8*Wo1*QN00*(v[0]*v[0]*v[0]*v[0]) - 4*S_eq*Wo1*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0]) - 8*Wo1*QN00*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0]) + S_eq*Wo1*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0]*v[0]*v[0]) + 2*Wo1*QN00*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0]*v[0]*v[0]) - 2*L2*Q_01_0*v[1] + Ls*Q_01_0*v[1] - 3*L3*QN01*Q_11_0*v[1] - 3*L1*Q_11_1*v[1] - 2*L2*Q_11_1*v[1] - 2*Ls*Q_11_1*v[1] - 3*L3*QN11*Q_11_1*v[1] - 3*L3*QN12*Q_11_2*v[1] - 2*L2*Q_12_2*v[1] + Ls*Q_12_2*v[1] + 4*Wo1*QN01*(v[0]*v[0]*v[0])*v[1] - 16*Wo1*QN01*(v[0]*v[0]*v[0]*v[0]*v[0])*v[1] + 4*Wo1*QN01*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0]*v[0])*v[1] - 4*S_eq*Wo1*v[0]*v[0]*v[1]*v[1] - 4*Wo1*QN00*v[0]*v[0]*v[1]*v[1] - 4*Wo1*QN11*v[0]*v[0]*v[1]*v[1] + 8*Wo1*QN00*(v[0]*v[0]*v[0]*v[0])*v[1]*v[1] - 8*Wo1*QN11*(v[0]*v[0]*v[0]*v[0])*v[1]*v[1] + S_eq*Wo1*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0])*v[1]*v[1] + 2*Wo1*QN11*(v[0]*v[0]*v[0]*v[0]*v[0]*v[0])*v[1]*v[1] - 20*Wo1*QN01*v[0]*(v[1]*v[1]*v[1]) + 16*Wo1*QN01*(v[0]*v[0]*v[0])*(v[1]*v[1]*v[1]) - 8*S_eq*Wo1*(v[1]*v[1]*v[1]*v[1]) - 16*Wo1*QN11*(v[1]*v[1]*v[1]*v[1]) + 12*S_eq*Wo1*v[0]*v[0]*(v[1]*v[1]*v[1]*v[1]) + 16*Wo1*QN00*v[0]*v[0]*(v[1]*v[1]*v[1]*v[1]) + 8*Wo1*QN11*v[0]*v[0]*(v[1]*v[1]*v[1]*v[1]) - 3*S_eq*Wo1*(v[0]*v[0]*v[0]*v[0])*(v[1]*v[1]*v[1]*v[1]) - 6*Wo1*QN00*(v[0]*v[0]*v[0]*v[0])*(v[1]*v[1]*v[1]*v[1]) + 32*Wo1*QN01*v[0]*(v[1]*v[1]*v[1]*v[1]*v[1]) - 12*Wo1*QN01*(v[0]*v[0]*v[0])*(v[1]*v[1]*v[1]*v[1]*v[1]) + 8*S_eq*Wo1*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) + 16*Wo1*QN11*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) - 5*S_eq*Wo1*v[0]*v[0]*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) - 4*Wo1*QN00*v[0]*v[0]*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) - 6*Wo1*QN11*v[0]*v[0]*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) - 8*Wo1*QN01*v[0]*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) - 2*S_eq*Wo1*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) - 4*Wo1*QN11*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) + (3*Lq_tilde*QN01 - (L2 + Ls)*Q_00_2 + (L2 + Ls)*Q_02_0 - 3*L3*QN02*Q_11_0 - 3*L3*QN12*Q_11_1 - 3*L1*Q_11_2 - (L2 + Ls - 3*L3*(QN00 + QN11))*Q_11_2 + (L2 - 2*Ls)*Q_12_1 + 4*Wo1*(QN12*v[1]*(v[0]*v[0] - 4*(v[0]*v[0]*v[0]*v[0]) + (v[0]*v[0]*v[0]*v[0]*v[0]*v[0]) + (-5 + 4*v[0]*v[0])*v[1]*v[1] + (8 - 3*v[0]*v[0])*(v[1]*v[1]*v[1]*v[1]) - 2*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1])) + QN02*v[0]*(-4*(v[0]*v[0]*v[0]*v[0]) + (v[0]*v[0]*v[0]*v[0]*v[0]*v[0]) + v[0]*v[0]*(4 + 4*v[1]*v[1] - 3*(v[1]*v[1]*v[1]*v[1])) - 2*(v[1]*v[1] - 4*(v[1]*v[1]*v[1]*v[1]) + (v[1]*v[1]*v[1]*v[1]*v[1]*v[1])))))*v[2] + Wo1*(S_eq*(-2 + v[0]*v[0] + v[1]*v[1])*(4*(v[0]*v[0]*v[0]*v[0]) + 2*v[1]*v[1] - 5*(v[1]*v[1]*v[1]*v[1]) - v[0]*v[0]*(4 + v[1]*v[1])) + 2*(QN11*(-pow(-2*v[0] + (v[0]*v[0]*v[0]),2) + 3*v[0]*v[0]*(-4 + v[0]*v[0])*v[1]*v[1] + (-4 + 3*v[0]*v[0])*(v[1]*v[1]*v[1]*v[1]) - (v[1]*v[1]*v[1]*v[1]*v[1]*v[1])) + 2*(QN01*v[0]*v[1]*(1 - 8*v[0]*v[0] + 3*(v[0]*v[0]*v[0]*v[0]) + 4*v[1]*v[1] - 3*(v[1]*v[1]*v[1]*v[1])) + QN00*(-2*(v[0]*v[0]*v[0]*v[0]) + (v[0]*v[0]*v[0]*v[0]*v[0]*v[0]) + v[1]*v[1] - 4*(v[1]*v[1]*v[1]*v[1]) + (v[1]*v[1]*v[1]*v[1]*v[1]*v[1])))))*v[2]*v[2] + 4*Wo1*(QN02*v[0]*(-2 + 3*v[0]*v[0] - 3*v[1]*v[1])*(-2 + v[0]*v[0] + v[1]*v[1]) + QN12*v[1]*(1 - 8*v[0]*v[0] + 3*(v[0]*v[0]*v[0]*v[0]) + 4*v[1]*v[1] - 3*(v[1]*v[1]*v[1]*v[1])))*(v[2]*v[2]*v[2]) + Wo1*(4*S_eq - 8*(QN00 + QN11) + 4*(-3*S_eq + 2*QN00 + 4*QN11)*v[0]*v[0] + 6*(S_eq - QN11)*(v[0]*v[0]*v[0]*v[0]) + 4*QN01*v[0]*(-4 + 3*v[0]*v[0])*v[1] + (-8*(QN00 + 2*QN11) + 3*(S_eq + 2*QN11)*v[0]*v[0])*v[1]*v[1] - 3*(S_eq - 2*(QN00 + QN11))*(v[1]*v[1]*v[1]*v[1]))*(v[2]*v[2]*v[2]*v[2]) + 4*Wo1*(-4 + 3*v[0]*v[0])*(QN02*v[0] + QN12*v[1])*(v[2]*v[2]*v[2]*v[2]*v[2]) + Wo1*(8*(QN00 + QN11) - 2*(2*QN00 + 3*QN11)*v[0]*v[0] + 4*QN01*v[0]*v[1] + 2*QN11*v[1]*v[1] + S_eq*(-4 + 4*v[0]*v[0] + v[1]*v[1]))*(v[2]*v[2]*v[2]*v[2]*v[2]*v[2]) + 4*Wo1*(QN02*v[0] + QN12*v[1])*(v[2]*v[2]*v[2]*v[2]*v[2]*v[2]*v[2]) + Wo1*(S_eq - 2*(QN00 + QN11))*(v[2]*v[2]*v[2]*v[2]*v[2]*v[2]*v[2]*v[2]))/T(3.));
}

template<class T> MARLICS_HD inline T fg_functional_derivative_12(const BulkParams& P, const BoundaryParams& B, const T QN[5], const T dQ[15], const T v[3])
{
  const T sigma=T(P.sigma), a=T(P.a), bb=T(P.bb), cc=T(P.cc), L1=T(P.L1), L2=T(P.L2), L3=T(P.L3), Lq=T(P.Lq), Ls=T(P.Ls), Lq_tilde=T(P.Lq_tilde), Lambda=T(P.Lambda), Lambda_s=T(P.Lambda_s), S_eq=T(P.S_eq), deltaepslon=T(P.deltaepslon), elecfieldx=T(P.elecfieldx), elecfieldy=T(P.elecfieldy), elecfieldz=T(P.elecfieldz);
  const T Wo1=T(B.Wo1);
 
 return Lambda_s*((Lq_tilde*(QN00 + 2*QN11)*v[0] - Ls*Q_01_2*v[0] - Ls*Q_02_1*v[0] - 2*L1*Q_12_0*v[0] - 2*L3*QN00*Q_12_0*v[0] - 2*L3*QN01*Q_12_1*v[0] - 2*L3*QN02*Q_12_2*v[0] - Lq_tilde*QN01*v[1] + L2*Q_00_2*v[1] - L2*Q_02_0*v[1] + L2*Q_11_2*v[1] - Ls*Q_11_2*v[1] - 2*L3*QN01*Q_12_0*v[1] - 2*L1*Q_12_1*v[1] - L2*Q_12_1*v[1] - Ls*Q_12_1*v[1] - 2*L3*QN11*Q_12_1*v[1] - 2*L3*QN12*Q_12_2*v[1] - 4*Wo1*QN02*(v[0]*v[0]*v[0])*v[1] - 4*Wo1*QN12*v[0]*v[0]*v[1]*v[1] - 4*Wo1*QN02*v[0]*(v[1]*v[1]*v[1]) - 4*Wo1*QN12*(v[1]*v[1]*v[1]*v[1]) + (Lq_tilde*QN02 + Ls*Q_00_1 + Ls*Q_11_1 - 2*L3*(QN02*Q_12_0 + QN12*Q_12_1) - 2*L1  *Q_12_2 + (-Ls + 2*L3*(QN00 + QN11))*Q_12_2 - L2*(Q_01_0 + Q_11_1 + Q_12_2) - 2*Wo1*v[0]*v[0]*(S_eq*v[0]*v[0]*(-4 + v[0]*v[0]) + 2*QN00*(1 - 4*v[0]*v[0] + (v[0]*v[0]*v[0]*v[0])))*v[1] - 4*Wo1*QN01*v[0]*(5 - 8*v[0]*v[0] + 2*(v[0]*v[0]*v[0]*v[0]))*v[1]*v[1] - 2*Wo1*(4*S_eq - 2*QN00 + 4*QN11 - 8*(S_eq + QN00 + QN11)*v[0]*v[0] + (3*S_eq + 4*QN00 + 2*QN11)*(v[0]*v[0]*v[0]*v[0]))*(v[1]*v[1]*v[1]) - 16*Wo1*QN01*v[0]*(-2 + v[0]*v[0])*(v[1]*v[1]*v[1]*v[1]) - 2*Wo1*(-4*(S_eq + 2*QN11) + (3*S_eq + 2*QN00 + 4*QN11)*v[0]*v[0])*(v[1]*v[1]*v[1]*v[1]*v[1]) - 8*Wo1*QN01*v[0]*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) - 2*Wo1*(S_eq + 2*QN11)*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1]*v[1]) - 4*Wo1*v[0]*v[0]*(QN01*v[0] + 2*S_eq*v[1]))*v[2] - 4*Wo1*(QN12*(v[0]*v[0] + 2*(3 - 4*v[0]*v[0] + (v[0]*v[0]*v[0]*v[0]))*v[1]*v[1] + 4*(-2 + v[0]*v[0])*(v[1]*v[1]*v[1]*v[1]) + 2*(v[1]*v[1]*v[1]*v[1]*v[1]*v[1])) + QN02*v[0]*v[1]*(5 - 8*v[1]*v[1] + 2*((v[0]*v[0]*v[0]*v[0]) + (v[1]*v[1]*v[1]*v[1]) + 2*v[0]*v[0]*(-2 + v[1]*v[1]))))*v[2]*v[2] - 2*Wo1*(2*QN01*v[0]*(1 + 4*v[1]*v[1]*(-2 + v[0]*v[0] + v[1]*v[1])) + v[1]*(S_eq*(-2 + v[0]*v[0] + v[1]*v[1])*(-2 + 3*v[0]*v[0] + 3*v[1]*v[1]) + 2*(QN00*(-3 + (v[0]*v[0]*v[0]*v[0]) + 4*v[1]*v[1] - (v[1]*v[1]*v[1]*v[1])) + QN11*(-2 + 4*v[0]*v[0] - (v[0]*v[0]*v[0]*v[0]) + (v[1]*v[1]*v[1]*v[1])))))*(v[2]*v[2]*v[2]) - 4*Wo1*(4*QN02*v[0]*v[1]*(-2 + v[0]*v[0] + v[1]*v[1]) + QN12*(1 + 4*v[1]*v[1]*(-2 + v[0]*v[0] + v[1]*v[1])))*(v[2]*v[2]*v[2]*v[2]) + 2*Wo1*v[1]*(S_eq*(4 - 3*v[0]*v[0] - 3*v[1]*v[1]) + 2*(-2*QN01*v[0]*v[1] + QN11*(-4 + 2*v[0]*v[0] + v[1]*v[1]) + QN00*(-4 + v[0]*v[0] + 2*v[1]*v[1])))*(v[2]*v[2]*v[2]*v[2]*v[2]) - 8*Wo1*v[1]*(QN02*v[0] + QN12*v[1])*(v[2]*v[2]*v[2]*v[2]*v[2]*v[2]) - 2*Wo1*(S_eq - 2*(QN00 + QN11))*v[1]*(v[2]*v[2]*v[2]*v[2]*v[2]*v[2]*v[2]))/T(2.));
}

template<class T> MARLICS_HD inline T fg_energy_calculation(const BulkParams& P, const BoundaryParams& B, const T QN[5], const T dQ[15], const T v[3])
{
  const T sigma=T(P.sigma), a=T(P.a), bb=T(P.bb), cc=T(P.cc), L1=T(P.L1), L2=T(P.L2), L3=T(P.L3), Lq=T(P.Lq), Ls=T(P.Ls), Lq_tilde=T(P.Lq_tilde), Lambda=T(P.Lambda), Lambda_s=T(P.Lambda_s), S_eq=T(P.S_eq), deltaepslon=T(P.deltaepslon), elecfieldx=T(P.elecfieldx), elecfieldy=T(P.elecfieldy), elecfieldz=T(P.elecfieldz);
  const T Wo1=T(B.Wo1);
 
 return 2*((Wo1*(v[1]*v[1]*pow(4*QN01*v[0]*(-1 + v[1]*v[1]) + v[1]*(S_eq*(-2 + v[0]*v[0] + v[1]*v[1]) + 2*(QN00*v[0]*v[0] + QN11*(-2 + v[1]*v[1]))) + 4*(QN02*v[0]*v[1] + QN12*(-1 + v[1]*v[1]))*v[2] + (S_eq - 2*(QN00 + QN11))*v[1]*v[2]*v[2],2) + 2*pow(QN01*(-2*v[1]*v[1] + v[0]*v[0]*(-2 + 4*v[1]*v[1])) + v[0]*v[1]*(S_eq*(-2 + v[0]*v[0] + v[1]*v[1]) + 2*(QN00*(-1 + v[0]*v[0]) + QN11*(-1 + v[1]*v[1]))) + 2*(QN02*(-1 + 2*v[0]*v[0])*v[1] + QN12*v[0]*(-1 + 2*v[1]*v[1]))*v[2] + (S_eq - 2*(QN00 + QN11))*v[0]*v[1]*v[2]*v[2],2) + pow((S_eq - 2*(QN00 + QN11))*(-1 + pow(v[0]*v[0] + v[1]*v[1],2)) - 4*(QN02*v[0] + QN12*v[1])*(v[0]*v[0] + v[1]*v[1])*v[2] + ((S_eq + 2*QN00)*v[0]*v[0] + 4*QN01*v[0]*v[1] + (S_eq + 2*QN11)*v[1]*v[1])*v[2]*v[2],2) + 2*pow((2*QN01*v[1] + v[0]*(S_eq - 2*(2*QN01*v[0]*v[1] + QN00*(-1 + 2*v[0]*v[0] + v[1]*v[1]) + QN11*(v[0]*v[0] + 2*v[1]*v[1]))))*v[2] + 2*QN12*v[0]*v[1]*(v[0]*v[0] + v[1]*v[1] - v[2]*v[2]) + 2*QN02*(1 + (v[0]*v[0]*v[0]*v[0]) - v[1]*v[1] + v[0]*v[0]*(-1 + v[1]*v[1] - v[2]*v[2])),2) + 2*pow((2*QN01*v[0]*(1 - 2*v[1]*v[1]) + v[1]*(S_eq - 2*(QN00*(2*v[0]*v[0] + v[1]*v[1]) + QN11*(-1 + v[0]*v[0] + 2*v[1]*v[1]))))*v[2] + 2*QN02*v[0]*v[1]*(v[0]*v[0] + v[1]*v[1] - v[2]*v[2]) + 2*QN12*(1 + (v[1]*v[1]*v[1]*v[1]) + v[0]*v[0]*(-1 + v[1]*v[1]) - v[1]*v[1]*(1 + v[2]*v[2])),2) + v[0]*v[0]*pow(S_eq*v[0]*(-2 + v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) + 2*(2*QN01*(-1 + v[0]*v[0])*v[1] - 2*QN02*v[2] + QN00*v[0]*(-2 + v[0]*v[0] - v[2]*v[2]) + v[0]*(2*(QN02*v[0] + QN12*v[1])*v[2] + QN11*(v[1]*v[1] - v[2]*v[2]))),2)))/T(8.));
}

#undef QN00
#undef QN01
#undef QN02
#undef QN11
#undef QN12
#undef Q_00_0
#undef Q_01_0
#undef Q_02_0
#undef Q_11_0
#undef Q_12_0
#undef Q_00_1
#undef Q_01_1
#undef Q_02_1
#undef Q_11_1
#undef Q_12_1
#undef Q_00_2
#undef Q_01_2
#undef Q_02_2
#undef Q_11_2
#undef Q_12_2
#undef Q_00_00
#undef Q_01_00
#undef Q_02_00
#undef Q_11_00
#undef Q_12_00
#undef Q_00_01
#undef Q_01_01
#undef Q_02_01
#undef Q_11_01
#undef Q_12_01
#undef Q_00_02
#undef Q_01_02
#undef Q_02_02
#undef Q_11_02
#undef Q_12_02
#undef Q_00_11
#undef Q_01_11
#undef Q_02_11
#undef Q_11_11
#undef Q_12_11
#undef Q_00_12
#undef Q_01_12
#undef Q_02_12
#undef Q_11_12
#undef Q_12_12
#undef Q_00_22
#undef Q_01_22
#undef Q_02_22
#undef Q_11_22
#undef Q_12_22

#endif
