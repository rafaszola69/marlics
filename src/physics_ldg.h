#ifndef PHYSICS_LDG_H
#define PHYSICS_LDG_H
// Bulk Landau-de Gennes functional derivatives and energy density,
// shared by the CPU and GPU code. Formulas moved unchanged from energy_ldg.cpp.
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

MARLICS_HD inline double ldg_functional_derivative_00(const BulkParams& P, const double QN[5], const double dQ[15], const double ddQ[30], const double v[3])
{
  const double sigma=P.sigma, a=P.a, bb=P.bb, cc=P.cc, L1=P.L1, L2=P.L2, L3=P.L3, Lq=P.Lq, Ls=P.Ls, Lq_tilde=P.Lq_tilde, Lambda=P.Lambda, Lambda_s=P.Lambda_s, S_eq=P.S_eq, deltaepslon=P.deltaepslon, elecfieldx=P.elecfieldx, elecfieldy=P.elecfieldy, elecfieldz=P.elecfieldz;
 
    return Lambda*((-deltaepslon*(-2*(elecfieldx*elecfieldx) + (elecfieldy*elecfieldy) + (elecfieldz*elecfieldz)) - 3*sigma*QN00 - 6*cc*QN00*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) + QN00*QN11 + (QN11*QN11) + (QN12*QN12)) - bb*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) - 2*QN00*QN11 - 2*((QN11*QN11) + (QN12*QN12))) + 6*Lq_tilde*(Q_01_2 - Q_02_1) + (3*L1 + 2*(L2 + Ls))*Q_00_00 + 3*L1*Q_00_11 + 3*L1*Q_00_22 + L2*Q_00_22 + Ls*Q_00_22 + L3*((Q_00_0*Q_00_0) + (Q_00_1*Q_00_1) - 2*(Q_00_2*Q_00_2) - 2*(Q_01_0*Q_01_0) + (Q_01_1*Q_01_1) + (Q_01_2*Q_01_2) - 2*(Q_02_0*Q_02_0) + (Q_02_1*Q_02_1) + (Q_02_2*Q_02_2) + Q_00_0*(3*(Q_01_1 + Q_02_2) - 2*Q_11_0) - 2*(Q_11_0*Q_11_0) + (Q_11_1*Q_11_1) + (Q_11_2*Q_11_2) - 2*(Q_12_0*Q_12_0) + (Q_12_1*Q_12_1) + Q_00_2*(3*Q_02_0 - 2*Q_11_2 + 3*Q_12_1) + (Q_12_2*Q_12_2) + Q_00_1*(3*Q_01_0 + 4*Q_11_1 + 3*Q_12_2) + 3*QN00*Q_00_00 + 6*QN01*Q_00_01 + 6*QN02*Q_00_02 + 3*QN11*Q_00_11 + 6*QN12*Q_00_12 - 3*(QN00 + QN11)*Q_00_22) + L2*Q_01_01 + Ls*Q_01_01 + L2*Q_02_02 + Ls*Q_02_02 - L2*Q_11_11 - Ls*Q_11_11 + L2*Q_11_22 + Ls*Q_11_22 - 2*(L2 + Ls)*Q_12_12)/3.);
}

MARLICS_HD inline double ldg_functional_derivative_01(const BulkParams& P, const double QN[5], const double dQ[15], const double ddQ[30], const double v[3])
{
  const double sigma=P.sigma, a=P.a, bb=P.bb, cc=P.cc, L1=P.L1, L2=P.L2, L3=P.L3, Lq=P.Lq, Ls=P.Ls, Lq_tilde=P.Lq_tilde, Lambda=P.Lambda, Lambda_s=P.Lambda_s, S_eq=P.S_eq, deltaepslon=P.deltaepslon, elecfieldx=P.elecfieldx, elecfieldy=P.elecfieldy, elecfieldz=P.elecfieldz;
 
    return Lambda*((-2*Q_00_2*(Lq_tilde + L3*Q_01_2) - L3*Q_00_0*(2*Q_00_1 - 2*Q_01_0 + Q_11_1) + 2*Lq_tilde*(Q_02_0 + Q_11_2 - Q_12_1) - 2*(-deltaepslon*elecfieldx*elecfieldy + sigma*QN01 + bb*QN01*(QN00 + QN11) + bb*QN02*QN12 + 2*cc*QN01*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) + QN00*QN11 + (QN11*QN11) + (QN12*QN12)) - L1*(Q_01_00 + Q_01_11 + Q_01_22)) + L3*(-2*Q_02_0*Q_02_1 + 2*Q_01_0*(Q_01_1 + Q_02_2) - Q_11_0*(Q_00_1 + 2*Q_11_1) + 2*Q_01_2*(Q_02_0 - Q_11_2 + Q_12_1) + 2*(-(Q_12_0*Q_12_1) + Q_01_1*(Q_11_1 + Q_12_2) + QN00*Q_01_00 + 2*QN01*Q_01_01 + 2*QN02*Q_01_02 + QN11*Q_01_11 + 2*QN12*Q_01_12 - (QN00 + QN11)*Q_01_22)) + (L2 + Ls)*(Q_00_01 + Q_01_00 + Q_01_11 + Q_02_12 + Q_11_01 + Q_12_02))/2.);
}

MARLICS_HD inline double ldg_functional_derivative_02(const BulkParams& P, const double QN[5], const double dQ[15], const double ddQ[30], const double v[3])
{
  const double sigma=P.sigma, a=P.a, bb=P.bb, cc=P.cc, L1=P.L1, L2=P.L2, L3=P.L3, Lq=P.Lq, Ls=P.Ls, Lq_tilde=P.Lq_tilde, Lambda=P.Lambda, Lambda_s=P.Lambda_s, S_eq=P.S_eq, deltaepslon=P.deltaepslon, elecfieldx=P.elecfieldx, elecfieldy=P.elecfieldy, elecfieldz=P.elecfieldz;
 
    return Lambda*((4*Lq_tilde*Q_00_1 + 2*Lq_tilde*(-Q_01_0 + Q_11_1 + Q_12_2) - 2*(-deltaepslon*elecfieldx*elecfieldz + sigma*QN02 - bb*QN02*QN11 + bb*QN01*QN12 + 2*cc*QN02*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) + QN00*QN11 + (QN11*QN11) + (QN12*QN12)) - L1*(Q_02_00 + Q_02_11 + Q_02_22)) + L3*(-(Q_00_2*(2*Q_02_2 + Q_11_0)) - Q_00_0*(2*Q_00_2 - 2*Q_02_0 + Q_11_2) + 2*(Q_01_1*Q_02_0 + Q_02_0*Q_02_2 - (Q_02_2 + Q_11_0)*Q_11_2 + Q_02_2*Q_12_1) + 2*(Q_01_0*(-Q_01_2 + Q_02_1) - Q_12_0*Q_12_2 + Q_02_1*(Q_11_1 + Q_12_2) + QN00*Q_02_00 + 2*QN01*Q_02_01 + 2*QN02*Q_02_02 + QN11*Q_02_11 + 2*QN12*Q_02_12 - (QN00 + QN11)*Q_02_22)) + (L2 + Ls)*(Q_01_12 + Q_02_00 + Q_02_22 - Q_11_02 + Q_12_01))/2.);
}

MARLICS_HD inline double ldg_functional_derivative_11(const BulkParams& P, const double QN[5], const double dQ[15], const double ddQ[30], const double v[3])
{
  const double sigma=P.sigma, a=P.a, bb=P.bb, cc=P.cc, L1=P.L1, L2=P.L2, L3=P.L3, Lq=P.Lq, Ls=P.Ls, Lq_tilde=P.Lq_tilde, Lambda=P.Lambda, Lambda_s=P.Lambda_s, S_eq=P.S_eq, deltaepslon=P.deltaepslon, elecfieldx=P.elecfieldx, elecfieldy=P.elecfieldy, elecfieldz=P.elecfieldz;
 
    return Lambda*((-deltaepslon*((elecfieldx*elecfieldx) - 2*(elecfieldy*elecfieldy) + (elecfieldz*elecfieldz)) + bb*(2*(QN00*QN00) - (QN01*QN01) + 2*(QN02*QN02) + 2*QN00*QN11 - (QN11*QN11) - (QN12*QN12)) - 3*QN11*(sigma + 2*cc*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) + QN00*QN11 + (QN11*QN11) + (QN12*QN12))) + 6*Lq_tilde*(-Q_01_2 + Q_12_0) - (L2 + Ls)*Q_00_00 + 3*L1*(Q_11_00 + Q_11_11 + Q_11_22) + L3*((Q_00_0*Q_00_0) - 2*(Q_00_1*Q_00_1) + (Q_00_2*Q_00_2) + (Q_01_0*Q_01_0) - 2*(Q_01_1*Q_01_1) + (Q_01_2*Q_01_2) + (Q_02_0*Q_02_0) - 2*(Q_02_1*Q_02_1) + (Q_02_2*Q_02_2) + (4*Q_00_0 + 3*(Q_01_1 + Q_02_2))*Q_11_0 + (Q_11_0*Q_11_0) - 2*Q_00_1*Q_11_1 + 3*Q_01_0*Q_11_1 + (Q_11_1*Q_11_1) - 2*Q_00_2*Q_11_2 + 3*Q_02_0*Q_11_2 - 2*(Q_11_2*Q_11_2) + (Q_12_0*Q_12_0) + 3*Q_11_2*Q_12_1 - 2*(Q_12_1*Q_12_1) + 3*Q_11_1*Q_12_2 + (Q_12_2*Q_12_2) + 3*QN00*Q_11_00 + 6*QN01*Q_11_01 + 6*QN02*Q_11_02 + 3*QN11*Q_11_11 + 6*QN12*Q_11_12 - 3*(QN00 + QN11)*Q_11_22) + (L2 + Ls)*(Q_00_22 + Q_01_01 - 2*Q_02_02 + 2*Q_11_11 + Q_11_22 + Q_12_12))/3.);
}

MARLICS_HD inline double ldg_functional_derivative_12(const BulkParams& P, const double QN[5], const double dQ[15], const double ddQ[30], const double v[3])
{
  const double sigma=P.sigma, a=P.a, bb=P.bb, cc=P.cc, L1=P.L1, L2=P.L2, L3=P.L3, Lq=P.Lq, Ls=P.Ls, Lq_tilde=P.Lq_tilde, Lambda=P.Lambda, Lambda_s=P.Lambda_s, S_eq=P.S_eq, deltaepslon=P.deltaepslon, elecfieldx=P.elecfieldx, elecfieldy=P.elecfieldy, elecfieldz=P.elecfieldz;
 
    return Lambda*((2*Lq_tilde*(Q_01_1 - Q_02_2 - 2*Q_11_0) - 2*Q_00_0*(Lq_tilde - L3*Q_12_0) - (L2 + Ls)*Q_00_12 + (L2 + Ls)*(Q_01_02 + Q_02_01 + Q_12_11 + Q_12_22) - 2*(-deltaepslon*elecfieldy*elecfieldz + bb*QN01*QN02 + (sigma - bb*QN00 + 2*cc*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) + QN00*QN11 + (QN11*QN11)))*QN12 + 2*cc*(QN12*QN12*QN12) - L1*(Q_12_00 + Q_12_11 + Q_12_22)) + L3*(-(Q_00_1*(2*Q_00_2 + Q_11_2)) + 2*Q_02_0*Q_12_2 - 2*(Q_11_2 - Q_12_1)*(Q_11_1 + Q_12_2) - Q_00_2*(Q_11_1 + 2*Q_12_2) + 2*(Q_01_1*(-Q_01_2 + Q_12_0) + Q_02_2*(-Q_02_1 + Q_12_0) + Q_01_0*Q_12_1 + QN00*Q_12_00 + 2*QN01*Q_12_01 + 2*QN02*Q_12_02 + QN11*Q_12_11 + 2*QN12*Q_12_12 - (QN00 + QN11)*Q_12_22)))/2.);
}

MARLICS_HD inline double ldg_energy_calculation(const BulkParams& P, const double QN[5], const double dQ[15], const double v[3])
{
  const double sigma=P.sigma, a=P.a, bb=P.bb, cc=P.cc, L1=P.L1, L2=P.L2, L3=P.L3, Lq=P.Lq, Ls=P.Ls, Lq_tilde=P.Lq_tilde, Lambda=P.Lambda, Lambda_s=P.Lambda_s, S_eq=P.S_eq, deltaepslon=P.deltaepslon, elecfieldx=P.elecfieldx, elecfieldy=P.elecfieldy, elecfieldz=P.elecfieldz;
 
    return ((2*deltaepslon*((elecfieldx*elecfieldx)*QN00 + 2*elecfieldx*(elecfieldy*QN01 + elecfieldz*QN02) + (elecfieldy*elecfieldy)*QN11 - (elecfieldz*elecfieldz)*(QN00 + QN11) + 2*elecfieldy*elecfieldz*QN12) + 2*sigma*((QN00*QN00) + (QN01*QN01) + (QN02*QN02) + QN00*QN11 + (QN11*QN11) + (QN12*QN12)) + 2*cc*pow((QN00*QN00) + (QN01*QN01) + (QN02*QN02) + QN00*QN11 + (QN11*QN11) + (QN12*QN12),2) - 2*bb*((QN00*QN00)*QN11 - (QN01*QN01)*QN11 + (QN02*QN02)*QN11 - 2*QN01*QN02*QN12 + QN00*(-(QN01*QN01) + (QN11*QN11) + (QN12*QN12))) + Ls*((Q_00_2*Q_00_2) + (Q_01_0*Q_01_0) + (Q_02_0*Q_02_0) + 2*Q_01_2*Q_02_1 + pow(Q_00_0 - Q_02_2,2) - 2*Q_02_2*Q_11_0 + Q_01_1*(Q_01_1 + 2*Q_11_0) + 2*Q_00_2*(Q_02_0 + Q_11_2) + 2*(Q_01_2 + Q_02_1)*Q_12_0 + pow(Q_11_2 + Q_12_1,2) + 2*Q_00_1*(Q_01_0 - Q_12_2) + pow(Q_11_1 - Q_12_2,2)) + 2*L1*((Q_00_0*Q_00_0) + (Q_00_1*Q_00_1) + (Q_00_2*Q_00_2) + (Q_01_0*Q_01_0) + (Q_01_1*Q_01_1) + (Q_01_2*Q_01_2) + (Q_02_0*Q_02_0) + (Q_02_1*Q_02_1) + (Q_02_2*Q_02_2) + Q_00_0*Q_11_0 + (Q_11_0*Q_11_0) + Q_00_1*Q_11_1 + (Q_11_1*Q_11_1) + Q_00_2*Q_11_2 + (Q_11_2*Q_11_2) + (Q_12_0*Q_12_0) + (Q_12_1*Q_12_1) + (Q_12_2*Q_12_2)) + 2*L3*(QN00*(Q_00_0*Q_00_0) + QN11*(Q_00_1*Q_00_1) - QN00*(Q_00_2*Q_00_2) - QN11*(Q_00_2*Q_00_2) + QN00*(Q_01_0*Q_01_0) + 2*QN01*Q_01_0*Q_01_1 + QN11*(Q_01_1*Q_01_1) + 2*QN02*Q_01_0*Q_01_2 + 2*QN12*Q_01_1*Q_01_2 - QN00*(Q_01_2*Q_01_2) - QN11*(Q_01_2*Q_01_2) + QN00*(Q_02_0*Q_02_0) + 2*QN01*Q_02_0*Q_02_1 + QN11*(Q_02_1*Q_02_1) + 2*QN02*Q_02_0*Q_02_2 + 2*QN12*Q_02_1*Q_02_2 - QN00*(Q_02_2*Q_02_2) - QN11*(Q_02_2*Q_02_2) + QN02*Q_00_2*Q_11_0 + QN00*(Q_11_0*Q_11_0) + QN12*Q_00_2*Q_11_1 + 2*QN01*Q_11_0*Q_11_1 + QN11*(Q_11_1*Q_11_1) - QN00*Q_00_2*Q_11_2 - QN11*Q_00_2*Q_11_2 + 2*QN02*Q_11_0*Q_11_2 + 2*QN12*Q_11_1*Q_11_2 - QN00*(Q_11_2*Q_11_2) - QN11*(Q_11_2*Q_11_2) + Q_00_0*(QN00*Q_11_0 + QN01*(2*Q_00_1 + Q_11_1) + QN02*(2*Q_00_2 + Q_11_2)) + Q_00_1*(QN01*Q_11_0 + QN11*Q_11_1 + QN12*(2*Q_00_2 + Q_11_2)) + QN00*(Q_12_0*Q_12_0) + 2*QN01*Q_12_0*Q_12_1 + QN11*(Q_12_1*Q_12_1) + 2*(QN02*Q_12_0 + QN12*Q_12_1)*Q_12_2 - (QN00 + QN11)*(Q_12_2*Q_12_2)) + 2*Lq_tilde*((-QN00 + QN11)*Q_01_2 + (2*QN00 + QN11)*Q_02_1 + QN12*(Q_00_0 - Q_01_1 + Q_02_2 + 2*Q_11_0) - (QN00 + 2*QN11)*Q_12_0 + QN01*(Q_00_2 - Q_02_0 - Q_11_2 + Q_12_1) - QN02*(2*Q_00_1 - Q_01_0 + Q_11_1 + Q_12_2)) + L2*(pow(Q_00_0 + Q_01_1 + Q_02_2,2) + pow(Q_00_2 - Q_02_0 + Q_11_2 - Q_12_1,2) + pow(Q_01_0 + Q_11_1 + Q_12_2,2)))/2.);
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
