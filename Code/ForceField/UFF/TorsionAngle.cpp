// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "TorsionAngle.h"
#include "Params.h"
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>

namespace ForceFields {
namespace UFF {
namespace Utils {
double calculateCosTorsion(const RDGeom::Point3D &p1, const RDGeom::Point3D &p2,
                           const RDGeom::Point3D &p3,
                           const RDGeom::Point3D &p4) {
  RDGeom::Point3D r1 = p1 - p2, r2 = p3 - p2, r3 = p2 - p3, r4 = p4 - p3;
  RDGeom::Point3D t1 = r1.crossProduct(r2);
  RDGeom::Point3D t2 = r3.crossProduct(r4);
  double d1 = t1.length(), d2 = t2.length();
  double cosPhi = t1.dotProduct(t2) / (d1 * d2);
  clipToOne(cosPhi);
  return cosPhi;
}

// used locally
bool isInGroup6(int num) {
  return (num == 8 || num == 16 || num == 34 || num == 52 || num == 84);
}

// used locally, implement equation 17 of the UFF paper.
double equation17(double bondOrder23, const AtomicParams *at2Params,
                  const AtomicParams *at3Params) {
  return 5. * sqrt(at2Params->U1 * at3Params->U1) *
         (1. + 4.18 * log(bondOrder23));
}

void calcTorsionGrad(RDGeom::Point3D *r, RDGeom::Point3D *t, double *d,
                     double **g, double &sinTerm, double &cosPhi) {
  // -------
  // dTheta/dx is trickier:
  double dCos_dT[6] = {1.0 / d[0] * (t[1].x - cosPhi * t[0].x),
                       1.0 / d[0] * (t[1].y - cosPhi * t[0].y),
                       1.0 / d[0] * (t[1].z - cosPhi * t[0].z),
                       1.0 / d[1] * (t[0].x - cosPhi * t[1].x),
                       1.0 / d[1] * (t[0].y - cosPhi * t[1].y),
                       1.0 / d[1] * (t[0].z - cosPhi * t[1].z)};

  g[0][0] += sinTerm * (dCos_dT[2] * r[1].y - dCos_dT[1] * r[1].z);
  g[0][1] += sinTerm * (dCos_dT[0] * r[1].z - dCos_dT[2] * r[1].x);
  g[0][2] += sinTerm * (dCos_dT[1] * r[1].x - dCos_dT[0] * r[1].y);

  g[1][0] += sinTerm *
             (dCos_dT[1] * (r[1].z - r[0].z) + dCos_dT[2] * (r[0].y - r[1].y) +
              dCos_dT[4] * (-r[3].z) + dCos_dT[5] * (r[3].y));
  g[1][1] += sinTerm *
             (dCos_dT[0] * (r[0].z - r[1].z) + dCos_dT[2] * (r[1].x - r[0].x) +
              dCos_dT[3] * (r[3].z) + dCos_dT[5] * (-r[3].x));
  g[1][2] += sinTerm *
             (dCos_dT[0] * (r[1].y - r[0].y) + dCos_dT[1] * (r[0].x - r[1].x) +
              dCos_dT[3] * (-r[3].y) + dCos_dT[4] * (r[3].x));

  g[2][0] += sinTerm *
             (dCos_dT[1] * (r[0].z) + dCos_dT[2] * (-r[0].y) +
              dCos_dT[4] * (r[3].z - r[2].z) + dCos_dT[5] * (r[2].y - r[3].y));
  g[2][1] += sinTerm *
             (dCos_dT[0] * (-r[0].z) + dCos_dT[2] * (r[0].x) +
              dCos_dT[3] * (r[2].z - r[3].z) + dCos_dT[5] * (r[3].x - r[2].x));
  g[2][2] += sinTerm *
             (dCos_dT[0] * (r[0].y) + dCos_dT[1] * (-r[0].x) +
              dCos_dT[3] * (r[3].y - r[2].y) + dCos_dT[4] * (r[2].x - r[3].x));

  g[3][0] += sinTerm * (dCos_dT[4] * r[2].z - dCos_dT[5] * r[2].y);
  g[3][1] += sinTerm * (dCos_dT[5] * r[2].x - dCos_dT[3] * r[2].z);
  g[3][2] += sinTerm * (dCos_dT[3] * r[2].y - dCos_dT[4] * r[2].x);
}
}  // namespace Utils

TorsionAngleContrib::TorsionAngleContrib(
    ForceField *owner, unsigned int idx1, unsigned int idx2, unsigned int idx3,
    unsigned int idx4, double bondOrder23, int atNum2, int atNum3,
    RDKit::Atom::HybridizationType hyb2, RDKit::Atom::HybridizationType hyb3,
    const AtomicParams *at2Params, const AtomicParams *at3Params,
    bool endAtomIsSP2) {
  PRECONDITION(owner, "bad owner");
  PRECONDITION(at2Params, "bad params pointer");
  PRECONDITION(at3Params, "bad params pointer");
  PRECONDITION((idx1 != idx2 && idx1 != idx3 && idx1 != idx4 && idx2 != idx3 &&
                idx2 != idx4 && idx3 != idx4),
               "degenerate points");
  URANGE_CHECK(idx1, owner->positions().size());
  URANGE_CHECK(idx2, owner->positions().size());
  URANGE_CHECK(idx3, owner->positions().size());
  URANGE_CHECK(idx4, owner->positions().size());

  dp_forceField = owner;
  d_at1Idx = idx1;
  d_at2Idx = idx2;
  d_at3Idx = idx3;
  d_at4Idx = idx4;

  calcTorsionParams(bondOrder23, atNum2, atNum3, hyb2, hyb3, at2Params,
                    at3Params, endAtomIsSP2);
}

void TorsionAngleContrib::calcTorsionParams(double bondOrder23, int atNum2,
                                            int atNum3,
                                            RDKit::Atom::HybridizationType hyb2,
                                            RDKit::Atom::HybridizationType hyb3,
                                            const AtomicParams *at2Params,
                                            const AtomicParams *at3Params,
                                            bool endAtomIsSP2) {
  PRECONDITION((hyb2 == RDKit::Atom::SP2 || hyb2 == RDKit::Atom::SP3) &&
                   (hyb3 == RDKit::Atom::SP2 || hyb3 == RDKit::Atom::SP3),
               "bad hybridizations");

  if (hyb2 == RDKit::Atom::SP3 && hyb3 == RDKit::Atom::SP3) {
    // general case:
    d_forceConstant = sqrt(at2Params->V1 * at3Params->V1);
    d_order = 3;
    d_cosTerm = -1;  // phi0=60

    // special case for single bonds between group 6 elements:
    if (bondOrder23 == 1.0 && Utils::isInGroup6(atNum2) &&
        Utils::isInGroup6(atNum3)) {
      double V2 = 6.8, V3 = 6.8;
      if (atNum2 == 8) {
        V2 = 2.0;
      }
      if (atNum3 == 8) {
        V3 = 2.0;
      }
      d_forceConstant = sqrt(V2 * V3);
      d_order = 2;
      d_cosTerm = -1;  // phi0=90
    }
  } else if (hyb2 == RDKit::Atom::SP2 && hyb3 == RDKit::Atom::SP2) {
    d_forceConstant = Utils::equation17(bondOrder23, at2Params, at3Params);
    d_order = 2;
    // FIX: is this angle term right?
    d_cosTerm = 1.0;  // phi0= 180
  } else {
    // SP2 - SP3,  this is, by default, independent of atom type in UFF:
    d_forceConstant = 1.0;
    d_order = 6;
    d_cosTerm = 1.0;  // phi0 = 0
    if (bondOrder23 == 1.0) {
      // special case between group 6 sp3 and non-group 6 sp2:
      if ((hyb2 == RDKit::Atom::SP3 && Utils::isInGroup6(atNum2) &&
           !Utils::isInGroup6(atNum3)) ||
          (hyb3 == RDKit::Atom::SP3 && Utils::isInGroup6(atNum3) &&
           !Utils::isInGroup6(atNum2))) {
        d_forceConstant = Utils::equation17(bondOrder23, at2Params, at3Params);
        d_order = 2;
        d_cosTerm = -1;  // phi0 = 90;
      }

      // special case for sp3 - sp2 - sp2
      // (i.e. the sp2 has another sp2 neighbor, like propene)
      else if (endAtomIsSP2) {
        d_forceConstant = 2.0;
        d_order = 3;
        d_cosTerm = -1;  // phi0 = 180;
      }
    }
  }
}
double TorsionAngleContrib::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(d_order == 2 || d_order == 3 || d_order == 6, "bad order");

  RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
                     pos[3 * d_at1Idx + 2]);
  RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
                     pos[3 * d_at2Idx + 2]);
  RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
                     pos[3 * d_at3Idx + 2]);
  RDGeom::Point3D p4(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
                     pos[3 * d_at4Idx + 2]);

  double cosPhi = Utils::calculateCosTorsion(p1, p2, p3, p4);
  double sinPhiSq = 1 - cosPhi * cosPhi;

  // E(phi) = V/2 * (1 - cos(n*phi_0)*cos(n*phi))
  double cosNPhi = 0.0;
  switch (d_order) {
    case 2:
      // cos(2x) = 1 - 2sin^2(x)
      cosNPhi = 1 - 2 * sinPhiSq;
      break;
    case 3:
      // cos(3x) = cos^3(x) - 3*cos(x)*sin^2(x) = 4cos^3(x) -3cos(x)
      cosNPhi = cosPhi * (cosPhi * cosPhi - 3. * sinPhiSq);
      break;
    case 6:
      // cos(6x) = 1 - 32*sin^6(x) + 48*sin^4(x) - 18*sin^2(x)
      cosNPhi =
          1 + sinPhiSq * (-32. * sinPhiSq * sinPhiSq + 48. * sinPhiSq - 18.);
      break;
  }
  double res = d_forceConstant / 2.0 * (1. - d_cosTerm * cosNPhi);
  // std::cout << " torsion(" << d_at1Idx << "," << d_at2Idx << "," << d_at3Idx
  // << "," << d_at4Idx << "): " << cosPhi << "(" << acos(cosPhi) << ")" << " ->
  // " << res << std::endl;
  // if(d_at2Idx==5&&d_at3Idx==6) std::cerr << " torsion(" << d_at1Idx << "," <<
  // d_at2Idx << "," << d_at3Idx << "," << d_at4Idx << "): " << cosPhi << "(" <<
  // acos(cosPhi) << ")" << " -> " << res << std::endl;
  return res;
}

void TorsionAngleContrib::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");

  double *g[4] = {&(grad[3 * d_at1Idx]), &(grad[3 * d_at2Idx]),
                  &(grad[3 * d_at3Idx]), &(grad[3 * d_at4Idx])};

  RDGeom::Point3D r[4];
  RDGeom::Point3D t[2];
  double d[2];
  double cosPhi;
  RDKit::ForceFieldsHelper::computeDihedral(
      pos, d_at1Idx, d_at2Idx, d_at3Idx, d_at4Idx, nullptr, &cosPhi, r, t, d);
  double sinPhiSq = 1.0 - cosPhi * cosPhi;
  double sinPhi = ((sinPhiSq > 0.0) ? sqrt(sinPhiSq) : 0.0);

  // dE/dPhi is independent of cartesians:
  double dE_dPhi = getThetaDeriv(cosPhi, sinPhi);
#if 0
      if(dE_dPhi!=dE_dPhi){
        std::cout << "\tNaN in Torsion("<<d_at1Idx<<","<<d_at2Idx<<","<<d_at3Idx<<","<<d_at4Idx<<")"<< std::endl;
        std::cout << "sin: " << sinPhi << std::endl;
        std::cout << "cos: " << cosPhi << std::endl;
      }

#endif

  double sinTerm =
      dE_dPhi * (isDoubleZero(sinPhi) ? (1.0 / cosPhi) : (1.0 / sinPhi));

  Utils::calcTorsionGrad(r, t, d, g, sinTerm, cosPhi);
}

double TorsionAngleContrib::getThetaDeriv(double cosTheta,
                                          double sinTheta) const {
  PRECONDITION(d_order == 2 || d_order == 3 || d_order == 6, "bad order");
  double sinThetaSq = sinTheta * sinTheta;
  // cos(6x) = 1 - 32*sin^6(x) + 48*sin^4(x) - 18*sin^2(x)

  double res = 0.0;
  switch (d_order) {
    case 2:
      res = 2 * sinTheta * cosTheta;
      break;
    case 3:
      // sin(3*x) = 3*sin(x) - 4*sin^3(x)
      res = sinTheta * (3 - 4 * sinThetaSq);
      break;
    case 6:
      // sin(6x) = cos(x) * [ 32*sin^5(x) - 32*sin^3(x) + 6*sin(x) ]
      res = cosTheta * sinTheta * (32 * sinThetaSq * (sinThetaSq - 1) + 6);
      break;
  }
  res *= d_forceConstant / 2.0 * d_cosTerm * -1 * d_order;

  return res;
}
}  // namespace UFF
}  // namespace ForceFields
