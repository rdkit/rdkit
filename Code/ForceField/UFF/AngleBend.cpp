// $Id$
//
//  Copyright (C) 2004-2013 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "AngleBend.h"
#include "BondStretch.h"
#include "Params.h"
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>

namespace ForceFields {
namespace UFF {

namespace Utils {
double calcAngleForceConstant(double theta0, double bondOrder12,
                              double bondOrder23, const AtomicParams *at1Params,
                              const AtomicParams *at2Params,
                              const AtomicParams *at3Params) {
  double cosTheta0 = cos(theta0);
  double r12 = calcBondRestLength(bondOrder12, at1Params, at2Params);
  double r23 = calcBondRestLength(bondOrder23, at2Params, at3Params);
  double r13 = sqrt(r12 * r12 + r23 * r23 - 2. * r12 * r23 * cosTheta0);
  double beta = 2. * Params::G / (r12 * r23);

  double preFactor = beta * at1Params->Z1 * at3Params->Z1 / int_pow<5>(r13);
  double rTerm = r12 * r23;
  double innerBit =
      3. * rTerm * (1. - cosTheta0 * cosTheta0) - r13 * r13 * cosTheta0;
  double res = preFactor * rTerm * innerBit;
  return res;
}

void calcAngleBendGrad(RDGeom::Point3D *r, double *dist, double **g,
                       double &dE_dTheta, double &cosTheta, double &sinTheta) {
  // -------
  // dTheta/dx is trickier:
  double dCos_dS[6] = {1.0 / dist[0] * (r[1].x - cosTheta * r[0].x),
                       1.0 / dist[0] * (r[1].y - cosTheta * r[0].y),
                       1.0 / dist[0] * (r[1].z - cosTheta * r[0].z),
                       1.0 / dist[1] * (r[0].x - cosTheta * r[1].x),
                       1.0 / dist[1] * (r[0].y - cosTheta * r[1].y),
                       1.0 / dist[1] * (r[0].z - cosTheta * r[1].z)};

  g[0][0] += dE_dTheta * dCos_dS[0] / (-sinTheta);
  g[0][1] += dE_dTheta * dCos_dS[1] / (-sinTheta);
  g[0][2] += dE_dTheta * dCos_dS[2] / (-sinTheta);

  g[1][0] += dE_dTheta * (-dCos_dS[0] - dCos_dS[3]) / (-sinTheta);
  g[1][1] += dE_dTheta * (-dCos_dS[1] - dCos_dS[4]) / (-sinTheta);
  g[1][2] += dE_dTheta * (-dCos_dS[2] - dCos_dS[5]) / (-sinTheta);

  g[2][0] += dE_dTheta * dCos_dS[3] / (-sinTheta);
  g[2][1] += dE_dTheta * dCos_dS[4] / (-sinTheta);
  g[2][2] += dE_dTheta * dCos_dS[5] / (-sinTheta);
}
}  // end of namespace Utils

AngleBendContrib::AngleBendContrib(ForceField *owner, unsigned int idx1,
                                   unsigned int idx2, unsigned int idx3,
                                   double bondOrder12, double bondOrder23,
                                   const AtomicParams *at1Params,
                                   const AtomicParams *at2Params,
                                   const AtomicParams *at3Params,
                                   unsigned int order) {
  PRECONDITION(owner, "bad owner");
  PRECONDITION(at1Params, "bad params pointer");
  PRECONDITION(at2Params, "bad params pointer");
  PRECONDITION(at3Params, "bad params pointer");
  PRECONDITION((idx1 != idx2 && idx2 != idx3 && idx1 != idx3),
               "degenerate points");
  URANGE_CHECK(idx1, owner->positions().size());
  URANGE_CHECK(idx2, owner->positions().size());
  URANGE_CHECK(idx3, owner->positions().size());
  // the following is a hack to get decent geometries
  // with 3- and 4-membered rings incorporating sp2 atoms
  double theta0 = at2Params->theta0;
  if (order >= 30) {
    switch (order) {
      case 30:
        theta0 = 150.0 / 180.0 * M_PI;
        break;
      case 35:
        theta0 = 60.0 / 180.0 * M_PI;
        break;
      case 40:
        theta0 = 135.0 / 180.0 * M_PI;
        break;
      case 45:
        theta0 = 90.0 / 180.0 * M_PI;
        break;
    }
    order = 0;
  }
  // end of the hack
  dp_forceField = owner;
  d_at1Idx = idx1;
  d_at2Idx = idx2;
  d_at3Idx = idx3;
  d_order = order;
  d_forceConstant = Utils::calcAngleForceConstant(
      theta0, bondOrder12, bondOrder23, at1Params, at2Params, at3Params);
  if (order == 0) {
    double sinTheta0 = sin(theta0);
    double cosTheta0 = cos(theta0);
    d_C2 = 1. / (4. * std::max(sinTheta0 * sinTheta0, 1e-8));
    d_C1 = -4. * d_C2 * cosTheta0;
    d_C0 = d_C2 * (2. * cosTheta0 * cosTheta0 + 1.);
  }
}

double AngleBendContrib::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  double dist1 = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
  double dist2 = dp_forceField->distance(d_at2Idx, d_at3Idx, pos);

  RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
                     pos[3 * d_at1Idx + 2]);
  RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
                     pos[3 * d_at2Idx + 2]);
  RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
                     pos[3 * d_at3Idx + 2]);
  RDGeom::Point3D p12 = p1 - p2;
  RDGeom::Point3D p32 = p3 - p2;
  double cosTheta = p12.dotProduct(p32) / (dist1 * dist2);
  clipToOne(cosTheta);
  // we need sin^2(theta) to get cos(2*theta), so compute that:
  double sinThetaSq = 1. - cosTheta * cosTheta;

  double angleTerm = getEnergyTerm(cosTheta, sinThetaSq);
  double res = d_forceConstant * angleTerm;

  return res;
}

void AngleBendContrib::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");

  double dist[2] = {dp_forceField->distance(d_at1Idx, d_at2Idx, pos),
                    dp_forceField->distance(d_at2Idx, d_at3Idx, pos)};

  RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
                     pos[3 * d_at1Idx + 2]);
  RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
                     pos[3 * d_at2Idx + 2]);
  RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
                     pos[3 * d_at3Idx + 2]);
  double *g[3] = {&(grad[3 * d_at1Idx]), &(grad[3 * d_at2Idx]),
                  &(grad[3 * d_at3Idx])};
  RDGeom::Point3D r[2] = {(p1 - p2) / dist[0], (p3 - p2) / dist[1]};
  double cosTheta = r[0].dotProduct(r[1]);
  clipToOne(cosTheta);
  double sinThetaSq = 1.0 - cosTheta * cosTheta;
  double sinTheta = std::max(sqrt(sinThetaSq), 1.0e-8);

  // std::cerr << "GRAD: " << cosTheta << " (" << acos(cosTheta)<< "), ";
  // std::cerr << sinTheta << " (" << asin(sinTheta)<< ")" << std::endl;

  // use the chain rule:
  // dE/dx = dE/dTheta * dTheta/dx

  // dE/dTheta is independent of cartesians:
  double dE_dTheta = getThetaDeriv(cosTheta, sinTheta);

  Utils::calcAngleBendGrad(r, dist, g, dE_dTheta, cosTheta, sinTheta);
}

double AngleBendContrib::getEnergyTerm(double cosTheta,
                                       double sinThetaSq) const {
  PRECONDITION(d_order == 0 || d_order == 1 || d_order == 2 || d_order == 3 ||
                   d_order == 4,
               "bad order");
  // cos(2x) = cos^2(x) - sin^2(x);
  double cos2Theta = cosTheta * cosTheta - sinThetaSq;

  double res = 0.0;
  if (d_order == 0) {
    res = d_C0 + d_C1 * cosTheta + d_C2 * cos2Theta;
  } else {
    switch (d_order) {
      case 1:
        res = -cosTheta;
        break;
      case 2:
        res = cos2Theta;
        break;
      case 3:
        // cos(3x) = cos^3(x) - 3*cos(x)*sin^2(x)
        res = cosTheta * (cosTheta * cosTheta - 3. * sinThetaSq);
        break;
      case 4:
        // cos(4x) = cos^4(x) - 6*cos^2(x)*sin^2(x)+sin^4(x)
        res = int_pow<4>(cosTheta) - 6. * cosTheta * cosTheta * sinThetaSq +
              sinThetaSq * sinThetaSq;
        break;
    }
    res = 1. - res;
    res /= (double)(d_order * d_order);
  }
  return res;
}

double AngleBendContrib::getThetaDeriv(double cosTheta, double sinTheta) const {
  PRECONDITION(d_order == 0 || d_order == 1 || d_order == 2 || d_order == 3 ||
                   d_order == 4,
               "bad order");

  double dE_dTheta = 0.0;
  double sin2Theta = 2. * sinTheta * cosTheta;

  if (d_order == 0) {
    dE_dTheta =
        -1. * d_forceConstant * (d_C1 * sinTheta + 2. * d_C2 * sin2Theta);
  } else {
    // E = k/n^2 [1-cos(n theta)]
    // dE = - k/n^2 * d cos(n theta)

    // these all use:
    // d cos(ax) = -a sin(ax)

    switch (d_order) {
      case 1:
        dE_dTheta = -sinTheta;
        break;
      case 2:
        // sin(2*x) = 2*cos(x)*sin(x)
        dE_dTheta = sin2Theta;
        break;
      case 3:
        // sin(3*x) = 3*sin(x) - 4*sin^3(x)
        dE_dTheta = sinTheta * (3. - 4. * sinTheta * sinTheta);
        break;
      case 4:
        // sin(4*x) = cos(x)*(4*sin(x) - 8*sin^3(x))
        dE_dTheta = cosTheta * sinTheta * (4. - 8. * sinTheta * sinTheta);
        break;
    }
    dE_dTheta *= d_forceConstant / (double)(d_order);
  }
  return dE_dTheta;
}
}  // namespace UFF
}  // namespace ForceFields
