//
//  Copyright (c) 2014-2024, Novartis Institutes for BioMedical Research and
//  other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MIFDescriptors.h"
#include <Geometry/point.h>
#include <Geometry/UniformRealValueGrid3D.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <ForceField/MMFF/Nonbonded.h>
#include <ForceField/UFF/Nonbonded.h>
#include <ForceField/UFF/Params.h>
#include <vector>
#include <iostream>
#include <fstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

namespace {
constexpr double M_55D = 0.95993108859688126730;
// angle C-O-H , C-N-H
constexpr double M_70_5D = 1.23045712265600235173;
constexpr double M_110D = 1.91986217719376253461;

constexpr double CUTOFF = 0.001;
constexpr double MIN_CUTOFF_VAL = CUTOFF * CUTOFF;
}  // namespace

namespace RDMIF {

std::unique_ptr<RDGeom::UniformRealValueGrid3D> constructGrid(
    const RDKit::ROMol &mol, int confId, double margin, double spacing) {
  PRECONDITION(mol.getNumConformers(), "No conformers available for molecule.");

  const std::vector<RDGeom::Point3D> &ptVect =
      mol.getConformer(confId).getPositions();
  RDGeom::Point3D originPt(0.0, 0.0, 0.0);
  RDGeom::Point3D marginPt(margin, margin, margin);
  const RDGeom::Point3D &firstPoint =
      (!ptVect.empty() ? ptVect.front() : originPt);

  auto minPt = firstPoint;
  auto maxPt = minPt;
  for (const auto &pt : ptVect) {
    for (auto i = 0u; i < pt.dimension(); ++i) {
      minPt[i] = std::min(minPt[i], pt[i]);
      maxPt[i] = std::max(maxPt[i], pt[i]);
    }
  }
  minPt -= marginPt;
  maxPt += marginPt;

  auto res = std::make_unique<RDGeom::UniformRealValueGrid3D>(
      maxPt.x - minPt.x, maxPt.y - minPt.y, maxPt.z - minPt.z, spacing, &minPt);

  return res;
}

VdWaals::VdWaals(const RDKit::ROMol &mol, int confId, double cutoff) {
  d_cutoff = std::max(cutoff * cutoff, MIN_CUTOFF_VAL);
  d_nAtoms = mol.getNumAtoms();
  d_pos.reserve(3 * d_nAtoms);
  d_R_star_ij.reserve(d_nAtoms);
  d_wellDepth.reserve(d_nAtoms);
  // this will throw a ConformerException if confId does not exist
  const RDKit::Conformer &conf = mol.getConformer(confId);
  d_mol.reset(new RDKit::ROMol(mol, false, conf.getId()));
}

void VdWaals::fillVectors() {
  const auto &conf = d_mol->getConformer();
  for (unsigned int i = 0; i < d_nAtoms; i++) {
    const RDGeom::Point3D &pt = conf.getAtomPos(i);
    d_pos.push_back(pt.x);
    d_pos.push_back(pt.y);
    d_pos.push_back(pt.z);
    fillVdwParamVectors(i);
  }
}

MMFFVdWaals::MMFFVdWaals(const RDKit::ROMol &mol, int confId,
                         unsigned int probeAtomType, bool scaling,
                         double cutoff)
    : VdWaals::VdWaals(mol, confId, cutoff), d_scaling(scaling) {
  d_props.reset(new RDKit::MMFF::MMFFMolProperties(*d_mol));
  if (!d_props->isValid()) {
    throw ValueErrorException(
        "No MMFF atom types available for at least one atom in molecule.");
  }
  d_mmffVdW = RDKit::MMFF::DefaultParameters::getMMFFVdW();
  d_probeParams = (*d_mmffVdW)(probeAtomType);
  fillVectors();
}

void MMFFVdWaals::fillVdwParamVectors(unsigned int atomIdx) {
  PRECONDITION(atomIdx < d_mol->getNumAtoms(), "atomIdx out of bounds");
  const auto iAtomType = d_props->getMMFFAtomType(atomIdx);
  auto params = (*d_mmffVdW)(iAtomType);
  auto rStarIJ = ForceFields::MMFF::Utils::calcUnscaledVdWMinimum(
      d_mmffVdW, params, d_probeParams);
  d_R_star_ij.push_back(rStarIJ);
  d_wellDepth.push_back(ForceFields::MMFF::Utils::calcUnscaledVdWWellDepth(
      rStarIJ, params, d_probeParams));
  // scaling for taking undirected H-Bonds into account
  if (d_scaling) {
    ForceFields::MMFF::Utils::scaleVdWParams(d_R_star_ij[atomIdx],
                                             d_wellDepth[atomIdx], d_mmffVdW,
                                             params, d_probeParams);
  }
}

UFFVdWaals::UFFVdWaals(const RDKit::ROMol &mol, int confId,
                       const std::string &probeAtomType, double cutoff)
    : VdWaals::VdWaals(mol, confId, cutoff) {
  d_uffParamColl = ForceFields::UFF::ParamCollection::getParams();
  d_probeParams = (*d_uffParamColl)(probeAtomType);
  const auto [params, haveParams] = RDKit::UFF::getAtomTypes(mol);
  if (!haveParams) {
    throw ValueErrorException(
        "No UFF atom types available for at least one atom in molecule.");
  }
  d_params = std::move(params);
  fillVectors();
}

void UFFVdWaals::fillVdwParamVectors(unsigned int atomIdx) {
  PRECONDITION(atomIdx < d_mol->getNumAtoms(), "atomIdx out of bounds");
  d_R_star_ij.push_back(d_probeParams->x1 * d_params[atomIdx]->x1);
  d_wellDepth.push_back(ForceFields::UFF::Utils::calcNonbondedDepth(
      d_probeParams, d_params[atomIdx]));
}

double VdWaals::operator()(double x, double y, double z, double thres) const {
  double res = 0.0;
  for (unsigned int i = 0, j = 0; i < d_nAtoms; ++i) {
    auto temp = x - d_pos[j++];
    auto dist2 = temp * temp;
    temp = y - d_pos[j++];
    dist2 += temp * temp;
    temp = z - d_pos[j++];
    dist2 += temp * temp;
    if (dist2 < thres) {
      dist2 = std::max(dist2, d_cutoff);
      res += calcEnergy(dist2, d_R_star_ij[i], d_wellDepth[i]);
    }
  }
  return res;
}

double UFFVdWaals::calcEnergy(double dist2, double x_ij,
                              double wellDepth) const {
  double r6 = x_ij / dist2;
  r6 *= r6 * r6;
  double r12 = r6 * r6;
  return wellDepth * (r12 - 2.0 * r6);
}

double MMFFVdWaals::calcEnergy(double dist2, double R_star_ij,
                               double wellDepth) const {
  return ForceFields::MMFF::Utils::calcVdWEnergy(sqrt(dist2), R_star_ij,
                                                 wellDepth);
}

namespace CoulombDetail {
constexpr double prefactor =
    1 / (4.0 * 3.141592 * 8.854188) * 1.602 * 1.602 * 6.02214129 * 10000;
}

Coulomb::Coulomb(const std::vector<double> &charges,
                 const std::vector<RDGeom::Point3D> &positions,
                 double probeCharge, bool absVal, double alpha, double cutoff)
    : d_nAtoms(charges.size()),
      d_absVal(absVal),
      d_cutoff(cutoff * cutoff),
      d_probe(probeCharge),
      d_alpha(alpha),
      d_charges(charges) {
  PRECONDITION(d_charges.size() == positions.size(),
               "Lengths of positions and charges vectors do not match.");
  d_pos.reserve(3 * d_nAtoms);
  for (const auto &position : positions) {
    d_pos.push_back(position.x);
    d_pos.push_back(position.y);
    d_pos.push_back(position.z);
  }
  if (fabs(d_alpha) < MIN_CUTOFF_VAL) {
    d_softcore = false;
    if (d_cutoff < MIN_CUTOFF_VAL) {
      d_cutoff = CUTOFF;
    }
  } else {
    d_softcore = true;
  }
}

Coulomb::Coulomb(const RDKit::ROMol &mol, int confId, double probeCharge,
                 bool absVal, const std::string &prop, double alpha,
                 double cutoff)
    : d_nAtoms(mol.getNumAtoms()),
      d_absVal(absVal),
      d_cutoff(cutoff * cutoff),
      d_probe(probeCharge),
      d_alpha(alpha) {
  d_charges.reserve(d_nAtoms);
  d_pos.reserve(3 * d_nAtoms);
  RDKit::Conformer conf = mol.getConformer(confId);
  for (unsigned int i = 0; i < d_nAtoms; ++i) {
    d_charges.push_back(mol.getAtomWithIdx(i)->getProp<double>(prop));
    const RDGeom::Point3D &pt = conf.getAtomPos(i);
    d_pos.push_back(pt.x);
    d_pos.push_back(pt.y);
    d_pos.push_back(pt.z);
  }
  if (fabs(d_alpha) < MIN_CUTOFF_VAL) {
    d_softcore = false;
    if (d_cutoff < MIN_CUTOFF_VAL) {
      d_cutoff = CUTOFF;
    }
  } else {
    d_softcore = true;
  }
}

double Coulomb::operator()(double x, double y, double z, double thres) const {
  double res = 0.0, dist2, temp;
  if (d_softcore) {
    for (unsigned int i = 0, j = 0; i < d_nAtoms; i++) {
      temp = x - d_pos[j++];
      dist2 = temp * temp;
      temp = y - d_pos[j++];
      dist2 += temp * temp;
      temp = z - d_pos[j++];
      dist2 += temp * temp;
      if (dist2 < thres) {
        res += d_charges[i] * (1.0 / sqrt(d_alpha + dist2));
      }
    }
  } else {
    for (unsigned int i = 0, j = 0; i < d_nAtoms; i++) {
      temp = x - d_pos[j++];
      dist2 = temp * temp;
      temp = y - d_pos[j++];
      dist2 += temp * temp;
      temp = z - d_pos[j++];
      dist2 += temp * temp;
      if (dist2 < thres) {
        dist2 = std::max(dist2, d_cutoff);
        res += d_charges[i] * (1.0 / sqrt(dist2));
      }
    }
  }
  res *= CoulombDetail::prefactor * d_probe;
  if (d_absVal) {
    res = -fabs(
        res);  // takes the negative absolute value of the interaction energy
  }
  return res;
}

CoulombDielectric::CoulombDielectric(
    const std::vector<double> &charges,
    const std::vector<RDGeom::Point3D> &positions, double probeCharge,
    bool absVal, double alpha, double cutoff, double epsilon, double xi)
    : d_nAtoms(charges.size()),
      d_absVal(absVal),
      d_cutoff(cutoff * cutoff),
      d_probe(probeCharge),
      d_epsilon(epsilon),
      d_xi(xi),
      d_alpha(alpha),
      d_charges(charges) {
  PRECONDITION(d_charges.size() == positions.size(),
               "Lengths of position and charge vectors do not match.");
  d_dielectric = (d_xi - d_epsilon) / (d_xi + d_epsilon);
  std::vector<unsigned int> neighbors(positions.size(), 0);

  d_dists.resize(d_nAtoms);
  d_sp.reserve(d_nAtoms);
  d_pos.reserve(3 * d_nAtoms);
  for (unsigned int i = 0; i < positions.size(); ++i) {
    d_pos.push_back(positions[i].x);
    d_pos.push_back(positions[i].y);
    d_pos.push_back(positions[i].z);
    for (unsigned int j = i + 1; j < positions.size(); ++j) {
      double dis = (positions[j] - positions[i]).length();
      if (dis < 4.0) {
        ++neighbors[i];
        ++neighbors[j];
      }
    }
  }
  for (unsigned int neighbor : neighbors) {
    switch (neighbor) {
      case 0:
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      case 6:
        d_sp.push_back(0.0);
        break;
      case 7:
        d_sp.push_back(.4);
        break;
      case 8:
        d_sp.push_back(.9);
        break;
      case 9:
        d_sp.push_back(1.4);
        break;
      case 10:
        d_sp.push_back(1.9);
        break;
      case 11:
        d_sp.push_back(2.6);
        break;
      default:
        d_sp.push_back(4.0);
    }
  }
  if (fabs(d_alpha) < MIN_CUTOFF_VAL) {
    d_softcore = false;
    if (d_cutoff < MIN_CUTOFF_VAL * MIN_CUTOFF_VAL) {
      d_cutoff = CUTOFF * CUTOFF;
    }
  } else {
    d_softcore = true;
  }
}

CoulombDielectric::CoulombDielectric(const RDKit::ROMol &mol, int confId,
                                     double probeCharge, bool absVal,
                                     const std::string &prop, double alpha,
                                     double cutoff, double epsilon, double xi)
    : d_nAtoms(mol.getNumAtoms()),
      d_absVal(absVal),
      d_cutoff(cutoff * cutoff),
      d_probe(probeCharge),
      d_epsilon(epsilon),
      d_xi(xi),
      d_alpha(alpha) {
  PRECONDITION(mol.getNumConformers() > 0, "No Conformers for Molecule");

  d_charges.reserve(d_nAtoms);
  d_sp.reserve(d_nAtoms);
  d_dists.resize(d_nAtoms);
  d_pos.reserve(3 * d_nAtoms);

  RDKit::Conformer conf = mol.getConformer(confId);
  for (unsigned int i = 0; i < d_nAtoms; ++i) {
    d_charges.push_back(mol.getAtomWithIdx(i)->getProp<double>(prop));
    const RDGeom::Point3D &pt = conf.getAtomPos(i);
    d_pos.push_back(pt.x);
    d_pos.push_back(pt.y);
    d_pos.push_back(pt.z);
  }
  d_dielectric = (d_xi - d_epsilon) / (d_xi + d_epsilon);

  std::vector<unsigned int> neighbors(d_charges.size(), 0);

  for (unsigned int i = 0; i < d_charges.size() - 1; ++i) {
    for (unsigned int j = i + 1; j < d_charges.size(); ++j) {
      double temp = d_pos[j * 3] - d_pos[i * 3];
      double dist2 = temp * temp;
      temp = d_pos[j * 3 + 1] - d_pos[i * 3 + 1];
      dist2 += temp * temp;
      temp = d_pos[j * 3 + 2] - d_pos[i * 3 + 2];
      dist2 += temp * temp;
      if (dist2 < 16.0) {
        neighbors[i]++;
        neighbors[j]++;
      }
    }
  }
  for (unsigned int neighbor : neighbors) {
    switch (neighbor) {
      case 0:
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      case 6:
        d_sp.push_back(0.0);
        break;
      case 7:
        d_sp.push_back(.4);
        break;
      case 8:
        d_sp.push_back(.9);
        break;
      case 9:
        d_sp.push_back(1.4);
        break;
      case 10:
        d_sp.push_back(1.9);
        break;
      case 11:
        d_sp.push_back(2.6);
        break;
      default:
        d_sp.push_back(4.0);
    }
  }
  if (fabs(d_alpha) < MIN_CUTOFF_VAL) {
    d_softcore = false;
    if (d_cutoff < MIN_CUTOFF_VAL * MIN_CUTOFF_VAL) {
      d_cutoff = CUTOFF * CUTOFF;
    }
  } else {
    d_softcore = true;
  }
}

double CoulombDielectric::operator()(double x, double y, double z,
                                     double thres) const {
  int neigh = 0;
  double res = 0.0, sq = 0.0;

  for (unsigned int i = 0, j = 0; i < d_nAtoms; ++i, j += 3) {
    double temp = x - d_pos[j];
    double dist2 = temp * temp;
    temp = y - d_pos[j + 1];
    dist2 += temp * temp;
    temp = z - d_pos[j + 2];
    dist2 += temp * temp;
    d_dists[i] = dist2;
    if (dist2 < 16.0) {
      neigh += 1;
    }
  }

  switch (neigh) {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
      sq = 0.0;
      break;
    case 7:
      sq = .4;
      break;
    case 8:
      sq = .9;
      break;
    case 9:
      sq = 1.4;
      break;
    case 10:
      sq = 1.9;
      break;
    case 11:
      sq = 2.6;
      break;
    default:
      sq = 4.0;
  }

  if (d_softcore) {
    for (unsigned int i = 0; i < d_nAtoms; i++) {
      if (d_dists[i] < thres) {
        res += d_charges[i] *
               (1 / sqrt(d_alpha + d_dists[i]) +
                d_dielectric / sqrt(d_alpha + d_dists[i] + 4.0 * d_sp[i] * sq));
      }
    }
  } else {
    for (unsigned int i = 0; i < d_nAtoms; i++) {
      if (d_dists[i] < thres) {
        double dist2 = std::max(d_dists[i], d_cutoff);
        res += d_charges[i] * (1 / sqrt(dist2) +
                               d_dielectric / sqrt(dist2 + 4.0 * d_sp[i] * sq));
      }
    }
  }
  res *= CoulombDetail::prefactor * (1 / d_xi) * d_probe;

  if (d_absVal) {
    res = -fabs(
        res);  // takes the negative absolute value of the interaction energy
  }
  return res;
}

namespace HBondDetail {
const double em[2][2] = {{-8.368, -11.715}, {-11.715, -16.736}};
const double rm[2][2] = {{3.2, 3.0}, {3.0, 2.8}};
const double K1 = 562.25380293;
const double K2 = 0.11697778;
const double bondlength[2] = {-0.972, -1.019};

double cos_2(double t, double, double) {
  if (t < M_PI_2) {
    double temp = cos(t);
    temp *= temp;
    return temp;
  } else {
    return 0.0;
  }
};

double cos_2_0(double t, double t_0 = 0.0, double t_i = 1.0) {
  double temp;
  if (t_i < M_55D) {
    temp = cos(t_0);
    temp *= temp;
    return temp;
  } else if (t_i < M_PI_2) {
    temp = cos(t);
    temp *= temp;
    return temp;
  } else {
    return 0.0;
  }
};

double cos_2_rot(double t, double, double) {
  t -= M_70_5D;
  if (t < M_PI_2) {
    double temp = cos(t);
    temp *= temp;
    return temp;
  } else {
    return 0.0;
  }
};

double cos_4(double t, double, double) {
  if (t < M_PI_2) {
    double temp = cos(t);
    temp *= temp;
    temp *= temp;
    return temp;
  } else {
    return 0.0;
  }
};

double cos_4_rot(double t, double, double) {
  t -= M_70_5D;
  if (t < M_PI_2) {
    double temp = cos(t);
    temp *= temp;
    temp *= temp;
    return temp;
  } else {
    return 0.0;
  }
};

double cos_6(double t, double, double) {
  if (t < M_PI_2) {
    double temp = cos(t);
    temp *= temp;
    temp *= temp * temp;
    return temp;
  } else {
    return 0.0;
  }
};

double cos_6_rot(double t, double, double) {
  t -= M_70_5D;
  if (t < M_PI_2) {
    double temp = cos(t);
    temp *= temp;
    temp *= temp * temp;
    return temp;
  } else {
    return 0.0;
  }
};

double cos_acc(double, double t_0, double t_i) {
  double temp;
  if (t_i < M_PI_2) {
    temp = cos(t_0) * (0.9 + 0.1 * sin(2 * t_i));
    return temp;
  } else if (t_i < M_110D) {
    temp = cos(t_i);
    temp *= temp;
    temp = K2 - temp;
    temp *= temp * temp;
    temp *= cos(t_0) * K1;
    return temp;
  } else {
    return 0.0;
  }
};

double no_dep(double, double, double) { return 1.0; };
}  // namespace HBondDetail

HBond::HBond(const RDKit::ROMol &mol, int confId,
             const std::string &probeAtomType, bool fixed, double cutoff)
    : d_cutoff(cutoff * cutoff) {
  if (d_cutoff < (MIN_CUTOFF_VAL * MIN_CUTOFF_VAL)) {
    d_cutoff = CUTOFF * CUTOFF;
  }

  if (probeAtomType == "O") {
    d_DAprop = 'D';
    d_probetype = O;
  } else if (probeAtomType == "OH") {
    d_DAprop = 'A';
    d_probetype = O;
  } else if (probeAtomType == "N") {
    d_DAprop = 'D';
    d_probetype = N;
  } else if (probeAtomType == "NH") {
    d_DAprop = 'A';
    d_probetype = N;
  } else {
    const std::string msg = "Probe atom type not supported: " + probeAtomType;
    BOOST_LOG(rdErrorLog) << msg << std::endl;
    throw ValueErrorException(msg);
  }

  d_nInteract = mol.getNumAtoms();  // number of atoms = highest possible number
                                    // of interactions

  std::vector<unsigned int> specialAtoms;
  findSpecials(mol, confId, fixed, specialAtoms);

  if (d_DAprop == 'A') {
    if (fixed) {
      findAcceptors(mol, confId, specialAtoms);
    } else {
      findAcceptorsUnfixed(mol, confId, specialAtoms);
    }
  } else if (d_DAprop == 'D') {
    if (fixed) {
      findDonors(mol, confId, specialAtoms);
    } else {
      findDonorsUnfixed(mol, confId, specialAtoms);
    }
  } else {  // this should never be the case
    BOOST_LOG(rdErrorLog) << "HBond: unknown target property d_DAprop: "
                          << d_DAprop << std::endl;
  }

  d_nInteract = d_targettypes.size();  // updated to number of interactions
  d_eneContrib.resize(d_nInteract, 0);
  d_vectTargetProbe.resize(d_nInteract * 3, 0);

  POSTCONDITION(
      d_nInteract * 3 == d_pos.size(),
      "Error in constructing H-Bond descriptor: Vector length mismatch (target atom types).");
  POSTCONDITION(
      d_nInteract * 3 == d_direction.size(),
      "Error in constructing H-Bond descriptor: Vector length mismatch (bond directions).");
  POSTCONDITION(
      d_nInteract == d_function.size(),
      "Error in constructing H-Bond descriptor: Vector length mismatch (angular functions).");
  POSTCONDITION(
      d_nInteract * 3 == d_plane.size(),
      "Error in constructing H-Bond descriptor: Vector length mismatch (lone pair planes).");
  POSTCONDITION(
      d_nInteract == d_lengths.size(),
      "Error in constructing H-Bond descriptor: Vector length mismatch (bondlengths).");
}

unsigned int HBond::findSpecials(const RDKit::ROMol &mol, int confId,
                                 bool fixed,
                                 std::vector<unsigned int> &specials) {
  using namespace HBondDetail;

  RDKit::MatchVectType matches;
  RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
  RDGeom::Point3D pos, dir, hbonddir, plane, bondDirection[12];
  unsigned int nbrs;
  unsigned int match = 0, nMatches = 0;

  const RDKit::Conformer &conf = mol.getConformer(confId);
  // RDKit::RWMol thr = *RDKit::SmilesToMol("C[C@H]([C@@H](C(=O))N)O");
  // //threonine, serine itself is a substructure of this
  static const auto ser = RDKit::v2::SmilesParse::MolFromSmiles(
      "[CH2]([C@@H](C(=O))N)O");  // serine
  // RDKit::RWMol his = *RDKit::SmilesToMol("Cc1cnc[nH]1"); //imidazole residue,
  // is correctly taken into account in 'normal' treatment

  match = RDKit::SubstructMatch(mol, *ser, matches);
  nMatches += match;

  for (auto &matche : matches) {
    const auto *atom = mol.getAtomWithIdx(matche.second);
    if (atom->getAtomicNum() == 8) {
      boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
      pos = conf.getAtomPos(matche.second);
      if (d_DAprop == 'A') {
        nbrs = 0;
        while (nbrIdx != endNbrs) {  // loop over atoms
          if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() == 1) {
            hbonddir = pos - conf.getAtomPos(*nbrIdx);
            hbonddir.normalize();
          } else {
            bondDirection[nbrs] = pos - conf.getAtomPos(*nbrIdx);
            bondDirection[nbrs].normalize();
            ++nbrs;
          }
          ++nbrIdx;
        }
        if (nbrs == 2) {
          dir = bondDirection[0] + hbonddir;
          plane =
              bondDirection[0].crossProduct(hbonddir);  // X-O-Y plane vector
          plane = plane.crossProduct(dir);              // lp plane vector
          if (fixed) {
            addVectElements(O, &cos_2_0, pos, dir, plane);
          } else {
            addVectElements(O, &cos_4_rot, pos, bondDirection[0]);
          }
          specials.push_back(matche.second);
        }
      } else if (d_DAprop == 'D') {
        if (fixed) {
          while (nbrIdx != endNbrs) {
            if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() ==
                1) {  // hydrogen neighbor necessary for H bond donation
              dir = conf.getAtomPos(*nbrIdx) - pos;
              addVectElements(O, &cos_6, pos, dir);
            }
            ++nbrIdx;
          }
        } else {
          while (nbrIdx != endNbrs) {
            if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() !=
                1) {  // search for non-hydrogen bond direction
              dir = pos - conf.getAtomPos(*nbrIdx);
              addVectElements(O, &cos_6_rot, pos, dir);
            }
            ++nbrIdx;
          }
        }
        specials.push_back(matche.second);
      } else {  // this should never be the case
        BOOST_LOG(rdErrorLog)
            << "HBond::operator(): unknown target property d_DAprop: "
            << d_DAprop << std::endl;
      }
    }
  }
  return nMatches;
}

/* General structure of findAcceptors, findAcceptorsUnfixed, findDonors,
 * findDonorsUnfixed functions: loop over all atoms:
 *   - check whether atom was already treated in specialAtoms
 *   - switch ( atomicNum ): find atoms which are able to donate/accept
 * hydrogen bonds if able:
 *   - switch ( number of neighbors ): find the right geometry
 *   - check for charge or aromatic neighborhood
 *   - calculate adequate hydrogen bond direction vector,
 * lone pair plane vector
 *   - choose correct angular function
 *   - add atomtype, angular function, position of atom,
 * hydrogen bond direction vector and lone pair plane vector to vectors for
 * calculating the interaction returns number of added interactions
 */
unsigned int HBond::findAcceptors(const RDKit::ROMol &mol, int confId,
                                  const std::vector<unsigned int> &specials) {
  using namespace HBondDetail;

  const RDKit::Conformer &conf =
      mol.getConformer(confId);  // get conformer of molecule
  RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
  RDKit::ROMol::ADJ_ITER secnbrIdx, secendNbrs;
  RDGeom::Point3D pos, dir, plane, bondDirection[12];
  unsigned int nbrs;          // no of neigbors
  unsigned int interact = 0;  // no of interactions
  bool aromaticnbr;           // aromatic neighbor atom?

  for (unsigned int i = 0; i < d_nInteract; i++) {  // loop over all atoms
    if (std::find(specials.begin(), specials.end(), i) !=
        specials.end()) {  // check whether atom was already treated specially
      interact++;
      continue;
    }

    const auto *atom = mol.getAtomWithIdx(i);  // get ptr to atom

    switch (atom->getAtomicNum()) {  // find atoms able to donate hydrogen bonds
      case 7:                        // Nitrogen
        boost::tie(nbrIdx, endNbrs) =
            mol.getAtomNeighbors(atom);  // get neighbors
        pos = conf.getAtomPos(i);        // get position of atom

        nbrs = 0;
        while (nbrIdx != endNbrs) {  // loop over all neigbors
          bondDirection[nbrs] =
              pos - conf.getAtomPos(*nbrIdx);  // store bond vectors in an array
          bondDirection[nbrs].normalize();
          ++nbrs;  // count neigbors
          ++nbrIdx;
        }

        switch (nbrs) {  // number of neigbors
          case 1:        // sp, eg. nitriles
            if (atom->getFormalCharge() <=
                0) {  // no positively charged nitrogens
              addVectElements(
                  N, &cos_2, pos,
                  bondDirection[0]);  // no differentiation between in-plane and
                                      // out-of-plane angle, plane not needed
              interact++;
            }
            break;
          case 2:  // eg. imines, heterocycles
            if (atom->getFormalCharge() <=
                0) {  // no positively charged nitrogen
              dir = bondDirection[0] +
                    bondDirection[1];  // get hydrogen bond direction
              plane = bondDirection[0].crossProduct(
                  bondDirection[1]);  // normal vector of lone pair plane =
                                      // plane of three atoms
              addVectElements(N, &cos_2, pos, dir, plane);
              interact++;
            }
            break;
          case 3:  // amines, iminium ions
            if (atom->getFormalCharge() <= 0 &&
                !(RDKit::MolOps::atomHasConjugatedBond(
                    atom))) {  // no positively charged nitrogen, no conjugated
                               // nitrogen (amide bonds!)
              dir = bondDirection[0] + bondDirection[1] +
                    bondDirection[2];  // get hydrogen bond direction
              addVectElements(N, &cos_2, pos,
                              dir);  // no differentiation between in-plane and
                                     // out-of-plane angle, plane not needed
              interact++;
            }
            break;
          case 0:  // unbound nitrogen atoms, no hydrogen bonding
          case 4:  // ammonium ions, no hydrogen bonding
            break;
          default:  // more than four bonds: not possible with nitrogen
            BOOST_LOG(rdErrorLog)
                << "HBond: Nitrogen atom bound to more than 4 neighbor atoms: Atom: "
                << i << std::endl;
        }
        break;
      case 8:  // Oxygen
        boost::tie(nbrIdx, endNbrs) =
            mol.getAtomNeighbors(atom);  // get neighbors
        pos = conf.getAtomPos(i);        // get position of atom

        nbrs = 0;
        aromaticnbr = false;
        while (nbrIdx != endNbrs) {  // loop over neighbors
          bondDirection[nbrs] =
              pos - conf.getAtomPos(*nbrIdx);  // store bond vectors in an array
          bondDirection[nbrs].normalize();     // normalization
          if (mol.getAtomWithIdx(*nbrIdx)
                  ->getIsAromatic()) {  // check whether neighbor is aromatic
                                        // (phenolic oxygens)
            aromaticnbr = true;
          }
          ++nbrs;  // count neighbors
          ++nbrIdx;
        }

        switch (nbrs) {  // no of neighbors
          case 1:        // carbonyl, carboxyl C=O, X=O (X=S,P,...), anions
                         // (alcoholates, carboxylates)
            --nbrIdx;
            boost::tie(secnbrIdx, secendNbrs) = mol.getAtomNeighbors(
                mol.getAtomWithIdx(*nbrIdx));  // get neighbors of neighbor atom
            while (secnbrIdx !=
                   secendNbrs) {  // loop over neighbors of neighbor atom
              if ((*secnbrIdx) !=
                  i) {  // second neighbor should not be the oxygen itself
                // bond direction of neighbor atom to second neighbor atom
                dir = conf.getAtomPos(*secnbrIdx) - conf.getAtomPos(*nbrIdx);
                break;  // we only need one bond vector
              }
              ++secnbrIdx;
            }
            plane = bondDirection[0].crossProduct(dir);  // lp plane vector
            if (atom->getFormalCharge() == 0) {  // carbonyl, carboxyl C=O, X=O
              addVectElements(O, &cos_acc, pos, bondDirection[0], plane);
              interact++;
            } else if (atom->getFormalCharge() ==
                       -1) {  // anion, eg. phenolate, carboxylic acid
              if (RDKit::MolOps::atomHasConjugatedBond(atom)) {
                // charged oxygen in conjungated system, eg phenolates or
                // carboxylates
                addVectElements(O, &cos_acc, pos, bondDirection[0], plane);
              } else {
                // non-conjugated anion, eg sp3-alcoholate
                addVectElements(O, &cos_2, pos, bondDirection[0], plane);
              }
              interact++;
            }
            break;
          case 2:  // alcohols, ethers, carboxyl OH
            dir = bondDirection[0] + bondDirection[1];  // get bond direction
            plane = bondDirection[0].crossProduct(
                bondDirection[1]);            // X-O-Y plane vector
            plane = plane.crossProduct(dir);  // lp plane vector
            if (aromaticnbr) {  // if oxygen bound to aromatic system, eg.
                                // phenol
              addVectElements(O, &cos_2, pos, dir, plane);
            } else {  // all other
              addVectElements(O, &cos_2_0, pos, dir, plane);
            }
            interact++;
            break;
          case 0:  // oxygen atoms, no hydrogen bonding
          case 3:  // only with positively charged atom possible, no hydrogen
                   // bonding
            break;
          default:  // more than 3 neighbors: not possible with oxygen
            BOOST_LOG(rdErrorLog)
                << "HBond: Oxygen atom bound to more than 3 neighbor atoms: Atom: "
                << i << std::endl;
        }
        break;
      // Halogens
      case 9:                                  // F
      case 17:                                 // Cl
      case 35:                                 // Br
      case 53:                                 // I
        pos = conf.getAtomPos(i);              // get position of atom
        dir = RDGeom::Point3D(1.0, 1.0, 1.0);  // no directionality needed
        addVectElements(
            N, &no_dep, pos,
            dir);  // type of halogens ~ nitrogen; no lp plane needed
        interact++;
        break;
      default:
        break;
    }
  }
  return interact;
}

unsigned int HBond::findAcceptorsUnfixed(
    const RDKit::ROMol &mol, int confId,
    const std::vector<unsigned int> &specials) {
  using namespace HBondDetail;

  const RDKit::Conformer &conf =
      mol.getConformer(confId);  // get conformer of molecule
  RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
  RDKit::ROMol::ADJ_ITER secnbrIdx, secendNbrs;
  RDGeom::Point3D pos, dir, plane, bondDirection[12];
  unsigned int nbrs,
      nonhnbrs;               // no of neighbors, no of nonhydrogen - neighbors
  unsigned int interact = 0;  // no of interactions
  bool aromaticnbr;           // aromatic neighbor atom?

  for (unsigned int i = 0; i < d_nInteract; i++) {  // loop over all atoms
    if (std::find(specials.begin(), specials.end(), i) !=
        specials.end()) {  // check whether atom was already treated specially
      interact++;
      continue;
    }

    const auto *atom = mol.getAtomWithIdx(i);  // get pointer to atom

    switch (atom->getAtomicNum()) {  // find atoms able to accept hydrogen bonds
                                     // (O, N, halogens)
      case 7:                        // Nitrogen
        boost::tie(nbrIdx, endNbrs) =
            mol.getAtomNeighbors(atom);  // get neighbors
        pos = conf.getAtomPos(i);        // get position of atom

        nbrs = 0;
        nonhnbrs = 0;
        while (nbrIdx != endNbrs) {  // loop over neighbors
          if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() !=
              1) {  // hydrogen bond directions are not included
            bondDirection[nonhnbrs] =
                pos -
                conf.getAtomPos(*nbrIdx);  // store bond direction in an array
            bondDirection[nonhnbrs].normalize();
            ++nonhnbrs;  // count non hydrogen neighbors
          }
          ++nbrs;  // count neighbors
          ++nbrIdx;
        }

        switch (nbrs) {  // number of neigbors
          case 1:        // sp, eg. nitriles,   no difference to fixed bonds
            if (atom->getFormalCharge() <=
                0) {  // no positively charged nitrogens
              addVectElements(
                  N, &cos_2, pos,
                  bondDirection[0]);  // no differentiation between in-plane and
                                      // out-of-plane angle, plane not needed
              interact++;
            }
            break;
          case 2:  // eg. imines, heterocycles
            if (atom->getFormalCharge() <=
                0) {                // no positively charged nitrogen
              if (nonhnbrs == 2) {  // secondary imines, heterocycles
                dir = bondDirection[0] +
                      bondDirection[1];  // get hydrogen bond direction
                plane = bondDirection[0].crossProduct(
                    bondDirection[1]);  // normal vector of lone pair plane =
                                        // plane of three atoms
                addVectElements(N, &cos_2, pos, dir, plane);
              } else if (nonhnbrs == 1) {  // primary imine, hydrogen is allowed
                                           // to swap places
                nbrIdx -= 2;
                nbrs = nonhnbrs;
                while (nbrIdx != endNbrs) {
                  if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() ==
                      1) {  // only bond directions to hydrogens are included
                    bondDirection[nbrs] = pos - conf.getAtomPos(*nbrIdx);
                    bondDirection[nbrs].normalize();
                    ++nbrs;
                  }
                  ++nbrIdx;
                }
                plane = bondDirection[0].crossProduct(
                    bondDirection[1]);  // normal vector of lone pair plane =
                                        // plane of three atoms
                addVectElements(N, &cos_acc, pos, bondDirection[0], plane);
              } else {  //[N-]H2 rotating
                addVectElements(N, &no_dep, pos,
                                RDGeom::Point3D(0.0, 0.0, 0.0));
              }
              interact++;
            }
            break;
          case 3:  // amines, iminium ions
            if (atom->getFormalCharge() <=
                0) {  // no iminium ions, no positively charged nitrogen
              if (nonhnbrs == 0) {  // ammonia
                addVectElements(N, &no_dep, pos,
                                RDGeom::Point3D(0.0, 0.0, 0.0));
              } else if (nonhnbrs == 1) {  // primary amines, rotation
                addVectElements(
                    N, &cos_2_rot, pos,
                    bondDirection[0]);  // no differentiation between in-plane
                                        // and out-of-plane angle, plane not
                                        // needed
              } else {  // secondary amines (not flexible) and tertiary amines,
                        // same as fixed hydrogen bonds
                if (atom->getFormalCharge() <= 0 &&
                    !(RDKit::MolOps::atomHasConjugatedBond(
                        atom))) {  // positively charged nitrogen, no conjugated
                                   // nitrogen (amide bonds!)
                  nbrIdx -= 3;
                  nbrs = nonhnbrs;
                  while (nbrIdx != endNbrs) {
                    if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() ==
                        1) {  // only bond directions to hydrogens are included
                      bondDirection[nbrs] = pos - conf.getAtomPos(*nbrIdx);
                      bondDirection[nbrs].normalize();
                      ++nbrs;
                    }
                    ++nbrIdx;
                  }
                  dir = bondDirection[0] + bondDirection[1] +
                        bondDirection[2];  // hydrogen bond direction
                  addVectElements(
                      N, &cos_2, pos,
                      dir);  // no differentiation between in-plane and
                             // out-of-plane angle, plane not needed
                }
              }
              interact++;
            }
            break;
          case 0:  // unbound nitrogen atoms, no hydrogen bonding
          case 4:  // ammonium ions, no hydrogen bonding
            break;
          default:  // more than four bonds: not possible with nitrogen
            BOOST_LOG(rdErrorLog)
                << "HBond: Nitrogen atom bound to more than 4 neighbor atoms: Atom: "
                << i << std::endl;
        }
        break;
      case 8:  // Oxygen
        boost::tie(nbrIdx, endNbrs) =
            mol.getAtomNeighbors(atom);  // get neighbors
        pos = conf.getAtomPos(i);        // get atom position

        nbrs = 0;
        nonhnbrs = 0;
        aromaticnbr = false;
        while (nbrIdx != endNbrs) {  // loop over neighbors
          if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() !=
              1) {  // bond directions to hydrogen are not included
            bondDirection[nonhnbrs] = pos - conf.getAtomPos(*nbrIdx);
            bondDirection[nonhnbrs].normalize();
            ++nonhnbrs;  // count non hydrogen neighbors
          }
          if (mol.getAtomWithIdx(*nbrIdx)
                  ->getIsAromatic()) {  // check whether aromatic neighbor
            aromaticnbr = true;
          }
          ++nbrs;  // count neighbors
          ++nbrIdx;
        }

        switch (nbrs) {  // no of neighbors
          case 1:        // carbonyl, carboxyl C=O, X=O (X=S,P,...), anions
                         // (alcoholates, carboxylates)
            --nbrIdx;
            boost::tie(secnbrIdx, secendNbrs) = mol.getAtomNeighbors(
                mol.getAtomWithIdx(*nbrIdx));  // get neighbors of neighbor atom
            while (secnbrIdx !=
                   secendNbrs) {  // loop over neighbors of neighbor atom
              if ((*secnbrIdx) !=
                  i) {  // second neighbor should not be the oxygen itself
                // bond direction of neighbor atom to second neighbor atom
                dir = conf.getAtomPos(*secnbrIdx) - conf.getAtomPos(*nbrIdx);
                break;  // we only need one vector
              }
              ++secnbrIdx;
            }
            plane = bondDirection[0].crossProduct(dir);  // lp plane vector
            if (atom->getFormalCharge() == 0) {  // carbonyl, carboxyl C=O, X=O
              addVectElements(O, &cos_acc, pos, bondDirection[0], plane);
              interact++;
            } else if (atom->getFormalCharge() ==
                       -1) {  // anion, eg. phenolate, carboxylic acid
              if (RDKit::MolOps::atomHasConjugatedBond(atom)) {
                // charged oxygen in conjungated system, eg phenolates or
                // carboxylates
                addVectElements(O, &cos_acc, pos, bondDirection[0], plane);
              } else {  // all other
                addVectElements(O, &cos_2, pos, bondDirection[0], plane);
              }
              interact++;
            }
            break;
          case 2:                 // alcohols, ethers, carboxyl OH
            if (nonhnbrs == 0) {  // water
              addVectElements(O, &no_dep, pos, RDGeom::Point3D(0.0, 0.0, 0.0));
            } else if (nonhnbrs == 1) {  // hydroxy groups
              if (aromaticnbr) {         // phenol
                nbrIdx -= 2;
                nbrs = nonhnbrs;
                while (nbrIdx != endNbrs) {  // loop over neighbors
                  if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() ==
                      1) {  // only bond directions to hydrogens are included
                    bondDirection[nbrs] =
                        pos - conf.getAtomPos(
                                  *nbrIdx);  // get O-H bond direction vector
                    bondDirection[nbrs].normalize();
                    ++nbrs;
                  }
                  ++nbrIdx;
                }
                plane = bondDirection[0].crossProduct(
                    bondDirection[1]);  // X-O-Y plane vector
                addVectElements(O, &cos_acc, pos, bondDirection[0], plane);
              } else {
                addVectElements(
                    O, &cos_2_rot, pos,
                    bondDirection[0]);  // no plane information needed
              }
            } else {  // ethers, same as fixed hydrogen bonds
              dir = bondDirection[0] +
                    bondDirection[1];  // get hydrogen bond direction
              plane = bondDirection[0].crossProduct(
                  bondDirection[1]);            // X-O-Y plane vector
              plane = plane.crossProduct(dir);  // lp plane vector
              addVectElements(O, &cos_2_0, pos, dir, plane);
            }
            interact++;
            break;
          case 0:  // oxygen atoms, no hydrogen bonding
          case 3:  // only with positively charged atom possible, no hydrogen
                   // bonding
            break;
          default:  // more than 3 neighbors: not possible with oxygen
            BOOST_LOG(rdErrorLog)
                << "HBond: Oxygen atom bound to more than 3 neighbor atoms: Atom: "
                << i << std::endl;
        }
        break;
      // Halogens
      case 9:                                  // F
      case 17:                                 // Cl
      case 35:                                 // Br
      case 53:                                 // I
        pos = conf.getAtomPos(i);              // get atoms position
        dir = RDGeom::Point3D(1.0, 1.0, 1.0);  // no directionality needed
        addVectElements(
            N, &no_dep, pos,
            dir);  // type of halogens ~ nitrogen; no lp plane needed
        interact++;
        break;
      default:
        break;
    }
  }
  return interact;
}

unsigned int HBond::findDonors(const RDKit::ROMol &mol, int confId,
                               const std::vector<unsigned int> &specials) {
  using namespace HBondDetail;

  const RDKit::Conformer &conf =
      mol.getConformer(confId);  // get conformer of molecule
  RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
  RDGeom::Point3D pos, dir;
  unsigned int interact = 0;

  for (unsigned int i = 0; i < d_nInteract; i++) {  // loop over all atoms
    if (std::find(specials.begin(), specials.end(), i) !=
        specials.end()) {  // check whether atom was already treated specially
      interact++;
      continue;
    }

    const auto *atom = mol.getAtomWithIdx(i);  // get ptr to atom

    switch (
        atom->getAtomicNum()) {  // find atoms able to donate hydrogen bonds (O,
                                 // N, of course only with attached hydrogen)
      case 7:                    // Nitrogen
        boost::tie(nbrIdx, endNbrs) =
            mol.getAtomNeighbors(atom);  // get neigbors
        pos = conf.getAtomPos(i);        // get position

        while (nbrIdx != endNbrs) {  // loop over neighbors
          if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() ==
              1) {  // hydrogen neighbor necessary for H bond donation
            dir = conf.getAtomPos(*nbrIdx) -
                  pos;  // bond direction vector, IMPORTANT: no normalization,
                        // because operator() needs length of vector in case of
                        // donors
            addVectElements(N, &cos_2, pos, dir);
            interact++;
          }
          ++nbrIdx;
        }
        break;
      case 8:  // Oxygen
        boost::tie(nbrIdx, endNbrs) =
            mol.getAtomNeighbors(atom);  // get neigbors
        pos = conf.getAtomPos(i);        // get position

        while (nbrIdx != endNbrs) {  // loop over neighbors
          if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() ==
              1) {  // hydrogen neighbor necessary for H bond donation
            dir = conf.getAtomPos(*nbrIdx) -
                  pos;  // bond direction vector, IMPORTANT: no normalization,
                        // because operator() needs length of vector in case of
                        // donors
            addVectElements(O, &cos_4, pos, dir);
            interact++;
          }
          ++nbrIdx;
        }
        break;
      default:
        break;
    }
  }
  return interact;
}

unsigned int HBond::findDonorsUnfixed(
    const RDKit::ROMol &mol, int confId,
    const std::vector<unsigned int> &specials) {
  using namespace HBondDetail;

  const auto &conf = mol.getConformer(confId);  // get conformer of molecule
  RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
  RDGeom::Point3D pos, hbonddir, dir, plane;
  unsigned int nbrs, nonhnbrs;  // no of neighbors, no of non hydrogen neighbors
  unsigned int interact = 0;
  bool aromaticnbr;

  for (unsigned int i = 0; i < d_nInteract; i++) {  // loop over all atoms
    if (std::find(specials.begin(), specials.end(), i) !=
        specials.end()) {  // check whether atom was already treated specially
      interact++;
      continue;  // skip loop for this atom
    }

    const auto *atom = mol.getAtomWithIdx(i);  // get ptr to atom

    switch (
        atom->getAtomicNum()) {  // find atoms able to donate hydrogen bonds (O,
                                 // N, of course only with attached hydrogen)
      case 7:                    // Nitrogen
        boost::tie(nbrIdx, endNbrs) =
            mol.getAtomNeighbors(atom);  // get neighbors
        pos = conf.getAtomPos(i);        // get position

        nbrs = 0;
        nonhnbrs = 0;
        while (nbrIdx != endNbrs) {  // loop over neighbors
          if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() != 1) {
            ++nonhnbrs;  // count non-hydrogen neighbors
          }
          ++nbrs;  // count neighbors
          ++nbrIdx;
        }
        nbrIdx -= nbrs;

        if (nonhnbrs != nbrs) {  // otherwise no hydrogens, no donation possible
          switch (nbrs) {        // number of neigbors
            case 2:              // eg. imines, heterocycles
              if (nonhnbrs == 0) {  //[N-]H2
                addVectElements(N, &no_dep, pos,
                                RDGeom::Point3D(0.0, 0.0, 0.0));
                interact++;
              } else {  // primary imine, swapping of hydrogen is allowed
                while (nbrIdx != endNbrs) {  // loop over neighbors
                  if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() != 1) {
                    dir = pos - conf.getAtomPos(*nbrIdx);  // the other bond
                    dir.normalize();

                  } else {
                    hbonddir = conf.getAtomPos(*nbrIdx) - pos;  // hydrogen bond
                  }
                  ++nbrIdx;
                }

                addVectElements(N, &cos_2, pos, hbonddir);  // first hbond
                interact++;

                // let's swap hydrogen and lp:
                plane = dir.crossProduct(hbonddir);  // lp plane vector
                plane =
                    plane.crossProduct(dir);  // plane through other bond
                                              // perpendicular to X-N-H plane
                dir = plane *
                      hbonddir.dotProduct(plane);  // projection of vector
                                                   // hbonddir on plane vector
                hbonddir -= dir * 2;  // mirroring of dir vector on plane
                addVectElements(
                    N, &cos_2, pos,
                    hbonddir);  // second hbond, hydrogen at other place
                interact++;
              }
              break;
            case 3:  // amines, iminium ions
              if (nonhnbrs ==
                  2) {  // sec amines, no rotation, same as fixed bonds
                while (nbrIdx != endNbrs) {  // loop over neighbors
                  if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() ==
                      1) {  // hydrogen neighbor necessary for H bond donation
                    hbonddir =
                        conf.getAtomPos(*nbrIdx) - pos;  // hbond direction
                  }
                  ++nbrIdx;
                }
                addVectElements(N, &cos_2, pos, hbonddir);
                interact++;
              } else if (nonhnbrs == 1) {    // primary amines, rotation
                while (nbrIdx != endNbrs) {  // loop over neighbors
                  if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() !=
                      1) {  // search for X-N bond
                    dir = pos - conf.getAtomPos(*nbrIdx);  // X-N bond direction
                    dir.normalize();
                  }
                  ++nbrIdx;
                }
                addVectElements(
                    N, &cos_2_rot, pos,
                    dir * (-bondlength[N]));  // vector dir has typcial N-H
                                              // bondlength, for approximate /
                                              // average calculation of angle p
                                              // (see operator())
                interact++;
              } else {  // ammonia
                addVectElements(N, &no_dep, pos,
                                RDGeom::Point3D(0.0, 0.0, 0.0));
                interact++;
              }
              break;
            case 4:
              if (nonhnbrs == 0) {  // ammonium ion
                addVectElements(N, &no_dep, pos,
                                RDGeom::Point3D(0.0, 0.0, 0.0));
                interact++;
              } else if (nonhnbrs == 1) {    // primary ammonium, rotation
                while (nbrIdx != endNbrs) {  // loop over neighbors
                  if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() !=
                      1) {  // search for X-N bond
                    dir = pos - conf.getAtomPos(*nbrIdx);  // X-N bond direction
                  }
                  ++nbrIdx;
                }
                addVectElements(
                    N, &cos_2_rot, pos,
                    dir * (-bondlength[N]));  // vector dir has typcial N-H
                                              // bondlength, for approximate /
                                              // average calculation of angle p
                                              // (see operator())
                interact++;
              } else {  // secondary or tertiary ammonium, no rotation, same as
                        // fixed bonds
                while (nbrIdx != endNbrs) {  // loop over neighbors
                  if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() ==
                      1) {  // hydrogen neighbor necessary for H bond donation
                    dir = pos -
                          conf.getAtomPos(*nbrIdx);  // hydrogen bond direction
                    addVectElements(N, &cos_2_rot, pos, dir);
                    interact++;
                  }
                  ++nbrIdx;
                }
              }
              break;
            default:  // more than four bonds: not possible with nitrogen
              BOOST_LOG(rdErrorLog)
                  << "HBond: Nitrogen atom bound to more than 4 neighbor atoms: Atom: "
                  << i << std::endl;
          }
        }
        break;
      case 8:  // Oxygen
        boost::tie(nbrIdx, endNbrs) =
            mol.getAtomNeighbors(atom);  // get neighbors
        pos = conf.getAtomPos(i);        // get position

        nbrs = 0;
        nonhnbrs = 0;
        aromaticnbr = false;
        while (nbrIdx != endNbrs) {  // loop over neighbors
          if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() != 1) {
            ++nonhnbrs;  // count non-hydrogen neighbors
          }
          if (mol.getAtomWithIdx(*nbrIdx)
                  ->getIsAromatic()) {  // check whether oxygen is bound to
                                        // aromatic system
            aromaticnbr = true;
          }
          ++nbrs;  // count neighbors
          ++nbrIdx;
        }
        nbrIdx -= nbrs;

        if (nonhnbrs != nbrs) {  // otherwise no hydrogen, no hydrogen bond
                                 // donation possible
          switch (nbrs) {        // no of neighbors
            case 1:              // hydroxyl
              addVectElements(O, &no_dep, pos, RDGeom::Point3D(0.0, 0.0, 0.0));
              interact++;
              break;
            case 2:
              if (nonhnbrs == 0) {  // water
                addVectElements(O, &no_dep, pos,
                                RDGeom::Point3D(0.0, 0.0, 0.0));
                interact++;
              } else {                       // OH groups
                while (nbrIdx != endNbrs) {  // loop over neighbors
                  if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() !=
                      1) {  // search for X-O bond
                    dir = pos - conf.getAtomPos(*nbrIdx);  // O-X bond direction
                    dir.normalize();
                  } else {
                    hbonddir =
                        conf.getAtomPos(*nbrIdx) - pos;  // O-H bond direction
                  }
                  ++nbrIdx;
                }
                if (aromaticnbr) {  // phenolic oxygen, allows h to swap places
                                    // in aromatic plane
                  addVectElements(O, &cos_4, pos, hbonddir);  // first hbond
                  interact++;

                  // let's swap hydrogen and lp:
                  plane = dir.crossProduct(hbonddir);  // lp plane vector
                  plane =
                      plane.crossProduct(dir);  // plane through other bond
                                                // perpendicular to X-N-H plane
                  dir = plane *
                        hbonddir.dotProduct(
                            plane);  // projection of vector dir on plane vector
                  hbonddir -= dir * 2;  // mirroring of dir vector on plane
                  addVectElements(
                      O, &cos_4, pos,
                      hbonddir);  // second hbond, hydrogen at other place
                  interact++;
                } else {  // all other oxygens, flexible
                  addVectElements(
                      O, &cos_4_rot, pos,
                      dir * (-bondlength[O]));  // vector dir has typcial O-H
                                                // bondlength, for approximate /
                                                // average calculation of angle
                                                // p (see operator())
                  interact++;
                }
              }
              break;
            case 3:                 // positively charged oxygen
              if (nonhnbrs == 0) {  // oxonium
                addVectElements(O, &no_dep, pos,
                                RDGeom::Point3D(0.0, 0.0, 0.0));
                interact++;
              } else if (nonhnbrs == 1) {    // R[O+]H2
                while (nbrIdx != endNbrs) {  // loop over neighbors
                  if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() !=
                      1) {  // search for X-O bond
                    dir = pos - conf.getAtomPos(*nbrIdx);  // X-O bond direction
                    dir.normalize();
                    addVectElements(
                        O, &cos_4_rot, pos,
                        dir * (-bondlength[O]));  // vector dir has typcial O-H
                                                  // bondlength, for approximate
                                                  // / average calculation of
                                                  // angle p (see operator())
                    interact++;
                  }
                  ++nbrIdx;
                }
              } else {                       // R1R2[O+]H, non-flexible
                while (nbrIdx != endNbrs) {  // loop over neighbors
                  if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() ==
                      1) {  // search for O-H bond
                    dir = conf.getAtomPos(*nbrIdx) - pos;  // O-H bond direction
                    addVectElements(O, &cos_4, pos, dir);
                    interact++;
                  }
                  ++nbrIdx;
                }
              }
              break;
            case 0:  // oxygen atom gas, no hydrogen bonding
              break;
            default:
              BOOST_LOG(rdErrorLog)
                  << "HBond: Oxygen atom bound to more than 3 neighbor atoms: Atom: "
                  << i << std::endl;
          }
        }
        break;
      default:
        break;
    }
  }
  return interact;
}

double HBond::operator()(double x, double y, double z, double thres) const {
  using namespace HBondDetail;

  if (d_nInteract < 1) {  // no interactions
    return 0.0;           // return 0.0
  }

  double res = 0.0;

  if (d_DAprop == 'A') {
    unsigned int minId = 0;
    double minEne = std::numeric_limits<double>::max();  // minimal energy
    double probeDirection[3];                            // direction of probe

    for (unsigned int i = 0, j = 0; i < d_nInteract;
         i++, j += 3) {  // calculation of energy contributions and searching
                         // for the favored probe direction ( direction of
                         // lowest energy contribution )
      d_vectTargetProbe[j] = x - d_pos[j];  // vector of interaction
      d_vectTargetProbe[j + 1] = y - d_pos[j + 1];
      d_vectTargetProbe[j + 2] = z - d_pos[j + 2];

      double dist2 =
          d_vectTargetProbe[j] * d_vectTargetProbe[j] +
          d_vectTargetProbe[j + 1] * d_vectTargetProbe[j + 1] +
          d_vectTargetProbe[j + 2] *
              d_vectTargetProbe[j +
                                2];  // calc of squared length of interaction
                                     //        std::cout << dist2 << std::endl;
      if (dist2 < thres) {
        double dis = sqrt(dist2);
        double distN[3] = {d_vectTargetProbe[j] / dis,
                           d_vectTargetProbe[j + 1] / dis,
                           d_vectTargetProbe[j + 2] / dis};

        dist2 = std::max(dist2, d_cutoff);

        double t = angle(
            distN[0], distN[1], distN[2], d_direction[j], d_direction[j + 1],
            d_direction[j + 2]);  // calc of angle between direction of hbond
                                  // and target-probe direction

        double eneTerm1 = rm[d_probetype][d_targettypes[i]];
        eneTerm1 *= eneTerm1;                   // squared rm
        eneTerm1 /= dist2;                      // division by squared distance
        double eneTerm2 = eneTerm1 * eneTerm1;  // fourth power of rm/distance
        eneTerm1 = eneTerm2 * eneTerm1;         // sixth power of rm/distance
        eneTerm2 *= eneTerm2;                   // eigth power of rm/distance
        eneTerm2 *= (4.0 * eneTerm1 - 3.0 * eneTerm2);
        eneTerm2 *=
            em[d_probetype][d_targettypes[i]];  // multiplication with em
        d_eneContrib[i] = eneTerm2;

        double t0, ti;
        if (d_function[i] == &cos_acc ||
            d_function[i] == &cos_2_0) {  // only if dependent of ti and t0
          double dotProd = d_vectTargetProbe[j] * d_plane[j] +
                           d_vectTargetProbe[j + 1] * d_plane[j + 1] +
                           d_vectTargetProbe[j + 2] * d_plane[j + 2];
          double vectTargetProbeInPlane[3] = {
              d_vectTargetProbe[j] -
                  d_plane[j] *
                      dotProd,  // projection of targetProbe vector on lp plane
              d_vectTargetProbe[j + 1] -
                  d_plane[j + 1] *
                      dotProd,  // projection of targetProbe vector on lp plane
              d_vectTargetProbe[j + 2] -
                  d_plane[j + 2] *
                      dotProd};  // projection of targetProbe vector on lp plane

          normalize(vectTargetProbeInPlane[0], vectTargetProbeInPlane[1],
                    vectTargetProbeInPlane[2]);

          // angles t_0 (out-of-lonepair-plane angle), t_i (in-lonepair-plane
          // angle),
          t0 = angle(distN[0], distN[1], distN[2], vectTargetProbeInPlane[0],
                     vectTargetProbeInPlane[1],
                     vectTargetProbeInPlane[2]);  // out of plane angle
          ti = angle(d_direction[j], d_direction[j + 1], d_direction[j + 2],
                     vectTargetProbeInPlane[0], vectTargetProbeInPlane[1],
                     vectTargetProbeInPlane[2]);  // in plane angle
        } else {
          t0 = 0.0;
          ti = 1.0;
        }
        d_eneContrib[i] *=
            (*(d_function[i]))(t, t0, ti);  // scaling of energy contribution

        if (d_eneContrib[i] <
            minEne) {  // check whether most favored interaction
          minEne = d_eneContrib[i];
          minId = i;
        }
      } else {
        d_eneContrib[i] = 0.0;
      }
    }
    minId *= 3;
    probeDirection[0] = -d_vectTargetProbe[minId];  // probe is directed to most
                                                    // favored interaction
    probeDirection[1] = -d_vectTargetProbe[++minId];
    probeDirection[2] = -d_vectTargetProbe[++minId];
    normalize(probeDirection[0], probeDirection[1], probeDirection[2]);

    for (unsigned int
             i = 0,
             j = 0;
         i < d_nInteract;
         i++,
             j +=
             3) {  // scaling to take probe direction into account and adding up
      double vectHydrogenTarget[3] = {
          probeDirection[0] * bondlength[d_probetype] - d_vectTargetProbe[j],
          probeDirection[1] * bondlength[d_probetype] -
              d_vectTargetProbe[j + 1],
          probeDirection[2] * bondlength[d_probetype] -
              d_vectTargetProbe[j + 2]};

      normalize(vectHydrogenTarget[0], vectHydrogenTarget[1],
                vectHydrogenTarget[2]);

      // p (angle between best hydrogen bond direction of probe and
      // hydrogen-acceptor vector)
      double p = angle(probeDirection[0], probeDirection[1], probeDirection[2],
                       vectHydrogenTarget[0], vectHydrogenTarget[1],
                       vectHydrogenTarget[2]);

      res +=
          d_eneContrib[i] *
          cos_2(p, 0.0, 1.0);  // scaling of energy contribution and adding up
    }
  } else {  // d_DAprop='D'
    unsigned int minId = 0;
    double minEne = std::numeric_limits<double>::max();  // minimal energy
    double probeDirection[3];                            // direction of prope

    for (unsigned int i = 0, j = 0; i < d_nInteract;
         i++, j += 3) {  // calculation of energy contributions and searching
                         // for the favored probe direction ( direction of
                         // lowest energy contribution )
      d_vectTargetProbe[j] = x - d_pos[j];  // vector of interaction
      d_vectTargetProbe[j + 1] = y - d_pos[j + 1];
      d_vectTargetProbe[j + 2] = z - d_pos[j + 2];

      double dist2 =
          d_vectTargetProbe[j] * d_vectTargetProbe[j] +
          d_vectTargetProbe[j + 1] * d_vectTargetProbe[j + 1] +
          d_vectTargetProbe[j + 2] *
              d_vectTargetProbe[j +
                                2];  // calc of squared length of interaction
                                     //    std::cout << dist2 << std::endl;
      if (dist2 < thres) {
        double dis = sqrt(dist2);
        double distN[3] = {d_vectTargetProbe[j] / dis,
                           d_vectTargetProbe[j + 1] / dis,
                           d_vectTargetProbe[j + 2] / dis};

        dist2 = std::max(dist2, d_cutoff);

        double t = angle(
            distN[0], distN[1], distN[2], d_direction[j], d_direction[j + 1],
            d_direction[j + 2]);  // calc of angle between direction of hbond
                                  // and target-probe direction

        double eneTerm1 = rm[d_probetype][d_targettypes[i]];
        eneTerm1 *= eneTerm1;                   // squared rm
        eneTerm1 /= dist2;                      // division by squared distance
        double eneTerm2 = eneTerm1 * eneTerm1;  // fourth power of rm/distance
        eneTerm1 = eneTerm2 * eneTerm1;         // sixth power of rm/distance
        eneTerm2 *= eneTerm2;                   // eigth power of rm/distance
        eneTerm2 *= (4.0 * eneTerm1 - 3.0 * eneTerm2);
        eneTerm2 *=
            em[d_probetype][d_targettypes[i]];  // multiplication with em
        d_eneContrib[i] = eneTerm2;

        d_eneContrib[i] *=
            (*(d_function[i]))(t, 0.0, 1.0);  // scaling of energy contribution

        if (d_eneContrib[i] <
            minEne) {  // check whether most favored interaction
          minEne = d_eneContrib[i];
          minId = i;
        }
      } else {
        d_eneContrib[i] = 0;
      }
    }
    minId *= 3;
    probeDirection[0] = -d_vectTargetProbe[minId];  // probe is directed to most
                                                    // favored interaction
    probeDirection[1] = -d_vectTargetProbe[++minId];
    probeDirection[2] = -d_vectTargetProbe[++minId];
    normalize(probeDirection[0], probeDirection[1], probeDirection[2]);

    for (unsigned int
             i = 0,
             j = 0;
         i < d_nInteract;
         i++,
             j +=
             3) {  // scaling to take probe direction into account and adding up
      double vectProbeHydrogen[3] = {
          d_direction[j] * d_lengths[i] - d_vectTargetProbe[j],
          d_direction[j + 1] * d_lengths[i] - d_vectTargetProbe[j + 1],
          d_direction[j + 2] * d_lengths[i] - d_vectTargetProbe[j + 2]};

      normalize(vectProbeHydrogen[0], vectProbeHydrogen[1],
                vectProbeHydrogen[2]);

      // p (angle between best hydrogen bond direction of probe and
      // hydrogen-acceptor vector)
      double p = angle(vectProbeHydrogen[0], vectProbeHydrogen[1],
                       vectProbeHydrogen[2], probeDirection[0],
                       probeDirection[1], probeDirection[2]);

      res +=
          d_eneContrib[i] *
          cos_2(p, 0.0, 1.0);  // scaling of energy contribution and adding up
    }
  }
  return res;
}

void HBond::addVectElements(atomtype type,
                            double (*funct)(double, double, double),
                            const RDGeom::Point3D &pos,
                            const RDGeom::Point3D &dir,
                            const RDGeom::Point3D &plane) {
  d_targettypes.push_back(type);
  d_function.push_back(funct);
  d_pos.push_back(pos.x);
  d_pos.push_back(pos.y);
  d_pos.push_back(pos.z);
  double len = dir.length();
  d_lengths.push_back(len);
  d_direction.push_back(dir.x / len);
  d_direction.push_back(dir.y / len);
  d_direction.push_back(dir.z / len);
  len = plane.length();
  if (len > 1.e-16) {
    d_plane.push_back(plane.x / len);
    d_plane.push_back(plane.y / len);
    d_plane.push_back(plane.z / len);
  } else {
    d_plane.push_back(0.0);
    d_plane.push_back(0.0);
    d_plane.push_back(0.0);
  }
}

void HBond::normalize(double &x, double &y, double &z) const {
  double temp = x * x + y * y + z * z;
  temp = sqrt(temp);
  x /= temp;
  y /= temp;
  z /= temp;
}

double HBond::angle(double x1, double y1, double z1, double x2, double y2,
                    double z2) const {
  double dotProd = x1 * x2 + y1 * y2 + z1 * z2;
  if (dotProd < -1.0) {
    dotProd = -1.0;
  } else if (dotProd > 1.0) {
    dotProd = 1.0;
  }
  return acos(dotProd);
}

Hydrophilic::Hydrophilic(const RDKit::ROMol &mol, int confId, bool fixed,
                         double cutoff) {
  d_hbondOH = HBond(mol, confId, "OH", fixed, cutoff);
  d_hbondO = HBond(mol, confId, "O", fixed, cutoff);
}

double Hydrophilic::operator()(double x, double y, double z,
                               double thres) const {
  double hbondO, hbondOH;
  hbondO = d_hbondO(x, y, z, thres);
  hbondOH = d_hbondOH(x, y, z, thres);
  return std::min(hbondO, hbondOH);
}

void writeToCubeStream(const RDGeom::UniformRealValueGrid3D &grd,
                       std::ostream &outStrm, const RDKit::ROMol *mol,
                       int confid) {
  PRECONDITION(outStrm, "bad stream");
  const double bohr = 0.529177249;
  int dimX = grd.getNumX();  //+2;
  int dimY = grd.getNumY();  //+2;
  int dimZ = grd.getNumZ();  //+2;
  auto spacing = grd.getSpacing() / bohr;
  auto offSet = grd.getOffset() / bohr;
  auto nAtoms = mol ? mol->getNumAtoms() : 0u;
  outStrm.setf(std::ios::left);
  outStrm
      << "Gaussian cube format generated by RDKit\n*************************\n";
  outStrm << std::setw(20) << std::setprecision(6) << nAtoms << std::setw(20)
          << std::setprecision(6) << offSet.x << std::setw(20)
          << std::setprecision(6) << offSet.y << std::setw(20)
          << std::setprecision(6) << offSet.z << std::endl;
  outStrm << std::setw(20) << std::setprecision(6) << dimX << std::setw(20)
          << std::setprecision(6) << spacing << std::setw(20)
          << std::setprecision(6) << 0 << std::setw(20) << std::setprecision(6)
          << 0 << std::endl
          << std::setw(20) << std::setprecision(6) << dimY << std::setw(20)
          << std::setprecision(6) << 0 << std::setw(20) << std::setprecision(6)
          << spacing << std::setw(20) << std::setprecision(6) << 0 << std::endl
          << std::setw(20) << std::setprecision(6) << dimZ << std::setw(20)
          << std::setprecision(6) << std::setw(20) << std::setprecision(6) << 0
          << std::setw(20) << std::setprecision(6) << 0 << std::setw(20)
          << std::setprecision(6) << spacing << std::endl;

  if (mol) {
    // this will throw a ConformerException if confId does not exist
    const auto &conf = mol->getConformer(confid);

    for (const auto &atom : mol->atoms()) {
      const auto &pt = conf.getAtomPos(atom->getIdx()) / bohr;
      outStrm << std::setw(20) << std::setprecision(6) << std::left
              << atom->getAtomicNum() << std::setw(20) << std::setprecision(6)
              << atom->getAtomicNum() << std::setw(20) << std::setprecision(6)
              << pt.x << std::setw(20) << std::setprecision(6) << std::setw(20)
              << std::setprecision(6) << pt.y << std::setw(20)
              << std::setprecision(6) << pt.z << std::endl;
    }
  }

  for (auto xi = 0u; xi < grd.getNumX(); ++xi) {
    for (auto yi = 0u; yi < grd.getNumY(); ++yi) {
      for (auto zi = 0u; zi < grd.getNumZ(); ++zi) {
        outStrm << std::setw(20) << std::setprecision(6) << std::left
                << static_cast<double>(
                       grd.getVal(grd.getGridIndex(xi, yi, zi)));
        // grd->d_numX-xi-1, grd->d_numY-yi-1, grd->d_numZ-zi-1
        if ((zi + 1) % 8 == 0) {
          outStrm << std::endl;
        }
      }
      outStrm << std::endl;
    }
    outStrm << std::endl;
  }
}

void writeToCubeFile(const RDGeom::UniformRealValueGrid3D &grd,
                     const std::string &filename, const RDKit::ROMol *mol,
                     int confid) {
  std::ofstream ofStrm(filename.c_str());
  writeToCubeStream(grd, ofStrm, mol, confid);
}

std::unique_ptr<RDKit::RWMol> readFromCubeStream(
    RDGeom::UniformRealValueGrid3D &grd, std::istream &inStrm) {
  PRECONDITION(inStrm, "bad stream");
  constexpr double bohr = 0.529177249;
  constexpr double spacingThreshold = 0.0001;
  if (inStrm.eof()) {
    return nullptr;
  }
  std::string string;
  int nAtoms;
  std::getline(inStrm, string);
  std::getline(inStrm, string);
  inStrm >> nAtoms;
  double x, y, z;
  inStrm >> x >> y >> z;
  const RDGeom::Point3D offSet(x * bohr, y * bohr, z * bohr);

  int dimX, dimY, dimZ;
  double spacingX, spacingY, spacingZ, temp1, temp2;
  inStrm >> dimX >> spacingX >> temp1 >> temp2;
  inStrm >> dimY >> temp1 >> spacingY >> temp2;
  inStrm >> dimZ >> temp1 >> temp2 >> spacingZ;

  if ((fabs(spacingX - spacingY) > spacingThreshold) ||
      (fabs(spacingX - spacingZ) > spacingThreshold)) {
    std::ostringstream errout;
    errout << "Same spacing in all directions needed";
    throw RDKit::FileParseException(errout.str());
  } else {
    spacingX *= bohr;
    grd = RDGeom::UniformRealValueGrid3D(spacingX * static_cast<double>(dimX),
                                         spacingX * static_cast<double>(dimY),
                                         spacingX * static_cast<double>(dimZ),
                                         spacingX, &offSet);
  }
  std::unique_ptr<RDKit::RWMol> molecule;
  if (nAtoms) {
    molecule.reset(new RDKit::RWMol());
    std::unique_ptr<RDKit::Conformer> conf(new RDKit::Conformer(nAtoms));
    int atomNum;
    for (auto i = 0; i < nAtoms; ++i) {
      inStrm >> atomNum >> temp1 >> x >> y >> z;
      RDKit::Atom atom(atomNum);
      molecule->addAtom(&atom, true, false);
      RDGeom::Point3D pos(x * bohr, y * bohr, z * bohr);
      conf->setAtomPos(i, pos);
    }
    molecule->addConformer(conf.release(), false);
  }
  for (auto xi = 0; xi < dimX; ++xi) {
    for (auto yi = 0; yi < dimY; ++yi) {
      for (auto zi = 0; zi < dimZ; ++zi) {
        double tempVal;
        inStrm >> tempVal;
        grd.setVal(grd.getGridIndex(xi, yi, zi), tempVal);
      }
    }
  }
  return molecule;
}

std::unique_ptr<RDKit::RWMol> readFromCubeFile(
    RDGeom::UniformRealValueGrid3D &grd, const std::string &filename) {
  std::ifstream ifStrm(filename.c_str());
  if (ifStrm.bad()) {
    std::ostringstream errout;
    errout << "Bad input file " << filename;
    throw RDKit::BadFileException(errout.str());
  };

  return readFromCubeStream(grd, ifStrm);
}
}  // namespace RDMIF
