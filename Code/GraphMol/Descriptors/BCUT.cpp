//
//  Copyright (C) 2020 Brian P. Kelley
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifdef RDK_HAS_EIGEN3
#include "BCUT.h"
#include "Crippen.h"
#include <Eigen/Dense>
#include <GraphMol/RDKitBase.h>
#include "GraphMol/PartialCharges/GasteigerCharges.h"
#include "GraphMol/PartialCharges/GasteigerParams.h"
#include <RDGeneral/types.h>

namespace RDKit {
namespace Descriptors {
constexpr int NUM_BCUTS = 8;

// diagonal elements are a property (atomic num, charge, etc)
// off diagonal are 1/sqrt(bond_order)
//  Original burden matrix was .1, .2, .3, .15 for single,double,triple or
//  aromatic all other elements are .001
namespace {
std::unique_ptr<Eigen::MatrixXd> make_burden(const ROMol &m) {
  auto num_atoms = m.getNumAtoms();
  std::unique_ptr<Eigen::MatrixXd> burden(
      new Eigen::MatrixXd(num_atoms, num_atoms));

  for (unsigned int i = 0; i < num_atoms; ++i) {
    for (unsigned int j = 0; j < num_atoms; ++j) {
      (*burden)(i, j) = (*burden)(j, i) = 0.001;
    }
  }

  for (auto &bond : m.bonds()) {
    unsigned int i = bond->getBeginAtomIdx();
    unsigned int j = bond->getEndAtomIdx();
    double score = 0.0;
    switch (bond->getBondType()) {
      case Bond::AROMATIC:
        // score = 0.15; orig burden
        score = 0.8164965809277261;  // 1/sqrt(1.5)
        break;
      case Bond::SINGLE:
        // score = 0.1;
        score = 1.0;  // 1/sqrt(1.0)
        break;
      case Bond::DOUBLE:
        // score = 0.2;
        score = 0.7071067811865475;  // 1/sqrt(2.0)
        break;
      case Bond::TRIPLE:
        // score = 0.3;
        score = 0.5773502691896258;  // 1/sqrt(3);
        break;
      default:
        CHECK_INVARIANT(
            0, "Bond order must be Single, Double, Triple or Aromatic");
    }
    (*burden)(i, j) = (*burden)(j, i) = score;
  }
  return burden;
}

std::pair<double, double> BCUT2D(std::unique_ptr<Eigen::MatrixXd> &burden,
                                 const std::vector<double> &atom_props) {
  if (atom_props.empty()) {
    return std::pair<double, double>(0.0, 0.0);
  }

  for (unsigned int i = 0; i < atom_props.size(); ++i) {
    (*burden)(i, i) = atom_props[i];
  }
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(*burden);
  auto eivals = es.eigenvalues();
  double lowest = eivals(0);
  double highest = eivals(atom_props.size() - 1);
  return std::pair<double, double>(highest, lowest);
}
}  // namespace

std::pair<double, double> BCUT2D(const ROMol &m,
                                 const std::vector<double> &atom_props) {
  unsigned int num_atoms = m.getNumAtoms();
  PRECONDITION(atom_props.size() == num_atoms,
               "Number of atom props not equal to number of atoms");

  if (num_atoms == 0) {
    return std::pair<double, double>(0, 0);
  }
  auto burden = make_burden(m);
  return BCUT2D(burden, atom_props);
}

std::pair<double, double> BCUT2D(const ROMol &m,
                                 const std::string &atom_double_prop) {
  std::vector<double> props;
  props.reserve(m.getNumAtoms());
  for (auto &atom : m.atoms()) {
    props.push_back(atom->getProp<double>(atom_double_prop));
  }
  return BCUT2D(m, props);
}

std::vector<double> BCUT2D(const ROMol &m) {
  if (!m.getNumAtoms()) {
    return std::vector<double>(NUM_BCUTS, 0.0);
  }
  std::unique_ptr<ROMol> mol(MolOps::removeAllHs(m));
  unsigned int numAtoms = mol->getNumAtoms();
  std::vector<double> masses;
  std::vector<double> charges;
  masses.reserve(numAtoms);
  charges.reserve(numAtoms);

  RDKit::computeGasteigerCharges(*mol, 12, true);
  for (auto &atom : mol->atoms()) {
    masses.push_back(atom->getMass());
    charges.push_back(
        atom->getProp<double>(common_properties::_GasteigerCharge));
  }

  std::vector<double> slogp(numAtoms, 0.0);
  std::vector<double> cmr(numAtoms, 0.0);
  getCrippenAtomContribs(*mol, slogp, cmr);
  // polarizability? - need model
  // slogp?  sasa?
  auto burden = make_burden(*mol);
  auto atom_bcut = BCUT2D(burden, masses);
  auto gasteiger = BCUT2D(burden, charges);
  auto logp = BCUT2D(burden, slogp);
  auto mr = BCUT2D(burden, cmr);
  std::vector<double> res = {
      atom_bcut.first, atom_bcut.second, gasteiger.first, gasteiger.second,
      logp.first,      logp.second,      mr.first,        mr.second};
  return res;
}
}  // namespace Descriptors
}  // namespace RDKit

#endif
