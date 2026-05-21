// $Id$
//
//  Copyright (C) 2003-2011 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "GasteigerCharges.h"
#include <RDGeneral/types.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include "GasteigerParams.h"

namespace Gasteiger {
using namespace RDKit;
/*! \brief split the formal charge across atoms of same type if we have a
 * conjugated system
 *
 * This function is called before the charge equivalization iteration is
 * started for the Gasteiger charges. If any of the atom involved in conjugated
 * system have formal charges set on them, this charges is equally distributed
 * across atoms of the same type in that conjugated system. So for example the
 * two nitrogens in the benzamidine system start the iteration with equal
 * charges of 0.5
 */
void splitChargeConjugated(const ROMol &mol, DOUBLE_VECT &charges) {
  INT_VECT marker;
  for (const auto at : mol.atoms()) {
    auto aix = at->getIdx();
    double formal = at->getFormalCharge();
    // std::cout << aix << " formal charges:" << formal << "\n";
    marker.resize(0);
    if ((fabs(formal) > EPS_DOUBLE) && (fabs(charges[aix]) < EPS_DOUBLE)) {
      marker.push_back(aix);
      for (const auto bnd1 : mol.atomBonds(at)) {
        if (bnd1->getIsConjugated()) {
          auto aax = bnd1->getOtherAtomIdx(aix);
          auto aat = mol.getAtomWithIdx(aax);
          for (const auto bnd2 : mol.atomBonds(aat)) {
            if (bnd1 != bnd2) {
              if (bnd2->getIsConjugated()) {
                auto yax = bnd2->getOtherAtomIdx(aax);
                auto yat = mol.getAtomWithIdx(yax);
                if (at->getAtomicNum() == yat->getAtomicNum()) {
                  formal += yat->getFormalCharge();
                  marker.push_back(yax);
                }
              }
            }
          }
        }
      }

      for (const auto mci : marker) {
        charges[mci] = (formal / marker.size());
      }
    }
  }
}
}  // end of namespace Gasteiger

namespace RDKit {
void computeGasteigerCharges(const ROMol *mol, int nIter,
                             bool throwOnParamFailure) {
  PRECONDITION(mol, "bad molecule");
  computeGasteigerCharges(*mol, nIter, throwOnParamFailure);
}
void computeGasteigerCharges(const ROMol &mol, int nIter,
                             bool throwOnParamFailure) {
  std::vector<double> chgs(mol.getNumAtoms());
  computeGasteigerCharges(mol, chgs, nIter, throwOnParamFailure);
}

/*! \brief compute the Gasteiger partial charges and return a new molecule with
 *the charges set
 *
 * Ref : J.Gasteiger, M. Marsili, "Iterative Equalization of Orbital
 * Electronegativity
 *  A Rapid Access to Atomic Charges", Tetrahedron Vol 36 p3219 1980
 */
void computeGasteigerCharges(const ROMol &mol, std::vector<double> &charges,
                             int nIter, bool throwOnParamFailure) {
  PRECONDITION(charges.size() >= mol.getNumAtoms(), "bad array size");

  const PeriodicTable *table = PeriodicTable::getTable();
  const GasteigerParams *params = GasteigerParams::getParams();

  const auto natms = mol.getNumAtoms();
  // space for parameters for each atom in the molecule
  std::vector<DOUBLE_VECT> atmPs;
  atmPs.reserve(natms);

  std::fill(charges.begin(), charges.end(), 0.0);

  DOUBLE_VECT hChrg(
      natms,
      0.0);  // total charge for the implicit hydrogens on each heavy atom

  DOUBLE_VECT ionX(natms, 0.0);
  DOUBLE_VECT energ(natms, 0.0);

  // deal with the conjugated system - distribute the formal charges on atoms of
  // same type in each conjugated system
  Gasteiger::splitChargeConjugated(mol, charges);

  for (const auto atom : mol.atoms()) {
    const std::string elem = table->getElementSymbol(atom->getAtomicNum());
    std::string mode;

    switch (atom->getHybridization()) {
      case Atom::SP3:
        mode = "sp3";
        break;
      case Atom::SP2:
        mode = "sp2";
        break;
      case Atom::SP:
        mode = "sp";
        break;
      default:
        if (atom->getAtomicNum() == 1) {
          // if it is hydrogen
          mode = "*";
        } else if (atom->getAtomicNum() == 16) {
          // we have a sulfur atom with no hybridization information
          // check how many oxygens we have on the sulfur
          int no = 0;
          for (const auto nbrAt : mol.atomNeighbors(atom)) {
            if (nbrAt->getAtomicNum() == 8) {
              no++;
            }
          }
          if (no == 2) {
            mode = "so2";
          } else if (no == 1) {
            mode = "so";
          } else {
            // some other sulfur state. Default to sp3
            mode = "sp3";
          }
        }
    }

    // if we get a unknown mode or element type the following will will throw an
    // exception
    atmPs.push_back(params->getParams(elem, mode, throwOnParamFailure));

    // set ionX parameters
    // if Hydrogen treat differently
    const auto idx = atom->getIdx();
    if (atom->getAtomicNum() == 1) {
      ionX[idx] = IONXH;
    } else {
      ionX[idx] = atmPs[idx][0] + atmPs[idx][1] + atmPs[idx][2];
    }
  }

  // parameters for hydrogen atoms for case where the hydrogen are not in the
  // graph (implicit hydrogens)
  DOUBLE_VECT hParams = params->getParams("H", "*", throwOnParamFailure);

  double damp = DAMP;
  for (auto itx = 0u; itx < static_cast<unsigned int>(nIter); itx++) {
    for (auto aix = 0u; aix < natms; aix++) {
      auto enr = atmPs[aix][0] +
                 charges[aix] * (atmPs[aix][1] + atmPs[aix][2] * charges[aix]);
      energ[aix] = enr;
    }

    for (const auto atom : mol.atoms()) {
      const auto aix = atom->getIdx();
      double dq = 0.0;
      for (const auto nbrAt : mol.atomNeighbors(atom)) {
        const auto nbr = nbrAt->getIdx();
        const auto dx = energ[nbr] - energ[aix];
        int sgn;
        if (dx < 0.0) {
          sgn = 0;
        } else {
          sgn = 1;
        }
        dq += dx / ((sgn * (ionX[aix] - ionX[nbr])) + ionX[nbr]);
      }
      // now loop over the implicit hydrogens and get their contributions
      // since hydrogens don't connect to anything else, update their charges at
      // the same time
      const auto niHs = atom->getTotalNumHs();
      if (niHs > 0) {
        const auto qHs = hChrg[aix] / niHs;
        const auto enr = hParams[0] + qHs * (hParams[1] + hParams[2] * qHs);
        const auto dx = enr - energ[aix];
        int sgn;
        if (dx < 0.0) {
          sgn = 0;
        } else {
          sgn = 1;
        }

        const auto dqH = dx / ((sgn * (ionX[aix] - IONXH)) + IONXH);
        dq += (niHs * dqH);

        // adjust the charges on the hydrogens simultaneously (possible because
        // each of the hydrogens have no other neighbors)
        hChrg[aix] -= (niHs * dqH * damp);
      }
      charges[aix] += (damp * dq);
    }

    damp *= DAMP_SCALE;
  }

  for (auto atom : mol.atoms()) {
    const auto aix = atom->getIdx();
    // set the charge on the heavy atom
    atom->setProp(common_properties::_GasteigerCharge, charges[aix], true);
    // set the implicit hydrogen charges
    atom->setProp(common_properties::_GasteigerHCharge, hChrg[aix], true);
  }
}
}  // namespace RDKit
