//  Copyright (C) 2014-2025 Greg Landrum and other RDKit contributors

//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <iostream>
#include <cmath>

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <ForceField/MMFF/Params.h>
#include <ForceField/MMFF/Contribs.h>
#include "AtomTyper.h"
#include "Builder.h"

namespace RDKit {
namespace MMFF {
using namespace ForceFields::MMFF;

namespace Tools {
// ------------------------------------------------------------------------
//
//
//
// ------------------------------------------------------------------------
void addBonds(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
              ForceFields::ForceField *field) {
  PRECONDITION(field, "bad ForceField");
  PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
  PRECONDITION(mmffMolProperties->isValid(),
               "missing atom types - invalid force-field");

  std::ostream &oStream = mmffMolProperties->getMMFFOStream();
  double totalBondStretchEnergy = 0.0;
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << "\n"
                 "B O N D   S T R E T C H I N G\n\n"
                 "------ATOMS------   ATOM TYPES   FF     BOND     IDEAL       "
                 "                 FORCE\n"
                 "  I        J          I    J   CLASS   LENGTH   LENGTH    "
                 "DIFF.    ENERGY   CONSTANT\n"
                 "-------------------------------------------------------------"
                 "------------------------"
              << std::endl;
    }
  }

  auto contrib = std::make_unique<BondStretchContrib>(field);
  bool hasContrib = false;
  for (ROMol::ConstBondIterator bi = mol.beginBonds(); bi != mol.endBonds();
       ++bi) {
    unsigned int idx1 = (*bi)->getBeginAtomIdx();
    unsigned int idx2 = (*bi)->getEndAtomIdx();
    unsigned int bondType;
    MMFFBond mmffBondParams;
    if (mmffMolProperties->getMMFFBondStretchParams(mol, idx1, idx2, bondType,
                                                    mmffBondParams)) {
      contrib->addTerm(idx1, idx2, &mmffBondParams);
      hasContrib = true;
      if (mmffMolProperties->getMMFFVerbosity()) {
        unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(idx1);
        unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(idx2);
        const Atom *iAtom = (*bi)->getBeginAtom();
        const Atom *jAtom = (*bi)->getEndAtom();
        const double dist = field->distance(idx1, idx2);
        const double bondStretchEnergy = MMFF::Utils::calcBondStretchEnergy(
            mmffBondParams.r0, mmffBondParams.kb, dist);
        if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
          oStream << std::left << std::setw(2) << iAtom->getSymbol() << " #"
                  << std::setw(5) << idx1 + 1 << std::setw(2)
                  << jAtom->getSymbol() << " #" << std::setw(5) << idx2 + 1
                  << std::right << std::setw(5) << iAtomType << std::setw(5)
                  << jAtomType << std::setw(6) << bondType << "  " << std::fixed
                  << std::setprecision(3) << std::setw(9) << dist
                  << std::setw(9) << mmffBondParams.r0 << std::setw(9)
                  << dist - mmffBondParams.r0 << std::setw(10)
                  << bondStretchEnergy << std::setw(10) << mmffBondParams.kb
                  << std::endl;
        }
        totalBondStretchEnergy += bondStretchEnergy;
      }
    }
  }
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << std::endl;
    }
    oStream << "TOTAL BOND STRETCH ENERGY      =" << std::right << std::setw(16)
            << std::fixed << std::setprecision(4) << totalBondStretchEnergy
            << std::endl;
  }
  if (hasContrib) {
    field->contribs().push_back(ForceFields::ContribPtr(contrib.release()));
  }
}

unsigned int twoBitCellPos(unsigned int nAtoms, int i, int j) {
  if (j < i) {
    std::swap(i, j);
  }

  return i * (nAtoms - 1) + i * (1 - i) / 2 + j;
}

void setTwoBitCell(boost::shared_array<std::uint8_t> &res, unsigned int pos,
                   std::uint8_t value) {
  unsigned int twoBitPos = pos / 4;
  unsigned int shift = 2 * (pos % 4);
  std::uint8_t twoBitMask = 3 << shift;
  res[twoBitPos] = ((res[twoBitPos] & (~twoBitMask)) | (value << shift));
}

std::uint8_t getTwoBitCell(boost::shared_array<std::uint8_t> &res,
                           unsigned int pos) {
  unsigned int twoBitPos = pos / 4;
  unsigned int shift = 2 * (pos % 4);
  std::uint8_t twoBitMask = 3 << shift;

  return ((res[twoBitPos] & twoBitMask) >> shift);
}

// ------------------------------------------------------------------------
//
// the two-bit matrix returned by this contains:
//   0: if atoms i and j are directly connected
//   1: if atoms i and j are connected via an atom
//   2: if atoms i and j are in a 1,4 relationship
//   3: otherwise
//
//  NOTE: the caller is responsible for calling delete []
//  on the result
//
// ------------------------------------------------------------------------
boost::shared_array<std::uint8_t> buildNeighborMatrix(const ROMol &mol) {
  const std::uint8_t RELATION_1_X_INIT = RELATION_1_X | (RELATION_1_X << 2) |
                                         (RELATION_1_X << 4) |
                                         (RELATION_1_X << 6);
  unsigned int nAtoms = mol.getNumAtoms();
  unsigned nTwoBitCells = (nAtoms * (nAtoms + 1) - 1) / 8 + 1;
  boost::shared_array<std::uint8_t> res(new std::uint8_t[nTwoBitCells]);
  std::memset(res.get(), RELATION_1_X_INIT, nTwoBitCells);
  for (ROMol::ConstBondIterator bondi = mol.beginBonds();
       bondi != mol.endBonds(); ++bondi) {
    setTwoBitCell(res,
                  twoBitCellPos(nAtoms, (*bondi)->getBeginAtomIdx(),
                                (*bondi)->getEndAtomIdx()),
                  RELATION_1_2);
    unsigned int bondiBeginAtomIdx = (*bondi)->getBeginAtomIdx();
    unsigned int bondiEndAtomIdx = (*bondi)->getEndAtomIdx();
    for (ROMol::ConstBondIterator bondj = bondi; ++bondj != mol.endBonds();) {
      int idx1 = -1;
      int idx3 = -1;
      unsigned int bondjBeginAtomIdx = (*bondj)->getBeginAtomIdx();
      unsigned int bondjEndAtomIdx = (*bondj)->getEndAtomIdx();
      if (bondiBeginAtomIdx == bondjBeginAtomIdx) {
        idx1 = bondiEndAtomIdx;
        idx3 = bondjEndAtomIdx;
      } else if (bondiBeginAtomIdx == bondjEndAtomIdx) {
        idx1 = bondiEndAtomIdx;
        idx3 = bondjBeginAtomIdx;
      } else if (bondiEndAtomIdx == bondjBeginAtomIdx) {
        idx1 = bondiBeginAtomIdx;
        idx3 = bondjEndAtomIdx;
      } else if (bondiEndAtomIdx == bondjEndAtomIdx) {
        idx1 = bondiBeginAtomIdx;
        idx3 = bondjBeginAtomIdx;
      } else {
        // check if atoms i and j are in a 1,4-relationship
        if ((mol.getBondBetweenAtoms(bondiBeginAtomIdx, bondjBeginAtomIdx)) &&
            (getTwoBitCell(res, twoBitCellPos(nAtoms, bondiEndAtomIdx,
                                              bondjEndAtomIdx)) ==
             RELATION_1_X)) {
          setTwoBitCell(res,
                        twoBitCellPos(nAtoms, bondiEndAtomIdx, bondjEndAtomIdx),
                        RELATION_1_4);
        } else if ((mol.getBondBetweenAtoms(bondiBeginAtomIdx,
                                            bondjEndAtomIdx)) &&
                   (getTwoBitCell(res, twoBitCellPos(nAtoms, bondiEndAtomIdx,
                                                     bondjBeginAtomIdx)) ==
                    RELATION_1_X)) {
          setTwoBitCell(
              res, twoBitCellPos(nAtoms, bondiEndAtomIdx, bondjBeginAtomIdx),
              RELATION_1_4);
        } else if ((mol.getBondBetweenAtoms(bondiEndAtomIdx,
                                            bondjBeginAtomIdx)) &&
                   (getTwoBitCell(res, twoBitCellPos(nAtoms, bondiBeginAtomIdx,
                                                     bondjEndAtomIdx)) ==
                    RELATION_1_X)) {
          setTwoBitCell(
              res, twoBitCellPos(nAtoms, bondiBeginAtomIdx, bondjEndAtomIdx),
              RELATION_1_4);
        } else if ((mol.getBondBetweenAtoms(bondiEndAtomIdx,
                                            bondjEndAtomIdx)) &&
                   (getTwoBitCell(res, twoBitCellPos(nAtoms, bondiBeginAtomIdx,
                                                     bondjBeginAtomIdx)) ==
                    RELATION_1_X)) {
          setTwoBitCell(
              res, twoBitCellPos(nAtoms, bondiBeginAtomIdx, bondjBeginAtomIdx),
              RELATION_1_4);
        }
      }
      if (idx1 > -1) {
        setTwoBitCell(res, twoBitCellPos(nAtoms, idx1, idx3), RELATION_1_3);
      }
    }
  }
  return res;
}

// ------------------------------------------------------------------------
//
//
//
// ------------------------------------------------------------------------
void addAngles(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
               ForceFields::ForceField *field) {
  PRECONDITION(field, "bad ForceField");
  PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
  PRECONDITION(mmffMolProperties->isValid(),
               "missing atom types - invalid force-field");

  std::ostream &oStream = mmffMolProperties->getMMFFOStream();
  unsigned int idx[3];
  const MMFFPropCollection *mmffProp = DefaultParameters::getMMFFProp();
  ROMol::ADJ_ITER nbr1Idx;
  ROMol::ADJ_ITER end1Nbrs;
  ROMol::ADJ_ITER nbr2Idx;
  ROMol::ADJ_ITER end2Nbrs;

  unsigned int nAtoms = mol.getNumAtoms();
  double totalAngleBendEnergy = 0.0;
  RDGeom::PointPtrVect points;
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << "\n"
                 "A N G L E   B E N D I N G\n\n"
                 "-----------ATOMS-----------    ATOM TYPES      FF    VALENCE "
                 "   IDEAL                          FORCE\n"
                 "  I        J        K          I    J    K   CLASS    ANGLE  "
                 "   ANGLE      DIFF.    ENERGY   CONSTANT\n"
                 "-------------------------------------------------------------"
                 "-----------------------------------------"
              << std::endl;
    }
    points = field->positions();
  }
  auto contrib = std::make_unique<AngleBendContrib>(field);
  bool hasContrib = false;
  for (idx[1] = 0; idx[1] < nAtoms; ++idx[1]) {
    const Atom *jAtom = mol.getAtomWithIdx(idx[1]);
    if (jAtom->getDegree() == 1) {
      continue;
    }
    unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(idx[1]);
    const MMFFProp *mmffPropParamsCentralAtom = (*mmffProp)(jAtomType);
    boost::tie(nbr1Idx, end1Nbrs) = mol.getAtomNeighbors(jAtom);
    for (; nbr1Idx != end1Nbrs; ++nbr1Idx) {
      const Atom *iAtom = mol[*nbr1Idx];
      idx[0] = iAtom->getIdx();
      boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(jAtom);
      for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
        if (nbr2Idx < (nbr1Idx + 1)) {
          continue;
        }
        const Atom *kAtom = mol[*nbr2Idx];
        idx[2] = kAtom->getIdx();
        unsigned int angleType;
        MMFFAngle mmffAngleParams;
        if (mmffMolProperties->getMMFFAngleBendParams(
                mol, idx[0], idx[1], idx[2], angleType, mmffAngleParams)) {
          hasContrib = true;
          contrib->addTerm(idx[0], idx[1], idx[2],
                           &mmffAngleParams, mmffPropParamsCentralAtom);
          if (mmffMolProperties->getMMFFVerbosity()) {
            unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(idx[0]);
            unsigned int kAtomType = mmffMolProperties->getMMFFAtomType(idx[2]);
            const RDGeom::Point3D p1((*(points[idx[0]]))[0],
                                     (*(points[idx[0]]))[1],
                                     (*(points[idx[0]]))[2]);
            const RDGeom::Point3D p2((*(points[idx[1]]))[0],
                                     (*(points[idx[1]]))[1],
                                     (*(points[idx[1]]))[2]);
            const RDGeom::Point3D p3((*(points[idx[2]]))[0],
                                     (*(points[idx[2]]))[1],
                                     (*(points[idx[2]]))[2]);
            const double cosTheta = MMFF::Utils::calcCosTheta(
                p1, p2, p3, field->distance(idx[0], idx[1]),
                field->distance(idx[1], idx[2]));
            const double theta = RAD2DEG * acos(cosTheta);
            const double angleBendEnergy = MMFF::Utils::calcAngleBendEnergy(
                mmffAngleParams.theta0, mmffAngleParams.ka,
                mmffPropParamsCentralAtom->linh, cosTheta);
            if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
              oStream << std::left << std::setw(2) << iAtom->getSymbol() << " #"
                      << std::setw(5) << idx[0] + 1 << std::setw(2)
                      << jAtom->getSymbol() << " #" << std::setw(5)
                      << idx[1] + 1 << std::setw(2) << kAtom->getSymbol()
                      << " #" << std::setw(5) << idx[2] + 1 << std::right
                      << std::setw(5) << iAtomType << std::setw(5) << jAtomType
                      << std::setw(5) << kAtomType << std::setw(6) << angleType
                      << "  " << std::fixed << std::setprecision(3)
                      << std::setw(10) << theta << std::setw(10)
                      << mmffAngleParams.theta0 << std::setw(10)
                      << theta - mmffAngleParams.theta0 << std::setw(10)
                      << angleBendEnergy << std::setw(10) << mmffAngleParams.ka
                      << std::endl;
            }
            totalAngleBendEnergy += angleBendEnergy;
          }
        }
      }
    }
  }
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << std::endl;
    }
    oStream << "TOTAL ANGLE BEND ENERGY        =" << std::right << std::setw(16)
            << std::fixed << std::setprecision(4) << totalAngleBendEnergy
            << std::endl;
  }
  if (hasContrib) {
    field->contribs().push_back(ForceFields::ContribPtr(contrib.release()));
  }
}

// ------------------------------------------------------------------------
//
//
//
// ------------------------------------------------------------------------
void addStretchBend(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
                    ForceFields::ForceField *field) {
  PRECONDITION(field, "bad ForceField");
  PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
  PRECONDITION(mmffMolProperties->isValid(),
               "missing atom types - invalid force-field");

  std::ostream &oStream = mmffMolProperties->getMMFFOStream();
  unsigned int idx[3];
  const MMFFPropCollection *mmffProp = DefaultParameters::getMMFFProp();
  ROMol::ADJ_ITER nbr1Idx;
  ROMol::ADJ_ITER end1Nbrs;
  ROMol::ADJ_ITER nbr2Idx;
  ROMol::ADJ_ITER end2Nbrs;

  unsigned int nAtoms = mol.getNumAtoms();
  double totalStretchBendEnergy = 0.0;
  RDGeom::PointPtrVect points;
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << "\n"
                 "S T R E T C H   B E N D I N G\n\n"
                 "-----------ATOMS-----------    ATOM TYPES      FF    VALENCE "
                 "    DELTA     DELTA     DELTA               F CON\n"
                 "  I        J        K          I    J    K   CLASS    ANGLE  "
                 "    ANGLE     R(I,J)    R(J,K)   ENERGY    I-J (J-K)\n"
                 "-------------------------------------------------------------"
                 "-----------------------------------------------------"
              << std::endl;
    }
    points = field->positions();
  }

  auto contrib = std::make_unique<StretchBendContrib>(field);
  bool contribAdded = false;

  for (idx[1] = 0; idx[1] < nAtoms; ++idx[1]) {
    const Atom *jAtom = mol.getAtomWithIdx(idx[1]);
    if (jAtom->getDegree() == 1) {
      continue;
    }
    unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(idx[1]);
    const MMFFProp *mmffPropParamsCentralAtom = (*mmffProp)(jAtomType);
    if (mmffPropParamsCentralAtom->linh) {
      continue;
    }
    boost::tie(nbr1Idx, end1Nbrs) = mol.getAtomNeighbors(jAtom);
    unsigned int i = 0;
    for (; nbr1Idx != end1Nbrs; ++nbr1Idx) {
      const Atom *iAtom = mol[*nbr1Idx];
      boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(jAtom);
      unsigned int j = 0;
      for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
        const Atom *kAtom = mol[*nbr2Idx];
        if (j < (i + 1)) {
          ++j;
          continue;
        }
        idx[0] = iAtom->getIdx();
        idx[2] = kAtom->getIdx();
        unsigned int stretchBendType;
        MMFFStbn mmffStbnParams;
        MMFFBond mmffBondParams[2];
        MMFFAngle mmffAngleParams;
        if (mmffMolProperties->getMMFFStretchBendParams(
                mol, idx[0], idx[1], idx[2], stretchBendType, mmffStbnParams,
                mmffBondParams, mmffAngleParams)) {
          contribAdded = true;
          contrib->addTerm(idx[0], idx[1], idx[2], &mmffStbnParams,
                           &mmffAngleParams, &mmffBondParams[0],
                           &mmffBondParams[1]);
          if (mmffMolProperties->getMMFFVerbosity()) {
            unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(idx[0]);
            unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(idx[1]);
            unsigned int kAtomType = mmffMolProperties->getMMFFAtomType(idx[2]);
            const double dist1 = field->distance(idx[0], idx[1]);
            const double dist2 = field->distance(idx[1], idx[2]);
            const RDGeom::Point3D p1((*(points[idx[0]]))[0],
                                     (*(points[idx[0]]))[1],
                                     (*(points[idx[0]]))[2]);
            const RDGeom::Point3D p2((*(points[idx[1]]))[0],
                                     (*(points[idx[1]]))[1],
                                     (*(points[idx[1]]))[2]);
            const RDGeom::Point3D p3((*(points[idx[2]]))[0],
                                     (*(points[idx[2]]))[1],
                                     (*(points[idx[2]]))[2]);
            const double cosTheta =
                MMFF::Utils::calcCosTheta(p1, p2, p3, dist1, dist2);
            const double theta = RAD2DEG * acos(cosTheta);
            const std::pair<double, double> forceConstants =
                MMFF::Utils::calcStbnForceConstants(&mmffStbnParams);
            const std::pair<double, double> stretchBendEnergies =
                MMFF::Utils::calcStretchBendEnergy(
                    dist1 - mmffBondParams[0].r0, dist2 - mmffBondParams[1].r0,
                    theta - mmffAngleParams.theta0, forceConstants);
            if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
              if (!isDoubleZero(forceConstants.first)) {
                oStream << std::left << std::setw(2) << iAtom->getSymbol()
                        << " #" << std::setw(5) << idx[0] + 1 << std::setw(2)
                        << jAtom->getSymbol() << " #" << std::setw(5)
                        << idx[1] + 1 << std::setw(2) << kAtom->getSymbol()
                        << " #" << std::setw(5) << idx[2] + 1 << std::right
                        << std::setw(5) << iAtomType << std::setw(5)
                        << jAtomType << std::setw(5) << kAtomType
                        << std::setw(6) << stretchBendType << "  " << std::fixed
                        << std::setprecision(3) << std::setw(10) << theta
                        << std::setw(10) << theta - mmffAngleParams.theta0
                        << std::setw(10) << dist1 - mmffBondParams[0].r0
                        << std::setw(10) << dist2 - mmffBondParams[1].r0
                        << std::setw(10) << stretchBendEnergies.first
                        << std::setw(10) << mmffStbnParams.kbaIJK << std::endl;
              }
              if (!isDoubleZero(forceConstants.second)) {
                oStream << std::left << std::setw(2) << kAtom->getSymbol()
                        << " #" << std::setw(5) << idx[2] + 1 << std::setw(2)
                        << jAtom->getSymbol() << " #" << std::setw(5)
                        << idx[1] + 1 << std::setw(2) << iAtom->getSymbol()
                        << " #" << std::setw(5) << idx[0] + 1 << std::right
                        << std::setw(5) << kAtomType << std::setw(5)
                        << jAtomType << std::setw(5) << iAtomType
                        << std::setw(6) << stretchBendType << "  " << std::fixed
                        << std::setprecision(3) << std::setw(10) << theta
                        << std::setw(10) << theta - mmffAngleParams.theta0
                        << std::setw(10) << dist1 - mmffBondParams[0].r0
                        << std::setw(10) << dist2 - mmffBondParams[1].r0
                        << std::setw(10) << stretchBendEnergies.second
                        << std::setw(10) << mmffStbnParams.kbaKJI << std::endl;
              }
            }
            totalStretchBendEnergy +=
                (stretchBendEnergies.first + stretchBendEnergies.second);
          }
        }
        ++j;
      }
      ++i;
    }
  }
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << std::endl;
    }
    oStream << "TOTAL STRETCH-BEND ENERGY      =" << std::right << std::setw(16)
            << std::fixed << std::setprecision(4) << totalStretchBendEnergy
            << std::endl;
  }
  if (contribAdded) {
    field->contribs().push_back(ForceFields::ContribPtr(contrib.release()));
  }
}

void addOop(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
            ForceFields::ForceField *field) {
  PRECONDITION(field, "bad ForceField");
  PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
  PRECONDITION(mmffMolProperties->isValid(),
               "missing atom types - invalid force-field");

  std::ostream &oStream = mmffMolProperties->getMMFFOStream();
  unsigned int idx[4];
  unsigned int atomType[4];
  unsigned int n[4];
  const Atom *atom[4];
  ROMol::ADJ_ITER nbrIdx;
  ROMol::ADJ_ITER endNbrs;

  double totalOopBendEnergy = 0.0;
  RDGeom::PointPtrVect points;
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << "\n"
                 "O U T - O F - P L A N E   B E N D I N G\n\n"
                 "--------------ATOMS---------------         ATOM TYPES        "
                 " OOP                FORCE\n"
                 "  I        J        K        L          I    J    K    L     "
                 "ANGLE    ENERGY   CONSTANT\n"
                 "-------------------------------------------------------------"
                 "-----------------------------"
              << std::endl;
    }
    points = field->positions();
  }
  bool hasContrib = false;
  auto contrib = std::make_unique<OopBendContrib>(field);
  for (idx[1] = 0; idx[1] < mol.getNumAtoms(); ++idx[1]) {
    atom[1] = mol.getAtomWithIdx(idx[1]);
    if (atom[1]->getDegree() != 3) {
      continue;
    }
    atomType[1] = mmffMolProperties->getMMFFAtomType(idx[1]);
    boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom[1]);
    unsigned int i = 0;
    for (; nbrIdx != endNbrs; ++nbrIdx) {
      atom[i] = mol[*nbrIdx];
      idx[i] = atom[i]->getIdx();
      atomType[i] = mmffMolProperties->getMMFFAtomType(idx[i]);
      if (!i) {
        ++i;
      }
      ++i;
    }
    MMFFOop mmffOopParams;
    // if no parameters could be found, we exclude this term (SURDOX02)
    if (!(mmffMolProperties->getMMFFOopBendParams(mol, idx[0], idx[1], idx[2],
                                                  idx[3], mmffOopParams))) {
      continue;
    }
    for (unsigned int i = 0; i < 3; ++i) {
      n[1] = 1;
      switch (i) {
        case 0:
          n[0] = 0;
          n[2] = 2;
          n[3] = 3;
          break;

        case 1:
          n[0] = 0;
          n[2] = 3;
          n[3] = 2;
          break;

        case 2:
          n[0] = 2;
          n[2] = 3;
          n[3] = 0;
          break;
      }
      contrib->addTerm(idx[n[0]], idx[n[1]], idx[n[2]], idx[n[3]], &mmffOopParams);
      hasContrib = true;
      if (mmffMolProperties->getMMFFVerbosity()) {
        const RDGeom::Point3D p1((*(points[idx[n[0]]]))[0],
                                 (*(points[idx[n[0]]]))[1],
                                 (*(points[idx[n[0]]]))[2]);
        const RDGeom::Point3D p2((*(points[idx[n[1]]]))[0],
                                 (*(points[idx[n[1]]]))[1],
                                 (*(points[idx[n[1]]]))[2]);
        const RDGeom::Point3D p3((*(points[idx[n[2]]]))[0],
                                 (*(points[idx[n[2]]]))[1],
                                 (*(points[idx[n[2]]]))[2]);
        const RDGeom::Point3D p4((*(points[idx[n[3]]]))[0],
                                 (*(points[idx[n[3]]]))[1],
                                 (*(points[idx[n[3]]]))[2]);
        const double chi = MMFF::Utils::calcOopChi(p1, p2, p3, p4);
        const double oopBendEnergy =
            MMFF::Utils::calcOopBendEnergy(chi, mmffOopParams.koop);
        if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
          oStream << std::left << std::setw(2) << atom[0]->getSymbol() << " #"
                  << std::setw(5) << idx[n[0]] + 1 << std::setw(2)
                  << atom[1]->getSymbol() << " #" << std::setw(5)
                  << idx[n[1]] + 1 << std::setw(2) << atom[2]->getSymbol()
                  << " #" << std::setw(5) << idx[n[2]] + 1 << std::setw(2)
                  << atom[3]->getSymbol() << " #" << std::setw(5)
                  << idx[n[3]] + 1 << std::right << std::setw(5)
                  << atomType[n[0]] << std::setw(5) << atomType[n[1]]
                  << std::setw(5) << atomType[n[2]] << std::setw(5)
                  << atomType[n[3]] << std::fixed << std::setprecision(3)
                  << std::setw(10) << chi << std::setw(10) << oopBendEnergy
                  << std::setw(10) << mmffOopParams.koop << std::endl;
        }
        totalOopBendEnergy += oopBendEnergy;
      }
    }
  }
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << std::endl;
    }
    oStream << "TOTAL OUT-OF-PLANE BEND ENERGY =" << std::right << std::setw(16)
            << std::fixed << std::setprecision(4) << totalOopBendEnergy
            << std::endl;
  }
  if (hasContrib) {
    field->contribs().push_back(ForceFields::ContribPtr(contrib.release()));
  }
}

// ------------------------------------------------------------------------
//
//
//
// ------------------------------------------------------------------------
const std::string DefaultTorsionBondSmarts::ds_string =
    "[!$(*#*)&!D1]~[!$(*#*)&!D1]";
boost::scoped_ptr<const ROMol> DefaultTorsionBondSmarts::ds_instance;

void DefaultTorsionBondSmarts::create() {
  ds_instance.reset(SmartsToMol(ds_string));
}

#ifdef RDK_BUILD_THREADSAFE_SSS
std::once_flag DefaultTorsionBondSmarts::ds_flag;
#endif
const ROMol *DefaultTorsionBondSmarts::query() {
#ifdef RDK_BUILD_THREADSAFE_SSS
  std::call_once(ds_flag, create);
#else
  static bool created = false;
  if (!created) {
    created = true;
    create();
  }
#endif
  return ds_instance.get();
}

void addTorsions(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
                 ForceFields::ForceField *field,
                 const std::string &torsionBondSmarts) {
  PRECONDITION(field, "bad ForceField");
  PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
  PRECONDITION(mmffMolProperties->isValid(),
               "missing atom types - invalid force-field");

  std::ostream &oStream = mmffMolProperties->getMMFFOStream();
  ROMol::ADJ_ITER nbr1Idx;
  ROMol::ADJ_ITER end1Nbrs;
  ROMol::ADJ_ITER nbr2Idx;
  ROMol::ADJ_ITER end2Nbrs;
  double totalTorsionEnergy = 0.0;
  RDGeom::PointPtrVect points;
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << "\n"
                 "T O R S I O N A L\n\n"
                 "--------------ATOMS---------------      ---ATOM TYPES---     "
                 "FF     TORSION              -----FORCE Params-----\n"
                 "  I        J        K        L          I    J    K    L   "
                 "CLASS     ANGLE    ENERGY       V1        V2        V3\n"
                 "-------------------------------------------------------------"
                 "-----------------------------------------------------"
              << std::endl;
    }
    points = field->positions();
  }
  std::vector<MatchVectType> matchVect;
  const ROMol *defaultQuery = DefaultTorsionBondSmarts::query();
  const ROMol *query = (torsionBondSmarts == DefaultTorsionBondSmarts::string())
                           ? defaultQuery
                           : SmartsToMol(torsionBondSmarts);
  TEST_ASSERT(query);
  unsigned int nHits = SubstructMatch(mol, *query, matchVect);
  if (query != defaultQuery) {
    delete query;
  }

  auto contrib = std::make_unique<TorsionAngleContrib>(field);
  bool contribAdded = false;

  for (unsigned int i = 0; i < nHits; ++i) {
    MatchVectType match = matchVect[i];
    TEST_ASSERT(match.size() == 2);
    int idx2 = match[0].second;
    int idx3 = match[1].second;
    const Bond *bond = mol.getBondBetweenAtoms(idx2, idx3);
    TEST_ASSERT(bond);
    const Atom *jAtom = mol.getAtomWithIdx(idx2);
    const Atom *kAtom = mol.getAtomWithIdx(idx3);
    if (((jAtom->getHybridization() == Atom::SP2) ||
         (jAtom->getHybridization() == Atom::SP3)) &&
        ((kAtom->getHybridization() == Atom::SP2) ||
         (kAtom->getHybridization() == Atom::SP3))) {
      ROMol::OEDGE_ITER beg1, end1;
      boost::tie(beg1, end1) = mol.getAtomBonds(jAtom);
      while (beg1 != end1) {
        const Bond *tBond1 = mol[*beg1];
        if (tBond1 != bond) {
          int idx1 = tBond1->getOtherAtomIdx(idx2);
          ROMol::OEDGE_ITER beg2, end2;
          boost::tie(beg2, end2) = mol.getAtomBonds(kAtom);
          while (beg2 != end2) {
            const Bond *tBond2 = mol[*beg2];
            if ((tBond2 != bond) && (tBond2 != tBond1)) {
              int idx4 = tBond2->getOtherAtomIdx(idx3);
              // make sure this isn't a three-membered ring:
              if (idx4 != idx1) {
                // we now have a torsion involving atoms (bonds):
                //  bIdx - (tBond1) - idx1 - (bond) - idx2 - (tBond2) - eIdx
                unsigned int torType;
                MMFFTor mmffTorParams;
                if (mmffMolProperties->getMMFFTorsionParams(
                        mol, idx1, idx2, idx3, idx4, torType, mmffTorParams)) {
                  contrib->addTerm(idx1, idx2, idx3, idx4, &mmffTorParams);
                  contribAdded = true;
                  if (mmffMolProperties->getMMFFVerbosity()) {
                    const Atom *iAtom = mol.getAtomWithIdx(idx1);
                    const Atom *lAtom = mol.getAtomWithIdx(idx4);
                    unsigned int iAtomType =
                        mmffMolProperties->getMMFFAtomType(idx1);
                    unsigned int jAtomType =
                        mmffMolProperties->getMMFFAtomType(idx2);
                    unsigned int kAtomType =
                        mmffMolProperties->getMMFFAtomType(idx3);
                    unsigned int lAtomType =
                        mmffMolProperties->getMMFFAtomType(idx4);
                    const RDGeom::Point3D p1((*(points[idx1]))[0],
                                             (*(points[idx1]))[1],
                                             (*(points[idx1]))[2]);
                    const RDGeom::Point3D p2((*(points[idx2]))[0],
                                             (*(points[idx2]))[1],
                                             (*(points[idx2]))[2]);
                    const RDGeom::Point3D p3((*(points[idx3]))[0],
                                             (*(points[idx3]))[1],
                                             (*(points[idx3]))[2]);
                    const RDGeom::Point3D p4((*(points[idx4]))[0],
                                             (*(points[idx4]))[1],
                                             (*(points[idx4]))[2]);
                    const double cosPhi =
                        MMFF::Utils::calcTorsionCosPhi(p1, p2, p3, p4);
                    const double torsionEnergy = MMFF::Utils::calcTorsionEnergy(
                        mmffTorParams.V1, mmffTorParams.V2, mmffTorParams.V3,
                        cosPhi);
                    if (mmffMolProperties->getMMFFVerbosity() ==
                        MMFF_VERBOSITY_HIGH) {
                      oStream
                          << std::left << std::setw(2) << iAtom->getSymbol()
                          << " #" << std::setw(5) << idx1 + 1 << std::setw(2)
                          << jAtom->getSymbol() << " #" << std::setw(5)
                          << idx2 + 1 << std::setw(2) << kAtom->getSymbol()
                          << " #" << std::setw(5) << idx3 + 1 << std::setw(2)
                          << lAtom->getSymbol() << " #" << std::setw(5)
                          << idx4 + 1 << std::right << std::setw(5) << iAtomType
                          << std::setw(5) << jAtomType << std::setw(5)
                          << kAtomType << std::setw(5) << lAtomType
                          << std::setw(6) << torType << "  " << std::fixed
                          << std::setprecision(3) << std::setw(10)
                          << RAD2DEG * acos(cosPhi) << std::setw(10)
                          << torsionEnergy << std::setw(10) << mmffTorParams.V1
                          << std::setw(10) << mmffTorParams.V2 << std::setw(10)
                          << mmffTorParams.V3 << std::endl;
                    }
                    totalTorsionEnergy += torsionEnergy;
                  }
                }
              }
            }
            beg2++;
          }
        }
        beg1++;
      }
    }
  }
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << std::endl;
    }
    oStream << "TOTAL TORSIONAL ENERGY         =" << std::right << std::setw(16)
            << std::fixed << std::setprecision(4) << totalTorsionEnergy
            << std::endl;
  }
  if (contribAdded) {
    field->contribs().push_back(ForceFields::ContribPtr(contrib.release()));
  }
}

// ------------------------------------------------------------------------
//
//
//
// ------------------------------------------------------------------------
void addVdW(const ROMol &mol, int confId, MMFFMolProperties *mmffMolProperties,
            ForceFields::ForceField *field,
            boost::shared_array<std::uint8_t> neighborMatrix,
            double nonBondedThresh, bool ignoreInterfragInteractions) {
  PRECONDITION(field, "bad ForceField");
  PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
  PRECONDITION(mmffMolProperties->isValid(),
               "missing atom types - invalid force-field");

  std::ostream &oStream = mmffMolProperties->getMMFFOStream();
  INT_VECT fragMapping;
  if (ignoreInterfragInteractions) {
    std::vector<ROMOL_SPTR> molFrags =
        MolOps::getMolFrags(mol, true, &fragMapping);
  }

  unsigned int nAtoms = mol.getNumAtoms();
  double totalVdWEnergy = 0.0;
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << "\n"
                 "V A N   D E R   W A A L S\n\n"
                 "------ATOMS------   ATOM TYPES                               "
                 "  WELL\n"
                 "  I        J          I    J    DISTANCE   ENERGY     R*     "
                 " DEPTH\n"
                 "-------------------------------------------------------------"
                 "-------"
              << std::endl;
    }
  }
  const Conformer &conf = mol.getConformer(confId);
  auto contrib = std::make_unique<VdWContrib>(field);
  bool hasContrib = false;
  for (unsigned int i = 0; i < nAtoms; ++i) {
    for (unsigned int j = i + 1; j < nAtoms; ++j) {
      if (ignoreInterfragInteractions && (fragMapping[i] != fragMapping[j])) {
        continue;
      }
      if (getTwoBitCell(neighborMatrix, twoBitCellPos(nAtoms, i, j)) >=
          RELATION_1_4) {
        double dist = (conf.getAtomPos(i) - conf.getAtomPos(j)).length();
        if (dist > nonBondedThresh) {
          continue;
        }
        MMFFVdWRijstarEps mmffVdWConstants;
        if (mmffMolProperties->getMMFFVdWParams(i, j, mmffVdWConstants)) {
          contrib->addTerm(i, j, &mmffVdWConstants);
          hasContrib = true;
          if (mmffMolProperties->getMMFFVerbosity()) {
            const Atom *iAtom = mol.getAtomWithIdx(i);
            const Atom *jAtom = mol.getAtomWithIdx(j);
            const double vdWEnergy = MMFF::Utils::calcVdWEnergy(
                dist, mmffVdWConstants.R_ij_star, mmffVdWConstants.epsilon);
            if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
              unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(i);
              unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(j);
              oStream << std::left << std::setw(2) << iAtom->getSymbol() << " #"
                      << std::setw(5) << i + 1 << std::setw(2)
                      << jAtom->getSymbol() << " #" << std::setw(5) << j + 1
                      << std::right << std::setw(5) << iAtomType << std::setw(5)
                      << jAtomType << "  " << std::fixed << std::setprecision(3)
                      << std::setw(9) << dist << std::setw(10) << vdWEnergy
                      << std::setw(9) << mmffVdWConstants.R_ij_star
                      << std::setw(9) << mmffVdWConstants.epsilon << std::endl;
            }
            totalVdWEnergy += vdWEnergy;
          }
        }
      }
    }
  }
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << std::endl;
    }
    oStream << "TOTAL VAN DER WAALS ENERGY     =" << std::right << std::setw(16)
            << std::fixed << std::setprecision(4) << totalVdWEnergy
            << std::endl;
  }
  if (hasContrib) {
    field->contribs().push_back(ForceFields::ContribPtr(contrib.release()));
  }
}

// ------------------------------------------------------------------------
//
//
//
// ------------------------------------------------------------------------
void addEle(const ROMol &mol, int confId, MMFFMolProperties *mmffMolProperties,
            ForceFields::ForceField *field,
            boost::shared_array<std::uint8_t> neighborMatrix,
            double nonBondedThresh, bool ignoreInterfragInteractions) {
  PRECONDITION(field, "bad ForceField");
  PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
  PRECONDITION(mmffMolProperties->isValid(),
               "missing atom types - invalid force-field");

  std::ostream &oStream = mmffMolProperties->getMMFFOStream();
  INT_VECT fragMapping;
  if (ignoreInterfragInteractions) {
    std::vector<ROMOL_SPTR> molFrags =
        MolOps::getMolFrags(mol, true, &fragMapping);
  }
  unsigned int nAtoms = mol.getNumAtoms();
  double totalEleEnergy = 0.0;
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << "\n"
                 "E L E C T R O S T A T I C\n\n"
                 "------ATOMS------   ATOM TYPES\n"
                 "  I        J          I    J    DISTANCE   ENERGY\n"
                 "--------------------------------------------------"
              << std::endl;
    }
  }
  const Conformer &conf = mol.getConformer(confId);
  double dielConst = mmffMolProperties->getMMFFDielectricConstant();
  std::uint8_t dielModel = mmffMolProperties->getMMFFDielectricModel();
  auto contrib = std::make_unique<EleContrib>(field);
  bool hasContrib = false;
  for (unsigned int i = 0; i < nAtoms; ++i) {
    for (unsigned int j = i + 1; j < nAtoms; ++j) {
      if (ignoreInterfragInteractions && (fragMapping[i] != fragMapping[j])) {
        continue;
      }
      std::uint8_t cell =
          getTwoBitCell(neighborMatrix, twoBitCellPos(nAtoms, i, j));
      bool is1_4 = (cell == RELATION_1_4);
      if (cell >= RELATION_1_4) {
        if (isDoubleZero(mmffMolProperties->getMMFFPartialCharge(i)) ||
            isDoubleZero(mmffMolProperties->getMMFFPartialCharge(j))) {
          continue;
        }
        double dist = (conf.getAtomPos(i) - conf.getAtomPos(j)).length();
        if (dist > nonBondedThresh) {
          continue;
        }
        double chargeTerm = mmffMolProperties->getMMFFPartialCharge(i) *
                            mmffMolProperties->getMMFFPartialCharge(j) /
                            dielConst;
        contrib->addTerm(i, j, chargeTerm, dielModel, is1_4);
        hasContrib = true;

        if (mmffMolProperties->getMMFFVerbosity()) {
          const unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(i);
          const unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(j);
          const Atom *iAtom = mol.getAtomWithIdx(i);
          const Atom *jAtom = mol.getAtomWithIdx(j);
          const double eleEnergy = MMFF::Utils::calcEleEnergy(
              i, j, dist, chargeTerm, dielModel, is1_4);
          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            oStream << std::left << std::setw(2) << iAtom->getSymbol() << " #"
                    << std::setw(5) << i + 1 << std::setw(2)
                    << jAtom->getSymbol() << " #" << std::setw(5) << j + 1
                    << std::right << std::setw(5) << iAtomType << std::setw(5)
                    << jAtomType << "  " << std::fixed << std::setprecision(3)
                    << std::setw(9) << dist << std::setw(10) << eleEnergy
                    << std::endl;
          }
          totalEleEnergy += eleEnergy;
        }
      }
    }
  }
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << std::endl;
    }
    oStream << "TOTAL ELECTROSTATIC ENERGY     =" << std::right << std::setw(16)
            << std::fixed << std::setprecision(4) << totalEleEnergy
            << std::endl;
  }
  if (hasContrib) {
    field->contribs().push_back(ForceFields::ContribPtr(contrib.release()));
  }
}

}  // end of namespace Tools

// ------------------------------------------------------------------------
//
//
//
// ------------------------------------------------------------------------
ForceFields::ForceField *constructForceField(ROMol &mol, double nonBondedThresh,
                                             int confId,
                                             bool ignoreInterfragInteractions) {
  MMFFMolProperties mmffMolProperties(mol);
  PRECONDITION(mmffMolProperties.isValid(),
               "missing atom types - invalid force-field");
  ForceFields::ForceField *res =
      constructForceField(mol, &mmffMolProperties, nonBondedThresh, confId,
                          ignoreInterfragInteractions);

  return res;
}

// ------------------------------------------------------------------------
//
//
//
// ------------------------------------------------------------------------
ForceFields::ForceField *constructForceField(
    ROMol &mol, MMFFMolProperties *mmffMolProperties, double nonBondedThresh,
    int confId, bool ignoreInterfragInteractions) {
  PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
  PRECONDITION(mmffMolProperties->isValid(),
               "missing atom types - invalid force-field");

  std::unique_ptr<ForceFields::ForceField> res(new ForceFields::ForceField());
  // add the atomic positions:
  Conformer &conf = mol.getConformer(confId);
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    res->positions().push_back(&(conf.getAtomPos(i)));
  }

  res->initialize();
  if (mmffMolProperties->getMMFFBondTerm()) {
    Tools::addBonds(mol, mmffMolProperties, res.get());
  }
  if (mmffMolProperties->getMMFFAngleTerm()) {
    Tools::addAngles(mol, mmffMolProperties, res.get());
  }
  if (mmffMolProperties->getMMFFStretchBendTerm()) {
    Tools::addStretchBend(mol, mmffMolProperties, res.get());
  }
  if (mmffMolProperties->getMMFFOopTerm()) {
    Tools::addOop(mol, mmffMolProperties, res.get());
  }
  if (mmffMolProperties->getMMFFTorsionTerm()) {
    Tools::addTorsions(mol, mmffMolProperties, res.get());
  }
  if (mmffMolProperties->getMMFFVdWTerm() ||
      mmffMolProperties->getMMFFEleTerm()) {
    boost::shared_array<std::uint8_t> neighborMat =
        Tools::buildNeighborMatrix(mol);
    if (mmffMolProperties->getMMFFVdWTerm()) {
      Tools::addVdW(mol, confId, mmffMolProperties, res.get(), neighborMat,
                    nonBondedThresh, ignoreInterfragInteractions);
    }
    if (mmffMolProperties->getMMFFEleTerm()) {
      Tools::addEle(mol, confId, mmffMolProperties, res.get(), neighborMat,
                    nonBondedThresh, ignoreInterfragInteractions);
    }
  }

  return res.release();
}
}  // namespace MMFF
}  // namespace RDKit
