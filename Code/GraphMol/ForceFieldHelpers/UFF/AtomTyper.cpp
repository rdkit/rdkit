//
//  Copyright (C) 2004-2021 Greg Landrun and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <iostream>
#include <GraphMol/RDKitBase.h>
#include <ForceField/UFF/Params.h>
#include <ForceField/UFF/BondStretch.h>
#include <ForceField/UFF/AngleBend.h>
#include <ForceField/UFF/TorsionAngle.h>
#include <ForceField/UFF/Inversion.h>
#include <ForceField/UFF/Nonbonded.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include "AtomTyper.h"

namespace RDKit {
namespace UFF {
using namespace ForceFields::UFF;

namespace Tools {
// ---------------------------------------------------------------
void addAtomChargeFlags(const Atom *atom, std::string &atomKey,
                        bool tolerateChargeMismatch) {
  PRECONDITION(atom, "bad atom");
  int totalValence = atom->getTotalValence();
  int fc = atom->getFormalCharge();
  // FIX: come up with some way of handling metals here
  switch (atom->getAtomicNum()) {
    // Atoms only +1 in default UFF params
    case 29:  // Cu
    case 47:  // Ag
      if (totalValence == 1 || fc == 1 || tolerateChargeMismatch) {
        atomKey += "+1";
      } else {
        BOOST_LOG(rdErrorLog)
            << "UFFTYPER: Unrecognized charge state for atom: "
            << atom->getIdx() << std::endl;
      }
      break;

    // Atoms only +2 in default UFF params
    case 4:   // Be
    case 20:  // Ca
    case 25:  // Mn
    case 26:  // Fe
    case 28:  // Ni
              //  case 30:  // Zn
    case 46:  // Pd
              //  case 48:  // Cd
    case 78:  // Pt
      if (totalValence == 2 || fc == 2 || tolerateChargeMismatch) {
        atomKey += "+2";
      } else {
        BOOST_LOG(rdErrorLog)
            << "UFFTYPER: Unrecognized charge state for atom: "
            << atom->getIdx() << std::endl;
      }
      break;

    // Atoms only +3 in default UFF params
    case 21:   // Sc
    case 24:   // Cr
    case 27:   // Co
               //  case 49:  // In
    case 79:   // Au
    case 89:   // Ac
    case 96:   // Cm
    case 97:   // Bk
    case 98:   // Cf
    case 99:   // Es
    case 100:  // Fm
    case 101:  // Md
    case 102:  // No
    case 103:  // Lr/Lw
      if (totalValence == 3 || fc == 3 || tolerateChargeMismatch) {
        atomKey += "+3";
      } else {
        BOOST_LOG(rdErrorLog)
            << "UFFTYPER: Unrecognized charge state for atom: "
            << atom->getIdx() << std::endl;
      }
      break;

    // Atoms only +4 in default UFF params
    case 2:   // He
    case 18:  // Ar
    case 22:  // Ti
    case 36:  // Kr
    case 54:  // Xe
    case 90:  // Th
    case 91:  // Pa
    case 92:  // U
    case 93:  // Np
    case 94:  // Pu
    case 95:  // Am
      if (totalValence == 4 || fc == 4 || tolerateChargeMismatch) {
        atomKey += "+4";
      } else {
        BOOST_LOG(rdErrorLog)
            << "UFFTYPER: Unrecognized charge state for atom: "
            << atom->getIdx() << std::endl;
      }
      break;

    // Atoms only +5 in default UFF params
    case 23:  // V
    case 41:  // Nb
    case 43:  // Tc
    case 73:  // Ta
      if (totalValence == 5 || fc == 5 || tolerateChargeMismatch) {
        atomKey += "+5";
      } else {
        BOOST_LOG(rdErrorLog)
            << "UFFTYPER: Unrecognized charge state for atom: "
            << atom->getIdx() << std::endl;
      }
      break;

    // Atoms only +6 in default UFF params
    case 42:  // Mo
      if (totalValence == 6 || fc == 6 || tolerateChargeMismatch) {
        atomKey += "+6";
      } else {
        BOOST_LOG(rdErrorLog)
            << "UFFTYPER: Unrecognized charge state for atom: "
            << atom->getIdx() << std::endl;
      }
      break;

    case 12:  // Mg
      switch (totalValence) {
        case 2:
          atomKey += "+2";
          break;
        default:
          if (tolerateChargeMismatch) {
            atomKey += "+2";
          }
          BOOST_LOG(rdErrorLog)
              << "UFFTYPER: Unrecognized charge state for atom: "
              << atom->getIdx() << std::endl;
      }
      break;
    case 13:  // Al
      if (totalValence != 3) {
        BOOST_LOG(rdErrorLog)
            << "UFFTYPER: Unrecognized charge state for atom: "
            << atom->getIdx() << std::endl;
      }
      break;
    case 14:  // Si
      if (totalValence != 4) {
        BOOST_LOG(rdErrorLog)
            << "UFFTYPER: Unrecognized charge state for atom: "
            << atom->getIdx() << std::endl;
      }
      break;
    case 15:  // P
      switch (totalValence) {
        case 3:
          atomKey += "+3";
          break;
        case 5:
          atomKey += "+5";
          break;
        default:
          if (tolerateChargeMismatch) {
            atomKey += "+5";
          }
          BOOST_LOG(rdErrorLog)
              << "UFFTYPER: Unrecognized charge state for atom: "
              << atom->getIdx() << std::endl;
      }
      break;
    case 16:  // S
      if (atom->getHybridization() != Atom::SP2) {
        switch (totalValence) {
          case 2:
            atomKey += "+2";
            break;
          case 4:
            atomKey += "+4";
            break;
          case 6:
            atomKey += "+6";
            break;
          default:
            if (tolerateChargeMismatch) {
              atomKey += "+6";
            }
            BOOST_LOG(rdErrorLog)
                << "UFFTYPER: Unrecognized charge state for atom: "
                << atom->getIdx() << std::endl;
        }
      }
      break;
    case 30:  // Zn
      switch (totalValence) {
        case 2:
          atomKey += "+2";
          break;
        default:
          if (tolerateChargeMismatch) {
            atomKey += "+2";
          }
          BOOST_LOG(rdErrorLog)
              << "UFFTYPER: Unrecognized charge state for atom: "
              << atom->getIdx() << std::endl;
      }
      break;
    case 31:  // Ga
      switch (totalValence) {
        case 3:
          atomKey += "+3";
          break;
        default:
          if (tolerateChargeMismatch) {
            atomKey += "+3";
          }
          BOOST_LOG(rdErrorLog)
              << "UFFTYPER: Unrecognized charge state for atom: "
              << atom->getIdx() << std::endl;
      }
      break;
    case 33:  // As
      switch (totalValence) {
        case 3:
          atomKey += "+3";
          break;
        default:
          if (tolerateChargeMismatch) {
            atomKey += "+3";
          }
          BOOST_LOG(rdErrorLog)
              << "UFFTYPER: Unrecognized charge state for atom: "
              << atom->getIdx() << std::endl;
      }
      break;
    case 34:  // Se
      switch (totalValence) {
        case 2:
          atomKey += "+2";
          break;
        default:
          if (tolerateChargeMismatch) {
            atomKey += "+2";
          }
          BOOST_LOG(rdErrorLog)
              << "UFFTYPER: Unrecognized charge state for atom: "
              << atom->getIdx() << std::endl;
      }
      break;
    case 48:  // Cd
      switch (totalValence) {
        case 2:
          atomKey += "+2";
          break;
        default:
          if (tolerateChargeMismatch) {
            atomKey += "+2";
          }
          BOOST_LOG(rdErrorLog)
              << "UFFTYPER: Unrecognized charge state for atom: "
              << atom->getIdx() << std::endl;
      }
      break;
    case 49:  // In
      switch (totalValence) {
        case 3:
          atomKey += "+3";
          break;
        default:
          if (tolerateChargeMismatch) {
            atomKey += "+3";
          }
          BOOST_LOG(rdErrorLog)
              << "UFFTYPER: Unrecognized charge state for atom: "
              << atom->getIdx() << std::endl;
      }
      break;
    case 51:  // Sb
      switch (totalValence) {
        case 3:
          atomKey += "+3";
          break;
        default:
          if (tolerateChargeMismatch) {
            atomKey += "+3";
          }
          BOOST_LOG(rdErrorLog)
              << "UFFTYPER: Unrecognized charge state for atom: "
              << atom->getIdx() << std::endl;
      }
      break;
    case 52:  // Te
      switch (totalValence) {
        case 2:
          atomKey += "+2";
          break;
        default:
          if (tolerateChargeMismatch) {
            atomKey += "+2";
          }
          BOOST_LOG(rdErrorLog)
              << "UFFTYPER: Unrecognized charge state for atom: "
              << atom->getIdx() << std::endl;
      }
      break;
    case 75:  // Re
      if (tolerateChargeMismatch) {
        if (atomKey == "Re6") {
          atomKey = "Re6+5";
        } else if (atomKey == "Re3") {
          atomKey = "Re3+7";
        }
      }
      BOOST_LOG(rdErrorLog)
          << "UFFTYPER: Unrecognized charge state for atom: " << atom->getIdx()
          << std::endl;
      break;
    case 80:  // Hg
      switch (totalValence) {
        case 2:
          atomKey += "+2";
          break;
        default:
          if (tolerateChargeMismatch) {
            atomKey += "+2";
          }
          BOOST_LOG(rdErrorLog)
              << "UFFTYPER: Unrecognized charge state for atom: "
              << atom->getIdx() << std::endl;
      }
      break;
    case 81:  // Tl
      switch (totalValence) {
        case 3:
          atomKey += "+3";
          break;
        default:
          if (tolerateChargeMismatch) {
            atomKey += "+3";
          }
          BOOST_LOG(rdErrorLog)
              << "UFFTYPER: Unrecognized charge state for atom: "
              << atom->getIdx() << std::endl;
      }
      break;
    case 82:  // Pb
      switch (totalValence) {
        case 3:
          atomKey += "+3";
          break;
        default:
          if (tolerateChargeMismatch) {
            atomKey += "+3";
          }
          BOOST_LOG(rdErrorLog)
              << "UFFTYPER: Unrecognized charge state for atom: "
              << atom->getIdx() << std::endl;
      }
      break;
    case 83:  // Bi
      switch (totalValence) {
        case 3:
          atomKey += "+3";
          break;
        default:
          if (tolerateChargeMismatch) {
            atomKey += "+3";
          }
          BOOST_LOG(rdErrorLog)
              << "UFFTYPER: Unrecognized charge state for atom: "
              << atom->getIdx() << std::endl;
      }
      break;
    case 84:  // Po
      switch (totalValence) {
        case 2:
          atomKey += "+2";
          break;
        default:
          if (tolerateChargeMismatch) {
            atomKey += "+2";
          }
          BOOST_LOG(rdErrorLog)
              << "UFFTYPER: Unrecognized charge state for atom: "
              << atom->getIdx() << std::endl;
      }
      break;
  }
  // lanthanides
  if (atom->getAtomicNum() >= 57 && atom->getAtomicNum() <= 71) {
    switch (totalValence) {
      case 6:
        atomKey += "+3";
        break;
      default:
        if (tolerateChargeMismatch) {
          atomKey += "+3";
        }
        BOOST_LOG(rdErrorLog)
            << "UFFTYPER: Unrecognized charge state for atom: "
            << atom->getIdx() << std::endl;
    }
  }
}

// ---------------------------------------------------------------
std::string getAtomLabel(const Atom *atom) {
  PRECONDITION(atom, "bad atom");
  int atNum = atom->getAtomicNum();
  std::string atomKey = atom->getSymbol();
  if (atomKey.size() == 1) {
    atomKey += '_';
  }
  PeriodicTable *table = PeriodicTable::getTable();

  // FIX: handle main group/organometallic cases better:
  if (atNum) {
    // do not do hybridization on alkali metals or halogens:
    if (table->getDefaultValence(atNum) == -1 ||
        (table->getNouterElecs(atNum) != 1 &&
         table->getNouterElecs(atNum) != 7)) {
      switch (atom->getAtomicNum()) {
        case 12:
        case 13:
        case 14:
        case 15:
        case 50:
        case 51:
        case 52:
        case 81:
        case 82:
        case 83:
        case 84:
          atomKey += '3';
          if (atom->getHybridization() != Atom::SP3) {
            BOOST_LOG(rdWarningLog)
                << "UFFTYPER: Warning: hybridization set to SP3 for atom "
                << atom->getIdx() << std::endl;
          }
          break;
        case 80:
          atomKey += '1';
          if (atom->getHybridization() != Atom::SP) {
            BOOST_LOG(rdWarningLog)
                << "UFFTYPER: Warning: hybridization set to SP for atom "
                << atom->getIdx() << std::endl;
          }
          break;
        default:
          switch (atom->getHybridization()) {
            case Atom::S:
              // don't need to do anything here
              break;
            case Atom::SP:
              atomKey += '1';
              break;

            case Atom::SP2:
              if ((atom->getIsAromatic() ||
                   MolOps::atomHasConjugatedBond(atom)) &&
                  (atNum == 6 || atNum == 7 || atNum == 8 || atNum == 16)) {
                atomKey += 'R';
              } else {
                atomKey += '2';
              }
              break;

            case Atom::SP3:
              atomKey += '3';
              break;

            case Atom::SP2D:
              atomKey += '4';
              break;

            case Atom::SP3D:
              atomKey += '5';
              break;

            case Atom::SP3D2:
              atomKey += '6';
              break;

            default:
              BOOST_LOG(rdErrorLog)
                  << "UFFTYPER: Unrecognized hybridization for atom: "
                  << atom->getIdx() << std::endl;
          }
      }
    }
  }
  // special cases by element type:
  addAtomChargeFlags(atom, atomKey);
  return atomKey;
}
}  // end of namespace Tools

// ---------------------------------------------------------------
std::pair<AtomicParamVect, bool> getAtomTypes(const ROMol &mol,
                                              const std::string &) {
  bool foundAll = true;
  auto params = ParamCollection::getParams();

  AtomicParamVect paramVect;
  paramVect.resize(mol.getNumAtoms());

  for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    const Atom *atom = mol.getAtomWithIdx(i);

    // construct the atom key:
    std::string atomKey = Tools::getAtomLabel(atom);

    // ok, we've got the atom key, now get the parameters:
    const AtomicParams *theParams = (*params)(atomKey);
    if (!theParams) {
      foundAll = false;
      BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized atom type: " << atomKey
                            << " (" << i << ")" << std::endl;
    }

    paramVect[i] = theParams;
  }

  return std::make_pair(paramVect, foundAll);
}

bool getUFFBondStretchParams(const ROMol &mol, unsigned int idx1,
                             unsigned int idx2, UFFBond &uffBondStretchParams) {
  auto params = ParamCollection::getParams();
  unsigned int idx[2] = {idx1, idx2};
  AtomicParamVect paramVect(2);
  unsigned int i;
  const Bond *bond = mol.getBondBetweenAtoms(idx1, idx2);
  bool res = bond != nullptr;
  for (i = 0; res && (i < 2); ++i) {
    const Atom *atom = mol.getAtomWithIdx(idx[i]);
    std::string atomKey = Tools::getAtomLabel(atom);
    paramVect[i] = (*params)(atomKey);
    res = paramVect[i] != nullptr;
  }
  if (res) {
    double bondOrder = bond->getBondTypeAsDouble();
    uffBondStretchParams.r0 =
        UFF::Utils::calcBondRestLength(bondOrder, paramVect[0], paramVect[1]);
    uffBondStretchParams.kb = UFF::Utils::calcBondForceConstant(
        uffBondStretchParams.r0, paramVect[0], paramVect[1]);
  }
  return res;
}

bool getUFFAngleBendParams(const ROMol &mol, unsigned int idx1,
                           unsigned int idx2, unsigned int idx3,
                           UFFAngle &uffAngleBendParams) {
  auto params = ParamCollection::getParams();
  unsigned int idx[3] = {idx1, idx2, idx3};
  AtomicParamVect paramVect(3);
  unsigned int i;
  const Bond *bond[2];
  bool res = true;
  for (i = 0; res && (i < 3); ++i) {
    if (i < 2) {
      bond[i] = mol.getBondBetweenAtoms(idx[i], idx[i + 1]);
      res = bond[i] != nullptr;
    }
    if (res) {
      const Atom *atom = mol.getAtomWithIdx(idx[i]);
      std::string atomKey = Tools::getAtomLabel(atom);
      paramVect[i] = (*params)(atomKey);
      res = paramVect[i] != nullptr;
    }
  }
  if (res) {
    double bondOrder12 = bond[0]->getBondTypeAsDouble();
    double bondOrder23 = bond[1]->getBondTypeAsDouble();
    uffAngleBendParams.theta0 = RAD2DEG * paramVect[1]->theta0;
    uffAngleBendParams.ka = UFF::Utils::calcAngleForceConstant(
        paramVect[1]->theta0, bondOrder12, bondOrder23, paramVect[0],
        paramVect[1], paramVect[2]);
  }
  return res;
}

bool getUFFTorsionParams(const ROMol &mol, unsigned int idx1, unsigned int idx2,
                         unsigned int idx3, unsigned int idx4,
                         UFFTor &uffTorsionParams) {
  auto params = ParamCollection::getParams();
  unsigned int idx[4] = {idx1, idx2, idx3, idx4};
  AtomicParamVect paramVect(2);
  unsigned int i;
  const Bond *bond = mol.getBondBetweenAtoms(idx2, idx3);
  int atNum[2];
  Atom::HybridizationType hyb[2];
  bool res = true;
  bool hasSP2 = false;
  for (i = 0; res && (i < 4); ++i) {
    if (i < 3) {
      res = mol.getBondBetweenAtoms(idx[i], idx[i + 1]) != nullptr;
    }
    const Atom *atom = mol.getAtomWithIdx(idx[i]);
    if ((i == 1) || (i == 2)) {
      unsigned int j = i - 1;
      atNum[j] = atom->getAtomicNum();
      hyb[j] = atom->getHybridization();
      std::string atomKey = Tools::getAtomLabel(atom);
      paramVect[j] = (*params)(atomKey);
      res = paramVect[j] != nullptr;
    } else if (atom->getHybridization() == Atom::SP2) {
      hasSP2 = true;
    }
  }
  if (res) {
    res = (((hyb[0] == RDKit::Atom::SP2) || (hyb[0] == RDKit::Atom::SP3)) &&
           ((hyb[1] == RDKit::Atom::SP2) || (hyb[1] == RDKit::Atom::SP3)));
  }
  if (res) {
    double bondOrder = bond->getBondTypeAsDouble();
    if ((hyb[0] == RDKit::Atom::SP3) && (hyb[1] == RDKit::Atom::SP3)) {
      // general case:
      uffTorsionParams.V = sqrt(paramVect[0]->V1 * paramVect[1]->V1);
      // special case for single bonds between group 6 elements:
      if (((int)(bondOrder * 10) == 10) && UFF::Utils::isInGroup6(atNum[0]) &&
          UFF::Utils::isInGroup6(atNum[1])) {
        double V2 = 6.8;
        double V3 = 6.8;
        if (atNum[0] == 8) {
          V2 = 2.0;
        }
        if (atNum[1] == 8) {
          V3 = 2.0;
        }
        uffTorsionParams.V = sqrt(V2 * V3);
      }
    } else if ((hyb[0] == RDKit::Atom::SP2) && (hyb[1] == RDKit::Atom::SP2)) {
      uffTorsionParams.V =
          UFF::Utils::equation17(bondOrder, paramVect[0], paramVect[1]);
    } else {
      // SP2 - SP3,  this is, by default, independent of atom type in UFF:
      uffTorsionParams.V = 1.0;
      if ((int)(bondOrder * 10) == 10) {
        // special case between group 6 sp3 and non-group 6 sp2:
        if (((hyb[0] == RDKit::Atom::SP3) && UFF::Utils::isInGroup6(atNum[0]) &&
             (!UFF::Utils::isInGroup6(atNum[1]))) ||
            ((hyb[1] == RDKit::Atom::SP3) && UFF::Utils::isInGroup6(atNum[1]) &&
             (!UFF::Utils::isInGroup6(atNum[0])))) {
          uffTorsionParams.V =
              UFF::Utils::equation17(bondOrder, paramVect[0], paramVect[1]);
        }
        // special case for sp3 - sp2 - sp2
        // (i.e. the sp2 has another sp2 neighbor, like propene)
        else if (hasSP2) {
          uffTorsionParams.V = 2.0;
        }
      }
    }
  }
  return res;
}

bool getUFFInversionParams(const ROMol &mol, unsigned int idx1,
                           unsigned int idx2, unsigned int idx3,
                           unsigned int idx4, UFFInv &uffInversionParams) {
  unsigned int idx[4] = {idx1, idx2, idx3, idx4};
  bool res = (mol.getBondBetweenAtoms(idx1, idx2) &&
              mol.getBondBetweenAtoms(idx2, idx3) &&
              mol.getBondBetweenAtoms(idx2, idx4));
  unsigned int i;
  // bool isAtom2C = false;
  bool isBoundToSP2O = false;
  unsigned int at2AtomicNum = 0;
  for (i = 0; res && (i < 4); ++i) {
    const Atom *atom = mol.getAtomWithIdx(idx[i]);
    if (i == 1) {
      at2AtomicNum = atom->getAtomicNum();
      if (res) {
        // if the central atom is not carbon, nitrogen, oxygen,
        // phosphorous, arsenic, antimonium or bismuth, skip it
        res = (!(((at2AtomicNum != 6) && (at2AtomicNum != 7) &&
                  (at2AtomicNum != 8) && (at2AtomicNum != 15) &&
                  (at2AtomicNum != 33) && (at2AtomicNum != 51) &&
                  (at2AtomicNum != 83)) ||
                 (atom->getDegree() != 3)));
      }
      if (res) {
        // if the central atom is carbon, nitrogen or oxygen
        // but hybridization is not sp2, skip it
        res = (!(((at2AtomicNum == 6) || (at2AtomicNum == 7) ||
                  (at2AtomicNum == 8)) &&
                 (atom->getHybridization() != Atom::SP2)));
      }
    } else if ((atom->getAtomicNum() == 8) &&
               (atom->getHybridization() == Atom::SP2)) {
      isBoundToSP2O = true;
    }
  }
  if (res) {
    isBoundToSP2O = (isBoundToSP2O && (at2AtomicNum == 6));
    boost::tuple<double, double, double, double> invCoeffForceCon =
        UFF::Utils::calcInversionCoefficientsAndForceConstant(at2AtomicNum,
                                                              isBoundToSP2O);
    uffInversionParams.K = boost::tuples::get<0>(invCoeffForceCon);
  }
  return res;
}

bool getUFFVdWParams(const ROMol &mol, unsigned int idx1, unsigned int idx2,
                     UFFVdW &uffVdWParams) {
  bool res = true;
  auto params = ParamCollection::getParams();
  unsigned int idx[2] = {idx1, idx2};
  AtomicParamVect paramVect(2);
  unsigned int i;
  for (i = 0; res && (i < 2); ++i) {
    const Atom *atom = mol.getAtomWithIdx(idx[i]);
    std::string atomKey = Tools::getAtomLabel(atom);
    paramVect[i] = (*params)(atomKey);
    res = paramVect[i] != nullptr;
  }
  if (res) {
    uffVdWParams.x_ij =
        UFF::Utils::calcNonbondedMinimum(paramVect[0], paramVect[1]);
    uffVdWParams.D_ij =
        UFF::Utils::calcNonbondedDepth(paramVect[0], paramVect[1]);
  }
  return res;
}
}  // namespace UFF
}  // namespace RDKit
