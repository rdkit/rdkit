//
//  Copyright (c) 2024, Glysade Inc
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// #include "node.h"
#include "utils.h"
#include "fragment.h"

namespace RDKit {
namespace ChemDraw {
bool parseBond(RWMol &mol, unsigned int fragmentId, CDXBond &bond,
               PageData &pagedata) {
  int bond_id = bond.GetObjectID();
  Atom *start_atom = pagedata.atomIds[bond.m_beginNodeID];
  Atom *end_atom = pagedata.atomIds[bond.m_endNodeID];
  if ((!start_atom || !end_atom)) {
    BOOST_LOG(rdErrorLog) << "Bad bond in CDXML skipping fragment "
                          << fragmentId << "..." << std::endl;
    return false;
  }
  Bond::BondType order = Bond::UNSPECIFIED;
  std::unique_ptr<QueryBond> qb;
  switch (bond.m_bondOrder) {
    case kCDXBondOrder_Single:
      order = Bond::BondType::SINGLE;
      break;
    case kCDXBondOrder_Double:
      order = Bond::BondType::DOUBLE;
      break;
    case kCDXBondOrder_Triple:
      order = Bond::BondType::TRIPLE;
      break;
    case kCDXBondOrder_Quadruple:
      order = Bond::BondType::QUADRUPLE;
      break;
    case kCDXBondOrder_Quintuple:
      order = Bond::BondType::QUINTUPLE;
      break;
    case kCDXBondOrder_Sextuple:
      order = Bond::BondType::HEXTUPLE;
      break;
    case kCDXBondOrder_OneHalf:
      order = Bond::BondType::AROMATIC;
      start_atom->setIsAromatic(true);
      end_atom->setIsAromatic(true);
      break;
    case kCDXBondOrder_TwoHalf:
      order = Bond::BondType::TWOANDAHALF;
      break;
    case kCDXBondOrder_ThreeHalf:
      order = Bond::BondType::THREEANDAHALF;
      break;
    case kCDXBondOrder_FourHalf:
      order = Bond::BondType::FOURANDAHALF;
      break;
    case kCDXBondOrder_FiveHalf:
      order = Bond::BondType::FIVEANDAHALF;
      break;
    case kCDXBondOrder_Dative:
      order = Bond::BondType::DATIVE;
      break;
    case kCDXBondOrder_Ionic:
      order = Bond::BondType::IONIC;
      break;
    case kCDXBondOrder_SingleOrDouble: {
      order = Bond::BondType::SINGLE;
      qb = std::make_unique<QueryBond>();
      qb->setQuery(makeSingleOrDoubleBondQuery());
      break;
    }
    case kCDXBondOrder_SingleOrAromatic: {
      order = Bond::BondType::SINGLE;
      qb = std::make_unique<QueryBond>();
      qb->setQuery(makeSingleOrAromaticBondQuery());
      break;
    }
    case kCDXBondOrder_DoubleOrAromatic: {
      order = Bond::BondType::DOUBLE;
      qb = std::make_unique<QueryBond>();
      qb->setQuery(makeDoubleOrAromaticBondQuery());
      break;
    }
    case kCDXBondOrder_Any: {
      qb = std::make_unique<QueryBond>();
      qb->setQuery(makeBondNullQuery());
      break;
    }
    case kCDXBondOrder_Hydrogen:
      BOOST_LOG(rdErrorLog)
          << "Unhandled bond order Hydrogen, skipping fragment" << std::endl;
      return false;
    case kCDXBondOrder_ThreeCenter:
      BOOST_LOG(rdErrorLog)
          << "Unhandled bond order ThreeCenter, skipping fragment" << std::endl;
      return false;
    case kCDXBondOrder_Half:
      BOOST_LOG(rdErrorLog)
          << "Unhandled bond order Half, skipping fragment" << std::endl;
      return false;
    default:
      BOOST_LOG(rdErrorLog) << "Bad bond, skipping fragment" << std::endl;
      return false;
  };

  // The RDKit only supports one direction for wedges so
  //  normalize it
  bool swap_bond_ends = false;
  switch (bond.m_display) {
    case kCDXBondDisplay_Solid:
      break;
    case kCDXBondDisplay_Dash:
      break;
    case kCDXBondDisplay_Hash:
      break;
    case kCDXBondDisplay_WedgedHashBegin:
      break;
    case kCDXBondDisplay_WedgedHashEnd:
      swap_bond_ends = true;
      break;
    case kCDXBondDisplay_Bold:
      break;
    case kCDXBondDisplay_WedgeBegin:
      break;
    case kCDXBondDisplay_WedgeEnd:
      swap_bond_ends = true;
      break;
    case kCDXBondDisplay_Wavy:
      break;
    case kCDXBondDisplay_HollowWedgeBegin:
      break;
    case kCDXBondDisplay_HollowWedgeEnd:
      break;
    case kCDXBondDisplay_WavyWedgeBegin:
      break;
    case kCDXBondDisplay_WavyWedgeEnd:
      break;
    case kCDXBondDisplay_Dot:
      break;
    case kCDXBondDisplay_DashDot:
      break;
    case kCDXBondDisplay_DottedHydrogen:
      break;
  }

  unsigned int bondIdx = 0;
  auto startIdx = start_atom->getIdx();
  auto endIdx = end_atom->getIdx();
  if (swap_bond_ends) {
    std::swap(startIdx, endIdx);
  }

  if (qb) {
    qb->setBeginAtomIdx(startIdx);
    qb->setEndAtomIdx(endIdx);
    bondIdx = mol.addBond(qb.release(), true) - 1;
  } else {
    bondIdx = mol.addBond(startIdx, endIdx, order) - 1;
  }

  Bond *bnd = mol.getBondWithIdx(bondIdx);
  if (order == Bond::BondType::AROMATIC) {
    bnd->setIsAromatic(true);
    bnd->getBeginAtom()->setIsAromatic(true);
    bnd->getEndAtom()->setIsAromatic(true);
  }
  bnd->setProp(CDX_BOND_ID, bond.GetObjectID());

  switch (bond.m_display) {
    case kCDXBondDisplay_WedgedHashBegin:
    case kCDXBondDisplay_WedgedHashEnd: {
      bnd->setBondDir(Bond::BondDir::BEGINDASH);
      bnd->setProp(common_properties::_MolFileBondCfg, 3);
    } break;
    case kCDXBondDisplay_WedgeBegin:
    case kCDXBondDisplay_WedgeEnd: {
      bnd->setBondDir(Bond::BondDir::BEGINWEDGE);
      bnd->setProp(common_properties::_MolFileBondCfg, 1);
    } break;
    case kCDXBondDisplay_Wavy: {
      switch (order) {
        case Bond::BondType::SINGLE:
          bnd->setBondDir(Bond::BondDir::UNKNOWN);
          bnd->setProp(common_properties::_MolFileBondCfg, 2);
          break;
        case Bond::BondType::DOUBLE:
          bnd->setBondDir(Bond::BondDir::EITHERDOUBLE);
          bnd->setStereo(Bond::STEREOANY);
          break;
        default:
          BOOST_LOG(rdWarningLog)
              << "ignoring Wavy bond set on a non double bond id: " << bond_id
              << std::endl;
      }
      break;

      default:
        break;
    }
  }
  return true;
}
}  // namespace ChemDraw
}  // namespace RDKit
