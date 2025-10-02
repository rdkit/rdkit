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
//
#include "chemdraw.h"
#include <GraphMol/Chirality.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/Depictor/RDDepictor.h>

#include "ChemDrawStartInclude.h"
#include "chemdraw/CDXStdObjects.h"
#include "ChemDrawEndInclude.h"


namespace RDKit {
namespace v2 {
const double DEFAULT_CDX_BOND_LENGTH = 14.4;

namespace {
// Do we need to set explicit hs in chemdraw, this uses basically the same
//  logic as SmilesWriter
bool needsExplicitHs(const Atom *atom) {
  auto num = atom->getAtomicNum();
  const INT_VECT &defaultVs = PeriodicTable::getTable()->getValenceList(num);
  int totalValence = atom->getTotalValence();
  bool nonStandard = false;

  if (atom->getNumRadicalElectrons()) {
    nonStandard = true;
  } else if ((num == 7 || num == 15) && atom->getIsAromatic() &&
             atom->getNumExplicitHs()) {
    // another type of "nonstandard" valence is an aromatic N or P with
    // explicit Hs indicated:
    nonStandard = true;
  } else {
    nonStandard = (totalValence != defaultVs.front() && atom->getTotalNumHs());
  }

  return nonStandard;
}
}  // namespace

std::string MolToChemDrawBlock(const ROMol &mol, CDXFormat format) {
  RWMol trmol(mol);
  MolOps::Kekulize(trmol);
  if (!trmol.getNumConformers()) {
    RDDepict::compute2DCoords(trmol);
  }

  CDXObjectID object_id = 1;
  CDXDocument document(object_id++);
  auto *page = new CDXPage(object_id++);
  document.m_bondLength = DEFAULT_CDX_BOND_LENGTH;
  document.m_flags |= CDXDocument::CDXDocumentProperty1::has_bondLength;
  auto *fragment = new CDXFragment(object_id++);
  page->AddChild(fragment);
  std::vector<CDXNode *> nodes;
  nodes.reserve(trmol.getNumAtoms());

  const Conformer *conf = nullptr;
  if (trmol.getNumConformers() == 0) {
    RDDepict::compute2DCoords(trmol);
  }
  conf = &trmol.getConformer(0);
  bool is3D = conf->is3D();

  // I REALLY don't know why this is 2*DEFAULT_CDX_BOND_LENGTH but it looks
  // right
  //   when loading the CDX into ChemDraw
  // We convert the average bond length into the target bond length here
  double target_bond_length = 2 * DEFAULT_CDX_BOND_LENGTH;
  double dist = 0.0;
  for (auto bond : trmol.bonds()) {
    auto pos1 = conf->getAtomPos(bond->getBeginAtomIdx());
    auto pos2 = conf->getAtomPos(bond->getEndAtomIdx());
    dist += (pos1 - pos2).length();
  }
  dist /= trmol.getNumBonds();
  double scale = is3D ? 1. : target_bond_length / dist;

  auto wedgeBonds = Chirality::pickBondsToWedge(trmol, nullptr, conf);

  for (auto &atom : trmol.atoms()) {
    auto *node = new CDXNode(object_id + atom->getIdx());
    auto pos = conf->getAtomPos(atom->getIdx());
    if (is3D) {
      node->Position3D(CDXPoint3D(CDXCoordinatefromPoints(pos.x),
                                  -CDXCoordinatefromPoints(pos.y),
                                  CDXCoordinatefromPoints(pos.z)));
    } else {
      node->Position(CDXPoint2D(CDXCoordinatefromPoints(scale * pos.x),
                                CDXCoordinatefromPoints(-scale * pos.y)));
    }
    node->m_nodeType = kCDXNodeType_Element;
    node->m_isotope = atom->getIsotope();
    node->m_elementNum = atom->getAtomicNum();
    // Use the same logic from the smiles writer needs brackets
    //
    // node->m_numHydrogens = atom->getNumExplicitHs() ?
    // atom->getNumExplicitHs()
    //                                                : kNumHydrogenUnspecified;
    node->m_numHydrogens =
        needsExplicitHs(atom) ? atom->getTotalNumHs() : kNumHydrogenUnspecified;
    node->m_charge = atom->getFormalCharge() * 0x1000000;
    if (atom->getFormalCharge() || atom->getNumRadicalElectrons() != 0) {
      node->m_numHydrogens =
          atom->getTotalNumHs();  // XXX is this right?  We seem to need to set
                                  // it with charges
    }
    if (atom->getNumRadicalElectrons()) {
      switch (atom->getNumRadicalElectrons()) {
        case 0:
          break;
        case 1:
          node->m_radical = kCDXRadical_Singlet;
          break;
        case 2:
          break;
        case 3:
          break;
      }
    }
    // this might be a bit slow, perhaps make into a map...
    unsigned int sgnum = 0;
    for (auto &sg : trmol.getStereoGroups()) {
      sgnum++;
      for (auto &sgatom : sg.getAtoms()) {
        if (atom->getIdx() == sgatom->getIdx()) {
          switch (sg.getGroupType()) {
            case StereoGroupType::STEREO_ABSOLUTE:
              node->m_enhancedStereoType = kCDXEnhancedStereo_Absolute;
              break;
            case StereoGroupType::STEREO_OR:
              node->m_enhancedStereoType = kCDXEnhancedStereo_Or;
              break;
            case StereoGroupType::STEREO_AND:
              node->m_enhancedStereoType = kCDXEnhancedStereo_And;
              break;
          }
          node->m_enhancedStereoGroupNum = sgnum;
        }
      }
    }
    nodes.push_back(node);
    fragment->AddChild(node);
  }

  for (auto &bond : trmol.bonds()) {
    auto *cdxbond =
        new CDXBond(object_id + mol.getNumAtoms() + bond->getIdx());

    int dirCode = 0;
    bool reverse = false;
    Chirality::GetMolFileBondStereoInfo(bond, wedgeBonds, conf, dirCode,
                                        reverse);

    switch (bond->getBondType()) {
      case Bond::BondType::SINGLE:
        cdxbond->m_bondOrder = kCDXBondOrder_Single;
        break;
      case Bond::DOUBLE:
        cdxbond->m_bondOrder = kCDXBondOrder_Double;
        break;
      case Bond::TRIPLE:
        cdxbond->m_bondOrder = kCDXBondOrder_Triple;
        break;
      case Bond::QUADRUPLE:
        cdxbond->m_bondOrder = kCDXBondOrder_Quadruple;
        break;
      case Bond::QUINTUPLE:
        cdxbond->m_bondOrder = kCDXBondOrder_Quintuple;
        break;
      case Bond::HEXTUPLE:
        cdxbond->m_bondOrder = kCDXBondOrder_Sextuple;
        break;
      case Bond::ONEANDAHALF:
        cdxbond->m_bondOrder = kCDXBondOrder_OneHalf;
        break;
      case Bond::TWOANDAHALF:
        cdxbond->m_bondOrder = kCDXBondOrder_TwoHalf;
        break;
      case Bond::THREEANDAHALF:
        cdxbond->m_bondOrder = kCDXBondOrder_ThreeHalf;
        break;
      case Bond::FOURANDAHALF:
        cdxbond->m_bondOrder = kCDXBondOrder_FourHalf;
        break;
      case Bond::FIVEANDAHALF:
        cdxbond->m_bondOrder = kCDXBondOrder_FiveHalf;
        break;
      case Bond::AROMATIC:
        cdxbond->m_bondOrder = kCDXBondOrder_OneHalf;
        break;
      case Bond::IONIC:
        cdxbond->m_bondOrder = kCDXBondOrder_Ionic;
        break;
      case Bond::HYDROGEN:
        cdxbond->m_bondOrder = kCDXBondOrder_Hydrogen;
        break;
      case Bond::THREECENTER:
        cdxbond->m_bondOrder = kCDXBondOrder_ThreeCenter;
        break;
      case Bond::DATIVE:
        cdxbond->m_bondOrder = kCDXBondOrder_Dative;
        break;
      case Bond::UNSPECIFIED: {
        auto query = describeQuery(bond);
        if (query == "DoubleOrAromaticBond 1 = val\n") {
          cdxbond->m_bondOrder = kCDXBondOrder_DoubleOrAromatic;
        } else if (query == "SingleOrAromaticBond 1 = val\n") {
          cdxbond->m_bondOrder = kCDXBondOrder_SingleOrAromatic;
        } else if (query == "SingleOrDoubleBond 1 = val\n") {
          cdxbond->m_bondOrder = kCDXBondOrder_SingleOrDouble;
        } else {
          cdxbond->m_bondOrder = kCDXBondOrder_Any;
        }
      } break;
      case Bond::DATIVEONE:
      case Bond::DATIVEL:
      case Bond::DATIVER:
      case Bond::OTHER:
      case Bond::ZERO:
        // unhandled
        break;
    }

    cdxbond->Connects(nodes[bond->getBeginAtomIdx()],
                      nodes[bond->getEndAtomIdx()]);

    switch (dirCode) {
      case 6:  // swap 1 and 6 due to swapped y
        cdxbond->m_display = reverse ? kCDXBondDisplay_WedgedHashEnd
                                     : kCDXBondDisplay_WedgedHashBegin;
        break;
      case 1:
        cdxbond->m_display =
            reverse ? kCDXBondDisplay_WedgeEnd : kCDXBondDisplay_WedgeBegin;
        break;
      case 3:
        cdxbond->m_display = kCDXBondDisplay_Wavy;
        break;
      default:
        break;
    }

    if (bond->getBondDir() == Bond::BondDir::EITHERDOUBLE ||
        bond->getBondDir() == Bond::BondDir::UNKNOWN) {
      cdxbond->m_display = kCDXBondDisplay_Wavy;
}

    fragment->AddChild(cdxbond);
  }

  document.AddChild(page);
  document.m_colorTable.m_colors
      .clear();  // if this isn't empty something fails.
  
  std::ostringstream os;
  if(format == CDXFormat::CDXML) {
    os << kCDXML_HeaderString;
    XMLDataSink ds(os);
    document.XMLWrite(ds);
  } else {
    CDXostream ds(os);
    CDXWriteDocToStorage(&document, ds);
  }
  return os.str();
}
}
}  // namespace RDKit
