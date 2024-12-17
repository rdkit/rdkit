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
#include <GraphMol/Depictor/RDDepictor.h>
#include <CDXStdObjects.h>

namespace RDKit {
std::string MolToChemDraw(const ROMol &mol, CDXFormat format) {
  RWMol trmol(mol);
  MolOps::Kekulize(trmol);
  if (!trmol.getNumConformers()) {
    RDDepict::compute2DCoords(trmol);
  }

  CDXObjectID object_id = 1;
  CDXDocument document(object_id++);
  CDXPage *page = new CDXPage(object_id++);
  CDXFragment *fragment = new CDXFragment(object_id++);
  page->AddChild(fragment);
  auto atom_start = object_id;
  std::vector<CDXNode *> nodes;
  nodes.reserve(trmol.getNumAtoms());

  const Conformer *conf = nullptr;
  if (trmol.getNumConformers() > 0) conf = &trmol.getConformer(0);

  auto wedgeBonds = Chirality::pickBondsToWedge(trmol, nullptr, conf);

  for (auto &atom : trmol.atoms()) {
    CDXNode *node = new CDXNode(object_id + atom->getIdx());
    auto pos = conf->getAtomPos(atom->getIdx());
    node->Position(CDXPoint2D(CDXCoordinatefromPoints(pos.x),
                              CDXCoordinatefromPoints(pos.y)));
    node->m_nodeType = kCDXNodeType_Element;
    node->m_isotope = atom->getIsotope();
    node->m_elementNum = atom->getAtomicNum();
    node->m_charge = atom->getFormalCharge();
    node->m_numHydrogens = atom->getNumExplicitHs() ? atom->getNumExplicitHs()
                                                    : kNumHydrogenUnspecified;
    nodes.push_back(node);
    fragment->AddChild(node);
  }

  for (auto &bond : trmol.bonds()) {
    CDXBond *cdxbond =
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
      case Bond::UNSPECIFIED:
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
      case 1:
        cdxbond->m_display = reverse ? kCDXBondDisplay_WedgedHashEnd
                                     : kCDXBondDisplay_WedgedHashBegin;
        break;
      case 6:
        cdxbond->m_display =
            reverse ? kCDXBondDisplay_WedgeEnd : kCDXBondDisplay_WedgeBegin;
        break;
      default:
        break;
    }

    if (bond->getBondDir() == Bond::BondDir::EITHERDOUBLE ||
        bond->getBondDir() == Bond::BondDir::UNKNOWN)
      cdxbond->m_display = kCDXBondDisplay_Wavy;

    fragment->AddChild(cdxbond);
  }

  document.AddChild(page);
  document.m_colorTable.m_colors
      .clear();  // if this isn't empty something fails.
  std::ostringstream os;
  os << kCDXML_HeaderString;
  XMLDataSink ds(os);
  document.XMLWrite(ds);
  return os.str();
}
}  // namespace RDKit
