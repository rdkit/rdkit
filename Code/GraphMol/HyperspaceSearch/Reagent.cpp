//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/HyperspaceSearch/Reagent.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace RDKit::HyperspaceSSSearch {

const std::unique_ptr<ROMol> &Reagent::mol() const {
  if (!d_mol) {
    d_mol = v2::SmilesParse::MolFromSmiles(d_smiles);
  }
  return d_mol;
}

const std::unique_ptr<ExplicitBitVect> &Reagent::pattFP() const {
  if (!d_pattFP) {
    d_pattFP.reset(PatternFingerprintMol(*mol(), 2048));
  }
  return d_pattFP;
}

const std::vector<std::shared_ptr<ROMol>> &Reagent::connRegions() const {
  if (d_connRegions.empty()) {
    auto cr = getConnRegion(*mol());
    std::vector<std::unique_ptr<ROMol>> tmpFrags;
    MolOps::getMolFrags(*cr, tmpFrags, false);
    for (auto &f : tmpFrags) {
      d_connRegions.push_back(std::shared_ptr<ROMol>(f.release()));
    }
  }
  return d_connRegions;
}

// Get the connector regions for this molecule, which are parts of the
// molecule within 3 bonds of an isotopically labelled dummy atom.  It's
// possible the molecule has ordinary wildcard dummy atoms that shouldn't
// have an isotopic label.  Returns an empty unique_ptr if there isn't a
// labelled dummy atom, which happens when the molecule is the whole original
// query before any splits were done.
std::unique_ptr<ROMol> getConnRegion(const ROMol &mol) {
  boost::dynamic_bitset<> inFrag(mol.getNumAtoms());
  for (const auto a : mol.atoms()) {
    if (!a->getAtomicNum() && a->getIsotope()) {
      inFrag[a->getIdx()] = true;
      for (const auto &n1 : mol.atomNeighbors(a)) {
        if (!inFrag[n1->getIdx()]) {
          inFrag[n1->getIdx()] = true;
          for (const auto &n2 : mol.atomNeighbors(n1)) {
            if (!inFrag[n2->getIdx()]) {
              inFrag[n2->getIdx()] = true;
              for (const auto &n3 : mol.atomNeighbors(n2)) {
                inFrag[n3->getIdx()] = true;
              }
            }
          }
        }
      }
    }
  }
  if (!inFrag.count()) {
    return std::unique_ptr<RWMol>();
  }

  std::unique_ptr<RWMol> molCp(new RWMol(mol));
  molCp->beginBatchEdit();
  for (auto &aCp : molCp->atoms()) {
    if (!inFrag[aCp->getIdx()]) {
      molCp->removeAtom(aCp);
    } else {
      if (!aCp->getAtomicNum()) {
        aCp->setIsotope(1);
        if (aCp->hasQuery()) {
          aCp->expandQuery(makeAtomIsotopeQuery(1), Queries::COMPOSITE_OR);
        }
      }
    }
  }
  molCp->commitBatchEdit();
  return molCp;
}

}  // namespace RDKit::HyperspaceSSSearch
