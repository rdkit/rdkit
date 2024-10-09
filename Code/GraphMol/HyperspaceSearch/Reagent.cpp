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
#include <GraphMol/MolPickler.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/HyperspaceSearch/Reagent.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace RDKit::HyperspaceSearch {

Reagent::Reagent(const RDKit::HyperspaceSearch::Reagent &other)
    : d_smiles(other.d_smiles), d_id(other.d_id) {}

Reagent::Reagent(RDKit::HyperspaceSearch::Reagent &&other)
    : d_smiles(std::move(other.d_smiles)), d_id(std::move(other.d_id)) {}

Reagent &Reagent::operator=(const RDKit::HyperspaceSearch::Reagent &other) {
  if (this == &other) {
    return *this;
  }
  d_smiles = other.d_smiles;
  d_id = other.d_id;
  if (other.d_mol) {
    d_mol.reset(new ROMol(*other.d_mol));
  } else {
    d_mol.reset();
  }
  if (other.d_pattFP) {
    d_pattFP.reset(new ExplicitBitVect(*other.d_pattFP));
  } else {
    d_pattFP.reset();
  }
  if (!other.d_connRegions.empty()) {
    d_connRegions.clear();
    std::transform(
        other.d_connRegions.begin(), other.d_connRegions.end(),
        std::back_inserter(d_connRegions),
        [&](const std::shared_ptr<ROMol> &m) -> std::shared_ptr<ROMol> {
          return std::shared_ptr<ROMol>(new ROMol(*m));
        });
  } else {
    d_connRegions.clear();
  }
  return *this;
}

Reagent &Reagent::operator=(RDKit::HyperspaceSearch::Reagent &&other) {
  if (this == &other) {
    return *this;
  }
  d_smiles = std::move(other.d_smiles);
  d_id = std::move(other.d_id);
  d_mol = std::exchange(other.d_mol, nullptr);
  d_pattFP = std::exchange(other.d_pattFP, nullptr);
  d_connRegions =
      std::vector<std::shared_ptr<ROMol>>(other.d_connRegions.size());
  for (size_t i = 0; i < d_connRegions.size(); ++i) {
    d_connRegions[i] = std::exchange(other.d_connRegions[i], nullptr);
  }
  return *this;
}

const std::unique_ptr<ROMol> &Reagent::mol() const {
  if (!d_mol) {
    d_mol = v2::SmilesParse::MolFromSmiles(d_smiles);
    d_mol->setProp<std::string>(common_properties::_Name, d_id);
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
    if (cr) {
      std::vector<std::unique_ptr<ROMol>> tmpFrags;
      MolOps::getMolFrags(*cr, tmpFrags, false);
      for (auto &f : tmpFrags) {
        d_connRegions.push_back(std::shared_ptr<ROMol>(f.release()));
      }
    }
  }
  return d_connRegions;
}

void Reagent::writeToDBStream(std::ostream &os) const {
  streamWrite(os, d_smiles);
  streamWrite(os, d_id);
  MolPickler::pickleMol(*mol(), os, PicklerOps::AllProps);
  auto pattFPstr = pattFP()->toString();
  streamWrite(os, pattFPstr);
  streamWrite(os, connRegions().size());
  for (const auto &cr : connRegions()) {
    MolPickler::pickleMol(*cr, os, PicklerOps::AllProps);
  }
}

void Reagent::readFromDBStream(std::istream &is) {
  streamRead(is, d_smiles, 0);
  streamRead(is, d_id, 0);
  d_mol.reset(new ROMol);
  MolPickler::molFromPickle(is, *d_mol);
  std::string pickle;
  streamRead(is, pickle, 0);
  d_pattFP.reset(new ExplicitBitVect(pickle));
  size_t numConnRegs;
  streamRead(is, numConnRegs);
  d_connRegions.resize(numConnRegs);
  for (size_t i = 0; i < numConnRegs; ++i) {
    d_connRegions[i].reset(new ROMol);
    MolPickler::molFromPickle(is, *d_connRegions[i]);
  }
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

}  // namespace RDKit::HyperspaceSearch
