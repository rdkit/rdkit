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
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/SynthonSpaceSearch/Synthon.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

namespace RDKit::SynthonSpaceSearch {

Synthon::Synthon(const std::string &smi, const std::string &id)
    : d_smiles(smi), d_id(id) {
  dp_mol = v2::SmilesParse::MolFromSmiles(d_smiles);
  if (!dp_mol) {
    // This should be rare, as it should be possible to assume that
    // the people who made the SynthonSpace know what they're doing.
    // Therefore, it's probably a corrupted or incorrect file, so
    // bring it all down.
    throw ValueErrorException("Unparsable synthon SMILES " + d_smiles +
                              " with ID " + d_id);
  }
  dp_mol->setProp<std::string>(common_properties::_Name, d_id);

  dp_pattFP.reset(PatternFingerprintMol(*getMol(), 2048));

  if (auto cr = getConnRegion(*getMol()); cr) {
    std::vector<std::unique_ptr<ROMol>> tmpFrags;
    MolOps::getMolFrags(*cr, tmpFrags, false);
    for (auto &f : tmpFrags) {
      d_connRegions.push_back(std::shared_ptr<ROMol>(f.release()));
    }
  }
}

Synthon::Synthon(const Synthon &other)
    : d_smiles(other.d_smiles),
      d_id(other.d_id),
      dp_mol(std::make_unique<ROMol>(*other.dp_mol)),
      dp_pattFP(std::make_unique<ExplicitBitVect>(*other.dp_pattFP)),
      d_connRegions(other.d_connRegions) {}

Synthon &Synthon::operator=(const Synthon &other) {
  if (this == &other) {
    return *this;
  }
  d_smiles = other.d_smiles;
  d_id = other.d_id;
  if (other.dp_mol) {
    dp_mol = std::make_unique<ROMol>(*other.dp_mol);
  } else {
    dp_mol.reset();
  }
  if (other.dp_pattFP) {
    dp_pattFP = std::make_unique<ExplicitBitVect>(*other.dp_pattFP);
  } else {
    dp_pattFP.reset();
  }
  if (!other.d_connRegions.empty()) {
    d_connRegions.clear();
    std::transform(
        other.d_connRegions.begin(), other.d_connRegions.end(),
        std::back_inserter(d_connRegions),
        [](const std::shared_ptr<ROMol> &m) -> std::shared_ptr<ROMol> {
          return std::make_shared<ROMol>(*m);
        });
  } else {
    d_connRegions.clear();
  }
  return *this;
}

const std::unique_ptr<ROMol> &Synthon::getMol() const { return dp_mol; }

const std::unique_ptr<ExplicitBitVect> &Synthon::getPattFP() const {
  return dp_pattFP;
}

const std::vector<std::shared_ptr<ROMol>> &Synthon::getConnRegions() const {
  return d_connRegions;
}

void Synthon::writeToDBStream(std::ostream &os) const {
  streamWrite(os, d_smiles);
  streamWrite(os, d_id);
  MolPickler::pickleMol(*getMol(), os, PicklerOps::AllProps);
  auto pattFPstr = getPattFP()->toString();
  streamWrite(os, pattFPstr);
  streamWrite(os, getConnRegions().size());
  for (const auto &cr : getConnRegions()) {
    MolPickler::pickleMol(*cr, os, PicklerOps::AllProps);
  }
}

void Synthon::readFromDBStream(std::istream &is) {
  streamRead(is, d_smiles, 0);
  streamRead(is, d_id, 0);
  dp_mol = std::make_unique<ROMol>();
  MolPickler::molFromPickle(is, *dp_mol);
  std::string pickle;
  streamRead(is, pickle, 0);
  dp_pattFP = std::make_unique<ExplicitBitVect>(pickle);
  size_t numConnRegs;
  streamRead(is, numConnRegs);
  d_connRegions.resize(numConnRegs);
  for (size_t i = 0; i < numConnRegs; ++i) {
    d_connRegions[i] = std::make_shared<ROMol>();
    MolPickler::molFromPickle(is, *d_connRegions[i]);
  }
}

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

}  // namespace RDKit::SynthonSpaceSearch
