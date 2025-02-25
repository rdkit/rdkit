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
#include <GraphMol/FileParsers/FileWriters.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/SynthonSpaceSearch/Synthon.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace RDKit::SynthonSpaceSearch {

Synthon::Synthon(const std::string &smi) : d_smiles(smi) {
  v2::SmilesParse::SmilesParserParams params;
  params.sanitize = false;
  dp_origMol = v2::SmilesParse::MolFromSmiles(d_smiles, params);
  if (!dp_origMol) {
    // This should be rare, as it should be possible to assume that
    // the people who made the SynthonSpace know what they're doing.
    // Therefore, it's probably a corrupted or incorrect file, so
    // bring it all down.
    throw ValueErrorException("Unparsable synthon SMILES " + d_smiles);
  }
}

Synthon::Synthon(const Synthon &other)
    : d_smiles(other.d_smiles),
      dp_origMol(std::make_unique<ROMol>(*other.dp_origMol)),
      dp_searchMol(std::make_unique<ROMol>(*other.dp_searchMol)),
      dp_pattFP(std::make_unique<ExplicitBitVect>(*other.dp_pattFP)),
      dp_FP(std::make_unique<ExplicitBitVect>(*other.dp_FP)),
      d_connRegions(other.d_connRegions) {}

Synthon &Synthon::operator=(const Synthon &other) {
  if (this == &other) {
    return *this;
  }
  d_smiles = other.d_smiles;
  if (other.dp_origMol) {
    dp_origMol = std::make_unique<ROMol>(*other.dp_origMol);
  } else {
    dp_origMol.reset();
  }
  if (other.dp_searchMol) {
    dp_searchMol = std::make_unique<ROMol>(*other.dp_searchMol);
  } else {
    dp_searchMol.reset();
  }
  if (other.dp_pattFP) {
    dp_pattFP = std::make_unique<ExplicitBitVect>(*other.dp_pattFP);
  } else {
    dp_pattFP.reset();
  }
  if (other.dp_FP) {
    dp_FP = std::make_unique<ExplicitBitVect>(*other.dp_FP);
  } else {
    dp_FP.reset();
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
const std::unique_ptr<ROMol> &Synthon::getOrigMol() const { return dp_origMol; }
const std::unique_ptr<ROMol> &Synthon::getSearchMol() const {
  return dp_searchMol;
}

const std::unique_ptr<ExplicitBitVect> &Synthon::getPattFP() const {
  return dp_pattFP;
}
const std::unique_ptr<ExplicitBitVect> &Synthon::getFP() const { return dp_FP; }

const std::vector<std::shared_ptr<ROMol>> &Synthon::getConnRegions() const {
  return d_connRegions;
}

void Synthon::setSearchMol(std::unique_ptr<ROMol> mol) {
  dp_searchMol = std::move(mol);
  // There are probably extraneous props on the atoms and bonds
  for (auto &atom : dp_searchMol->atoms()) {
    atom->clearProp("molNum");
    atom->clearProp("idx");
  }
  for (auto &bond : dp_searchMol->bonds()) {
    bond->clearProp("molNum");
    bond->clearProp("idx");
  }
  finishInitialization();
}

void Synthon::setFP(std::unique_ptr<ExplicitBitVect> fp) {
  dp_FP = std::move(fp);
}

void Synthon::writeToDBStream(std::ostream &os) const {
  streamWrite(os, d_smiles);
  MolPickler::pickleMol(*dp_origMol, os, PicklerOps::AllProps);
  MolPickler::pickleMol(*dp_searchMol, os, PicklerOps::AllProps);
  const auto pattFPstr = getPattFP()->toString();
  streamWrite(os, pattFPstr);
  streamWrite(os, static_cast<std::uint64_t>(getConnRegions().size()));
  for (const auto &cr : getConnRegions()) {
    MolPickler::pickleMol(*cr, os, PicklerOps::AllProps);
  }
  if (dp_FP) {
    streamWrite(os, true);
    streamWrite(os, dp_FP->toString());
  } else {
    streamWrite(os, false);
  }
}

void Synthon::readFromDBStream(std::istream &is) {
  streamRead(is, d_smiles, 0);
  dp_origMol = std::make_unique<ROMol>();
  MolPickler::molFromPickle(is, *dp_origMol);
  dp_searchMol = std::make_unique<ROMol>();
  MolPickler::molFromPickle(is, *dp_searchMol);
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
  bool haveFP = false;
  streamRead(is, haveFP);
  if (haveFP) {
    std::string pickle;
    streamRead(is, pickle, 0);
    dp_FP = std::make_unique<ExplicitBitVect>(pickle);
  }
}

void Synthon::finishInitialization() {
  dp_pattFP.reset(PatternFingerprintMol(*dp_searchMol, PATT_FP_NUM_BITS));
  d_connRegions.clear();
  if (const auto cr = details::buildConnRegion(*dp_searchMol); cr) {
    std::vector<std::unique_ptr<ROMol>> tmpFrags;
    MolOps::getMolFrags(*cr, tmpFrags, false);
    for (auto &f : tmpFrags) {
      d_connRegions.push_back(std::shared_ptr<ROMol>(f.release()));
    }
  }
}

}  // namespace RDKit::SynthonSpaceSearch
