//
// Copyright (C) 2025 David Cosgrove and other RDKit contributors.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <cmath>
#include <limits>

#include <boost/numeric/conversion/cast.hpp>

#include <RDGeneral/export.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/EnumerateStereoisomers/EnumerateStereoisomers.h>
#include <GraphMol/EnumerateStereoisomers/Flippers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace RDKit {
namespace EnumerateStereoisomers {
StereoisomerEnumerator::StereoisomerEnumerator(
    const ROMol &mol, const StereoEnumerationOptions &options, bool verbose)
    : d_mol(mol), d_options(options), d_verbose(verbose) {
  if (d_mol.getNumConformers()) {
    Chirality::wedgeMolBonds(d_mol, &d_mol.getConformer());
  }
  buildFlippers();
  // Clear unhelpful stuff out
  for (auto atom : d_mol.atoms()) {
    atom->clearProp("_CIPCode");
  }
  for (auto bond : d_mol.bonds()) {
    if (bond->getBondDir() == Bond::BondDir::EITHERDOUBLE ||
        bond->getBondDir() == Bond::BondDir::UNKNOWN) {
      bond->setBondDir(Bond::BondDir::NONE);
    }
  }
  if (d_flippers.empty()) {
    d_totalPoss = 1;
  } else {
    d_mol.setProp<std::string>("_MolFileChiralFlag", "1");
    try {
      d_totalPoss = boost::numeric_cast<unsigned long>(std::pow(
          std::uint64_t(2), static_cast<std::uint64_t>(d_flippers.size())));
    } catch (boost::numeric::positive_overflow &e) {
      d_totalPoss = std::numeric_limits<std::uint64_t>::max();
    }
  }

  if (d_options.maxIsomers) {
    d_numToReturn = std::min(getStereoisomerCount(), d_options.maxIsomers);
  } else {
    d_numToReturn = getStereoisomerCount();
  }
  if (d_options.randomSeed == -1) {
    d_randGen.reset(new std::mt19937(std::random_device()()));
  } else {
    d_randGen.reset(new std::mt19937(d_options.randomSeed));
  }
}

std::uint64_t StereoisomerEnumerator::getStereoisomerCount() const {
  return d_totalPoss;
}

std::unique_ptr<ROMol> StereoisomerEnumerator::next() {
  if (d_numReturned == d_numToReturn) {
    return std::unique_ptr<ROMol>();
  }
  if (d_flippers.empty()) {
    ++d_numReturned;
    return std::unique_ptr<ROMol>(new ROMol(d_mol));
  }
  else {
    auto isomer = generateRandomIsomer();
    ++d_numReturned;
    return isomer;
  }
}

void StereoisomerEnumerator::buildFlippers() {
  auto sis = Chirality::findPotentialStereo(d_mol, true, true);
  for (const auto &si : sis) {
    if (d_options.onlyUnassigned &&
        si.specified != Chirality::StereoSpecified::Unknown &&
        si.specified != Chirality::StereoSpecified::Unspecified) {
      continue;
    }
    if (si.type == Chirality::StereoType::Atom_Tetrahedral) {
      d_flippers.push_back(std::unique_ptr<details::Flipper>(
          new details::AtomFlipper(d_mol, si)));
    } else if (si.type == Chirality::StereoType::Bond_Double) {
      std::unique_ptr<details::BondFlipper> newFlipper(
          new details::BondFlipper(d_mol, si));
      if (newFlipper->dp_bond) {
        d_flippers.push_back(std::move(newFlipper));
      }
    } else if (si.type == Chirality::StereoType::Bond_Atropisomer) {
      std::unique_ptr<details::AtropisomerFlipper> newFlipper(
          new details::AtropisomerFlipper(d_mol, si));
      d_flippers.push_back(std::move(newFlipper));
    }
  }

  if (d_options.onlyUnassigned) {
    // otherwise these will be counted twice
    for (const auto &group : d_mol.getStereoGroups()) {
      if (group.getGroupType() != StereoGroupType::STEREO_ABSOLUTE) {
        d_flippers.push_back(std::unique_ptr<details::Flipper>(
            new details::StereoGroupFlipper(group)));
      }
    }
  }
}

std::unique_ptr<ROMol> StereoisomerEnumerator::generateRandomIsomer() {
  boost::dynamic_bitset<> nextConfig{d_flippers.size()};
  while (d_seen.size() < d_totalPoss) {
    for (size_t i = 0; i < d_flippers.size(); i++) {
      bool config = d_randDis(*d_randGen);
      nextConfig[i] = config;
    }
    if (d_seen.find(nextConfig) == d_seen.end()) {
      d_seen.insert(nextConfig);
      for (size_t i = 0; i < d_flippers.size(); i++) {
        d_flippers[i]->flip(nextConfig[i]);
      }
      // We don't need StereoGroups any more so remove them.
      std::unique_ptr<ROMol> isomer;
      if (!d_mol.getStereoGroups().empty()) {
        isomer.reset(new RWMol(d_mol));
        isomer->setStereoGroups(std::vector<StereoGroup>());
      } else {
        isomer.reset(new ROMol(d_mol));
      }
      MolOps::setDoubleBondNeighborDirections(*isomer);
      isomer->clearComputedProps(false);
      MolOps::assignStereochemistry(*isomer, true, true, true);
      if (d_options.unique) {
        auto smi =
            MolToCXSmiles(*isomer, SmilesWriteParams(),
                          SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS);
        if (d_generatedIsomers.find(smi) != d_generatedIsomers.end()) {
          continue;
        }
        d_generatedIsomers.insert(smi);
      }

      if (d_options.tryEmbedding) {
        if (embeddable(*isomer)) {
          return isomer;
        }
        if (d_verbose) {
          BOOST_LOG(rdInfoLog)
              << MolToSmiles(*isomer) << "     failed to embed." << std::endl;
        }
      } else {
        return isomer;
      }
    }
  }
  return std::unique_ptr<ROMol>();
}

bool StereoisomerEnumerator::embeddable(ROMol &isomer) {
  std::unique_ptr<ROMol> isomerWithHs(MolOps::addHs(isomer));
  auto cid = DGeomHelpers::EmbedMolecule(*isomerWithHs);
  if (cid >= 0) {
    Conformer *conf = new Conformer(isomer.getNumAtoms());
    for (unsigned int i = 0; i < isomer.getNumAtoms(); i++) {
      conf->setAtomPos(i, isomerWithHs->getConformer().getAtomPos(i));
    }
    isomer.addConformer(conf);
  }
  return cid >= 0;
}

}  // namespace EnumerateStereoisomers
}  // namespace RDKit
