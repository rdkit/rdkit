//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <cmath>

#include <GraphMol/MolOps.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/EnumerateStereoisomers/EnumerateStereoisomers.h>
#include <GraphMol/EnumerateStereoisomers/Flippers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace RDKit::EnumerateStereoisomers {
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
    if (bond->getBondDir() == Bond::BondDir::EITHERDOUBLE) {
      bond->setBondDir(Bond::BondDir::NONE);
    }
  }
  if (d_flippers.empty()) {
    d_totalPoss = 1;
  } else {
    d_mol.setProp<std::string>("_MolFileChiralFlag", "1");
    d_totalPoss = std::pow(2u, d_flippers.size());
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

unsigned int StereoisomerEnumerator::getStereoisomerCount() const {
  return d_totalPoss;
}

std::unique_ptr<ROMol> StereoisomerEnumerator::next() {
  if (d_numReturned == d_numToReturn) {
    return std::unique_ptr<ROMol>();
  }
  auto isomer = generateRandomIsomer();
  ++d_numReturned;
  return isomer;
}

void StereoisomerEnumerator::buildFlippers() {
  auto sis = Chirality::findPotentialStereo(d_mol, true, true);
  for (const auto &si : sis) {
    if (d_options.onlyUnassigned &&
        si.specified != Chirality::StereoSpecified::Unknown &&
        si.specified != Chirality::StereoSpecified::Unspecified) {
      continue;
    }
    std::cout << si.type << " : " << si.centeredOn;
    if (si.controllingAtoms.empty()) {
      std::cout << std::endl;
    } else {
      std::cout << " :: ";
      for (auto ca : si.controllingAtoms) {
        std::cout << ca << " ";
      }
      std::cout << std::endl;
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
    for (auto &group : d_mol.getStereoGroups()) {
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
    if (auto it = d_seen.find(nextConfig); it == d_seen.end()) {
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
            MolToCXSmiles(d_mol, SmilesWriteParams(),
                          SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS);
        if (auto it = d_generatedIsomers.find(smi);
            it != d_generatedIsomers.end()) {
          continue;
        }
        d_generatedIsomers.insert(smi);
      }

      if (d_options.tryEmbedding) {
        if (embeddable(*isomer)) {
          return isomer;
        }
        if (d_verbose) {
          std::cout << MolToSmiles(*isomer) << "     failed to embed."
                    << std::endl;
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

}  // namespace RDKit::EnumerateStereoisomers
