//
//  Copyright (C) 2023 Novartis Biomedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <cmath>
#include <regex>
#include <sstream>
#include "Pipeline.h"
#include "Validate.h"
#include "Metal.h"
#include "Normalize.h"
#include "Charge.h"
#include "Fragment.h"
#include <RDGeneral/FileParseException.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Chirality.h>

namespace RDKit {
namespace MolStandardize {

void PipelineResult::append(PipelineStatus newStatus, const std::string &info) {
  status = static_cast<PipelineStatus>(status | newStatus);
  log.push_back({newStatus, info});
}

PipelineResult Pipeline::run(const std::string &molblock) const {
  PipelineResult result;
  result.status = NO_EVENT;
  result.inputMolData = molblock;

  // parse the molblock into an RWMol instance
  result.stage = static_cast<uint32_t>(PipelineStage::PARSING_INPUT);
  RWMOL_SPTR mol = parse(molblock, result, options);
  if (!mol || ((result.status & PIPELINE_ERROR) != NO_EVENT &&
               !options.reportAllFailures)) {
    return result;
  }

  RWMOL_SPTR_PAIR output;

  if (mol->getNumAtoms() == 0 && options.allowEmptyMolecules) {
    output = {mol, mol};
  } else {
    // we try sanitization and validation on a copy, because we want to preserve
    // the original input molecule for later
    RWMOL_SPTR molCopy{new RWMol(*mol)};
    for (const auto &[stage, operation] : validationSteps) {
      result.stage = stage;
      molCopy = operation(molCopy, result, options);
      if (!molCopy || ((result.status & PIPELINE_ERROR) != NO_EVENT &&
                       !options.reportAllFailures)) {
        return result;
      }
    }

    for (const auto &[stage, operation] : standardizationSteps) {
      result.stage = stage;
      mol = operation(mol, result, options);
      if (!mol || ((result.status & PIPELINE_ERROR) != NO_EVENT &&
                   !options.reportAllFailures)) {
        return result;
      }
    }
    if (makeParent) {
      result.stage = static_cast<uint32_t>(PipelineStage::MAKE_PARENT);
      output = makeParent(mol, result, options);
      if (!output.first || !output.second ||
          ((result.status & PIPELINE_ERROR) != NO_EVENT &&
           !options.reportAllFailures)) {
        return result;
      }
    } else {
      output = {mol, mol};
    }
  }

  // serialize as MolBlocks
  result.stage = static_cast<uint32_t>(PipelineStage::SERIALIZING_OUTPUT);
  serialize(output, result, options);
  if ((result.status & PIPELINE_ERROR) != NO_EVENT &&
      !options.reportAllFailures) {
    return result;
  }

  result.stage = static_cast<uint32_t>(PipelineStage::COMPLETED);

  return result;
}

namespace Operations {
RWMOL_SPTR parse(const std::string &molblock, PipelineResult &result,
                 const PipelineOptions &options) {
  v2::FileParsers::MolFileParserParams params;
  // we don't want to sanitize the molecule at this stage
  params.sanitize = false;
  // Hs wouldn't be anyway removed if the mol is not sanitized
  params.removeHs = false;
  // strict parsing is configurable via the pipeline options
  params.strictParsing = options.strictParsing;

  RWMOL_SPTR mol{};

  try {
    mol.reset(v2::FileParsers::MolFromMolBlock(molblock, params).release());
  } catch (FileParseException &e) {
    result.append(INPUT_ERROR, e.what());
  }

  if (!mol) {
    result.append(INPUT_ERROR,
                  "Could not instantiate a valid molecule from input");
  }

  return mol;
}

void serialize(RWMOL_SPTR_PAIR output, PipelineResult &result,
               const PipelineOptions &options) {
  const ROMol &outputMol = *output.first;
  const ROMol &parentMol = *output.second;

  try {
    if (!options.outputV2000) {
      result.outputMolData = MolToV3KMolBlock(outputMol);
      result.parentMolData = MolToV3KMolBlock(parentMol);
    } else {
      try {
        result.outputMolData = MolToV2KMolBlock(outputMol);
        result.parentMolData = MolToV2KMolBlock(parentMol);
      } catch (ValueErrorException &e) {
        result.append(OUTPUT_ERROR,
                      "Can't write molecule to V2000 output format: " +
                          std::string(e.what()));
      }
    }
  } catch (const std::exception &e) {
    result.append(OUTPUT_ERROR, "Can't write molecule to output format: " +
                                    std::string(e.what()));
  } catch (...) {
    result.append(
        OUTPUT_ERROR,
        "An unexpected error occurred while serializing the output structures.");
  }
}
RWMOL_SPTR prepareForValidation(RWMOL_SPTR mol, PipelineResult &result,
                                const PipelineOptions &) {
  // Prepare the mol for validation.

  try {
    // The general intention is about validating the original input, and
    // therefore limit the sanitization to the minimum, but it's not very useful
    // to record a valence validation error for issues like a badly drawn nitro
    // group that would be later fixed during by the normalization step.
    //
    // Some sanitization also needs to be performed in order to assign the
    // stereochemistry (which needs to happen prior to reapplying the wedging,
    // see below), and we need to find radicals, in order to support the
    // corresponding validation criterion.
    constexpr unsigned int sanitizeOps =
        (MolOps::SANITIZE_CLEANUP | MolOps::SANITIZE_SYMMRINGS |
         MolOps::SANITIZE_CLEANUP_ORGANOMETALLICS |
         MolOps::SANITIZE_FINDRADICALS);
    unsigned int failedOp = 0;
    MolOps::sanitizeMol(*mol, failedOp, sanitizeOps);

    // We want to restore the original MolBlock wedging, but this step may in
    // some cases overwrite the ENDDOWNRIGHT/ENDUPRIGHT info that describes the
    // configuration of double bonds adjacent to stereocenters. We therefore
    // first assign the stereochemistry, and then restore the wedging.
    constexpr bool cleanIt = true;
    constexpr bool force = true;
    constexpr bool flagPossible = true;
    MolOps::assignStereochemistry(*mol, cleanIt, force, flagPossible);
    Chirality::reapplyMolBlockWedging(*mol);
  } catch (MolSanitizeException &) {
    result.append(
        PREPARE_FOR_VALIDATION_ERROR,
        "An error occurred while preparing the molecule for validation.");
  }

  return mol;
}

namespace {
// The error messages from the ValidationMethod classes include some metadata
// in a string prefix that are not particularly useful within the context of
// this Pipeline. The function below removes that prefix.
static const std::regex prefix("^(ERROR|INFO): \\[.+\\] ");
std::string removeErrorPrefix(const std::string &message) {
  return std::regex_replace(message, prefix, "");
}
}  // namespace

RWMOL_SPTR validate(RWMOL_SPTR mol, PipelineResult &result,
                    const PipelineOptions &options) {
  auto applyValidation = [&mol, &result, &options](
                             const ValidationMethod &v,
                             PipelineStatus status) -> bool {
    auto errors = v.validate(*mol, options.reportAllFailures);
    for (const auto &error : errors) {
      result.append(status, removeErrorPrefix(error));
    }
    return errors.empty();
  };

  // check for undesired features in the input molecule (e.g., query
  // atoms/bonds)
  FeaturesValidation featuresValidation(options.allowEnhancedStereo,
                                        options.allowAromaticBondType,
                                        options.allowDativeBondType);
  if (!applyValidation(featuresValidation, FEATURES_VALIDATION_ERROR) &&
      !options.reportAllFailures) {
    return mol;
  }

  // check the number of atoms and valence status
  RDKitValidation rdkitValidation;
  if (!applyValidation(rdkitValidation, BASIC_VALIDATION_ERROR) &&
      !options.reportAllFailures) {
    return mol;
  }

  // disallow radicals
  DisallowedRadicalValidation radicalValidation;
  if (!applyValidation(radicalValidation, BASIC_VALIDATION_ERROR) &&
      !options.reportAllFailures) {
    return mol;
  }

  // validate the isotopic numbers (if any are specified)
  IsotopeValidation isotopeValidation(true);
  if (!applyValidation(isotopeValidation, BASIC_VALIDATION_ERROR) &&
      !options.reportAllFailures) {
    return mol;
  }

  // verify that the input is a 2D structure
  Is2DValidation is2DValidation(options.is2DZeroThreshold);
  if (!applyValidation(is2DValidation, IS2D_VALIDATION_ERROR) &&
      !options.reportAllFailures) {
    return mol;
  }

  // validate the 2D layout (check for clashing atoms and abnormally long bonds)
  Layout2DValidation layout2DValidation(
      options.atomClashLimit, options.bondLengthLimit,
      options.allowLongBondsInRings, options.allowAtomBondClashExemption,
      options.minMedianBondLength);
  if (!applyValidation(layout2DValidation, LAYOUT2D_VALIDATION_ERROR) &&
      !options.reportAllFailures) {
    return mol;
  }

  // verify that the specified stereochemistry is formally correct
  StereoValidation stereoValidation;
  if (!applyValidation(stereoValidation, STEREO_VALIDATION_ERROR) &&
      !options.reportAllFailures) {
    return mol;
  }

  return mol;
}

RWMOL_SPTR prepareForStandardization(RWMOL_SPTR mol, PipelineResult &result,
                                     const PipelineOptions &) {
  // Prepare the mol for standardization.

  try {
    MolOps::sanitizeMol(*mol);
  } catch (MolSanitizeException &) {
    result.append(
        PREPARE_FOR_STANDARDIZATION_ERROR,
        "An error occurred while preparing the molecule for standardization.");
  }

  return mol;
}

RWMOL_SPTR standardize(RWMOL_SPTR mol, PipelineResult &result,
                       const PipelineOptions &options) {
  auto smiles = MolToSmiles(*mol);
  auto reference = smiles;

  // bonding to metals
  try {
    MetalDisconnectorOptions mdOpts;
    MetalDisconnector metalDisconnector(mdOpts);
    std::unique_ptr<ROMol> metalNof{SmartsToMol(options.metalNof)};
    metalDisconnector.setMetalNof(*metalNof);
    std::unique_ptr<ROMol> metalNon{SmartsToMol(options.metalNon)};
    metalDisconnector.setMetalNon(*metalNon);
    metalDisconnector.disconnectInPlace(*mol);
  } catch (...) {
    result.append(
        METAL_STANDARDIZATION_ERROR,
        "An error occurred while processing the bonding of metal species.");
    return mol;
  }

  smiles = MolToSmiles(*mol);
  if (smiles != reference) {
    result.append(METALS_DISCONNECTED,
                  "One or more metal atoms were disconnected.");
  }
  reference = smiles;

  // functional groups
  try {
    std::unique_ptr<Normalizer> normalizer{};
    if (options.normalizerData.empty()) {
      normalizer.reset(new Normalizer);
    } else {
      std::istringstream sstr(options.normalizerData);
      normalizer.reset(new Normalizer(sstr, options.normalizerMaxRestarts));
    }
    // normalizeInPlace() may return an ill-formed molecule if
    // the sanitization of a transformed structure failed
    // => use normalize() instead (also see GitHub #7189)
    mol.reset(static_cast<RWMol *>(normalizer->normalize(*mol)));
    mol->updatePropertyCache(false);
  } catch (...) {
    result.append(
        NORMALIZER_STANDARDIZATION_ERROR,
        "An error occurred while normalizing the representation of some functional groups");
    return mol;
  }

  smiles = MolToSmiles(*mol);
  if (smiles != reference) {
    result.append(NORMALIZATION_APPLIED,
                  "The representation of some functional groups was adjusted.");
  }
  reference = smiles;

  // keep the largest fragment
  try {
    LargestFragmentChooser fragmentChooser;
    fragmentChooser.chooseInPlace(*mol);
  } catch (...) {
    result.append(
        FRAGMENT_STANDARDIZATION_ERROR,
        "An error occurred while removing the disconnected fragments");
    return mol;
  }

  smiles = MolToSmiles(*mol);
  if (smiles != reference) {
    result.append(
        FRAGMENTS_REMOVED,
        "One or more disconnected fragments (e.g., counterions) were removed.");
  }

  // The stereochemistry is not assigned until after we are done modifying the
  // molecular graph:
  constexpr bool cleanIt = true;
  constexpr bool force = true;
  constexpr bool flagPossible = true;
  MolOps::assignStereochemistry(*mol, cleanIt, force, flagPossible);

  return mol;
}

RWMOL_SPTR reapplyWedging(RWMOL_SPTR mol, PipelineResult &result,
                          const PipelineOptions &) {
  // in general, we want to restore the bond wedging from the input molblock,
  // but we prefer to not use any wavy bonds, because of their ambiguity
  // in some configurations.

  // we therefore proceed in two steps, we first reapply the molblock wedging
  // and then revert the changes related to double bonds with undefined/unknown
  // stereochemistry and change single bonds with "unknown" direction into plain
  // single bonds.

  // in order to do so, we need to keep track of the current bond configuration
  // settings.
  using BondInfo = std::tuple<Bond::BondType, Bond::BondDir, Bond::BondStereo>;
  std::map<unsigned int, BondInfo> oldBonds;
  for (auto bond : mol->bonds()) {
    oldBonds[bond->getIdx()] = {bond->getBondType(), bond->getBondDir(),
                                bond->getStereo()};
  }

  // 1) restore the original wedging from the input MolBlock
  Chirality::reapplyMolBlockWedging(*mol);

  // 2) revert the changes related to double bonds with stereo type "either":
  //    restore the STEREOANY direction of double bonds that have a substituent
  //    with direction UNKNOWN and are now STEREONONE
  for (auto bond : mol->bonds()) {
    if (bond->getBondType() != Bond::DOUBLE) {
      continue;
    }
    Bond::BondStereo oldStereo = std::get<2>(oldBonds[bond->getIdx()]);
    Bond::BondStereo newStereo = bond->getStereo();
    bool hasAdjacentWavy{false};
    for (auto atom : {bond->getBeginAtom(), bond->getEndAtom()}) {
      for (auto adjacentBond : mol->atomBonds(atom)) {
        if (adjacentBond == bond) {
          continue;
        }
        if (adjacentBond->getBondDir() == Bond::UNKNOWN) {
          hasAdjacentWavy = true;
        }
      }
    }
    if (hasAdjacentWavy && oldStereo == Bond::STEREOANY &&
        newStereo == Bond::STEREONONE) {
      bond->setStereo(Bond::STEREOANY);
      result.append(
          NORMALIZATION_APPLIED,
          "Double bond " + std::to_string(bond->getIdx()) +
              " was assigned an undefined/unknown stereochemical configuration");
    }
  }

  // 3) set the bond direction to NONE for bonds with direction UNKNOWN
  for (auto bond : mol->bonds()) {
    if (bond->getBondDir() != Bond::UNKNOWN) {
      continue;
    }
    bond->setBondDir(Bond::NONE);
    result.append(NORMALIZATION_APPLIED, "The \"wavy\" style of bond " +
                                             std::to_string(bond->getIdx()) +
                                             " was removed");
  }

  return mol;
}

RWMOL_SPTR cleanup2D(RWMOL_SPTR mol, PipelineResult & /*result*/,
                     const PipelineOptions &options) {
  // scale the atoms coordinates
  // and make sure that z coords are set to 0 (some z coords may be non-null
  // albeit smaller than the validation threshold - these noisy coords may in
  // some cases also interfere with the perception of stereochemistry by some
  // tools e.g., inchi)
  if (options.scaledMedianBondLength > 0. && mol->getNumConformers()) {
    auto &conf = mol->getConformer();
    double medianBondLength =
        sqrt(Layout2DValidation::squaredMedianBondLength(*mol, conf));
    if (medianBondLength > options.minMedianBondLength) {
      double scaleFactor = options.scaledMedianBondLength / medianBondLength;
      unsigned int natoms = conf.getNumAtoms();
      for (unsigned int i = 0; i < natoms; ++i) {
        auto pos = conf.getAtomPos(i) * scaleFactor;
        pos.z = 0.;
        conf.setAtomPos(i, pos);
      }
    }
  }

  return mol;
}

namespace {
void replaceDativeBonds(RWMOL_SPTR mol) {
  bool modified{false};
  for (auto bond : mol->bonds()) {
    if (bond->getBondType() != Bond::BondType::DATIVE) {
      continue;
    }
    auto donor = bond->getBeginAtom();
    donor->setFormalCharge(donor->getFormalCharge() + 1);
    auto acceptor = bond->getEndAtom();
    acceptor->setFormalCharge(acceptor->getFormalCharge() - 1);
    bond->setBondType(Bond::BondType::SINGLE);
    modified = true;
  }
  if (modified) {
    mol->updatePropertyCache(false);
  }
}

void removeHsAtProtonatedSites(RWMOL_SPTR mol) {
  boost::dynamic_bitset<> protons{mol->getNumAtoms(), 0};
  for (auto atom : mol->atoms()) {
    if (atom->getAtomicNum() != 1 || atom->getDegree() != 1) {
      continue;
    }
    for (auto neighbor : mol->atomNeighbors(atom)) {
      if (neighbor->getFormalCharge() > 0) {
        protons.set(atom->getIdx());
      }
    }
  }
  if (protons.any()) {
    for (int idx = mol->getNumAtoms() - 1; idx >= 0; --idx) {
      if (!protons[idx]) {
        continue;
      }
      auto atom = mol->getAtomWithIdx(idx);
      for (auto bond : mol->atomBonds(atom)) {
        auto neighbor = bond->getOtherAtom(atom);
        neighbor->setNumExplicitHs(neighbor->getNumExplicitHs() + 1);
        break;  // there are no other bonds anyways
      }
      mol->removeAtom(atom);
    }
    mol->updatePropertyCache(false);
  }
}
}  // namespace

RWMOL_SPTR_PAIR makeParent(RWMOL_SPTR mol, PipelineResult &result,
                           const PipelineOptions &) {
  auto reference = MolToSmiles(*mol);

  RWMOL_SPTR parent{new RWMol(*mol)};

  // A "parent" structure is constructed here, in order to provide a
  // representation of the original input that may be more suitable for
  // identification purposes even though it may not reflect the most stable
  // physical state or nicest representation for the compound.
  //
  // The two steps that are currently implemented for this procedure consist in
  // normalizing the overall charge status and replacing any explicit dative
  // bonds.
  //
  // If the input was submitted in an unsuitable protonation status, the
  // neutralized parent structure may become the actual output from the
  // standardization.

  // overall charge status
  try {
    // The Uncharger implementation wouldn't identify the positively
    // charged sites with adjacent explicit Hs correctly (it's a quite
    // unlikely configuration, but potentially possible considering that
    // the pipeline operates on unsanitized input).
    //
    // If present, these Hs are therefore removed from the molecular graph
    // prior to neutralization.
    removeHsAtProtonatedSites(parent);

    static const bool canonicalOrdering = false;
    static const bool force = true;
    static const bool protonationOnly = true;
    Uncharger uncharger(canonicalOrdering, force, protonationOnly);
    uncharger.unchargeInPlace(*parent);
  } catch (...) {
    result.append(
        CHARGE_STANDARDIZATION_ERROR,
        "An error occurred while normalizing the compound's charge status");
    return {{}, {}};
  }

  // Check if `mol` was submitted in a suitable ionization state
  int parentCharge{};
  for (auto atom : parent->atoms()) {
    parentCharge += atom->getFormalCharge();
  }

  int molCharge{};
  for (auto atom : mol->atoms()) {
    molCharge += atom->getFormalCharge();
  }

  // If mol is neutral or in a protonation state that partially or fully
  // balances the non-neutralizable charged sites in the parent structure,
  // then mol is accepted. Otherwise, it is replaced by its parent.
  if ((molCharge > 0 && molCharge > parentCharge) ||
      (molCharge < 0 && molCharge < parentCharge)) {
    mol = parent;
  }

  auto smiles = MolToSmiles(*mol);
  if (smiles != reference) {
    result.append(PROTONATION_CHANGED, "The protonation state was adjusted.");
  }
  reference = smiles;

  // normalize the dative bonds
  replaceDativeBonds(parent);

  return {mol, parent};
}
}  // namespace Operations

}  // namespace MolStandardize
}  // namespace RDKit
