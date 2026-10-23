//
//  Copyright (C) 2026 Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//
#include "MacroMolTemplateLibrary.h"

#include <GraphMol/Chirality.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <cstdlib>
#include <fstream>
#include <iterator>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace RDKit {
namespace detail {
// Defined in DefaultMonomers.cpp, generated from Data/DefaultMonomers.json
std::string_view getEmbeddedDefaultMonomerDefinitions();
}  // namespace detail

namespace {

constexpr const char *monomerLibraryEnvVar = "RDKIT_MONOMER_LIBRARY";

struct MonomerDefinitions {
  std::string json;
  std::string source;  // describes where the definitions came from
};

// Uses the JSON file named by $RDKIT_MONOMER_LIBRARY if it is set, otherwise
// the default definitions embedded at build time.
MonomerDefinitions readMonomerDefinitions() {
  const auto *fileName = std::getenv(monomerLibraryEnvVar);
  if (!fileName || !*fileName) {
    return {std::string(detail::getEmbeddedDefaultMonomerDefinitions()),
            "embedded default monomer definitions"};
  }
  std::ifstream input(fileName);
  if (!input) {
    throw std::runtime_error(std::string("Could not open monomer library '") +
                             fileName + "' specified by " +
                             monomerLibraryEnvVar);
  }
  return {std::string(std::istreambuf_iterator<char>(input),
                      std::istreambuf_iterator<char>()),
          std::string("monomer library '") + fileName + "' specified by " +
              monomerLibraryEnvVar};
}

std::string getString(const boost::property_tree::ptree &object,
                      const std::string &key) {
  // ptree::get<std::string>() would return an empty string for an object or
  // array, so check that the value has no children. An empty object can't be
  // distinguished from an empty string, so empty values are rejected too.
  const auto value = object.get_child_optional(key);
  if (!value || !value->empty() || value->data().empty()) {
    throw std::runtime_error("Missing or invalid '" + key +
                             "' in monomer definition");
  }
  return value->data();
}

MonomerClass monomerClassFromPolymerType(const std::string &polymerType) {
  if (polymerType == "PEPTIDE") {
    return MonomerClass::AminoAcid;
  }
  if (polymerType == "RNA") {
    return MonomerClass::NucleicAcid;
  }
  if (polymerType == "CHEM") {
    return MonomerClass::Chemical;
  }
  throw std::runtime_error("Unsupported polymer type in monomer definition: " +
                           polymerType);
}

std::unique_ptr<RWMol> parseCXSmiles(const std::string &smiles,
                                     const std::string &what,
                                     const std::string &name) {
  v2::SmilesParse::SmilesParserParams params;
  params.removeHs = false;
  params.parseName = false;
  auto mol = v2::SmilesParse::MolFromSmiles(smiles, params);
  if (!mol) {
    throw std::runtime_error("Could not parse " + what + " for monomer " +
                             name);
  }
  // E/Z labels depend on CIP ranks, which change when dummy atoms are
  // replaced, so express double bond stereo relative to the stereo atoms
  for (const auto bond : mol->bonds()) {
    if ((bond->getStereo() == Bond::STEREOE ||
         bond->getStereo() == Bond::STEREOZ) &&
        bond->getStereoAtoms().size() == 2) {
      bond->setStereo(Chirality::translateEZLabelToCisTrans(bond->getStereo()));
    }
  }
  return mol;
}

// Returns n if the atom is a dummy atom labelled _Rn, otherwise 0.
int getRGroupNumber(const Atom *atom) {
  std::string label;
  if (atom->getAtomicNum() != 0 ||
      !atom->getPropIfPresent(common_properties::atomLabel, label) ||
      label.size() < 3 || label.compare(0, 2, "_R") != 0 ||
      label.find_first_not_of("0123456789", 2) != std::string::npos ||
      label[2] == '0') {
    return 0;
  }
  return std::stoi(label.substr(2));
}

// Replaces the core dummy atom at coreDummyIdx with the R-group, and returns
// the resulting leaving group. The R-group atom bonded to the R-group's dummy
// atom takes the place of the core dummy atom, keeping its index and bonds so
// that the stereochemistry of the attachment atom is preserved. Any further
// R-group atoms are appended to the molecule.
MacroMolLeavingGroup attachRGroup(RWMol &mol, unsigned int coreDummyIdx,
                                  int rGroupNumber, const ROMol &rGroup,
                                  const std::string &name) {
  const auto rGroupLabel = "R" + std::to_string(rGroupNumber);
  const Atom *rGroupDummy = nullptr;
  for (const auto atom : rGroup.atoms()) {
    if (atom->getAtomicNum() == 0) {
      if (rGroupDummy || getRGroupNumber(atom) != rGroupNumber) {
        throw std::runtime_error(
            "R-group " + rGroupLabel + " for monomer " + name +
            " must contain exactly one dummy atom, labelled _" + rGroupLabel);
      }
      rGroupDummy = atom;
    }
  }
  if (!rGroupDummy || rGroupDummy->getDegree() != 1) {
    throw std::runtime_error("The dummy atom in R-group " + rGroupLabel +
                             " for monomer " + name +
                             " must have exactly one neighbor");
  }
  const auto rGroupLeavingAtom = *rGroup.atomNeighbors(rGroupDummy).begin();

  const auto coreDummy = mol.getAtomWithIdx(coreDummyIdx);
  if (coreDummy->getDegree() != 1) {
    throw std::runtime_error("Attachment point _" + rGroupLabel +
                             " in monomer " + name +
                             " must have exactly one neighbor");
  }
  const auto attachAtomIdx = (*mol.atomNeighbors(coreDummy).begin())->getIdx();

  std::vector<int> newIdxs(rGroup.getNumAtoms(), -1);
  newIdxs[rGroupLeavingAtom->getIdx()] = coreDummyIdx;
  mol.replaceAtom(coreDummyIdx, rGroupLeavingAtom);
  std::vector<unsigned int> atomIdxs{coreDummyIdx};
  for (const auto atom : rGroup.atoms()) {
    if (atom != rGroupDummy && atom != rGroupLeavingAtom) {
      newIdxs[atom->getIdx()] = mol.addAtom(atom->copy(), false, true);
      atomIdxs.push_back(newIdxs[atom->getIdx()]);
    }
  }
  // the R-group's dummy atom corresponds to the core attachment atom
  const auto mapIdx = [&](unsigned int rGroupIdx) {
    return rGroupIdx == rGroupDummy->getIdx() ? attachAtomIdx
                                              : newIdxs[rGroupIdx];
  };
  for (const auto bond : rGroup.bonds()) {
    if (bond->getBeginAtom() == rGroupDummy ||
        bond->getEndAtom() == rGroupDummy) {
      continue;  // replaced by the core bond to the attachment atom
    }
    auto newBond = bond->copy();
    newBond->setBeginAtomIdx(mapIdx(bond->getBeginAtomIdx()));
    newBond->setEndAtomIdx(mapIdx(bond->getEndAtomIdx()));
    for (auto &stereoAtomIdx : newBond->getStereoAtoms()) {
      stereoAtomIdx = mapIdx(stereoAtomIdx);
    }
    mol.addBond(newBond, true);
  }

  // The leaving atom's chiral tag refers to the order of its bonds in the
  // R-group, which differs from the order in the combined molecule (where the
  // bond to the attachment atom comes first), so adjust it for the permutation.
  const auto leavingAtom = mol.getAtomWithIdx(coreDummyIdx);
  const auto chiralTag = leavingAtom->getChiralTag();
  if (chiralTag != Atom::CHI_UNSPECIFIED) {
    if (chiralTag != Atom::CHI_TETRAHEDRAL_CW &&
        chiralTag != Atom::CHI_TETRAHEDRAL_CCW) {
      throw std::runtime_error(
          "Unsupported chirality on the atom bonded to "
          "the attachment point of R-group " +
          rGroupLabel + " for monomer " + name);
    }
    Chirality::restoreTetrahedralParity(*rGroupLeavingAtom, *leavingAtom,
                                        mapIdx);
  }
  return {std::move(atomIdxs), attachAtomIdx, static_cast<int>(coreDummyIdx),
          rGroupNumber};
}

// Builds a template from CORE_SMILES, a CXSMILES in which the attachment points
// are dummy atoms labelled _R1, _R2, etc., and RGROUPS, an object mapping the
// labels R1, R2, etc. to CXSMILES for the leaving groups. Each leaving group
// contains a dummy atom with the same label marking its attachment point.
std::unique_ptr<MacroMolTemplate> makeTemplate(
    const boost::property_tree::ptree &definition) {
  const auto symbol = getString(definition, "SYMBOL");
  const auto name = getString(definition, "PDBCODE");
  const auto coreSmiles = getString(definition, "CORE_SMILES");
  const auto monomerClass =
      monomerClassFromPolymerType(getString(definition, "POLYMER_TYPE"));
  const auto subclass = getString(definition, "MONOMER_TYPE");
  const auto rGroupDefinitions = definition.get_child_optional("RGROUPS");
  if (!rGroupDefinitions) {
    throw std::runtime_error("Missing 'RGROUPS' in definition of monomer " +
                             name);
  }

  auto mol = parseCXSmiles(coreSmiles, "CORE_SMILES", name);
  std::vector<std::pair<unsigned int, int>> coreDummies;
  for (const auto atom : mol->atoms()) {
    if (const auto rGroupNumber = getRGroupNumber(atom)) {
      coreDummies.emplace_back(atom->getIdx(), rGroupNumber);
    }
  }
  if (coreDummies.size() != rGroupDefinitions->size()) {
    throw std::runtime_error(
        "The attachment points in CORE_SMILES of monomer " + name +
        " do not match its RGROUPS");
  }

  std::vector<MacroMolLeavingGroup> leavingGroups;
  for (const auto &[coreDummyIdx, rGroupNumber] : coreDummies) {
    const auto rGroupLabel = "R" + std::to_string(rGroupNumber);
    const auto rGroupSmiles =
        rGroupDefinitions->get_child_optional(rGroupLabel);
    if (!rGroupSmiles || !rGroupSmiles->empty()) {
      throw std::runtime_error("Missing or invalid R-group " + rGroupLabel +
                               " in RGROUPS of monomer " + name);
    }
    const auto rGroup =
        parseCXSmiles(rGroupSmiles->data(), "R-group " + rGroupLabel, name);
    leavingGroups.push_back(
        attachRGroup(*mol, coreDummyIdx, rGroupNumber, *rGroup, name));
  }
  // Double bond stereo is fully described by the bond stereo and stereo atoms.
  // The bond directions from the SMILES are incomplete where dummy atoms were
  // replaced, so regenerate them before perceiving the stereochemistry (which
  // uses them with the legacy stereo perception), as parsing SMILES would.
  for (const auto bond : mol->bonds()) {
    if (bond->getBondDir() == Bond::ENDUPRIGHT ||
        bond->getBondDir() == Bond::ENDDOWNRIGHT) {
      bond->setBondDir(Bond::NONE);
    }
  }
  MolOps::sanitizeMol(*mol);
  MolOps::setDoubleBondNeighborDirections(*mol);
  const bool cleanIt = true;
  const bool force = true;
  MolOps::assignStereochemistry(*mol, cleanIt, force);

  std::vector<bool> isLeavingAtom(mol->getNumAtoms(), false);
  for (const auto &leavingGroup : leavingGroups) {
    for (const auto atomIdx : leavingGroup.atomIdxs) {
      isLeavingAtom[atomIdx] = true;
    }
  }
  std::vector<unsigned int> mainAtomIdxs;
  for (unsigned int atomIdx = 0; atomIdx < mol->getNumAtoms(); ++atomIdx) {
    if (!isLeavingAtom[atomIdx]) {
      mainAtomIdxs.push_back(atomIdx);
    }
  }

  // The original data is the monomer's JSON definition. Note that the boost
  // JSON writer escapes '/' as \/, which is valid JSON.
  std::ostringstream originalData;
  boost::property_tree::write_json(originalData, definition, false);
  auto originalDataStr = originalData.str();
  if (!originalDataStr.empty() && originalDataStr.back() == '\n') {
    originalDataStr.pop_back();
  }

  MacroMolTemplateBuilder builder(*mol, monomerClass, name, symbol,
                                  std::move(originalDataStr));
  builder.setMainGroup(std::move(mainAtomIdxs));
  builder.setSubclass(subclass);
  for (auto &leavingGroup : leavingGroups) {
    builder.addLeavingGroup(std::move(leavingGroup));
  }
  return builder.build();
}

// Builds a library from the monomer definitions. Any error, including
// duplicate monomers, is reported along with the source of the definitions.
MacroMolTemplateLibrary loadMacroMolTemplates() {
  const auto definitions = readMonomerDefinitions();
  try {
    std::istringstream input(definitions.json);
    boost::property_tree::ptree tree;
    boost::property_tree::read_json(input, tree);
    MacroMolTemplateLibrary result;
    for (const auto &definition : tree) {
      result.addTemplate(makeTemplate(definition.second));
    }
    return result;
  } catch (const std::exception &e) {
    throw std::runtime_error("Error loading " + definitions.source + ": " +
                             e.what());
  }
}

}  // namespace

void addBuiltinMacroMolTemplates(MacroMolTemplateLibrary &templates) {
  for (const auto *builtinTemplate :
       getGlobalMacroMolTemplateLibrary().getTemplates()) {
    templates.addTemplate(std::make_unique<MacroMolTemplate>(*builtinTemplate));
  }
}

const MacroMolTemplateLibrary &getGlobalMacroMolTemplateLibrary() {
  // If loading throws, the static is left uninitialized and the next call
  // retries from scratch, so a partially populated library is never exposed.
  static const MacroMolTemplateLibrary templates = loadMacroMolTemplates();
  return templates;
}

}  // namespace RDKit
