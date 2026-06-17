//
//  Copyright (c) 2007-2026 Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/wstring.h>
#include <nanobind/stl/tuple.h>

#include <GraphMol/MolPickler.h>
#include <GraphMol/nbWrap/props.hpp>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionRunner.h>
#include <GraphMol/ChemReactions/PreprocessRxn.h>
#include <GraphMol/ChemReactions/SanitizeRxn.h>
#include <GraphMol/MarvinParse/MarvinParser.h>
#include <GraphMol/Depictor/DepictUtils.h>
#include <GraphMol/FilterCatalog/FunctionalGroupHierarchy.h>

#include <RDBoost/Wrap_nb.h>
#include <nanobind/stl/shared_ptr.h>

#include <RDGeneral/Exceptions.h>
#include <GraphMol/SanitException.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/ChemReactions/ReactionFingerprints.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>

namespace nb = nanobind;
using namespace nb::literals;

void wrap_enumeration(nb::module_ &m);

namespace RDKit {

std::string pyObjectToString(nb::object input) {
  if (nb::isinstance<nb::str>(input)) {
    return nb::cast<std::string>(input);
  }
  std::wstring ws = nb::cast<std::wstring>(input);
  return std::string(ws.begin(), ws.end());
}

nb::bytes ReactionToBinaryWithProps(const ChemicalReaction &self,
                                    unsigned int props) {
  std::string res;
  ReactionPickler::pickleReaction(self, res, props);
  return nb::bytes(res.c_str(), res.length());
}

nb::bytes ReactionToBinary(const ChemicalReaction &self) {
  return ReactionToBinaryWithProps(self,
                                   MolPickler::getDefaultPickleProperties());
}
std::string ReactionToBinaryString(const ChemicalReaction &self) {
  std::string res;
  ReactionPickler::pickleReaction(self, res,
                                  MolPickler::getDefaultPickleProperties());
  return res;
}

nb::tuple RunReactants(ChemicalReaction *self, nb::object reactants,
                       unsigned int maxProducts) {
  if (!self->isInitialized()) {
    NOGIL gil;
    self->initReactantMatchers();
  }
  MOL_SPTR_VECT reacts;
  unsigned int len1 = nb::len(reactants);
  reacts.resize(len1);
  for (unsigned int i = 0; i < len1; ++i) {
    nb::object mol_obj = reactants[i];
    if (mol_obj.is_none()) {
      throw nb::value_error("reaction called with None reactants");
    }
    reacts[i] = ROMOL_SPTR(&nb::cast<ROMol &>(mol_obj), [](ROMol *) {});
  }
  std::vector<MOL_SPTR_VECT> mols;
  {
    NOGIL gil;
    mols = self->runReactants(reacts, maxProducts);
  }
  nb::list res;
  for (const auto &mol_vec : mols) {
    nb::list inner;
    for (const auto &mol : mol_vec) {
      inner.append(toStd(mol));
    }
    res.append(nb::tuple(inner));
  }
  return nb::tuple(res);
}

nb::tuple RunReactant(ChemicalReaction *self, ROMol &reactant,
                      unsigned int reactionIdx) {
  ROMOL_SPTR react(&reactant, [](ROMol *) {});
  std::vector<MOL_SPTR_VECT> mols;
  {
    NOGIL gil;
    if (!self->isInitialized()) {
      self->initReactantMatchers();
    }
    mols = self->runReactant(react, reactionIdx);
  }
  nb::list res;
  for (const auto &mol_vec : mols) {
    nb::list inner;
    for (const auto &mol : mol_vec) {
      inner.append(toStd(mol));
    }
    res.append(nb::tuple(inner));
  }
  return nb::tuple(res);
}

bool RunReactantInPlace(ChemicalReaction *self, ROMol &reactant,
                        bool removeUnmatchedAtoms) {
  auto react = static_cast<RWMol *>(&reactant);
  bool res = false;
  {
    NOGIL gil;
    if (!self->isInitialized()) {
      self->initReactantMatchers();
    }
    res = self->runReactant(*react, removeUnmatchedAtoms);
  }
  return res;
}

std::tuple<unsigned int, unsigned int> ValidateReaction(
    const ChemicalReaction *self, bool silent) {
  unsigned int numWarn, numError;
  self->validate(numWarn, numError, silent);
  return std::make_tuple(numWarn, numError);
}

ROMol *GetProductTemplate(const ChemicalReaction *self, unsigned int which) {
  if (which >= self->getNumProductTemplates()) {
    throw nb::value_error("requested template index too high");
  }
  auto iter = self->beginProductTemplates();
  iter += which;
  return const_cast<ROMol *>(iter->get());
}

ROMol *GetReactantTemplate(const ChemicalReaction *self, unsigned int which) {
  if (which >= self->getNumReactantTemplates()) {
    throw nb::value_error("requested template index too high");
  }
  auto iter = self->beginReactantTemplates();
  iter += which;
  return const_cast<ROMol *>(iter->get());
}

ROMol *GetAgentTemplate(const ChemicalReaction *self, unsigned int which) {
  if (which >= self->getNumAgentTemplates()) {
    throw nb::value_error("requested template index too high");
  }
  auto iter = self->beginAgentTemplates();
  iter += which;
  return const_cast<ROMol *>(iter->get());
}

void RemoveUnmappedReactantTemplates(ChemicalReaction *self,
                                     double thresholdUnmappedAtoms,
                                     bool moveToAgentTemplates,
                                     nb::object targetList) {
  if (targetList.is_none()) {
    self->removeUnmappedReactantTemplates(thresholdUnmappedAtoms,
                                          moveToAgentTemplates);
  } else {
    MOL_SPTR_VECT tmp;
    self->removeUnmappedReactantTemplates(thresholdUnmappedAtoms,
                                          moveToAgentTemplates, &tmp);
    nb::list molList = nb::cast<nb::list>(targetList);
    for (auto &mol : tmp) {
      molList.append(toStd(mol));
    }
  }
}

void RemoveUnmappedProductTemplates(ChemicalReaction *self,
                                    double thresholdUnmappedAtoms,
                                    bool moveToAgentTemplates,
                                    nb::object targetList) {
  if (targetList.is_none()) {
    self->removeUnmappedProductTemplates(thresholdUnmappedAtoms,
                                         moveToAgentTemplates);
  } else {
    MOL_SPTR_VECT tmp;
    self->removeUnmappedProductTemplates(thresholdUnmappedAtoms,
                                         moveToAgentTemplates, &tmp);
    nb::list molList = nb::cast<nb::list>(targetList);
    for (auto &mol : tmp) {
      molList.append(toStd(mol));
    }
  }
}

void RemoveAgentTemplates(ChemicalReaction &self, nb::object targetList) {
  if (targetList.is_none()) {
    self.removeAgentTemplates();
  } else {
    MOL_SPTR_VECT tmp;
    self.removeAgentTemplates(&tmp);
    nb::list molList = nb::cast<nb::list>(targetList);
    for (auto &mol : tmp) {
      molList.append(toStd(mol));
    }
  }
}

void Compute2DCoordsForReaction(RDKit::ChemicalReaction &rxn,
                                double spacing = 2.0, bool updateProps = true,
                                bool canonOrient = false,
                                unsigned int nFlipsPerSample = 0,
                                unsigned int nSamples = 0, int sampleSeed = 0,
                                bool permuteDeg4Nodes = false,
                                double bondLength = -1) {
  double oBondLen = RDDepict::BOND_LEN;
  if (bondLength > 0) {
    RDDepict::BOND_LEN = bondLength;
  }
  RDDepict::compute2DCoordsForReaction(rxn, spacing, updateProps, canonOrient,
                                       nFlipsPerSample, nSamples, sampleSeed,
                                       permuteDeg4Nodes);
  if (bondLength > 0) {
    RDDepict::BOND_LEN = oBondLen;
  }
}

bool IsMoleculeReactantOfReaction(const ChemicalReaction &rxn,
                                  const ROMol &mol) {
  unsigned int which;
  return isMoleculeReactantOfReaction(rxn, mol, which);
}

bool IsMoleculeProductOfReaction(const ChemicalReaction &rxn,
                                 const ROMol &mol) {
  unsigned int which;
  return isMoleculeProductOfReaction(rxn, mol, which);
}

bool IsMoleculeAgentOfReaction(const ChemicalReaction &rxn, const ROMol &mol) {
  unsigned int which;
  return isMoleculeAgentOfReaction(rxn, mol, which);
}

ChemicalReaction *ReactionFromSmarts(const char *smarts, nb::dict replDict,
                                     bool useSmiles) {
  PRECONDITION(smarts, "null SMARTS string");
  std::map<std::string, std::string> replacements;
  for (auto [k, v] : replDict) {
    replacements[nb::cast<std::string>(k)] = nb::cast<std::string>(v);
  }
  return RxnSmartsToChemicalReaction(smarts, &replacements, useSmiles);
}

ChemicalReaction *ReactionFromSmiles(const char *smiles, nb::dict replDict) {
  return ReactionFromSmarts(smiles, replDict, true);
}

ChemicalReaction *ReactionFromMrvFile(const char *rxnFilename, bool sanitize,
                                      bool removeHs) {
  ChemicalReaction *newR = nullptr;
  try {
    newR = MrvFileToChemicalReaction(rxnFilename, sanitize, removeHs);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw nb::python_error();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return newR;
}

ChemicalReaction *ReactionFromMrvBlock(nb::object imolBlock, bool sanitize,
                                       bool removeHs) {
  std::istringstream inStream(pyObjectToString(imolBlock));
  ChemicalReaction *newR = nullptr;
  try {
    newR = MrvDataStreamToChemicalReaction(inStream, sanitize, removeHs);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return newR;
}

nb::tuple ReactionsFromCDXMLFile(const char *filename, bool sanitize,
                                 bool removeHs) {
  std::vector<std::unique_ptr<ChemicalReaction>> rxns;
  try {
    rxns = CDXMLFileToChemicalReactions(filename, sanitize, removeHs);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw nb::python_error();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  nb::list res;
  for (auto &rxn : rxns) {
    res.append(nb::cast(rxn.release(), nb::rv_policy::take_ownership));
  }
  return nb::tuple(res);
}

nb::tuple ReactionsFromCDXMLBlock(nb::object imolBlock, bool sanitize,
                                  bool removeHs) {
  std::istringstream inStream(pyObjectToString(imolBlock));
  std::vector<std::unique_ptr<ChemicalReaction>> rxns;
  try {
    rxns = CDXMLDataStreamToChemicalReactions(inStream, sanitize, removeHs);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  nb::list res;
  for (auto &rxn : rxns) {
    res.append(nb::cast(rxn.release(), nb::rv_policy::take_ownership));
  }
  return nb::tuple(res);
}

nb::tuple GetReactingAtoms(const ChemicalReaction &self, bool mappedAtomsOnly) {
  nb::list res;
  VECT_INT_VECT rAs = getReactingAtoms(self, mappedAtomsOnly);
  for (auto &rA : rAs) {
    nb::list inner;
    for (int atom : rA) {
      inner.append(atom);
    }
    res.append(nb::tuple(inner));
  }
  return nb::tuple(res);
}

nb::object AddRecursiveQueriesToReaction(ChemicalReaction &self,
                                         nb::dict queryDict,
                                         std::string propName, bool getLabels) {
  std::map<std::string, ROMOL_SPTR> queries;
  for (auto [k, v] : queryDict) {
    ROMol *m = nb::cast<ROMol *>(v);
    ROMOL_SPTR nm(new ROMol(*m));
    queries[nb::cast<std::string>(k)] = nm;
  }

  if (getLabels) {
    std::vector<std::vector<std::pair<unsigned int, std::string>>> labels;
    addRecursiveQueriesToReaction(self, queries, propName, &labels);

    nb::list reactantLabels;
    for (auto &label : labels) {
      nb::list tmpLabels;
      for (auto &j : label) {
        nb::list tmpPair;
        tmpPair.append(j.first);
        tmpPair.append(j.second);
        tmpLabels.append(nb::tuple(tmpPair));
      }
      reactantLabels.append(nb::tuple(tmpLabels));
    }
    return nb::tuple(reactantLabels);
  } else {
    addRecursiveQueriesToReaction(self, queries, propName);
    return nb::none();
  }
}

nb::object PreprocessReaction(ChemicalReaction &reaction, nb::dict queryDict,
                              std::string propName) {
  std::map<std::string, ROMOL_SPTR> queries;
  unsigned int size = nb::len(queryDict);
  if (!size) {
    const bool normalized = true;
    queries = GetFlattenedFunctionalGroupHierarchy(normalized);
  } else {
    for (auto [k, v] : queryDict) {
      ROMol *m = nb::cast<ROMol *>(v);
      ROMOL_SPTR nm(new ROMol(*m));
      queries[nb::cast<std::string>(k)] = nm;
    }
  }

  unsigned int nReactants = reaction.getNumReactantTemplates();
  unsigned int nProducts = reaction.getNumProductTemplates();
  unsigned int nWarn, nError;
  reaction.validate(nWarn, nError);
  std::vector<std::vector<std::pair<unsigned int, std::string>>> labels;

  if (!nError) {
    preprocessReaction(reaction, nWarn, nError, labels, queries, propName);
  }

  nb::list reactantLabels;
  for (auto &label : labels) {
    nb::list tmpLabels;
    for (auto &j : label) {
      nb::list tmpPair;
      tmpPair.append(j.first);
      tmpPair.append(j.second);
      tmpLabels.append(nb::tuple(tmpPair));
    }
    reactantLabels.append(nb::tuple(tmpLabels));
  }
  return nb::make_tuple(nWarn, nError, nReactants, nProducts,
                        nb::tuple(reactantLabels));
}

RxnOps::SanitizeRxnFlags sanitizeReaction(
    ChemicalReaction &rxn, uint64_t sanitizeOps,
    const MolOps::AdjustQueryParameters &params, bool catchErrors) {
  unsigned int operationsThatFailed = 0;
  try {
    RxnOps::sanitizeRxn(rxn, operationsThatFailed, sanitizeOps, params);
  } catch (...) {
    if (!catchErrors) {
      throw;
    }
  }
  return static_cast<RxnOps::SanitizeRxnFlags>(operationsThatFailed);
}

nb::bytes addReactionToPNGStringHelper(const ChemicalReaction &rxn,
                                       nb::bytes png, bool includePkl,
                                       bool includeSmiles, bool includeSmarts,
                                       bool includeRxn) {
  std::string cstr(static_cast<const char *>(png.data()), png.size());
  auto res = addChemicalReactionToPNGString(
      rxn, cstr, includePkl, includeSmiles, includeSmarts, includeRxn);
  return nb::bytes(res.c_str(), res.length());
}

nb::bytes addReactionToPNGFileHelper(const ChemicalReaction &rxn,
                                     nb::object fname, bool includePkl,
                                     bool includeSmiles, bool includeSmarts,
                                     bool includeRxn) {
  std::string cstr = nb::cast<std::string>(fname);
  auto res = addChemicalReactionToPNGFile(rxn, cstr, includePkl, includeSmiles,
                                          includeSmarts, includeRxn);
  return nb::bytes(res.c_str(), res.length());
}

}  // namespace RDKit

NB_MODULE(rdChemReactions, m) {
  m.doc() =
      "Module containing classes and functions for working with chemical "
      "reactions.";

  nb::exception<RDKit::ChemicalReactionParserException>(
      m, "ChemicalReactionParserException", PyExc_ValueError);
  nb::exception<RDKit::ChemicalReactionException>(
      m, "ChemicalReactionException", PyExc_ValueError);

  nb::enum_<RDKit::FingerprintType>(m, "FingerprintType")
      .value("AtomPairFP", RDKit::FingerprintType::AtomPairFP)
      .value("TopologicalTorsion", RDKit::FingerprintType::TopologicalTorsionFP)
      .value("MorganFP", RDKit::FingerprintType::MorganFP)
      .value("RDKitFP", RDKit::FingerprintType::RDKitFP)
      .value("PatternFP", RDKit::FingerprintType::PatternFP);

  nb::class_<RDKit::ReactionFingerprintParams>(
      m, "ReactionFingerprintParams",
      R"DOC(A class for storing parameters to manipulate the calculation of
fingerprints of chemical reactions.)DOC")
      .def(nb::init<>(), "Constructor, takes no arguments")
      .def(nb::init<bool, double, unsigned int, int, unsigned int,
                    RDKit::FingerprintType>(),
           "includeAgents"_a, "bitRatioAgents"_a, "nonAgentWeight"_a,
           "agentWeight"_a, "fpSize"_a, "fpType"_a)
      .def_rw("fpSize", &RDKit::ReactionFingerprintParams::fpSize)
      .def_rw("fpType", &RDKit::ReactionFingerprintParams::fpType)
      .def_rw("bitRatioAgents",
              &RDKit::ReactionFingerprintParams::bitRatioAgents)
      .def_rw("nonAgentWeight",
              &RDKit::ReactionFingerprintParams::nonAgentWeight)
      .def_rw("agentWeight", &RDKit::ReactionFingerprintParams::agentWeight)
      .def_rw("includeAgents", &RDKit::ReactionFingerprintParams::includeAgents)
      .def("__setattr__", &safeSetattr);

  nb::class_<RDKit::ChemicalReaction>(
      m, "ChemicalReaction",
      R"DOC(A class for storing and applying chemical reactions.

Sample Usage:
  >>> from rdkit import Chem
  >>> from rdkit.Chem import rdChemReactions
  >>> rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')
  >>> reacts = (Chem.MolFromSmiles('C(=O)O'),Chem.MolFromSmiles('CNC'))
  >>> products = rxn.RunReactants(reacts)
  >>> len(products)
  1
  >>> len(products[0])
  1
  >>> Chem.MolToSmiles(products[0][0])
  'CN(C)C=O'
)DOC")
      .def(nb::new_([]() { return new RDKit::ChemicalReaction(); }),
           "Constructor, takes no arguments")
      .def(nb::new_([](nb::bytes b) {
             return new RDKit::ChemicalReaction(
                 std::string(static_cast<const char *>(b.data()), b.size()));
           }),
           "binStr"_a)
      .def(nb::new_([](const RDKit::ChemicalReaction &other) {
             return new RDKit::ChemicalReaction(other);
           }),
           "other"_a)
      .def("GetNumReactantTemplates",
           &RDKit::ChemicalReaction::getNumReactantTemplates,
           "returns the number of reactants this reaction expects")
      .def("GetNumProductTemplates",
           &RDKit::ChemicalReaction::getNumProductTemplates,
           "returns the number of products this reaction generates")
      .def("GetNumAgentTemplates",
           &RDKit::ChemicalReaction::getNumAgentTemplates,
           "returns the number of agents this reaction expects")
      .def(
          "AddReactantTemplate",
          [](RDKit::ChemicalReaction &rxn, RDKit::ROMol &mol) -> unsigned int {
            return rxn.addReactantTemplate(
                RDKit::ROMOL_SPTR(&mol, [](RDKit::ROMol *) {}));
          },
          nb::keep_alive<1, 2>(), "mol"_a,
          "adds a reactant (a Molecule) to the reaction")
      .def(
          "AddProductTemplate",
          [](RDKit::ChemicalReaction &rxn, RDKit::ROMol &mol) -> unsigned int {
            return rxn.addProductTemplate(
                RDKit::ROMOL_SPTR(&mol, [](RDKit::ROMol *) {}));
          },
          nb::keep_alive<1, 2>(), "mol"_a, "adds a product (a Molecule)")
      .def(
          "AddAgentTemplate",
          [](RDKit::ChemicalReaction &rxn, RDKit::ROMol &mol) -> unsigned int {
            return rxn.addAgentTemplate(
                RDKit::ROMOL_SPTR(&mol, [](RDKit::ROMol *) {}));
          },
          nb::keep_alive<1, 2>(), "mol"_a, "adds a agent (a Molecule)")
      .def("RemoveUnmappedReactantTemplates",
           RDKit::RemoveUnmappedReactantTemplates,
           "thresholdUnmappedAtoms"_a = 0.2, "moveToAgentTemplates"_a = true,
           "targetList"_a = nb::none(),
           R"DOC(Removes molecules with an atom mapping ratio below
thresholdUnmappedAtoms from reactant templates to the agent
templates or to a given targetList)DOC")
      .def("RemoveUnmappedProductTemplates",
           RDKit::RemoveUnmappedProductTemplates,
           "thresholdUnmappedAtoms"_a = 0.2, "moveToAgentTemplates"_a = true,
           "targetList"_a = nb::none(),
           R"DOC(Removes molecules with an atom mapping ratio below
thresholdUnmappedAtoms from product templates to the agent
templates or to a given targetList)DOC")
      .def(
          "RemoveAgentTemplates", RDKit::RemoveAgentTemplates,
          "targetList"_a = nb::none(),
          R"DOC(Removes agents from reaction. If targetList is provide the agents
will be transferred to that list.)DOC")
      .def(
          "RunReactants", RDKit::RunReactants, "reactants"_a,
          "maxProducts"_a = 1000,
          R"DOC(apply the reaction to a sequence of reactant molecules and return
the products as a tuple of tuples.  If maxProducts is not zero,
 stop the reaction when maxProducts have been generated [default=1000])DOC")
      .def("RunReactant", RDKit::RunReactant, "reactant"_a, "reactionIdx"_a,
           "apply the reaction to a single reactant")
      .def("RunReactantInPlace", RDKit::RunReactantInPlace, "reactant"_a,
           "removeUnmatchedAtoms"_a = true,
           R"DOC(apply the reaction to a single reactant in place. The reactant
itself is modified. This can only be used for single reactant -
single product reactions.)DOC")
      .def("Initialize", &RDKit::ChemicalReaction::initReactantMatchers,
           "silent"_a = false,
           "initializes the reaction so that it can be used")
      .def("IsInitialized", &RDKit::ChemicalReaction::isInitialized,
           "checks if the reaction is ready for use")
      .def("Validate", RDKit::ValidateReaction, "silent"_a = false,
           "checks the reaction for potential problems, returns "
           "(numWarnings,numErrors)")
      .def("GetProductTemplate", RDKit::GetProductTemplate, "which"_a,
           nb::rv_policy::reference_internal,
           "returns one of our product templates")
      .def("GetReactantTemplate", RDKit::GetReactantTemplate, "which"_a,
           nb::rv_policy::reference_internal,
           "returns one of our reactant templates")
      .def("GetAgentTemplate", RDKit::GetAgentTemplate, "which"_a,
           nb::rv_policy::reference_internal,
           "returns one of our agent templates")
      .def("_setImplicitPropertiesFlag",
           &RDKit::ChemicalReaction::setImplicitPropertiesFlag, "val"_a,
           "EXPERT USER: indicates that the reaction can have implicit "
           "properties")
      .def("_getImplicitPropertiesFlag",
           &RDKit::ChemicalReaction::getImplicitPropertiesFlag,
           "EXPERT USER: returns whether or not the reaction can have implicit "
           "properties")
      .def("ToBinary", RDKit::ReactionToBinary,
           "Returns a binary string representation of the reaction.")
      .def(
          "ToBinary",
          [](const RDKit::ChemicalReaction &self, nb::object propertyFlags) {
            unsigned int val;
            if (nb::isinstance<nb::int_>(propertyFlags)) {
              val = nb::cast<unsigned int>(propertyFlags);
            } else {
              val = nb::cast<unsigned int>(propertyFlags.attr("value"));
            }
            return RDKit::ReactionToBinaryWithProps(self, val);
          },
          "propertyFlags"_a,
          "Returns a binary string representation of the reaction.")
      .def("IsMoleculeReactant", RDKit::IsMoleculeReactantOfReaction, "mol"_a,
           "returns whether or not the molecule has a substructure match to "
           "one of the reactants.")
      .def("IsMoleculeProduct", RDKit::IsMoleculeProductOfReaction, "mol"_a,
           "returns whether or not the molecule has a substructure match to "
           "one of the products.")
      .def("IsMoleculeAgent", RDKit::IsMoleculeAgentOfReaction, "mol"_a,
           "returns whether or not the molecule has a substructure match to "
           "one of the agents.")
      .def("GetReactingAtoms", RDKit::GetReactingAtoms,
           "mappedAtomsOnly"_a = false,
           "returns a sequence of sequences with the atoms that change in the "
           "reaction")
      .def("AddRecursiveQueriesToReaction",
           RDKit::AddRecursiveQueriesToReaction, "queries"_a = nb::dict(),
           "propName"_a = "molFileValue", "getLabels"_a = false,
           "adds recursive queries and returns reactant labels")
      .def(
          "GetReactants",
          [](const RDKit::ChemicalReaction &rxn) {
            nb::list res;
            for (const auto &mol : rxn.getReactants()) {
              res.append(toStd(mol));
            }
            return res;
          },
          "get the reactant templates")
      .def(
          "GetProducts",
          [](const RDKit::ChemicalReaction &rxn) {
            nb::list res;
            for (const auto &mol : rxn.getProducts()) {
              res.append(toStd(mol));
            }
            return res;
          },
          "get the product templates")
      .def(
          "GetAgents",
          [](const RDKit::ChemicalReaction &rxn) {
            nb::list res;
            for (const auto &mol : rxn.getAgents()) {
              res.append(toStd(mol));
            }
            return res;
          },
          "get the agent templates")
      .def(
          "GetSubstructParams",
          [](RDKit::ChemicalReaction &rxn)
              -> RDKit::SubstructMatchParameters & {
            return rxn.getSubstructParams();
          },
          nb::rv_policy::reference_internal,
          "get the parameter object controlling the substructure matching")
      // properties
      .def("SetProp", RDKit::MolSetProp<RDKit::ChemicalReaction, std::string>,
           "key"_a, "val"_a, "computed"_a = false,
           R"DOC(Sets a molecular property

  ARGUMENTS:
    - key: the name of the property to be set (a string).
    - value: the property value (a string).
    - computed: (optional) marks the property as being computed.
                Defaults to False.)DOC")
      .def("SetDoubleProp", RDKit::MolSetProp<RDKit::ChemicalReaction, double>,
           "key"_a, "val"_a, "computed"_a = false,
           R"DOC(Sets a double valued molecular property

  ARGUMENTS:
    - key: the name of the property to be set (a string).
    - value: the property value as a double.
    - computed: (optional) marks the property as being computed.
                Defaults to 0.)DOC")
      .def("SetIntProp", RDKit::MolSetProp<RDKit::ChemicalReaction, int>,
           "key"_a, "val"_a, "computed"_a = false,
           R"DOC(Sets an integer valued molecular property

  ARGUMENTS:
    - key: the name of the property to be set (an unsigned number).
    - value: the property value as an integer.
    - computed: (optional) marks the property as being computed.
                Defaults to False.)DOC")
      .def("SetUnsignedProp",
           RDKit::MolSetProp<RDKit::ChemicalReaction, unsigned int>, "key"_a,
           "val"_a, "computed"_a = false,
           R"DOC(Sets an unsigned integer valued molecular property

  ARGUMENTS:
    - key: the name of the property to be set (a string).
    - value: the property value as an unsigned integer.
    - computed: (optional) marks the property as being computed.
                Defaults to False.)DOC")
      .def("SetBoolProp", RDKit::MolSetProp<RDKit::ChemicalReaction, bool>,
           "key"_a, "val"_a, "computed"_a = false,
           R"DOC(Sets a boolean valued molecular property

  ARGUMENTS:
    - key: the name of the property to be set (a string).
    - value: the property value as a bool.
    - computed: (optional) marks the property as being computed.
                Defaults to False.)DOC")
      .def(
          "HasProp", RDKit::MolHasProp<RDKit::ChemicalReaction>, "key"_a,
          R"DOC(Queries a molecule to see if a particular property has been assigned.

  ARGUMENTS:
    - key: the name of the property to check for (a string).)DOC")
      .def("GetProp", RDKit::GetProp<RDKit::ChemicalReaction, std::string>,
           "key"_a,
           R"DOC(Returns the value of the property.

  ARGUMENTS:
    - key: the name of the property to return (a string).

  RETURNS: a string

  NOTE:
    - If the property has not been set, a KeyError exception will be raised.)DOC")
      .def("GetDoubleProp", RDKit::GetProp<RDKit::ChemicalReaction, double>,
           "key"_a,
           R"DOC(Returns the double value of the property if possible.

  ARGUMENTS:
    - key: the name of the property to return (a string).

  RETURNS: a double

  NOTE:
    - If the property has not been set, a KeyError exception will be raised.)DOC")
      .def("GetIntProp", RDKit::GetProp<RDKit::ChemicalReaction, int>, "key"_a,
           R"DOC(Returns the integer value of the property if possible.

  ARGUMENTS:
    - key: the name of the property to return (a string).

  RETURNS: an integer

  NOTE:
    - If the property has not been set, a KeyError exception will be raised.)DOC")
      .def("GetUnsignedProp",
           RDKit::GetProp<RDKit::ChemicalReaction, unsigned int>, "key"_a,
           R"DOC(Returns the unsigned int value of the property if possible.

  ARGUMENTS:
    - key: the name of the property to return (a string).

  RETURNS: an unsigned integer

  NOTE:
    - If the property has not been set, a KeyError exception will be raised.)DOC")
      .def("GetBoolProp", RDKit::GetProp<RDKit::ChemicalReaction, bool>,
           "key"_a,
           R"DOC(Returns the Bool value of the property if possible.

  ARGUMENTS:
    - key: the name of the property to return (a string).

  RETURNS: a bool

  NOTE:
    - If the property has not been set, a KeyError exception will be raised.)DOC")
      .def("ClearProp", RDKit::MolClearProp<RDKit::ChemicalReaction>, "key"_a,
           R"DOC(Removes a property from the reaction.

  ARGUMENTS:
    - key: the name of the property to clear (a string).)DOC")
      .def("ClearComputedProps",
           RDKit::MolClearComputedProps<RDKit::ChemicalReaction>,
           "Removes all computed properties from the reaction.")
      .def("GetPropNames", &RDKit::ChemicalReaction::getPropList,
           "includePrivate"_a = false, "includeComputed"_a = false,
           R"DOC(Returns a tuple with all property names for this reaction.

  ARGUMENTS:
    - includePrivate: (optional) toggles inclusion of private properties in the result set.
                      Defaults to 0.
    - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                      Defaults to 0.

  RETURNS: a tuple of strings)DOC")
      .def("GetPropsAsDict", RDKit::GetPropsAsDict<RDKit::ChemicalReaction>,
           "includePrivate"_a = false, "includeComputed"_a = false,
           "autoConvertStrings"_a = true,
           R"DOC(Returns a dictionary populated with the reaction's properties.
 n.b. Some properties are not able to be converted to python types.

  ARGUMENTS:
    - includePrivate: (optional) toggles inclusion of private properties in the result set.
                      Defaults to False.
    - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                      Defaults to False.

  RETURNS: a dictionary)DOC")
      // pickle support
      .def("__setstate__", setObjectState<RDKit::ChemicalReaction>)
      .def("__getstate__", getObjectState<RDKit::ChemicalReaction,
                                          RDKit::ReactionToBinaryString>);

  m.def("ReactionFromSmarts", RDKit::ReactionFromSmarts, "SMARTS"_a,
        "replacements"_a = nb::dict(), "useSmiles"_a = false,
        R"DOC(construct a ChemicalReaction from a reaction SMARTS string.
see the documentation for rdkit.Chem.MolFromSmiles for an explanation
of the replacements argument.)DOC",
        nb::rv_policy::take_ownership);

  m.def("ReactionToSmarts",
        (std::string(*)(
            const RDKit::ChemicalReaction &))RDKit::ChemicalReactionToRxnSmarts,
        "reaction"_a,
        "construct a reaction SMARTS string for a ChemicalReaction");

  m.def("ReactionFromSmiles", RDKit::ReactionFromSmiles, "SMILES"_a,
        "replacements"_a = nb::dict(),
        R"DOC(construct a ChemicalReaction from a reaction SMILES string.
see the documentation for rdkit.Chem.MolFromSmiles for an explanation
of the replacements argument.)DOC",
        nb::rv_policy::take_ownership);

  m.def("ReactionToSmiles",
        (std::string(*)(const RDKit::ChemicalReaction &,
                        bool))RDKit::ChemicalReactionToRxnSmiles,
        "reaction"_a, "canonical"_a = true,
        "construct a reaction SMILES string for a ChemicalReaction");

  m.def("ReactionToSmarts",
        (std::string(*)(const RDKit::ChemicalReaction &,
                        const RDKit::SmilesWriteParams &))
            RDKit::ChemicalReactionToRxnSmarts,
        "reaction"_a, "params"_a,
        "construct a reaction SMARTS string for a ChemicalReaction");

  m.def("ReactionToSmiles",
        (std::string(*)(const RDKit::ChemicalReaction &,
                        const RDKit::SmilesWriteParams &))
            RDKit::ChemicalReactionToRxnSmiles,
        "reaction"_a, "params"_a,
        "construct a reaction SMILES string for a ChemicalReaction");

  m.def("ReactionToCXSmarts",
        (std::string(*)(const RDKit::ChemicalReaction &))
            RDKit::ChemicalReactionToCXRxnSmarts,
        "reaction"_a,
        "construct a reaction SMARTS string for a ChemicalReaction");

  m.def("ReactionToCXSmiles",
        (std::string(*)(const RDKit::ChemicalReaction &,
                        bool))RDKit::ChemicalReactionToCXRxnSmiles,
        "reaction"_a, "canonical"_a = true,
        "construct a reaction SMILES string for a ChemicalReaction");

  m.def("ReactionToCXSmarts",
        (std::string(*)(const RDKit::ChemicalReaction &,
                        const RDKit::SmilesWriteParams &,
                        std::uint32_t))RDKit::ChemicalReactionToCXRxnSmarts,
        "reaction"_a, "params"_a,
        "flags"_a = RDKit::SmilesWrite::CXSmilesFields::CX_ALL,
        "construct a reaction CXSMARTS string for a ChemicalReaction");

  m.def("ReactionToCXSmiles",
        (std::string(*)(const RDKit::ChemicalReaction &,
                        const RDKit::SmilesWriteParams &,
                        std::uint32_t))RDKit::ChemicalReactionToCXRxnSmiles,
        "reaction"_a, "params"_a,
        "flags"_a = RDKit::SmilesWrite::CXSmilesFields::CX_ALL,
        "construct a reaction CXSMILES string for a ChemicalReaction");

  m.def("ReactionFromRxnFile", RDKit::RxnFileToChemicalReaction, "filename"_a,
        "sanitize"_a = false, "removeHs"_a = false, "strictParsing"_a = true,
        "construct a ChemicalReaction from an MDL rxn file",
        nb::rv_policy::take_ownership);

  m.def("ReactionFromRxnBlock", RDKit::RxnBlockToChemicalReaction, "rxnblock"_a,
        "sanitize"_a = false, "removeHs"_a = false, "strictParsing"_a = true,
        "construct a ChemicalReaction from a string in MDL rxn format",
        nb::rv_policy::take_ownership);

  m.def("ReactionFromMrvFile", RDKit::ReactionFromMrvFile, "filename"_a,
        "sanitize"_a = false, "removeHs"_a = false,
        "construct a ChemicalReaction from an Marvin (mrv) rxn file",
        nb::rv_policy::take_ownership);

  m.def("ReactionFromMrvBlock", RDKit::ReactionFromMrvBlock, "rxnblock"_a,
        "sanitize"_a = false, "removeHs"_a = false,
        "construct a ChemicalReaction from a string in Marvin (mrv) format",
        nb::rv_policy::take_ownership);

  m.def("MrvFileIsReaction", RDKit::MrvFileIsReaction, "filename"_a,
        "returns whether or not an MRV file contains reaction data");

  m.def("MrvBlockIsReaction", RDKit::MrvBlockIsReaction, "mrvData"_a,
        "returns whether or not an MRV block contains reaction data");

  m.def("ReactionsFromCDXMLFile", RDKit::ReactionsFromCDXMLFile, "filename"_a,
        "sanitize"_a = false, "removeHs"_a = false,
        "construct a tuple of ChemicalReactions from a CDXML rxn file");

  m.def("ReactionsFromCDXMLBlock", RDKit::ReactionsFromCDXMLBlock, "rxnblock"_a,
        "sanitize"_a = false, "removeHs"_a = false,
        "construct a tuple of ChemicalReactions from a string in CDXML format");

  m.def("ReactionToRxnBlock", RDKit::ChemicalReactionToRxnBlock, "reaction"_a,
        "separateAgents"_a = false, "forceV3000"_a = false,
        "construct a string in MDL rxn format for a ChemicalReaction");

  m.def("ReactionToMrvBlock", RDKit::ChemicalReactionToMrvBlock, "reaction"_a,
        "prettyPrint"_a = false,
        "construct a string in Marvin (MRV) rxn format for a ChemicalReaction");

  m.def("ReactionToMrvFile", RDKit::ChemicalReactionToMrvFile, "reaction"_a,
        "filename"_a, "prettyPrint"_a = false,
        "write a Marvin (MRV) rxn file for a ChemicalReaction");

  m.def("ReactionToV3KRxnBlock", RDKit::ChemicalReactionToV3KRxnBlock,
        "reaction"_a, "separateAgents"_a = false,
        "construct a string in MDL v3000 rxn format for a ChemicalReaction");

#ifdef RDK_USE_BOOST_IOSTREAMS
  m.def("ReactionFromPNGFile", RDKit::PNGFileToChemicalReaction, "fname"_a,
        "construct a ChemicalReaction from metadata in a PNG file",
        nb::rv_policy::take_ownership);
  m.def(
      "ReactionFromPNGString",
      [](nb::bytes data) {
        return RDKit::PNGStringToChemicalReaction(
            std::string(static_cast<const char *>(data.data()), data.size()));
      },
      "data"_a, "construct a ChemicalReaction from an string with PNG data",
      nb::rv_policy::take_ownership);
  m.def(
      "ReactionMetadataToPNGFile", RDKit::addReactionToPNGFileHelper, "mol"_a,
      "filename"_a, "includePkl"_a = true, "includeSmiles"_a = true,
      "includeSmarts"_a = false, "includeMol"_a = false,
      R"DOC(Reads the contents of a PNG file and adds metadata about a reaction to
it. The modified file contents are returned.)DOC");
  m.def("ReactionMetadataToPNGString", RDKit::addReactionToPNGStringHelper,
        "mol"_a, "pngdata"_a, "includePkl"_a = true, "includeSmiles"_a = true,
        "includeSmarts"_a = false, "includeRxn"_a = false,
        R"DOC(Adds metadata about a reaction to the PNG string passed in.
The modified string is returned.)DOC");
#endif

  m.def("ReactionFromMolecule", RDKit::RxnMolToChemicalReaction, "mol"_a,
        "construct a ChemicalReaction from an molecule if the RXN role "
        "property of the molecule is set",
        nb::rv_policy::take_ownership);

  m.def(
      "ReactionToMolecule", RDKit::ChemicalReactionToRxnMol, "reaction"_a,
      "construct a molecule for a ChemicalReaction with RXN role property set",
      nb::rv_policy::take_ownership);

  m.def("Compute2DCoordsForReaction", RDKit::Compute2DCoordsForReaction,
        "reaction"_a, "spacing"_a = 1.0, "updateProps"_a = true,
        "canonOrient"_a = true, "nFlipsPerSample"_a = 0, "nSample"_a = 0,
        "sampleSeed"_a = 0, "permuteDeg4Nodes"_a = false, "bondLength"_a = -1.0,
        R"DOC(Compute 2D coordinates for a reaction.
  ARGUMENTS:
     - reaction - the reaction of interest
     - spacing - the amount of space left between components of the reaction
     - canonOrient - orient the reactants and products in a canonical way
     - updateProps - if set, properties such as conjugation and
        hybridization will be calculated for the reactant and product
        templates before generating coordinates. This should result in
        better depictions, but can lead to errors in some cases.
     - nFlipsPerSample - number of rotatable bonds that are
                flipped at random at a time.
     - nSample - Number of random samplings of rotatable bonds.
     - sampleSeed - seed for the random sampling process.
     - permuteDeg4Nodes - allow permutation of bonds at a degree 4
                 node during the sampling process
     - bondLength - change the default bond length for depiction
)DOC");

  m.def("CreateDifferenceFingerprintForReaction",
        RDKit::DifferenceFingerprintChemReaction, "reaction"_a,
        "ReactionFingerPrintParams"_a = RDKit::DefaultDifferenceFPParams,
        R"DOC(construct a difference fingerprint for a ChemicalReaction by
subtracting the reactant fingerprint from the product fingerprint)DOC",
        nb::rv_policy::take_ownership);

  m.def("CreateStructuralFingerprintForReaction",
        RDKit::StructuralFingerprintChemReaction, "reaction"_a,
        "ReactionFingerPrintParams"_a = RDKit::DefaultStructuralFPParams,
        R"DOC(construct a structural fingerprint for a ChemicalReaction by
concatenating the reactant fingerprint and the product fingerprint)DOC",
        nb::rv_policy::take_ownership);

  m.def("IsReactionTemplateMoleculeAgent",
        RDKit::isReactionTemplateMoleculeAgent, "molecule"_a,
        "agentThreshold"_a,
        "tests if a molecule can be classified as an agent depending on "
        "the ratio of mapped atoms and a give threshold");

  m.def("HasReactionAtomMapping", RDKit::hasReactionAtomMapping, "rxn"_a,
        "tests if a reaction obtains any atom mapping");

  m.def("HasReactionSubstructMatch", RDKit::hasReactionSubstructMatch,
        "reaction"_a, "queryReaction"_a, "includeAgents"_a = false,
        "tests if the queryReaction is a substructure of a reaction");

  m.def("HasAgentTemplateSubstructMatch", RDKit::hasAgentTemplateSubstructMatch,
        "reaction"_a, "queryReaction"_a,
        "tests if the agents of a queryReaction are the same as those of "
        "a reaction");

  m.def("HasProductTemplateSubstructMatch",
        RDKit::hasProductTemplateSubstructMatch, "reaction"_a,
        "queryReaction"_a,
        "tests if the products of a queryReaction are substructures of "
        "the products of a reaction");

  m.def("HasReactantTemplateSubstructMatch",
        RDKit::hasReactantTemplateSubstructMatch, "reaction"_a,
        "queryReaction"_a,
        "tests if the reactants of a queryReaction are substructures of "
        "the reactants of a reaction");

  m.def(
      "UpdateProductsStereochemistry", RDKit::updateProductsStereochem,
      "reaction"_a,
      R"DOC(Caution: This is an expert-user function which will change a property (molInversionFlag) of your products.
          This function is called by default using the RXN or SMARTS parser for reactions and should really only be called if reactions have been constructed some other way.
          The function updates the stereochemistry of the product by considering 4 different cases: inversion, retention, removal, and introduction)DOC");

  m.def(
      "ReduceProductToSideChains",
      [](RDKit::ROMol &product, bool addDummyAtoms) {
        RDKit::ROMOL_SPTR sptr(&product, [](RDKit::ROMol *) {});
        return RDKit::reduceProductToSideChains(sptr, addDummyAtoms);
      },
      "product"_a, "addDummyAtoms"_a = true,
      R"DOC(reduce the product of a reaction to the side chains added by the reaction.
              The output is a molecule with attached wildcards indicating where the product was attached.
              The dummy atom has the same reaction-map number as the product atom (if available).)DOC",
      nb::rv_policy::take_ownership);

  m.def("RemoveMappingNumbersFromReactions",
        RDKit::removeMappingNumbersFromReactions, "reaction"_a,
        "Removes the mapping numbers from the molecules of a reaction");

  m.def("PreprocessReaction", RDKit::PreprocessReaction, "reaction"_a,
        "queries"_a = nb::dict(), "propName"_a = "molFileValue",
        R"DOC(A function for preprocessing reactions with more specific queries.
Queries are indicated by labels on atoms (molFileAlias property by default)
When these labels are found, more specific queries are placed on the atoms.
By default, the available quieries come from
  FilterCatalog.GetFlattenedFunctionalGroupHierarchy(True)n
Sample Usage:
  >>> from rdkit import Chem, RDConfig
  >>> from rdkit.Chem import MolFromSmiles, AllChem
  >>> from rdkit.Chem.rdChemReactions import PreprocessReaction
  >>> import os
  >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','boronic1.rxn')
  >>> rxn = AllChem.ReactionFromRxnFile(testFile)
  >>> rxn.Initialize()
  >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
  >>> nWarn
  0
  >>> nError
  0
  >>> nReacts
  2
  >>> nProds
  1
  >>> reactantLabels
  (((0, 'halogen.bromine.aromatic'),), ((1, 'boronicacid'),))

If there are functional group labels in the input reaction (via atoms with molFileValue properties),
the corresponding atoms will have queries added to them so that they only match such things. We can
see this here:
  >>> rxn = AllChem.ReactionFromRxnFile(testFile)
  >>> rxn.Initialize()
  >>> r1 = rxn.GetReactantTemplate(0)
  >>> m1 = Chem.MolFromSmiles('CCBr')
  >>> m2 = Chem.MolFromSmiles('c1ccccc1Br')

These both match because the reaction file itself just has R1-Br:
  >>> m1.HasSubstructMatch(r1)
  True
  >>> m2.HasSubstructMatch(r1)
  True

After preprocessing, we only match the aromatic Br:
  >>> d = PreprocessReaction(rxn)
  >>> m1.HasSubstructMatch(r1)
  False
  >>> m2.HasSubstructMatch(r1)
  True

We also support or queries in the values field (separated by commas):
  >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','azide_reaction.rxn')
  >>> rxn = AllChem.ReactionFromRxnFile(testFile)
  >>> rxn.Initialize()
  >>> reactantLabels = PreprocessReaction(rxn)[-1]
  >>> reactantLabels
  (((1, 'azide'),), ((1, 'carboxylicacid,acidchloride'),))
  >>> m1 = Chem.MolFromSmiles('CC(=O)O')
  >>> m2 = Chem.MolFromSmiles('CC(=O)Cl')
  >>> m3 = Chem.MolFromSmiles('CC(=O)N')
  >>> r2 = rxn.GetReactantTemplate(1)
  >>> m1.HasSubstructMatch(r2)
  True
  >>> m2.HasSubstructMatch(r2)
  True
  >>> m3.HasSubstructMatch(r2)
  False

unrecognized final group types are returned as None:
  >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','bad_value1.rxn')
  >>> rxn = AllChem.ReactionFromRxnFile(testFile)
  >>> rxn.Initialize()
  >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
  Traceback (most recent call last):
    ...
  KeyError: 'boromicacid'

One unrecognized group type in a comma-separated list makes the whole thing fail:
  >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','bad_value2.rxn')
  >>> rxn = AllChem.ReactionFromRxnFile(testFile)
  >>> rxn.Initialize()
  >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
  Traceback (most recent call last):
    ...
  KeyError: 'carboxylicacid,acidchlroide'
  >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','bad_value3.rxn')
  >>> rxn = AllChem.ReactionFromRxnFile(testFile)
  >>> rxn.Initialize()
  >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
  Traceback (most recent call last):
    ...
  KeyError: 'carboxyliccaid,acidchloride'
  >>> rxn = rdChemReactions.ChemicalReaction()
  >>> rxn.Initialize()
  >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
  >>> reactantLabels
  ()
  >>> reactantLabels == ()
  True
)DOC");

  nb::enum_<RDKit::RxnOps::SanitizeRxnFlags>(m, "SanitizeFlags")
      .value("SANITIZE_NONE", RDKit::RxnOps::SANITIZE_NONE)
      .value("SANITIZE_ATOM_MAPS", RDKit::RxnOps::SANITIZE_ATOM_MAPS)
      .value("SANITIZE_RGROUP_NAMES", RDKit::RxnOps::SANITIZE_RGROUP_NAMES)
      .value("SANITIZE_ADJUST_REACTANTS",
             RDKit::RxnOps::SANITIZE_ADJUST_REACTANTS)
      .value("SANITIZE_MERGEHS", RDKit::RxnOps::SANITIZE_MERGEHS)
      .value("SANITIZE_ALL", RDKit::RxnOps::SANITIZE_ALL)
      .export_values();

  m.def("GetDefaultAdjustParams", RDKit::RxnOps::DefaultRxnAdjustParams,
        "Returns the default adjustment parameters for reactant templates");

  m.def("GetChemDrawRxnAdjustParams", RDKit::RxnOps::ChemDrawRxnAdjustParams,
        "(deprecated, see MatchOnlyAtRgroupsAdjustParams)\n\tReturns the "
        "chemdraw style adjustment parameters for reactant templates");

  m.def(
      "MatchOnlyAtRgroupsAdjustParams",
      RDKit::RxnOps::MatchOnlyAtRgroupsAdjustParams,
      "Only match at the specified rgroup locations in the reactant templates");

  m.def(
      "SanitizeRxn",
      [](RDKit::ChemicalReaction &rxn, unsigned int sanitizeOps,
         nb::object params, bool catchErrors) {
        if (params.is_none()) {
          return RDKit::sanitizeReaction(
              rxn, sanitizeOps, RDKit::RxnOps::DefaultRxnAdjustParams(),
              catchErrors);
        }
        return RDKit::sanitizeReaction(
            rxn, sanitizeOps,
            nb::cast<const RDKit::MolOps::AdjustQueryParameters &>(params),
            catchErrors);
      },
      "rxn"_a,
      "sanitizeOps"_a = static_cast<unsigned int>(RDKit::RxnOps::SANITIZE_ALL),
      "params"_a = nb::none(), "catchErrors"_a = false,
      R"DOC(Does some sanitization of the reactant and product templates of a reaction.

    - The reaction is modified in place.
    - If sanitization fails, an exception will be thrown unless catchErrors is set

  ARGUMENTS:

    - rxn: the reaction to be modified
    - sanitizeOps: (optional) reaction sanitization operations to be carried out
      these should be constructed by or'ing together the
      operations in rdkit.Chem.rdChemReactions.SanitizeFlags
    - optional adjustment parameters for changing the meaning of the substructure
      matching done in the templates.  The default is
      rdkit.Chem.rdChemReactions.DefaultRxnAdjustParams which aromatizes
      kekule structures if possible.
    - catchErrors: (optional) if provided, instead of raising an exception
      when sanitization fails (the default behavior), the
      first operation that failed (as defined in rdkit.Chem.rdChemReactions.SanitizeFlags)
      is returned. Zero is returned on success.

  The operations carried out by default are:
    1) fixRGroups(): sets R group labels on mapped dummy atoms when possible
    2) fixAtomMaps(): attempts to set atom maps on unmapped R groups
    3) adjustTemplate(): calls adjustQueryProperties() on all reactant templates
    4) fixHs(): merges explicit Hs in the reactant templates that don't map to heavy atoms
)DOC");

  m.def(
      "SanitizeRxnAsMols", RDKit::RxnOps::sanitizeRxnAsMols, "rxn"_a,
      "sanitizeOps"_a = static_cast<unsigned int>(RDKit::MolOps::SANITIZE_ALL),
      "Does the usual molecular sanitization on each reactant, agent, and product of the reaction");

  wrap_enumeration(m);
}
