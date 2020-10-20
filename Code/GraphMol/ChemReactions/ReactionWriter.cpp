//
//  Copyright (c) 2010-2018, Novartis Institutes for BioMedical Research Inc.
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

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/MolOps.h>
#include <sstream>

namespace {

void setRXNRoleOfAllMoleculeAtoms(RDKit::ROMol &mol, int role) {
  for(auto *oAtom : mol.atoms()) {
    oAtom->setProp(RDKit::common_properties::molRxnRole, role);
  }
}

std::string molToString(RDKit::ROMol &mol, bool toSmiles) {
  std::string res = "";
  if (toSmiles) {
    res = MolToSmiles(mol, true);
  } else {
    res = MolToSmarts(mol, true);
  }
  std::vector<int> mapping;
  if (RDKit::MolOps::getMolFrags(mol, mapping) > 1) {
    res = "(" + res + ")";
  }
  return res;
}

std::string chemicalReactionTemplatesToString(
    const RDKit::ChemicalReaction &rxn, RDKit::ReactionMoleculeType type,
    bool toSmiles, bool canonical) {
  std::string res = "";
  std::vector<std::string> vfragsmi;
  auto begin = getStartIterator(rxn, type);
  auto end = getEndIterator(rxn, type);
  for (; begin != end; ++begin) {
    vfragsmi.push_back(molToString(**begin, toSmiles));
  }
  if (canonical) {
    std::sort(vfragsmi.begin(), vfragsmi.end());
  }
  for (unsigned i = 0; i < vfragsmi.size(); ++i) {
    res += vfragsmi[i];
    if (i < vfragsmi.size() - 1) {
      res += ".";
    }
  }
  return res;
}

std::string chemicalReactionToRxnToString(const RDKit::ChemicalReaction &rxn,
                                          bool toSmiles, bool canonical) {
  std::string res = "";
  res += chemicalReactionTemplatesToString(rxn, RDKit::Reactant, toSmiles,
                                           canonical);
  res += ">";
  res +=
      chemicalReactionTemplatesToString(rxn, RDKit::Agent, toSmiles, canonical);
  res += ">";
  res += chemicalReactionTemplatesToString(rxn, RDKit::Product, toSmiles,
                                           canonical);
  return res;
}
}  // namespace

namespace RDKit {

//! returns the reaction SMARTS for a reaction
std::string ChemicalReactionToRxnSmarts(const ChemicalReaction &rxn) {
  return chemicalReactionToRxnToString(rxn, false, false);
};

//! returns the reaction SMILES for a reaction
std::string ChemicalReactionToRxnSmiles(const ChemicalReaction &rxn,
                                        bool canonical) {
  return chemicalReactionToRxnToString(rxn, true, canonical);
};

#if 1
//! returns an RXN block for a reaction
std::string ChemicalReactionToRxnBlock(const ChemicalReaction &rxn,
                                       bool separateAgents) {
  std::ostringstream res;
  res << "$RXN\n\n      RDKit\n\n";
  if (separateAgents) {
    res << std::setw(3) << rxn.getNumReactantTemplates() << std::setw(3)
        << rxn.getNumProductTemplates() << std::setw(3)
        << rxn.getNumAgentTemplates() << "\n";
  } else {
    res << std::setw(3)
        << (rxn.getNumReactantTemplates() + rxn.getNumAgentTemplates())
        << std::setw(3) << rxn.getNumProductTemplates() << "\n";
  }

  for (auto iter = rxn.beginReactantTemplates();
       iter != rxn.endReactantTemplates(); ++iter) {
    // to write the mol block, we need ring information:
    MolOps::findSSSR(**iter);
    res << "$MOL\n";
    res << MolToMolBlock(**iter, true, -1, false);
  }
  if (!separateAgents) {
    for (auto iter = rxn.beginAgentTemplates(); iter != rxn.endAgentTemplates();
         ++iter) {
      // to write the mol block, we need ring information:
      MolOps::findSSSR(**iter);
      res << "$MOL\n";
      res << MolToMolBlock(**iter, true, -1, false);
    }
  }
  for (auto iter = rxn.beginProductTemplates();
       iter != rxn.endProductTemplates(); ++iter) {
    // to write the mol block, we need ring information:
    MolOps::findSSSR(**iter);
    res << "$MOL\n";
    res << MolToMolBlock(**iter, true, -1, false);
  }
  if (separateAgents) {
    for (auto iter = rxn.beginAgentTemplates(); iter != rxn.endAgentTemplates();
         ++iter) {
      // to write the mol block, we need ring information:
      MolOps::findSSSR(**iter);
      res << "$MOL\n";
      res << MolToMolBlock(**iter, true, -1, false);
    }
  }
  return res.str();
};
#endif

//! returns a ROMol with RXNMolRole used for a reaction
ROMol *ChemicalReactionToRxnMol(const ChemicalReaction &rxn) {
  auto *res = new RWMol();

  for (auto iter = rxn.beginReactantTemplates();
       iter != rxn.endReactantTemplates(); ++iter) {
    setRXNRoleOfAllMoleculeAtoms(*iter->get(), 1);
    res->insertMol(*iter->get());
  }
  for (auto iter = rxn.beginProductTemplates();
       iter != rxn.endProductTemplates(); ++iter) {
    setRXNRoleOfAllMoleculeAtoms(*iter->get(), 2);
    res->insertMol(*iter->get());
  }
  for (auto iter = rxn.beginAgentTemplates(); iter != rxn.endAgentTemplates();
       ++iter) {
    setRXNRoleOfAllMoleculeAtoms(*iter->get(), 3);
    res->insertMol(*iter->get());
  }
  return (ROMol *)res;
}
}  // namespace RDKit
