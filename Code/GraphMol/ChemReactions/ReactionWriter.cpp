//
//  Copyright (c) 2010-2024, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
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
#include <GraphMol/FileParsers/FileParserUtils.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>

#include <sstream>

namespace {

void setRXNRoleOfAllMoleculeAtoms(RDKit::ROMol &mol, int role) {
  RDKit::ROMol::ATOM_ITER_PAIR atItP = mol.getVertices();
  while (atItP.first != atItP.second) {
    RDKit::Atom *oAtom = mol[*(atItP.first++)];
    oAtom->setProp(RDKit::common_properties::molRxnRole, role);
  }
}

std::string molToString(RDKit::ROMol &mol, bool toSmiles,
                        const RDKit::SmilesWriteParams &params) {
  std::string res = "";
  if (toSmiles) {
    res = MolToSmiles(mol, params);
  } else {
    res = MolToSmarts(mol, params);
  }
  std::vector<int> mapping;
  if (RDKit::MolOps::getMolFrags(mol, mapping) > 1) {
    res = "(" + res + ")";
  }
  return res;
}

std::string chemicalReactionTemplatesToString(
    const RDKit::ChemicalReaction &rxn, RDKit::ReactionMoleculeType type,
    bool toSmiles, const RDKit::SmilesWriteParams &params) {
  std::string res = "";
  std::vector<std::string> vfragsmi;
  auto begin = getStartIterator(rxn, type);
  auto end = getEndIterator(rxn, type);
  for (; begin != end; ++begin) {
    vfragsmi.push_back(molToString(**begin, toSmiles, params));
  }
  if (params.canonical) {
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

/*
* next step: figuring out how CX extensions work -- what properties do these extensions set? I have some examples in the description of the PR

* printprops as dict...
*
*/
void insertTemplates(const RDKit::ChemicalReaction &rxn,
                    RDKit::RWMol &rwmol, std::vector<unsigned int> &atomOrdering,
                    std::vector<unsigned int> &bondOrdering, RDKit::ReactionMoleculeType type) {
  auto at_count = rwmol.getNumAtoms();
  auto bnd_count = rwmol.getNumBonds();
  std::vector<unsigned int> prevAtomOrdering;
  std::vector<unsigned int> prevBondOrdering;

  for (auto begin = getStartIterator(rxn, type);
        begin != getEndIterator(rxn, type); ++begin) {
    rwmol.insertMol(**begin);

    // getting the smiles atom and bond output ordering
    std::vector<unsigned int> prevAtomOrdering;
    std::vector<unsigned int> prevBondOrdering;
    (*begin)->getProp(RDKit::common_properties::_smilesAtomOutputOrder, prevAtomOrdering);
    (*begin)->getProp(RDKit::common_properties::_smilesBondOutputOrder, prevBondOrdering);
    for (auto i : prevAtomOrdering) {
      atomOrdering.push_back(i + at_count);
    }
    for (auto i : prevBondOrdering) {
      bondOrdering.push_back(i + bnd_count);
    }
    at_count = rwmol.getNumAtoms();
    bnd_count = rwmol.getNumBonds();
  }
}

std::string chemicalReactionToRxnToString(
    const RDKit::ChemicalReaction &rxn, bool toSmiles,
    const RDKit::SmilesWriteParams &params, bool includeCX) {
  std::cout << "here1\n";
  std::string res = "";
  res +=
      chemicalReactionTemplatesToString(rxn, RDKit::Reactant, toSmiles, params);
  std::cout << "here2\n";
  res += ">";
  res += chemicalReactionTemplatesToString(rxn, RDKit::Agent, toSmiles, params);
  std::cout << "here3\n";
  res += ">";
  res +=
      chemicalReactionTemplatesToString(rxn, RDKit::Product, toSmiles, params);
  std::cout << "here4\n";

  std::cout << "in here\n";
  std::cout << "includeCX: " << params.includeCX << std::endl;


  if (includeCX) {
    // I think this will work? combine reactants, agents and products into a
    // single molecule and get the cxextension for that
    // blank rwmol
    RDKit::RWMol rwmol;
    std::vector<unsigned int> atomOrdering;
    std::vector<unsigned int> bondOrdering;

    // TODO
    /*
     * what I'm doing right now is making a 'super molecule' out of all the reactants, agents
     * and products, and getting the total atom & bond ordering, then using that to get the extension
     * and add that to the reaction cxsmarts
     * 
     * what we want to do instead is send a vector of all the reactants, agents & products to getCXExtensions
     * and have that return the CX extension directly 
     *
    */
    insertTemplates(rxn, rwmol, atomOrdering, bondOrdering, RDKit::Reactant);
    insertTemplates(rxn, rwmol, atomOrdering, bondOrdering, RDKit::Agent);
    insertTemplates(rxn, rwmol, atomOrdering, bondOrdering, RDKit::Product);

    rwmol.setProp(RDKit::common_properties::_smilesAtomOutputOrder, atomOrdering, true); // TODO: Need to get since the SMILES atom order might not be the same as the index order for our RDKit molecule...
    rwmol.setProp(RDKit::common_properties::_smilesBondOutputOrder, bondOrdering, true); // in python can get _smilesAtomOutputOrder property and see how it works...

    // ignore atom properties -- really this is to avoid the atom map numbers, maybe we just want to remove
    // those

    // make vector of all reactants, agents and product
    // MOL_SPTR_VECT all_templates; 

    // insert reactants, agents, and products into all_templates

    auto flags = RDKit::SmilesWrite::CXSmilesFields::CX_ATOM_PROPS ^ RDKit::SmilesWrite::CXSmilesFields::CX_ALL;
    auto ext = RDKit::SmilesWrite::getCXExtensions(rwmol, flags); // this should take vect<ROMol> of reactants, agents, and products
    if (!ext.empty()) {
      res += " ";
      res += ext;
    }
  }

  return res;
}

void write_template(std::ostringstream &res, RDKit::ROMol &tpl) {
  RDKit::RWMol trwmol(tpl);

  if (trwmol.needsUpdatePropertyCache()) {
    trwmol.updatePropertyCache(false);
  }
  RDKit::FileParserUtils::moveAdditionalPropertiesToSGroups(trwmol);

  res << RDKit::FileParserUtils::getV3000CTAB(trwmol, -1);
}

}  // namespace

namespace RDKit {

//! returns the reaction SMARTS for a reaction
std::string ChemicalReactionToRxnSmarts(const ChemicalReaction &rxn,
                                        const SmilesWriteParams &params) {
  return chemicalReactionToRxnToString(rxn, false, params, false);
};

//! returns the reaction SMILES for a reaction
std::string ChemicalReactionToRxnSmiles(const ChemicalReaction &rxn,
                                        const SmilesWriteParams &params) {
  return chemicalReactionToRxnToString(rxn, true, params, false);
};


//! returns the reaction SMARTS for a reaction with CX extension
std::string ChemicalReactionToRxnCXSmarts(const ChemicalReaction &rxn,
                                        const SmilesWriteParams &params) {
  std::cout << "before calling\n";
  return chemicalReactionToRxnToString(rxn, false, params, true);
};

//! returns the reaction SMILES for a reaction with CX extension
std::string ChemicalReactionToRxnCXSmiles(const ChemicalReaction &rxn,
                                        const SmilesWriteParams &params) {
  return chemicalReactionToRxnToString(rxn, true, params, true);
};

//! returns an RXN block for a reaction
std::string ChemicalReactionToV3KRxnBlock(const ChemicalReaction &rxn,
                                          bool separateAgents) {
  std::ostringstream res;
  res << "$RXN V3000\n\n      RDKit\n\n";

  if (separateAgents) {
    res << "M  V30 COUNTS " << rxn.getNumReactantTemplates() << " "
        << rxn.getNumProductTemplates() << " " << rxn.getNumAgentTemplates();
  } else {
    res << "M  V30 COUNTS "
        << rxn.getNumReactantTemplates() + rxn.getNumAgentTemplates() << " "
        << rxn.getNumProductTemplates();
  }
  res << "\n";

  res << "M  V30 BEGIN REACTANT\n";
  for (const auto &rt : boost::make_iterator_range(
           rxn.beginReactantTemplates(), rxn.endReactantTemplates())) {
    write_template(res, *rt);
  }
  if (!separateAgents) {
    for (const auto &at : boost::make_iterator_range(rxn.beginAgentTemplates(),
                                                     rxn.endAgentTemplates())) {
      write_template(res, *at);
    }
  }
  res << "M  V30 END REACTANT\n";
  res << "M  V30 BEGIN PRODUCT\n";
  for (const auto &pt : boost::make_iterator_range(rxn.beginProductTemplates(),
                                                   rxn.endProductTemplates())) {
    write_template(res, *pt);
  }
  res << "M  V30 END PRODUCT\n";
  if (separateAgents) {
    res << "M  V30 BEGIN AGENT\n";
    for (const auto &rt : boost::make_iterator_range(rxn.beginAgentTemplates(),
                                                     rxn.endAgentTemplates())) {
      write_template(res, *rt);
    }
    res << "M  V30 END AGENT\n";
  }

  res << "M  END\n";
  return res.str();
}

//! returns an RXN block for a reaction
std::string ChemicalReactionToRxnBlock(const ChemicalReaction &rxn,
                                       bool separateAgents, bool forceV3000) {
  if (forceV3000) {
    return ChemicalReactionToV3KRxnBlock(rxn, separateAgents);
  }
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
    res << "$MOL\n";
    res << MolToMolBlock(**iter, true, -1, false);
  }
  if (!separateAgents) {
    for (auto iter = rxn.beginAgentTemplates(); iter != rxn.endAgentTemplates();
         ++iter) {
      res << "$MOL\n";
      res << MolToMolBlock(**iter, true, -1, false);
    }
  }
  for (auto iter = rxn.beginProductTemplates();
       iter != rxn.endProductTemplates(); ++iter) {
    res << "$MOL\n";
    res << MolToMolBlock(**iter, true, -1, false);
  }
  if (separateAgents) {
    for (auto iter = rxn.beginAgentTemplates(); iter != rxn.endAgentTemplates();
         ++iter) {
      res << "$MOL\n";
      res << MolToMolBlock(**iter, true, -1, false);
    }
  }
  return res.str();
};

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
