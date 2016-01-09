//
//  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
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
//       products derived from this software without specific prior written permission.
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
#include "SanitizeRxn.h"
#include <GraphMol/RDKitBase.h>

namespace RDKit {
namespace RxnOps {

namespace {
unsigned int getMaxProp(ChemicalReaction &rxn, const std::string &prop) {
  unsigned int max_atom;
  for(MOL_SPTR_VECT::iterator it = rxn.beginReactantTemplates();
      it != rxn.endReactantTemplates();
      ++it) {
    for(unsigned int idx=0;idx<(*it)->getNumAtoms(); ++idx) {
      Atom *atom = (*it)->getAtomWithIdx(idx);
      unsigned int map;
      if (atom->getPropIfPresent(prop, map)) {
        if (map > max_atom)
          max_atom = map;
      }
    }
  }

  for(MOL_SPTR_VECT::iterator it = rxn.beginAgentTemplates();
      it != rxn.endAgentTemplates();
      ++it) {
    for(unsigned int idx=0;idx<(*it)->getNumAtoms(); ++idx) {
      Atom *atom = (*it)->getAtomWithIdx(idx);
      unsigned int map;
      if (atom->getPropIfPresent(prop, map)) {
        if (map > max_atom)
          max_atom = map;
      }
    }
  }

  for(MOL_SPTR_VECT::iterator it = rxn.beginProductTemplates();
      it != rxn.endProductTemplates();
      ++it) {
    for(unsigned int idx=0;idx<(*it)->getNumAtoms(); ++idx) {
      Atom *atom = (*it)->getAtomWithIdx(idx);
      unsigned int map;
      if (atom->getPropIfPresent(prop, map)) {
        if (map > max_atom)
          max_atom = map;
      }
    }
  }
  
  return max_atom;
}
}

void fixAtomMaps(ChemicalReaction &rxn) {
  unsigned int max_atom_map = getMaxProp(rxn,common_properties::molAtomMapNumber);
  std::map<unsigned int, unsigned int> potential_mappings;
  
  for(MOL_SPTR_VECT::iterator it = rxn.beginReactantTemplates();
      it != rxn.endReactantTemplates();
      ++it) {
    for(unsigned int idx=0;idx<(*it)->getNumAtoms(); ++idx) {
      Atom *atom = (*it)->getAtomWithIdx(idx);
      unsigned int rlabel;
      // find any rlabels with no atom maps, and try to set atom maps
      //  Do I need to worry about the Agents?  I hope not!
      if (atom->getPropIfPresent(common_properties::_MolFileRLabel, rlabel)) {
        if(!atom->hasProp(common_properties::molAtomMapNumber)) {
          if(potential_mappings.find(rlabel) != potential_mappings.end()) {
            throw RxnSanitizeException(std::string("Duplicated RLabels"));
          }
          potential_mappings[rlabel] = rlabel+max_atom_map;
          atom->setProp(common_properties::molAtomMapNumber, rlabel+max_atom_map);
        }
      }
    }
  }
  
  if (!potential_mappings.size())
    return; // everything is ok!

  for(MOL_SPTR_VECT::iterator it = rxn.beginProductTemplates();
      it != rxn.endProductTemplates();
      ++it) {
    for(unsigned int idx=0;idx<(*it)->getNumAtoms(); ++idx) {
      Atom *atom = (*it)->getAtomWithIdx(idx);
      unsigned int rlabel;
      if (atom->getPropIfPresent(common_properties::_MolFileRLabel, rlabel)) {
        if(!atom->hasProp(common_properties::molAtomMapNumber)) {
          atom->setProp(common_properties::molAtomMapNumber,
                        potential_mappings[rlabel]);
        } else {
          if (atom->getProp<unsigned int>(common_properties::molAtomMapNumber) !=
              potential_mappings[rlabel]) {
            throw RxnSanitizeException("Unmapped Rlabel on product template, "
                                       "could not resolve atom map");
          }
        }
      }
    }
  }
}

// if we have query atoms without rlabels, make proper rlabels
void fixRGroups(ChemicalReaction &rxn) {
  unsigned int maxrlabel = getMaxProp(rxn, common_properties::_MolFileRLabel);
  std::map<unsigned int, unsigned int> remapped;
  for(MOL_SPTR_VECT::iterator it = rxn.beginReactantTemplates();
      it != rxn.endReactantTemplates();
      ++it) {
    for(unsigned int idx=0;idx<(*it)->getNumAtoms(); ++idx) {
          Atom *atom = (*it)->getAtomWithIdx(idx);
          unsigned int map;
          if(atom->getPropIfPresent(common_properties::molAtomMapNumber, map)) {
            if(atom->getAtomicNum() == 0) {
              if(!atom->hasProp(common_properties::_MolFileRLabel)) {
                atom->setProp(common_properties::_MolFileRLabel, maxrlabel+map);
                remapped[map] = maxrlabel+map;
              }
            }
          }
    }
  }
  
  if (!remapped.size())
    return;

  for(MOL_SPTR_VECT::iterator it = rxn.beginProductTemplates();
      it != rxn.endProductTemplates();
      ++it) {
    for(unsigned int idx=0;idx<(*it)->getNumAtoms(); ++idx) {
      Atom *atom = (*it)->getAtomWithIdx(idx);
      unsigned int map;
      if(atom->getPropIfPresent(common_properties::molAtomMapNumber, map)) {
        if(atom->getAtomicNum() == 0) {
          if(!atom->hasProp(common_properties::_MolFileRLabel)) {
            atom->setProp(common_properties::_MolFileRLabel, maxrlabel+map);
          } else {
            if (atom->getProp<unsigned int>(common_properties::_MolFileRLabel) !=
                remapped[map]) {
              // this might or might not actually be illegal :)
              throw RxnSanitizeException("Mismatched mappings for rgroups");
            }
          }
        }
      }
    }
  }
  return;
}

// might throw mol sanitization exception??? wrap in RxnSanitize?
void fixReactantTemplateAromaticity(ChemicalReaction &rxn) {
  unsigned int ops;
  for(MOL_SPTR_VECT::iterator it = rxn.beginReactantTemplates();
      it != rxn.endReactantTemplates();
      ++it) {
    // Cheat here, we know that this came from a RWMol
    RWMol * rw = dynamic_cast<RWMol*>(it->get());
    if (rw)
      sanitizeMol(*rw, ops, MolOps::SANITIZE_SETAROMATICITY);
    else
      PRECONDITION(rw, "Oops, not really a RWMol?");
  }
}

void fixHs(ChemicalReaction &rxn) {
  const bool mergeUnmappedOnly = true;
  for(MOL_SPTR_VECT::iterator it = rxn.beginReactantTemplates();
      it != rxn.endReactantTemplates();
      ++it) {
    RWMol * rw = dynamic_cast<RWMol*>(it->get());
    if (rw)
      MolOps::mergeQueryHs(*rw, mergeUnmappedOnly);
    else
      PRECONDITION(rw, "Oops, not really a RWMol?");    
  }
}

void sanitizeRxn(ChemicalReaction &rxn,
                 unsigned int &operationsThatFailed,
                 unsigned int ops)
{
  operationsThatFailed = SANITIZE_NONE;

  if (ops & SANITIZE_ATOM_MAPS) {
    operationsThatFailed = SANITIZE_ATOM_MAPS;
    fixAtomMaps(rxn);
  }

  if (ops & SANITIZE_RGROUP_NAMES) {
    operationsThatFailed = SANITIZE_RGROUP_NAMES;
    fixRGroups(rxn);
  }

  if (ops & SANITIZE_REAGENT_AROMATICITY) {
    operationsThatFailed = SANITIZE_REAGENT_AROMATICITY;
    fixReactantTemplateAromaticity(rxn);
  }

  if (ops & SANITIZE_MERGEHS) {
    operationsThatFailed = SANITIZE_MERGEHS;
    fixHs(rxn);
  }

  operationsThatFailed = SANITIZE_NONE;
}

void sanitizeRxn(ChemicalReaction &rxn) {
  unsigned int ops = 0;
  return sanitizeRxn(rxn, ops, SANITIZE_ALL);
}


}
}
