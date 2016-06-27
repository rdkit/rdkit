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
#include "RGroupTemplate.h"
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include "../../ReactionRunner.h"
#include <algorithm>

namespace RDKit {
// if we have more than 50 ring closures in a smiles string, this technique fails.
//   We may need to add better bookkeeping, but for now we support up to
//   15 mapped atoms (actualy max atom map number is 15)
const int arbitrary_ring_closure_offset = 40;
const int MAX_MAPNO = 15;

namespace {
std::string LabelSidechain(unsigned int mapno, unsigned int rgroupCount = 0) {
  PRECONDITION(rgroupCount < 4, "Too many rgroups attached to a single atom");
  unsigned int closure = arbitrary_ring_closure_offset + mapno*4 + rgroupCount;
  PRECONDITION(closure < 100, "Too many atom maps in Reaction for Template Enumeration");  
  return std::string("%") + boost::lexical_cast<std::string>(closure);
}

// Default label for a scaffold closure
std::string LabelScaffold(unsigned int mapno) {
  unsigned int closure = arbitrary_ring_closure_offset + mapno*4;
  PRECONDITION(closure < 100, "Too many atom maps in Reaction for Template Enumeration");  
  return std::string("%") + boost::lexical_cast<std::string>(closure);
}

// Expand the rgroups for this scaffold label
//  Example
// "%51" => "%51%52%53" to support three attachment points
std::string ExpandScaffoldLabel(unsigned int mapno, int numRgroups) {
  std::string res;
  for (int i=0; i<numRgroups; ++i) {
    res += LabelSidechain(mapno, i);
  }
  return res;
}

namespace {
std::string GetSmilesBond(int bondtype) {
  if      (bondtype == static_cast<int>(Bond::DOUBLE))   return "=";
  else if (bondtype == static_cast<int>(Bond::TRIPLE))   return "#";
  else if (bondtype == static_cast<int>(Bond::AROMATIC)) return ":";
  return "";
}
}

bool AddRGroupClosures(ROMol *mol,
                       RGroupTemplate &rtemplate) {
  for(ROMol::AtomIterator it=mol->beginAtoms(); it != mol->endAtoms(); ++it) {
    std::vector<int> atomRgroups;
    std::vector<int> bondRgroups;
    if ((*it)->getPropIfPresent(common_properties::_rgroupAtomMaps,
                                atomRgroups)) {
      const bool hasbonds = (*it)->getPropIfPresent(
          common_properties::_rgroupBonds, bondRgroups);

      for(size_t i=0;i<atomRgroups.size();++i) {
        const int &mapno = atomRgroups[i];
        const int bondtype = (hasbonds && bondRgroups.size() >= i) ?
            bondRgroups[i] : -1;

        rtemplate.addReactionMapping(mapno);
        
        std::string rgroupClosure;
        (*it)->getPropIfPresent(common_properties::_supplementalSmilesLabel,
                                rgroupClosure);
        rgroupClosure += GetSmilesBond(bondtype) + LabelSidechain(mapno);

        (*it)->setProp(common_properties::_supplementalSmilesLabel, rgroupClosure);
      }
    }
  }
  return true;
}

bool GetScaffoldFromRxn(const ChemicalReaction &rxn,
                               ReactantTemplates &templates) {
  PRECONDITION(rxn.getProducts().size() == 1,
               "Can only deal with one product in a reaction for now");
  templates.m_mappings.clear();
  templates.m_scaffoldSmiles = "";
  
  boost::shared_ptr<RWMol> scaffold(new RWMol(*rxn.getProducts()[0]));
  unsigned int sanitizeOps;
  MolOps::sanitizeMol(*scaffold.get(), sanitizeOps, MolOps::SANITIZE_FINDRADICALS);
      
  std::vector<Atom*> atomsToRemove;

  std::vector<unsigned int> &mappings = templates.m_mappings;
  boost::array<int, 100> counts;
  std::fill(counts.begin(), counts.end(), 0);

  for(ROMol::AtomIterator ai=scaffold->beginAtoms(); ai != scaffold->endAtoms();
      ++ai) {
    if ((*ai)->getAtomicNum() == 0) { // proper group
      atomsToRemove.push_back((*ai));
      int mapno=-1;
      if ((*ai)->getPropIfPresent(common_properties::molAtomMapNumber, mapno)) {
        if (mapno > MAX_MAPNO) {
          BOOST_LOG(rdErrorLog)<<"Mapped atom number too high in reaction\n" <<
              "\tGot " << mapno << " maximum is " << MAX_MAPNO << "\n";
          return false;
        }
        
        ROMol::ADJ_ITER nbrIdx,endNbrs;
        boost::tie(nbrIdx,endNbrs) = scaffold->getAtomNeighbors(*ai);
        int count = 0;
        while(nbrIdx != endNbrs) {
          Atom *nbr = scaffold->getAtomWithIdx(*nbrIdx);
          ++count;
          std::string rgroup;
          nbr->getPropIfPresent(common_properties::_supplementalSmilesLabel,
                                rgroup);
          rgroup += LabelScaffold(mapno);//, counts[mapno]);
          counts[mapno]++;
          mappings.push_back(mapno);
          nbr->setProp(common_properties::_supplementalSmilesLabel, rgroup);

          ++nbrIdx;
        }
        PRECONDITION(count==1, "Too many bonds from an rgroup");
      }
    } else {
      if ((*ai)->getImplicitValence()) { // Can we make any bonds?
        int mapno=-1;
        if ((*ai)->getPropIfPresent(common_properties::molAtomMapNumber, mapno)) {
          if (mapno > MAX_MAPNO) {
            BOOST_LOG(rdErrorLog)<<"Mapped atom number too high in reaction\n" <<
                "\tGot " << mapno << " maximum is " << MAX_MAPNO << "\n";
            return false;
          }
          
          if (!counts[mapno]) {
            mappings.push_back(mapno);
          }

          (*ai)->setProp(common_properties::_supplementalSmilesLabel,
                         LabelScaffold(mapno));//, counts[mapno]));
          counts[mapno]++;
        }
      }
      if ((*ai)->hasProp(common_properties::molAtomMapNumber))
        (*ai)->clearProp(common_properties::molAtomMapNumber);
    }
  }
  for (size_t i=0;i<atomsToRemove.size();++i)
    scaffold->removeAtom(atomsToRemove[i]);
      
  MolOps::sanitizeMol(*scaffold.get(), sanitizeOps, MolOps::SANITIZE_ADJUSTHS);
  templates.m_scaffoldSmiles = MolToSmiles(*scaffold.get(), true);
  return templates.m_mappings.size() > 0;
}
}

int RGroupTemplate::addReactionMapping(unsigned int mapno) {
  for(size_t i=0;i<mappings.size();++i) {
    if (mappings[i] == rdcast<int>(mapno)) {
      counts[i]++;
      return counts[i];
    }
  }
  mappings.push_back(mapno);
  counts.push_back(1);
  return 1;
}

bool RGroupTemplate::initialize(ROMOL_SPTR &product,
                                unsigned int /*reactantIdx*/,
                                unsigned int /*bbIdx*/) {
  const bool addDummyAtoms = false;
  ROMOL_SPTR sidechain(reduceProductToSideChains(product, addDummyAtoms));
  sideChainMol.reset();
  mappings.clear();
  counts.clear();
  smiles = "";
  if (!AddRGroupClosures(sidechain.get(), *this)) {
    mappings.clear();
    return false;
  } else {
    coreAndSideChainMol = product;
    sideChainMol = sidechain;
    
    std::vector<int> temp;
    const bool isomeric=true;
    smiles = MolToSmiles(*sideChainMol, isomeric);
    return true;
  }
}
  
bool RGroupTemplate::isValid() const {
  // can't handle multiple bonds to one reactant atom yet...
  return sideChainMol.get();// && rgroups.size() == 1; 
}

namespace {
bool replace(std::string& str, const std::string& from, const std::string& to) {
  size_t start_pos = str.find(from);
  if(start_pos == std::string::npos) {
    PRECONDITION(0, "Missing closure in template scaffold");
    return false;
  }
  str.replace(start_pos, from.length(), to);
  return true;
}
}

std::string ReactantTemplates::smiles(const RGROUPS &reactantIds) const {
  std::string res;
  // keep track of the mappings used
  boost::array<int, 100> counts;
  std::fill(counts.begin(), counts.end(), 0);
  
  for (size_t i=0; i<reactantIds.size(); ++i) {
    const RGroupTemplate &rgroup = m_templates[i][reactantIds[i]];
    res += "." + rgroup.smiles;
    for (size_t mapidx=0; mapidx<rgroup.mappings.size(); ++mapidx) {
      unsigned int mapno = rgroup.mappings[mapidx];
      counts[mapno] += rgroup.counts[mapidx];
    }
  }

  // Adjusts the scaffold to match the closures
  std::string scaffold = m_scaffoldSmiles;
  for (size_t i=0; i<m_mappings.size(); ++i) {
    unsigned int mapno = m_mappings[i];
    if (!counts[mapno]) { // remove unused closure
      replace(scaffold, LabelScaffold(mapno), "");
    } else if (counts[i] > 1) { // expands closures
      replace(scaffold,
              LabelScaffold(mapno),
              ExpandScaffoldLabel(mapno, counts[i]));
    }
  }
  
  return scaffold + res;
}

// extract the bbs into the reaction templates
bool ReactantsToTemplates(ReactantTemplates &templates,
                          const ChemicalReaction &rxn,
                          const BBS &bbs) {
  // get the scaffold smiles from the BBs -- note that this will ONLY WORK
  //  with proper RGroups, not mapped atom with no rgroups -- working on that
  if (!GetScaffoldFromRxn(rxn, templates))
    return false;
  
  {
    VectVectRGroupTemplate vect(bbs.size()) ;
    templates.m_templates.swap( vect );
  }
    
  for(size_t reactantIdx=0;reactantIdx<bbs.size();++reactantIdx) {
    std::set<std::string> known_templates;
      
    for(size_t bbIdx=0; bbIdx<bbs[reactantIdx].size(); ++bbIdx) {
      std::vector<MOL_SPTR_VECT> res = run_Reactant(rxn,
                                                    bbs[reactantIdx][bbIdx],
                                                    reactantIdx);
      templates.m_templates.reserve(bbs[reactantIdx].size());
        
      for(size_t j=0; j<res.size(); ++j) {
        for(size_t prodidx=0; prodidx<res[j].size(); ++prodidx) {
          boost::dynamic_bitset<> bv(res[j][prodidx]->getNumAtoms());
          RGroupTemplate temp;
          temp.initialize(res[j][prodidx], reactantIdx, bbIdx);
          
          if (temp.isValid()) {
            if (known_templates.find(temp.smiles) == known_templates.end()) {
              templates.m_templates[reactantIdx].push_back( temp );
              known_templates.insert(temp.smiles);
            }
          } else {
            /// skip???
            PRECONDITION(0, "Invalid reactant");
          }
        }
      }
    }
  }
  return true;
}

}
