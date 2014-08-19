// $Id$
//
//  Copyright (c) 2007-2014, Novartis Institutes for BioMedical Research Inc.
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

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <boost/algorithm/string.hpp>
#include <vector>
#include <string>

namespace RDKit {
  namespace DaylightParserUtils {
    void updateProductsStereochem(ChemicalReaction *rxn){
      // EFF: this isn't the speediest code the world has ever seen,
      //   but we hopefully aren't going to be calling this a lot.
      for(MOL_SPTR_VECT::const_iterator prodIt=rxn->beginProductTemplates();
          prodIt!=rxn->endProductTemplates();++prodIt){
        for(ROMol::AtomIterator prodAtomIt=(*prodIt)->beginAtoms();
            prodAtomIt!=(*prodIt)->endAtoms();++prodAtomIt){
          if(!(*prodAtomIt)->hasProp("molAtomMapNumber")){
            // if we have stereochemistry specified, it's automatically creating stereochem:
            (*prodAtomIt)->setProp("molInversionFlag",4);
            continue;
          } 
          int mapNum;
          (*prodAtomIt)->getProp("molAtomMapNumber",mapNum);
          for(MOL_SPTR_VECT::const_iterator reactIt=rxn->beginReactantTemplates();
              reactIt!=rxn->endReactantTemplates();++reactIt){
            for(ROMol::AtomIterator reactAtomIt=(*reactIt)->beginAtoms();
                reactAtomIt!=(*reactIt)->endAtoms();++reactAtomIt){
              if(!(*reactAtomIt)->hasProp("molAtomMapNumber")) continue;
              int reactMapNum;
              (*reactAtomIt)->getProp("molAtomMapNumber",reactMapNum);
              if(reactMapNum==mapNum){
                // finally, in the bowels of the nesting, we get to some actual
                // work:
                if((*prodAtomIt)->getChiralTag()!=Atom::CHI_UNSPECIFIED &&
                   (*prodAtomIt)->getChiralTag()!=Atom::CHI_OTHER) {
                  if((*reactAtomIt)->getChiralTag()!=Atom::CHI_UNSPECIFIED &&
                     (*reactAtomIt)->getChiralTag()!=Atom::CHI_OTHER){
                    // both have stereochem specified:
                    if((*reactAtomIt)->getChiralTag()==(*prodAtomIt)->getChiralTag()){
                      (*prodAtomIt)->setProp("molInversionFlag",2);
                    } else {
                      // FIX: this is technically fragile: it should be checking
                      // if the atoms both have tetrahedral chirality. However,
                      // at the moment that's the only chirality available, so there's
                      // no need to go monkeying around.
                      (*prodAtomIt)->setProp("molInversionFlag",1);
                    }
                  } else {
                    // stereochem in the product, but not in the reactant
                    (*prodAtomIt)->setProp("molInversionFlag",4);
                  }
                } else if((*reactAtomIt)->getChiralTag()!=Atom::CHI_UNSPECIFIED &&
                          (*reactAtomIt)->getChiralTag()!=Atom::CHI_OTHER){
                  // stereochem in the reactant, but not the product:
                  (*prodAtomIt)->setProp("molInversionFlag",3);                  
                }
              }
            }
          }
        }
      }
    }
    
    std::vector<std::string> splitSmartsIntoComponents(const std::string &reactText){
      std::vector<std::string> res;
      unsigned int pos=0;
      unsigned int blockStart=0;
      unsigned int level=0;
      unsigned int inBlock=0;
      while(pos<reactText.size()){
        if(reactText[pos]=='('){
          if(pos==blockStart){
            inBlock=1;
          }
          ++level;
        } else if(reactText[pos]==')'){
          if(level==1 && inBlock){
            // this closes a block
            inBlock=2;
          }
          --level;
        } else if(level==0 && reactText[pos]=='.'){
          if(inBlock==2){
            std::string element=reactText.substr(blockStart+1,pos-blockStart-2);
            res.push_back(element);
          } else {
            std::string element=reactText.substr(blockStart,pos-blockStart);
            res.push_back(element);
          }
          blockStart = pos+1;
          inBlock=0;
        }
        ++pos;
      }
      if(blockStart<pos){
        if(inBlock==2){
          std::string element=reactText.substr(blockStart+1,pos-blockStart-2);
          res.push_back(element);
        } else {
          std::string element=reactText.substr(blockStart,pos-blockStart);
          res.push_back(element);
        }
      }
      return res;
    }
  } // end of namespace DaylightParserUtils
  
  ChemicalReaction * RxnSmartsToChemicalReaction(const std::string &text,
                                                 std::map<std::string,std::string> *replacements,
                                                 bool useSmiles) {
    std::size_t pos1=text.find(">");
    std::size_t pos2=text.rfind(">");
    if(pos1==std::string::npos){
      throw ChemicalReactionParserException("a reaction requires at least one reactant and one product");
    }
    if(text.find(">",pos1+1)!=pos2){
      throw ChemicalReactionParserException("multi-step reactions not supported");
    }
    
    std::string reactText=text.substr(0,pos1);
    std::string agentText="";
    if(pos2!=pos1+1){
      agentText=text.substr(pos1+1,(pos2-pos1)-1);
    }
    std::string productText=text.substr(pos2+1);

    // recognize changes within the same molecules, e.g., intra molecular bond formation
    // therefore we need to correctly interpret parenthesis and dots in the reaction smarts
    std::vector<std::string> reactSmarts=DaylightParserUtils::splitSmartsIntoComponents(reactText);
    std::vector<std::string> productSmarts=DaylightParserUtils::splitSmartsIntoComponents(productText);

    ChemicalReaction *rxn=new ChemicalReaction();
   
    for(std::vector<std::string>::const_iterator txtIt=reactSmarts.begin();
        txtIt!=reactSmarts.end();++txtIt){
      ROMol *mol;
      if(!useSmiles){
        mol=SmartsToMol(*txtIt,0,false,replacements);
      }
      else{
        mol=SmilesToMol(*txtIt,0,false,replacements);
      }
      if(!mol){
        std::string errMsg="Problems constructing reactant from SMARTS: ";
        errMsg += *txtIt;
        throw ChemicalReactionParserException(errMsg);
      }
      rxn->addReactantTemplate(ROMOL_SPTR(mol));        
    }

    for(std::vector<std::string>::const_iterator txtIt=productSmarts.begin();
        txtIt!=productSmarts.end();++txtIt){
      ROMol *mol;
      if(!useSmiles){
        mol=SmartsToMol(*txtIt,0,false,replacements);
      }
      else{
        mol=SmilesToMol(*txtIt,0,false,replacements);
      }
      if(!mol){
        std::string errMsg="Problems constructing product from SMARTS: ";
        errMsg += *txtIt;
        throw ChemicalReactionParserException(errMsg);
      }
      rxn->addProductTemplate(ROMOL_SPTR(mol));
    }
    DaylightParserUtils::updateProductsStereochem(rxn);

    ROMol *agentMol;
    //allow a reaction template to have no agent specified
    if(agentText.size() != 0){
      if(!useSmiles) agentMol=SmartsToMol(agentText,0,false,replacements);
      else agentMol=SmilesToMol(agentText,0,false,replacements);

      if(!agentMol){
        std::string errMsg="Problems constructing agent from SMARTS: ";
        errMsg += agentText;
        throw ChemicalReactionParserException(errMsg);
      }
      std::vector<ROMOL_SPTR> agents=MolOps::getMolFrags(*agentMol,false);
      delete agentMol;
      for(std::vector<ROMOL_SPTR>::iterator aIt=agents.begin();
        aIt!=agents.end();++aIt){
        rxn->addAgentTemplate(*aIt);
      }
    }
    
    // "SMARTS"-based reactions have implicit properties
    rxn->setImplicitPropertiesFlag(true);

    return rxn;    
  }

}
