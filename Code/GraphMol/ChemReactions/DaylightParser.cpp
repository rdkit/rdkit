// $Id$
//
//  Copyright (c) 2007, Novartis Institutes for BioMedical Research Inc.
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
          if(!(*prodAtomIt)->hasProp("molAtomMapNumber")) continue;
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
                                                 std::map<std::string,std::string> *replacements) {
    std::size_t pos=text.find(">>");
    if(pos==std::string::npos){
      throw ChemicalReactionParserException("a reaction requires at least one reactant and one product");
    }
    if(text.rfind(">>")!=pos){
      throw ChemicalReactionParserException("multi-step reactions not supported");
    }
    
    std::string reactText=text.substr(0,pos);
    std::vector<std::string> reactSmarts=DaylightParserUtils::splitSmartsIntoComponents(reactText);

    std::string productText=text.substr(pos+2);
    
    ChemicalReaction *rxn=new ChemicalReaction();
   
    for(std::vector<std::string>::const_iterator txtIt=reactSmarts.begin();
        txtIt!=reactSmarts.end();++txtIt){
      ROMol *mol=SmartsToMol(*txtIt,0,false,replacements);
      if(!mol){
        std::string errMsg="Problems constructing reactant from SMARTS: ";
        errMsg += *txtIt;
        throw ChemicalReactionParserException(errMsg);
      }
      rxn->addReactantTemplate(ROMOL_SPTR(mol));        
    }
    //std::cerr << " ---------------------------------------------------------" << std::endl;    
    ROMol *prodMol=SmartsToMol(productText,0,false,replacements);
    if(!prodMol){
      std::string errMsg="Problems constructing product from SMARTS: ";
      errMsg += productText;
      throw ChemicalReactionParserException(errMsg);
    }
    std::vector<ROMOL_SPTR> prods=MolOps::getMolFrags(*prodMol,false);
    delete prodMol;
    for(std::vector<ROMOL_SPTR>::iterator pIt=prods.begin();
        pIt!=prods.end();++pIt){
      rxn->addProductTemplate(*pIt);        
    }
    //std::cerr << " ---------------------------------------------------------" << std::endl;    
    DaylightParserUtils::updateProductsStereochem(rxn);
    
    // "SMARTS"-based reactions have implicit properties
    rxn->setImplicitPropertiesFlag(true);
    
    return rxn;    
  }

}
