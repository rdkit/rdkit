// $Id$
//
//  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
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
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <sstream>

namespace RDKit {

  //! returns the reaction SMARTS for a reaction
  std::string ChemicalReactionToRxnSmarts(ChemicalReaction &rxn){
    std::string res="";
    for(MOL_SPTR_VECT::const_iterator iter=rxn.beginReactantTemplates();
	iter != rxn.endReactantTemplates();++iter){
      if(iter!=rxn.beginReactantTemplates()) res +=".";
      res += MolToSmarts(**iter,true);
    }
    res += ">>";
    for(MOL_SPTR_VECT::const_iterator iter=rxn.beginProductTemplates();
	iter != rxn.endProductTemplates();++iter){
      if(iter!=rxn.beginProductTemplates()) res +=".";
      res += MolToSmarts(**iter,true);
    }
    return res;
  };
#if 1
  //! returns an RXN block for a reaction
  std::string ChemicalReactionToRxnBlock(const ChemicalReaction &rxn){
    std::ostringstream res;
    res<<"$RXN\n\n      RDKit\n\n";
    res<<std::setw(3)<<rxn.getNumReactantTemplates()<<std::setw(3)<<rxn.getNumProductTemplates()<<"\n";

    
    for(MOL_SPTR_VECT::const_iterator iter=rxn.beginReactantTemplates();
	iter != rxn.endReactantTemplates();++iter){
      // to write the mol block, we need ring information:
      MolOps::findSSSR(**iter);
      res<<"$MOL\n";
      res << MolToMolBlock(**iter,true,-1,false);
    }

    for(MOL_SPTR_VECT::const_iterator iter=rxn.beginProductTemplates();
	iter != rxn.endProductTemplates();++iter){
      // to write the mol block, we need ring information:
      MolOps::findSSSR(**iter);
      res<<"$MOL\n";
      res << MolToMolBlock(**iter,true,-1,false);
    }
    return res.str();
  };
#endif
}  
