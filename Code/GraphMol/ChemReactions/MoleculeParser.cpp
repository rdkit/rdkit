// $Id$
//
//  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
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
#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ChemReactions/ReactionParser.h>

namespace{

bool testForSameRXNRoleOfAllMoleculeAtoms(const RDKit::ROMol &mol, int role)
{
  RDKit::ROMol::ATOM_ITER_PAIR atItP = mol.getVertices();
  while(atItP.first != atItP.second ){
    const RDKit::Atom *oAtom=mol[*(atItP.first++)].get();
    if(oAtom->hasProp("molRxnRole") && oAtom->getProp<int>("molRxnRole")!=role){
  	  return false;
    }
  }
  return true;
}

int getRXNRoleOfMolecule(const RDKit::ROMol &mol)
{
  RDKit::ROMol::ATOM_ITER_PAIR atItP = mol.getVertices();
  while(atItP.first != atItP.second ){
    const RDKit::Atom *oAtom=mol[*(atItP.first++)].get();
    if(oAtom->hasProp("molRxnRole")){
  	  return oAtom->getProp<int>("molRxnRole");
    }
  }
  return -1;
}

}

namespace RDKit {
  
  ChemicalReaction * RxnMolToChemicalReaction(const ROMol &mol)
  {
	ChemicalReaction *rxn=new ChemicalReaction();

	MOL_SPTR_VECT fragments = MolOps::getMolFrags(mol);

	unsigned countFragments = 0;
	for(MOL_SPTR_VECT::iterator iter = fragments.begin();
		    iter != fragments.end(); ++iter, countFragments++){
	  int role = getRXNRoleOfMolecule(*iter->get());
      if(!testForSameRXNRoleOfAllMoleculeAtoms(*iter->get(), role)){
    	BOOST_LOG(rdWarningLog)<<">> Atoms within one molecule have different RXN roles.\n";
        continue;
      }
      switch(role){
      case 1:
    	rxn->addReactantTemplate(*iter);
        break;
      case 2:
    	rxn->addProductTemplate(*iter);
        break;
      case 3:
    	rxn->addAgentTemplate(*iter);
        break;
      default:
    	BOOST_LOG(rdWarningLog)<<">> Fragment "<< countFragments <<" not included in the reaction, atoms do not have a correct RXN role.\n";
      }
	}
    return rxn;    
  }

}
