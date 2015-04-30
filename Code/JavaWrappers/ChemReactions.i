/* 
* $Id$
*
*  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
*  All rights reserved.
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met: 
*
*     * Redistributions of source code must retain the above copyright 
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above
*       copyright notice, this list of conditions and the following 
*       disclaimer in the documentation and/or other materials provided 
*       with the distribution.
*     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
*       nor the names of its contributors may be used to endorse or promote 
*       products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

%{
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <vector>

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
%}

%include <GraphMol/ChemReactions/Reaction.h>
%include <GraphMol/ChemReactions/ReactionParser.h>
%ignore RDKit::ChemicalReaction::validate(unsigned int &,unsigned int &,bool);
%ignore RDKit::ChemicalReaction::validate(unsigned int &,unsigned int &);
%ignore RDKit::isMoleculeReactantOfReaction(const ChemicalReaction &r,const ROMol &,
                                            unsigned int &);
%ignore RDKit::isMoleculeProductOfReaction(const ChemicalReaction &r,const ROMol &,
                                            unsigned int &);

%extend RDKit::ChemicalReaction {
static RDKit::ChemicalReaction *ReactionFromSmarts(std::string sma,bool useSmiles=false){
  RDKit::ChemicalReaction *res=RDKit::RxnSmartsToChemicalReaction(sma,0,useSmiles);
  if(res) res->initReactantMatchers();
  return res;
};
static RDKit::ChemicalReaction *ReactionFromRxnBlock(std::string block){
  RDKit::ChemicalReaction *res=RDKit::RxnBlockToChemicalReaction(block);
  if(res) res->initReactantMatchers();
  return res;
};
static RDKit::ChemicalReaction *ReactionFromRxnFile(std::string filename){
  RDKit::ChemicalReaction *res=RDKit::RxnFileToChemicalReaction(filename);
  if(res) res->initReactantMatchers();
  return res;
};
static std::string ReactionToSmarts(ChemicalReaction &rxn) {
  return RDKit::ChemicalReactionToRxnSmarts(rxn);
}
static std::string ReactionToRxnBlock (const ChemicalReaction &rxn) {
  return RDKit::ChemicalReactionToRxnBlock(rxn);
};
void compute2DCoordsForReaction(double spacing,
                                  bool updateProps,
                                  bool canonOrient,
                                  unsigned int nFlipsPerSample,
                                  unsigned int nSamples,
                                  int sampleSeed,
                                  bool permuteDeg4Nodes){

       RDDepict::compute2DCoordsForReaction(*($self),
                                  spacing,
                                  updateProps,
                                  canonOrient,
                                  nFlipsPerSample,
                                  nSamples,
                                  sampleSeed,
                                  permuteDeg4Nodes);
};
std::vector<int> ToBinary(){
  std::string sres;
  RDKit::ReactionPickler::pickleReaction(*($self),sres);
  std::vector<int> res(sres.length());
  std::copy(sres.begin(),sres.end(),res.begin());
  return res;
};
static RDKit::ChemicalReaction *RxnFromBinary(std::vector<int> pkl){
  std::string sres;
  sres.resize(pkl.size());
  std::copy(pkl.begin(),pkl.end(),sres.begin());
  RDKit::ChemicalReaction *res=new RDKit::ChemicalReaction(sres);
  return res;
};

/* A Java-accessible validation function */
    std::pair<int,int> *validateReaction(bool silent=false) {
	std::pair<int,int> *res = new std::pair<int,int>();
	// Use some local unsigned ints so that we don't create a new and useless type at the Java
	// level.
	unsigned int first;
	unsigned int second;
	($self)->validate(first, second, silent);
	res->first = (int) first;
	res->second = (int) second;
        return res;
    };

  bool validate() {
    unsigned int nErr=0,nWarn=0;
    bool res=$self->validate(nErr,nWarn);
    return res;
  };

}
