/* 
* $Id: Bond.i 2519 2013-05-17 03:01:18Z glandrum $
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

%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"
%include "std_pair.i"

%{
#include <RDGeneral/types.h>
#include <GraphMol/Bond.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
%}

%ignore RDKit::Bond::getValenceContrib(const Atom *) const;
%ignore RDKit::Bond::Match(const Bond *) const;
%ignore RDKit::Bond::setBeginAtom(Atom *at);
%ignore RDKit::Bond::setEndAtom(Atom *at);
%ignore RDKit::getTwiceBondType(const RDKit::Bond &b);
%ignore RDKit::Bond::setQuery;
%ignore RDKit::Bond::expandQuery;

%include <GraphMol/Bond.h>

%extend RDKit::Bond {
  std::string getProp(const std::string key){
    std::string res;
    ($self)->getProp(key, res);
    return res;
  }

  /* Methods from MolFileStereoChem.h */
  Bond::BondDir DetermineBondWedgeState(const std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>> &wedgeBonds,
                                        const RDKit::Conformer *conf) {
    return RDKit::DetermineBondWedgeState(($self), wedgeBonds, conf);
  }
  
  /* Based on corresponding methods in Atom.i */
   bool IsInRing(){
    if(!($self)->getOwningMol().getRingInfo()->isInitialized()){
      RDKit::MolOps::findSSSR(($self)->getOwningMol());
    }
    return ($self)->getOwningMol().getRingInfo()->numBondRings(($self)->getIdx())!=0;
  }

  bool IsInRingSize(int size){
    if(!($self)->getOwningMol().getRingInfo()->isInitialized()){
      RDKit::MolOps::findSSSR(($self)->getOwningMol());
    }
    return ($self)->getOwningMol().getRingInfo()->isBondInRingOfSize(($self)->getIdx(),size);
  }

}

