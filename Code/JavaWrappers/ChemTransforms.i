/* 
* $Id$
*
*  Copyright (c) 2011, Novartis Institutes for BioMedical Research Inc.
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
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/Bond.h>
// Fixes annoying compilation namespace issue
typedef RDKit::MatchVectType MatchVectType;

// Partial binding to one overload of RDKit::fragmentOnBonds
RDKit::ROMol *fragmentMolOnBonds(
    const RDKit::ROMol &mol, const std::vector<int> &bondIndices,
    bool addDummies = true, 
    const std::vector<std::pair<unsigned int, unsigned int>> *dummyLabels = nullptr,
    std::vector<int> *nCutsPerAtom = nullptr) {
        // The c# wrapper is unable to wrap vector<unsigned int> so pass in vector<int> and convert
        std::vector<unsigned int> uBondIndices(bondIndices.begin(), bondIndices.end());
        std::vector<unsigned int> uNCutsPerAtom;
        if (nCutsPerAtom) {
            uNCutsPerAtom.insert(uNCutsPerAtom.begin(), nCutsPerAtom->begin(), nCutsPerAtom->end());
        }
        auto cutsPerAtomPtr = nCutsPerAtom ? &uNCutsPerAtom : nullptr;
        // Can't handle the 5th argument of vector<BondType> as swig will not wrap vector of enum
        return RDKit::MolFragmenter::fragmentOnBonds(mol, uBondIndices, addDummies, dummyLabels, nullptr, cutsPerAtomPtr);
 }

// Fix std::unique_ptr issue
RDKit::ROMol * new_molzip(
	       const RDKit::ROMol &a, const RDKit::ROMol &b,
	       const RDKit::MolzipParams &params=RDKit::MolzipParams()) {
  return RDKit::molzip(a, b, params).release();;
}

RDKit::ROMol * new_molzip(
	       const RDKit::ROMol &a,
	       const RDKit::MolzipParams &params=RDKit::MolzipParams()) {
  return RDKit::molzip(a, params).release();
}

RDKit::ROMol * new_molzip(
               std::vector<RDKit::ROMOL_SPTR> &mols,
               const RDKit::MolzipParams &params=RDKit::MolzipParams()) {
  return RDKit::molzip(mols, params).release();
}
%}

%include "std_string.i"
%include "std_vector.i"

%newobject deleteSubstructs;
%newobject replaceSidechains;
%newobject replaceCores;
%newobject MurckoDecompose;
%template(StringMolMap) std::map<std::string,boost::shared_ptr<RDKit::ROMol> >;
%include <GraphMol/Bond.h>
%include <GraphMol/ChemTransforms/ChemTransforms.h>

%newobject fragmentMolOnBonds;
%newobject new_molzip;
%ignore fragmentOnBonds;
%ignore molzip;
%rename("fragmentOnBonds") fragmentMolOnBonds;
%rename("molzip") new_molzip;

RDKit::ROMol *fragmentMolOnBonds(
    const RDKit::ROMol &mol, const std::vector<int> &bondIndices,
    bool addDummies = true, 
    const std::vector<std::pair<unsigned int, unsigned int>> *dummyLabels = nullptr,
    std::vector<int> *nCutsPerAtom = nullptr);

RDKit::ROMol * new_molzip(
			  const RDKit::ROMol &a, const RDKit::ROMol &b,
			  const RDKit::MolzipParams &params=RDKit::MolzipParams());

RDKit::ROMol * new_molzip(
			  const RDKit::ROMol &a,
			  const RDKit::MolzipParams &params=RDKit::MolzipParams());


RDKit::ROMol * new_molzip(
                          std::vector<RDKit::ROMOL_SPTR> &mols,
                          const RDKit::MolzipParams &params=RDKit::MolzipParams());

%ignore fragmentOnSomeBonds;
%ignore constructFragmenterAtomTypes;
%ignore constructBRICSAtomTypes;
%ignore constructFragmenterBondTypes;
%ignore constructBRICSBondTypes;

%newobject fragmentOnBRICSBonds;
%template(UIntMolMap) std::map<unsigned int,boost::shared_ptr<RDKit::ROMol> >;
%include <GraphMol/ChemTransforms/MolFragmenter.h>
