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
#include <GraphMol/QueryOps.h>
%}

%ignore RDKit::RecursiveStructureQuery::d_mutex;
%ignore RDKit::HasPropWithValueQuery::getPair;
%ignore RDKit::HasPropWithValueQueryBase::getPair;
%ignore RDKit::queryIsAtomBridgeheadInternal;
%ignore RDKit::convertComplexNameToQuery2;
%ignore RDKit::queryAtomRingMembership2;
%ignore RDKit::RecursiveStructureQuery2;
%ignore RDKit::queryAtomAromatic2;
%ignore RDKit::queryAtomAliphatic2;
%ignore RDKit::queryAtomExplicitDegree2;
%ignore RDKit::queryAtomTotalDegree2;
%ignore RDKit::queryAtomNonHydrogenDegree2;
%ignore RDKit::queryAtomHeavyAtomDegree2;
%ignore RDKit::queryAtomHCount2;
%ignore RDKit::queryAtomImplicitHCount2;
%ignore RDKit::queryAtomHasImplicitH2;
%ignore RDKit::queryAtomImplicitValence2;
%ignore RDKit::queryAtomExplicitValence2;
%ignore RDKit::queryAtomTotalValence2;
%ignore RDKit::queryAtomUnsaturated2;
%ignore RDKit::queryAtomNum2;
%ignore RDKit::queryAtomType2;
%ignore RDKit::queryAtomMass2;
%ignore RDKit::queryAtomIsotope2;
%ignore RDKit::queryAtomFormalCharge2;
%ignore RDKit::queryAtomNegativeFormalCharge2;
%ignore RDKit::queryAtomHybridization2;
%ignore RDKit::queryAtomNumRadicalElectrons2;
%ignore RDKit::queryAtomHasChiralTag2;
%ignore RDKit::queryAtomMissingChiralTag2;
%ignore RDKit::queryAtomHasHeteroatomNbrs2;
%ignore RDKit::queryAtomNumHeteroatomNbrs2;
%ignore RDKit::queryAtomHasAliphaticHeteroatomNbrs2;
%ignore RDKit::queryAtomNumAliphaticHeteroatomNbrs2;
%ignore RDKit::queryAtomBondProduct2;
%ignore RDKit::queryAtomAllBondProduct2;
%ignore RDKit::queryBondOrder2;
%ignore RDKit::queryBondIsSingleOrAromatic2;
%ignore RDKit::queryBondIsDoubleOrAromatic2;
%ignore RDKit::queryBondIsSingleOrDouble2;
%ignore RDKit::queryBondIsSingleOrDoubleOrAromatic2;
%ignore RDKit::queryBondDir2;
%ignore RDKit::queryIsBondInNRings2;
%ignore RDKit::queryBondHasStereo2;
%ignore RDKit::queryIsAtomInNRings2;
%ignore RDKit::queryIsAtomInRing2;
%ignore RDKit::queryAtomHasRingBond2;
%ignore RDKit::queryAtomMinRingSize2;
%ignore RDKit::queryBondMinRingSize2;
%ignore RDKit::queryAtomRingBondCount2;
%ignore RDKit::queryAtomIsInRingOfSize2;
%ignore RDKit::queryBondIsInRingOfSize2;
%ignore RDKit::makeAtomRangeQuery2;
%ignore RDKit::makeAtomNumQuery2;
%ignore RDKit::makeAtomTypeQuery2;
%ignore RDKit::makeAtomImplicitValenceQuery2;
%ignore RDKit::makeAtomExplicitValenceQuery2;
%ignore RDKit::makeAtomTotalValenceQuery2;
%ignore RDKit::makeAtomExplicitDegreeQuery2;
%ignore RDKit::makeAtomTotalDegreeQuery2;
%ignore RDKit::makeAtomHeavyAtomDegreeQuery2;
%ignore RDKit::makeAtomHCountQuery2;
%ignore RDKit::makeAtomHasImplicitHQuery2;
%ignore RDKit::makeAtomImplicitHCountQuery2;
%ignore RDKit::makeAtomAromaticQuery2;
%ignore RDKit::makeAtomAliphaticQuery2;
%ignore RDKit::makeAtomMassQuery2;
%ignore RDKit::makeAtomIsotopeQuery2;
%ignore RDKit::makeAtomFormalChargeQuery2;
%ignore RDKit::makeAtomNegativeFormalChargeQuery2;
%ignore RDKit::makeAtomHybridizationQuery2;
%ignore RDKit::makeAtomNumRadicalElectronsQuery2;
%ignore RDKit::makeAtomHasChiralTagQuery2;
%ignore RDKit::makeAtomMissingChiralTagQuery2;
%ignore RDKit::makeAtomUnsaturatedQuery2;
%ignore RDKit::makeAtomInRingQuery2;
%ignore RDKit::makeAtomInNRingsQuery2;
%ignore RDKit::makeAtomInRingOfSizeQuery2;
%ignore RDKit::makeAtomMinRingSizeQuery2;
%ignore RDKit::makeAtomRingBondCountQuery2;
%ignore RDKit::makeAAtomQuery2;
%ignore RDKit::makeAHAtomQuery2;
%ignore RDKit::makeQAtomQuery2;
%ignore RDKit::makeQHAtomQuery2;
%ignore RDKit::makeXAtomQuery2;
%ignore RDKit::makeXHAtomQuery2;
%ignore RDKit::makeMAtomQuery2;
%ignore RDKit::makeMHAtomQuery2;
%ignore RDKit::convertComplexNameToQuery2;
%ignore RDKit::makeAtomHasRingBondQuery2;
%ignore RDKit::makeAtomNumHeteroatomNbrsQuery2;
%ignore RDKit::makeAtomHasHeteroatomNbrsQuery2;
%ignore RDKit::makeAtomNumAliphaticHeteroatomNbrsQuery2;
%ignore RDKit::makeAtomHasAliphaticHeteroatomNbrsQuery2;
%ignore RDKit::makeAtomNonHydrogenDegreeQuery2;
%ignore RDKit::makeBondOrderEqualsQuery2;
%ignore RDKit::makeBondDirEqualsQuery2;
%ignore RDKit::QueryOps::replaceAtomWithQueryAtom(RDMol *,atomindex_t);
// %ignore RDKit::;
// %ignore RDKit::;
// %ignore RDKit::;
// %ignore RDKit::;
// %ignore RDKit::;
// %ignore RDKit::;
// %ignore RDKit::;
// %ignore RDKit::;
// %ignore RDKit::;
// %ignore RDKit::;
// %ignore RDKit::;
// %ignore RDKit::;
// %ignore RDKit::;
// %ignore RDKit::;
// %ignore RDKit::;
// %ignore RDKit::;
// %ignore RDKit::;





%include <Query/QueryObjects.h>
%include <Query/EqualityQuery.h>
%include <GraphMol/QueryOps.h>
