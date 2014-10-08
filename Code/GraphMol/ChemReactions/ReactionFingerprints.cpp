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

#include <GraphMol/RDKitBase.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include "GraphMol/ChemReactions/ReactionFingerprints.h"
#include <GraphMol/Fingerprints/MACCS.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <DataStructs/SparseIntVect.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>

namespace{

ExplicitBitVect *concatenateBitsets(const ExplicitBitVect& first, const ExplicitBitVect& second)
{
	ExplicitBitVect *res = new ExplicitBitVect(first.size() + second.size());

	for(unsigned i = 0; i < first.size(); i++){
		if(first[i]){
			res->setBit(i);
		}
	}
	for(unsigned j = 0; j < second.size(); j++){
		if(second[j]){
			res->setBit(first.size()+j);
		}
	}
    return res;
}

RDKit::SparseIntVect<boost::uint32_t> *generateFingerprint(RDKit::ROMol &mol,
		unsigned int fpSize,
		RDKit::FingerprintType t)
{
	mol.updatePropertyCache(false);
	RDKit::SparseIntVect<boost::uint32_t> *res;
	switch(t){
	case RDKit::AtomPairFP:
	{
		RDKit::SparseIntVect<boost::int32_t> *tmp1 =
				RDKit::AtomPairs::getHashedAtomPairFingerprint(mol, fpSize);
		res = new RDKit::SparseIntVect<boost::uint32_t>(fpSize);
        BOOST_FOREACH(RDKit::SparseIntVect<boost::int32_t>::StorageType::value_type val,tmp1->getNonzeroElements()){
        	res->setVal(static_cast<boost::uint32_t>(val.first),val.second);
        }
        delete tmp1;
	}
		break;
	case RDKit::TopologicalTorsion:
	{
		RDKit::SparseIntVect<boost::int64_t> *tmp2 =
				RDKit::AtomPairs::getHashedTopologicalTorsionFingerprint(mol, fpSize);
		res = new RDKit::SparseIntVect<boost::uint32_t>(fpSize);
		BOOST_FOREACH(RDKit::SparseIntVect<boost::int64_t>::StorageType::value_type val,tmp2->getNonzeroElements()){
		  res->setVal(static_cast<boost::uint32_t>(val.first),val.second);
		}
		delete tmp2;
	}
		break;
	case RDKit::MorganFP:
	{
		if(!mol.getRingInfo()->isInitialized()){
          mol.updatePropertyCache();
          RDKit::MolOps::findSSSR(mol);
		}
		res = RDKit::MorganFingerprints::getHashedFingerprint(mol, 2,fpSize);
	}
		break;
	default:
		std::stringstream err;
		err << ">> unsupported fingerprint type" << std::endl;
		throw RDKit::ChemicalReactionException(err.str());
	}
	return res;
}

ExplicitBitVect *generateFingerprintAsBitVect(RDKit::ROMol &mol,
		unsigned int fpSize,
		RDKit::FingerprintType t)
{
	mol.updatePropertyCache(false);
	ExplicitBitVect *res;
	switch(t){
	case RDKit::AtomPairFP:
		res = RDKit::AtomPairs::getHashedAtomPairFingerprintAsBitVect(mol, fpSize);
		break;
	case RDKit::RDKitFP:
		res = RDKit::RDKFingerprintMol(mol,1, 7, fpSize);
		break;
	case RDKit::TopologicalTorsion:
		res = RDKit::AtomPairs::getHashedTopologicalTorsionFingerprintAsBitVect(mol);
		break;
	case RDKit::MorganFP:
	{
		if(!mol.getRingInfo()->isInitialized()){
          mol.updatePropertyCache();
          RDKit::MolOps::findSSSR(mol);
		}
		res = RDKit::MorganFingerprints::getFingerprintAsBitVect(mol, 2, fpSize);
	}
		break;
	case RDKit::PatternFP:
		res = RDKit::PatternFingerprintMol(mol, fpSize);
		break;
	default:
		std::stringstream err;
		err << ">> unsupported fingerprint type" << std::endl;
		throw RDKit::ChemicalReactionException(err.str());
	}
	return res;
}

}

namespace RDKit{

SparseIntVect<boost::uint32_t> *generateFingerprintChemReactionAsCountVect(const ChemicalReaction &rxn,
		unsigned int fpSize,
		FingerprintType t,
		ReactionMoleculeType mt)
{
	PRECONDITION(fpSize!=0,"fpSize==0");

	SparseIntVect<boost::uint32_t> *result = new SparseIntVect<boost::uint32_t>(fpSize);
    RDKit::MOL_SPTR_VECT::const_iterator begin = getStartIterator(rxn, mt);
    RDKit::MOL_SPTR_VECT::const_iterator end = getEndIterator(rxn, mt);
    for(; begin != end; ++begin){
		SparseIntVect<boost::uint32_t> *tmp = generateFingerprint(**begin, fpSize, t);
		(*result) += *tmp;
		delete tmp;
	}
	return result;
}

ExplicitBitVect *generateFingerprintChemReactionAsBitVect(const ChemicalReaction &rxn,
		unsigned int fpSize,
		FingerprintType t,
		ReactionMoleculeType mt)
{
	PRECONDITION(fpSize!=0,"fpSize==0");

    ExplicitBitVect *result = new ExplicitBitVect(fpSize);
    RDKit::MOL_SPTR_VECT::const_iterator begin = getStartIterator(rxn, mt);
    RDKit::MOL_SPTR_VECT::const_iterator end = getEndIterator(rxn, mt);
    for(; begin != end; ++begin){
    	ExplicitBitVect *tmp = generateFingerprintAsBitVect(**begin, fpSize, t);
    	(*result) |= *tmp;
    	delete tmp;
    }
	return result;
}

// caller owns the result, it must be deleted
ExplicitBitVect *StructuralFingerprintChemReaction(const ChemicalReaction &rxn,
		const ReactionFingerprintParams &params)
{
	PRECONDITION(params.fpSize!=0,"fpSize==0");

	unsigned int fpSize_final = params.fpSize/2;
	if(params.includeAgents){
	  unsigned agent_fp_size = params.fpSize/3;
      if(params.bitRatioAgents < 1.0){
	    agent_fp_size = int(ceil(static_cast<double>(params.fpSize) * params.bitRatioAgents));
      }
	  unsigned scaling = !(agent_fp_size%2) ? agent_fp_size : agent_fp_size - 1;
	  fpSize_final = (params.fpSize - scaling)/2;
	}
	unsigned int fpSize_agent = params.fpSize - 2*fpSize_final;

    ExplicitBitVect *reactantFP = generateFingerprintChemReactionAsBitVect(rxn, fpSize_final, params.fpType, Reactant);
    ExplicitBitVect *productFP = generateFingerprintChemReactionAsBitVect(rxn, fpSize_final, params.fpType, Product);
    ExplicitBitVect *res = new ExplicitBitVect;
    /* concatenate the two bitvectors */
    (*res) = *reactantFP + *productFP;
    if(fpSize_agent > 0){
      ExplicitBitVect *agentFP = generateFingerprintChemReactionAsBitVect(rxn, fpSize_agent, params.fpType, Agent);
      (*res) += *agentFP;
      delete agentFP;
    }

    delete reactantFP;
    delete productFP;
	return res;
}

SparseIntVect<boost::uint32_t> *DifferenceFingerprintChemReaction(const ChemicalReaction &rxn,
		const ReactionFingerprintParams &params)
{
	PRECONDITION(params.fpSize!=0,"fpSize==0");
	PRECONDITION(params.fpType>0 & params.fpType<4 ,"Fingerprinttype not supported");

	SparseIntVect<boost::uint32_t> *reactantFP = generateFingerprintChemReactionAsCountVect(rxn, params.fpSize, params.fpType, Reactant);
	SparseIntVect<boost::uint32_t> *productFP = generateFingerprintChemReactionAsCountVect(rxn, params.fpSize, params.fpType, Product);
	SparseIntVect<boost::uint32_t>* res = new SparseIntVect<boost::uint32_t>;
	if(params.includeAgents){
	  SparseIntVect<boost::uint32_t> *agentFP = generateFingerprintChemReactionAsCountVect(rxn, params.fpSize, params.fpType, Agent);
	  (*productFP) *= params.nonAgentWeight;
	  (*reactantFP) *= params.nonAgentWeight;
      (*agentFP) *= params.agentWeight;
  	  (*res) = *productFP - *reactantFP + *agentFP;
      delete agentFP;
	}
    /* subtract product FP from reactant FP */
	else{
	  (*res) = *productFP - *reactantFP;
	}

    delete reactantFP;
    delete productFP;
	return res;
}

}
