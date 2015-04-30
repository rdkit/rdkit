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
#ifndef RD_REACTIONFINGERPRINTS_H
#define RD_REACTIONFINGERPRINTS_H

#include <DataStructs/SparseIntVect.h>
#include <DataStructs/ExplicitBitVect.h>

namespace RDKit{

  class ChemicalReaction;

  enum FingerprintType{
	  AtomPairFP = 1,
	  TopologicalTorsion,
	  MorganFP,
	  RDKitFP,
	  PatternFP
  };

  //! A struct for storing parameters to manipulate
  //! the calculation of fingerprints of chemical reactions
  /*!
     Different parameters can be chosen to influence the generation
     of chemical reaction fingerprints. Generally different setting
     should be used for structural or difference fingerprints.

     \param includeAgents        include the agents of a reaction for fingerprint generation
     \param bitRatioAgents       in structural fingerprints it determines the ratio of bits of
                                 the agents in the fingerprint
     \param nonAgentWeight       in difference fingerprints weight factor for reactants and products
                                 compared to agents
     \param agentWeight          if agents are included, agents could be weighted compared to reactants
                                 and products in difference fingerprints
     \param fpSize               number of bits of the fingerprint
     \param fpType               kind of fingerprint used, e.g AtompairFP. Be aware that only AtompairFP,
                                 TopologicalTorsion and MorganFP were supported in the difference fingerprint.
   */
  struct ReactionFingerprintParams{

	  ReactionFingerprintParams():
		includeAgents(false),
		bitRatioAgents(0.2),
		nonAgentWeight(10),
		agentWeight(1),
		fpSize(2048),
		fpType(AtomPairFP){}

	  ReactionFingerprintParams(bool includeAgents, double bitRatioAgents,
		  unsigned int nonAgentWeight, int agentWeight,
	      unsigned int fpSize,FingerprintType fpType):
    	    includeAgents(includeAgents),
    	    bitRatioAgents(bitRatioAgents),
    		nonAgentWeight(nonAgentWeight),
    		agentWeight(agentWeight),
			fpSize(fpSize),
			fpType(fpType){}

    bool includeAgents;
    double bitRatioAgents;
    unsigned int nonAgentWeight;
    int agentWeight;
    unsigned int fpSize;
    FingerprintType fpType;
  };

  const ReactionFingerprintParams DefaultStructuralFPParams(true, 0.2, 1, 1, 4096, PatternFP);
  const ReactionFingerprintParams DefaultDifferenceFPParams(true, 0.0, 10, 1, 2048, AtomPairFP);

  //! Generates a structural fingerprint for a reaction
  //! to use in screening
  /*!
     A structural fingerprint is generated as an ExplicitBitVect to use for searching
     e.g. substructure in reactions. By default the fingerprint is generated as 4096 BitVect
     using a PatternFP for reactants and products and tentatively agents which
     were finally  concatenated

    \param rxn:          the reaction to be fingerprinted
    \param params:       specific settings to manipulate fingerprint generation

    \return the reaction fingerprint, as an ExplicitBitVect

    <b>Notes:</b>
      - the caller is responsible for <tt>delete</tt>ing the result
  */
  ExplicitBitVect *StructuralFingerprintChemReaction(const ChemicalReaction &rxn,
		  const ReactionFingerprintParams &params = DefaultStructuralFPParams);


  //! Generates a difference fingerprint for a reaction
  //! to use in similarity search of reactions
  /*!
     A difference fingerprint is generated as a SparseIntVect to use for
     similarity search of reactions. By default the fingerprint is generated as 2048 bit
     hashed fingerprint subtracting AtompairFP of the reactants from the products' AtompairFP
     and tentatively the agent AtompairFP is added

    \param rxn:          the reaction to be fingerprinted
    \param params:       specific settings to manipulate fingerprint generation

    \return the reaction fingerprint, as an SparseIntVec

    <b>Notes:</b>
      - the caller is responsible for <tt>delete</tt>ing the result
  */
  SparseIntVect<boost::uint32_t> *DifferenceFingerprintChemReaction(const ChemicalReaction &rxn,
		  const ReactionFingerprintParams &params = DefaultDifferenceFPParams);

}

#endif
