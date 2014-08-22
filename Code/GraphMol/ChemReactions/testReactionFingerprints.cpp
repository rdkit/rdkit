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
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
#include <string>
#include <iostream>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include "GraphMol/ChemReactions/ReactionFingerprints.h"
#include <DataStructs/BitOps.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>
#include <GraphMol/Fingerprints/AtomPairs.h>

using namespace RDKit;

void testStructuralFingerprintsReaction(){

	BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
	BOOST_LOG(rdErrorLog) << "    Test Reaction StructuralFingerprint" << std::endl;
	{
		std::string reaction, reactionq;
		ChemicalReaction* rxn, *rxnq;

		reaction = "C1CCCCC1>>C1CCNCC1";
		reactionq = "C1CCCCC1>>C1CCNCC1";

		rxn = RxnSmartsToChemicalReaction(reaction,0,true);
		rxnq = RxnSmartsToChemicalReaction(reactionq,0,true);
		TEST_ASSERT(rxn);
		TEST_ASSERT(rxnq);
		ReactionFingerprintParams params;
		params.fpType = PatternFP;
		params.fpSize = 4096;

		ExplicitBitVect *rxnFP = StructuralFingerprintChemReaction(*rxn, params);
		ExplicitBitVect *rxnqFP = StructuralFingerprintChemReaction(*rxnq, params);
		TEST_ASSERT(AllProbeBitsMatch(*rxnqFP,*rxnFP));

		TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
		TEST_ASSERT(rxn->getNumProductTemplates() == 1);
		TEST_ASSERT(rxnq->getNumReactantTemplates() == 1);
		TEST_ASSERT(rxnq->getNumProductTemplates() == 1);

		MOL_SPTR_VECT::const_iterator reacts_iter = rxn->beginReactantTemplates();
		MOL_SPTR_VECT::const_iterator reacts_iterq = rxnq->beginReactantTemplates();
		MOL_SPTR_VECT::const_iterator products_iter = rxn->beginProductTemplates();
		MOL_SPTR_VECT::const_iterator products_iterq = rxnq->beginProductTemplates();

	    MatchVectType mv;
		TEST_ASSERT(SubstructMatch(*reacts_iter->get(),*reacts_iterq->get(),mv))
	    TEST_ASSERT(SubstructMatch(*products_iter->get(),*products_iterq->get(),mv))

		delete rxn;
		delete rxnq;
		delete rxnFP;
		delete rxnqFP;
	}
	{
		std::string reaction, reactionq;
		ChemicalReaction* rxn, *rxnq;

		reaction = "C1CCCCC1>>C1CCNCC1";
		reactionq = "C1CCCCC1>>C1CCOCC1";

		rxn = RxnSmartsToChemicalReaction(reaction,0,true);
		rxnq = RxnSmartsToChemicalReaction(reactionq,0,true);
		TEST_ASSERT(rxn);
		TEST_ASSERT(rxnq);

		ReactionFingerprintParams params;
		params.fpType = PatternFP;
		params.fpSize = 4096;

		ExplicitBitVect *rxnFP = StructuralFingerprintChemReaction(*rxn, params);
		ExplicitBitVect *rxnqFP = StructuralFingerprintChemReaction(*rxnq, params);

		TEST_ASSERT(!AllProbeBitsMatch(*rxnqFP,*rxnFP));

		TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
		TEST_ASSERT(rxn->getNumProductTemplates() == 1);
		TEST_ASSERT(rxnq->getNumReactantTemplates() == 1);
		TEST_ASSERT(rxnq->getNumProductTemplates() == 1);

		MOL_SPTR_VECT::const_iterator reacts_iter = rxn->beginReactantTemplates();
		MOL_SPTR_VECT::const_iterator reacts_iterq = rxnq->beginReactantTemplates();
		MOL_SPTR_VECT::const_iterator products_iter = rxn->beginProductTemplates();
		MOL_SPTR_VECT::const_iterator products_iterq = rxnq->beginProductTemplates();

	    MatchVectType mv;
		TEST_ASSERT(SubstructMatch(*reacts_iter->get(),*reacts_iterq->get(),mv))
	    TEST_ASSERT(!SubstructMatch(*products_iter->get(),*products_iterq->get(),mv))

		delete rxn;
		delete rxnq;
		delete rxnFP;
		delete rxnqFP;
	}
	{
		std::string reaction, reactionq, reactionq2;
		ChemicalReaction* rxn, *rxnq, *rxnq2;

		reaction = "C1CCCCC1>>C1CCNCC1";
		reactionq = ">>C1CCNCC1";
		reactionq2 = ">>C1CCOCC1";

		rxn = RxnSmartsToChemicalReaction(reaction,0,true);
		rxnq = RxnSmartsToChemicalReaction(reactionq,0,true);
		rxnq2 = RxnSmartsToChemicalReaction(reactionq2,0,true);
		TEST_ASSERT(rxn);
		TEST_ASSERT(rxnq);
		TEST_ASSERT(rxnq2);

		ReactionFingerprintParams params;
		params.fpType = PatternFP;
		params.fpSize = 4096;
		ExplicitBitVect *rxnFP = StructuralFingerprintChemReaction(*rxn, params);
		ExplicitBitVect *rxnqFP = StructuralFingerprintChemReaction(*rxnq, params);
		ExplicitBitVect *rxnq2FP = StructuralFingerprintChemReaction(*rxnq2, params);

		TEST_ASSERT(AllProbeBitsMatch(*rxnqFP,*rxnFP));
		TEST_ASSERT(!AllProbeBitsMatch(*rxnq2FP,*rxnFP));

		TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
		TEST_ASSERT(rxn->getNumProductTemplates() == 1);
		TEST_ASSERT(rxnq->getNumReactantTemplates() == 0);
		TEST_ASSERT(rxnq->getNumProductTemplates() == 1);
		TEST_ASSERT(rxnq2->getNumReactantTemplates() == 0);
		TEST_ASSERT(rxnq2->getNumProductTemplates() == 1);

		MOL_SPTR_VECT::const_iterator products_iter = rxn->beginProductTemplates();
		MOL_SPTR_VECT::const_iterator products_iterq = rxnq->beginProductTemplates();
		MOL_SPTR_VECT::const_iterator products_iterq2 = rxnq2->beginProductTemplates();

	    MatchVectType mv;
	    TEST_ASSERT(SubstructMatch(*products_iter->get(),*products_iterq->get(),mv))
	    TEST_ASSERT(!SubstructMatch(*products_iter->get(),*products_iterq2->get(),mv))

		delete rxn;
		delete rxnq;
		delete rxnq2;
		delete rxnFP;
		delete rxnqFP;
		delete rxnq2FP;
	}
	{
		std::string reaction, reactionq, reactionq2;
		ChemicalReaction* rxn, *rxnq, *rxnq2;

		reaction = "C1CCCCC1>>C1CCNCC1";
		reactionq = "C1CCCCC1>>";
		reactionq2 = "C1CCOCC1>>";

		rxn = RxnSmartsToChemicalReaction(reaction,0,true);
		rxnq = RxnSmartsToChemicalReaction(reactionq,0,true);
		rxnq2 = RxnSmartsToChemicalReaction(reactionq2,0,true);
		TEST_ASSERT(rxn);
		TEST_ASSERT(rxnq);
		TEST_ASSERT(rxnq2);

		ReactionFingerprintParams params;
		params.fpType = PatternFP;
		params.fpSize = 4096;

		ExplicitBitVect *rxnFP = StructuralFingerprintChemReaction(*rxn, params);
		ExplicitBitVect *rxnqFP = StructuralFingerprintChemReaction(*rxnq, params);
		ExplicitBitVect *rxnq2FP = StructuralFingerprintChemReaction(*rxnq2, params);

		TEST_ASSERT(AllProbeBitsMatch(*rxnqFP,*rxnFP));
		TEST_ASSERT(!AllProbeBitsMatch(*rxnq2FP,*rxnFP));

		TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
		TEST_ASSERT(rxn->getNumProductTemplates() == 1);
		TEST_ASSERT(rxnq->getNumReactantTemplates() == 1);
		TEST_ASSERT(rxnq->getNumProductTemplates() == 0);
		TEST_ASSERT(rxnq2->getNumReactantTemplates() == 1);
		TEST_ASSERT(rxnq2->getNumProductTemplates() == 0);

		MOL_SPTR_VECT::const_iterator react_iter = rxn->beginReactantTemplates();
		MOL_SPTR_VECT::const_iterator react_iterq = rxnq->beginReactantTemplates();
		MOL_SPTR_VECT::const_iterator react_iterq2 = rxnq2->beginReactantTemplates();

	    MatchVectType mv;
	    TEST_ASSERT(SubstructMatch(*react_iter->get(),*react_iterq->get(),mv))
	    TEST_ASSERT(!SubstructMatch(*react_iter->get(),*react_iterq2->get(),mv))

		delete rxn;
		delete rxnq;
		delete rxnq2;
		delete rxnFP;
		delete rxnqFP;
		delete rxnq2FP;
	}
	{
		std::string reaction, reactionq, reactionq2;
		ChemicalReaction* rxn, *rxnq, *rxnq2;

		reaction = "C1CCCCC1>>C1CCNCC1";
		reactionq = "CCC>>CNC";
		reactionq2 = "CCCCC>>CCCCN";

		rxn = RxnSmartsToChemicalReaction(reaction,0,true);
		rxnq = RxnSmartsToChemicalReaction(reactionq,0,true);
		rxnq2 = RxnSmartsToChemicalReaction(reactionq2,0,true);
		TEST_ASSERT(rxn);
		TEST_ASSERT(rxnq);
		TEST_ASSERT(rxnq2);

		ReactionFingerprintParams params;
		params.fpType = PatternFP;
		params.fpSize = 4096;

		ExplicitBitVect *rxnFP = StructuralFingerprintChemReaction(*rxn, params);
		ExplicitBitVect *rxnqFP = StructuralFingerprintChemReaction(*rxnq, params);
		ExplicitBitVect *rxnq2FP = StructuralFingerprintChemReaction(*rxnq2, params);

		TEST_ASSERT(AllProbeBitsMatch(*rxnqFP,*rxnFP));
		TEST_ASSERT(AllProbeBitsMatch(*rxnq2FP,*rxnFP));

		TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
		TEST_ASSERT(rxn->getNumProductTemplates() == 1);
		TEST_ASSERT(rxnq->getNumReactantTemplates() == 1);
		TEST_ASSERT(rxnq->getNumProductTemplates() == 1);
		TEST_ASSERT(rxnq2->getNumReactantTemplates() == 1);
		TEST_ASSERT(rxnq2->getNumProductTemplates() == 1);

		MOL_SPTR_VECT::const_iterator react_iter = rxn->beginReactantTemplates();
		MOL_SPTR_VECT::const_iterator react_iterq = rxnq->beginReactantTemplates();
		MOL_SPTR_VECT::const_iterator react_iterq2 = rxnq2->beginReactantTemplates();
		MOL_SPTR_VECT::const_iterator products_iter = rxn->beginProductTemplates();
		MOL_SPTR_VECT::const_iterator products_iterq = rxnq->beginProductTemplates();
		MOL_SPTR_VECT::const_iterator products_iterq2 = rxnq2->beginProductTemplates();

	    MatchVectType mv;
	    TEST_ASSERT(SubstructMatch(*react_iter->get(),*react_iterq->get(),mv))
	    TEST_ASSERT(SubstructMatch(*react_iter->get(),*react_iterq2->get(),mv))
	    TEST_ASSERT(SubstructMatch(*products_iter->get(),*products_iterq->get(),mv))
	    TEST_ASSERT(SubstructMatch(*products_iter->get(),*products_iterq2->get(),mv))

		delete rxn;
		delete rxnq;
		delete rxnq2;
		delete rxnFP;
		delete rxnqFP;
		delete rxnq2FP;
	}
	{
		std::string reaction, reactionq, reactionq2;
		ChemicalReaction* rxn, *rxnq, *rxnq2;
	    unsigned int nWarn,nError,which;
		reaction = "C1CCCCC1>>C1CCNCC1";
		reactionq = "CCC>>CNC";
		reactionq2 = "CCCCC>>CCCCN";

		rxn = RxnSmartsToChemicalReaction(reaction,0,true);
		rxnq = RxnSmartsToChemicalReaction(reactionq,0,true);
		rxnq2 = RxnSmartsToChemicalReaction(reactionq2,0,true);
		TEST_ASSERT(rxn);
		TEST_ASSERT(rxnq);
		TEST_ASSERT(rxnq2);

		ReactionFingerprintParams params;
		params.fpType = PatternFP;
		params.fpSize = 4096;

		ExplicitBitVect *rxnFP = StructuralFingerprintChemReaction(*rxn, params);
		ExplicitBitVect *rxnqFP = StructuralFingerprintChemReaction(*rxnq, params);
		ExplicitBitVect *rxnq2FP = StructuralFingerprintChemReaction(*rxnq2, params);

		TEST_ASSERT(AllProbeBitsMatch(*rxnqFP,*rxnFP));
		TEST_ASSERT(AllProbeBitsMatch(*rxnq2FP,*rxnFP));

	    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *rxnq));
        TEST_ASSERT(hasReactionSubstructMatch(*rxn, *rxnq2));

		delete rxn;
		delete rxnq;
		delete rxnq2;
		delete rxnFP;
		delete rxnqFP;
		delete rxnq2FP;
	}
	{
		std::string reaction, reactionq, reactionq2;
		ChemicalReaction* rxn, *rxnq, *rxnq2;

		reaction = "C1CCCCC1>>C1CCNCC1";
		reactionq = "C1CCCCC1>>C1CCNCC1";
		reactionq2 = "C1CCCCC1>>C1CCOCC1";

		rxn = RxnSmartsToChemicalReaction(reaction,0,true);
		rxnq = RxnSmartsToChemicalReaction(reactionq,0,true);
		rxnq2 = RxnSmartsToChemicalReaction(reactionq2,0,true);
		TEST_ASSERT(rxn);
		TEST_ASSERT(rxnq);
		TEST_ASSERT(rxnq2);

		ReactionFingerprintParams params;
		params.fpType = MorganFP;
		params.fpSize = 4096;

		ExplicitBitVect *rxnFP = StructuralFingerprintChemReaction(*rxn, params);
		ExplicitBitVect *rxnqFP = StructuralFingerprintChemReaction(*rxnq, params);
		ExplicitBitVect *rxnq2FP = StructuralFingerprintChemReaction(*rxnq2, params);

		TEST_ASSERT(AllProbeBitsMatch(*rxnqFP,*rxnFP));
		TEST_ASSERT(!AllProbeBitsMatch(*rxnq2FP,*rxnFP));

		delete rxn;
		delete rxnq;
		delete rxnq2;
		delete rxnFP;
		delete rxnqFP;
		delete rxnq2FP;
	}
	{
		std::string reaction, reactionq, reactionq2;
		ChemicalReaction* rxn, *rxnq, *rxnq2;

		reaction = "C1CCCCC1>>C1CCNCC1";
		reactionq = "C1CCCCC1>>C1CCNCC1";
		reactionq2 = "C1CCCCC1>>C1CCOCC1";

		rxn = RxnSmartsToChemicalReaction(reaction,0,true);
		rxnq = RxnSmartsToChemicalReaction(reactionq,0,true);
		rxnq2 = RxnSmartsToChemicalReaction(reactionq2,0,true);
		TEST_ASSERT(rxn);
		TEST_ASSERT(rxnq);
		TEST_ASSERT(rxnq2);

		ReactionFingerprintParams params;
		params.fpSize = 4096;

		ExplicitBitVect *rxnFP = StructuralFingerprintChemReaction(*rxn, params);
		ExplicitBitVect *rxnqFP = StructuralFingerprintChemReaction(*rxnq, params);
		ExplicitBitVect *rxnq2FP = StructuralFingerprintChemReaction(*rxnq2, params);

		TEST_ASSERT(AllProbeBitsMatch(*rxnqFP,*rxnFP));
		TEST_ASSERT(!AllProbeBitsMatch(*rxnq2FP,*rxnFP));

		delete rxn;
		delete rxnq;
		delete rxnq2;
		delete rxnFP;
		delete rxnqFP;
		delete rxnq2FP;
	}
	{
		std::string reaction, reactionq, reactionq2;
		ChemicalReaction* rxn, *rxnq, *rxnq2;

		reaction = "C1CCCCC1>>C1CCNCC1";
		reactionq = "C1CCCCC1>>C1CCNCC1";
		reactionq2 = "C1CCCCC1>>C1CCOCC1";

		rxn = RxnSmartsToChemicalReaction(reaction,0,true);
		rxnq = RxnSmartsToChemicalReaction(reactionq,0,true);
		rxnq2 = RxnSmartsToChemicalReaction(reactionq2,0,true);
		TEST_ASSERT(rxn);
		TEST_ASSERT(rxnq);
		TEST_ASSERT(rxnq2);

		ReactionFingerprintParams params;
		params.fpType = TopologicalTorsion;
		params.fpSize = 4096;

		ExplicitBitVect *rxnFP = StructuralFingerprintChemReaction(*rxn, params);
		ExplicitBitVect *rxnqFP = StructuralFingerprintChemReaction(*rxnq, params);
		ExplicitBitVect *rxnq2FP = StructuralFingerprintChemReaction(*rxnq2, params);

		TEST_ASSERT(AllProbeBitsMatch(*rxnqFP,*rxnFP));
		TEST_ASSERT(!AllProbeBitsMatch(*rxnq2FP,*rxnFP));

		delete rxn;
		delete rxnq;
		delete rxnq2;
		delete rxnFP;
		delete rxnqFP;
		delete rxnq2FP;
	}
	{
		std::string reaction, reactionq, reactionq2;
		ChemicalReaction* rxn, *rxnq, *rxnq2;

		reaction = "C1CCCCC1>>C1CCNCC1";
		reactionq = "C1CCCCC1>>C1CCNCC1";
		reactionq2 = "C1CCCCC1>>C1CCOCC1";

		rxn = RxnSmartsToChemicalReaction(reaction,0,true);
		rxnq = RxnSmartsToChemicalReaction(reactionq,0,true);
		rxnq2 = RxnSmartsToChemicalReaction(reactionq2,0,true);
		TEST_ASSERT(rxn);
		TEST_ASSERT(rxnq);
		TEST_ASSERT(rxnq2);

		ReactionFingerprintParams params;
		params.fpType = RDKitFP;
		params.fpSize = 4096;

		ExplicitBitVect *rxnFP = StructuralFingerprintChemReaction(*rxn, params);
		ExplicitBitVect *rxnqFP = StructuralFingerprintChemReaction(*rxnq, params);
		ExplicitBitVect *rxnq2FP = StructuralFingerprintChemReaction(*rxnq2, params);

		TEST_ASSERT(AllProbeBitsMatch(*rxnqFP,*rxnFP));
		TEST_ASSERT(!AllProbeBitsMatch(*rxnq2FP,*rxnFP));

		delete rxn;
		delete rxnq;
		delete rxnq2;
		delete rxnFP;
		delete rxnqFP;
		delete rxnq2FP;
	}
	{
		std::string reaction, reactionq, reactionq2;
		ChemicalReaction* rxn, *rxnq, *rxnq2;

		reaction = "C1CCCCC1>>C1CCNCC1";
		reactionq = "C1CCCCC1>>C1CCNCC1";
		reactionq2 = "C1CCCCC1>>C1CCOCC1";

		rxn = RxnSmartsToChemicalReaction(reaction,0,true);
		rxnq = RxnSmartsToChemicalReaction(reactionq,0,true);
		rxnq2 = RxnSmartsToChemicalReaction(reactionq2,0,true);
		TEST_ASSERT(rxn);
		TEST_ASSERT(rxnq);
		TEST_ASSERT(rxnq2);

		ReactionFingerprintParams params;
		params.fpType = RDKitFP;
		params.fpSize = 4096;

		ExplicitBitVect *rxnFP = StructuralFingerprintChemReaction(*rxn, params);
		ExplicitBitVect *rxnqFP = StructuralFingerprintChemReaction(*rxnq, params);
		ExplicitBitVect *rxnq2FP = StructuralFingerprintChemReaction(*rxnq2, params);

		TEST_ASSERT(AllProbeBitsMatch(*rxnqFP,*rxnFP));
		TEST_ASSERT(!AllProbeBitsMatch(*rxnq2FP,*rxnFP));

		delete rxn;
		delete rxnq;
		delete rxnq2;
		delete rxnFP;
		delete rxnqFP;
		delete rxnq2FP;
	}
	{
		std::string reaction, reactionq, reactionq2;
		ChemicalReaction* rxn, *rxnq, *rxnq2;

		reaction = "C1CCCCC1>C(=O)O.[Na]>C1CCNCC1";
		reactionq = "C1CCCCC1>C(=O)O.[Na]>C1CCNCC1";
		reactionq2 = "C1CCCCC1>>C1CCNCC1";

		rxn = RxnSmartsToChemicalReaction(reaction,0,true);
		rxnq = RxnSmartsToChemicalReaction(reactionq,0,true);
		rxnq2 = RxnSmartsToChemicalReaction(reactionq2,0,true);
		TEST_ASSERT(rxn);
		TEST_ASSERT(rxnq);
		TEST_ASSERT(rxnq2);

		ExplicitBitVect *rxnFP = StructuralFingerprintChemReaction(*rxn, DefaultStructuralFPParams);
		ExplicitBitVect *rxnqFP = StructuralFingerprintChemReaction(*rxnq, DefaultStructuralFPParams);
		ExplicitBitVect *rxnq2FP = StructuralFingerprintChemReaction(*rxnq2, DefaultStructuralFPParams);

		TEST_ASSERT(AllProbeBitsMatch(*rxnqFP,*rxnFP));
		TEST_ASSERT(AllProbeBitsMatch(*rxnq2FP,*rxnFP));

		delete rxn;
		delete rxnq;
		delete rxnq2;
		delete rxnFP;
		delete rxnqFP;
		delete rxnq2FP;
	}
	{
		std::string reaction, reactionq, reactionq2;
		ChemicalReaction* rxn, *rxnq, *rxnq2;

		reaction = "C1CCCCC1>[Na]>C1CCNCC1";
		reactionq = "C1CCCCC1>C(=O)O.[Na]>C1CCNCC1";
		reactionq2 = "C1CCCCC1>>C1CCNCC1";

		rxn = RxnSmartsToChemicalReaction(reaction,0,true);
		rxnq = RxnSmartsToChemicalReaction(reactionq,0,true);
		rxnq2 = RxnSmartsToChemicalReaction(reactionq2,0,true);
		TEST_ASSERT(rxn);
		TEST_ASSERT(rxnq);
		TEST_ASSERT(rxnq2);

		ExplicitBitVect *rxnFP = StructuralFingerprintChemReaction(*rxn, DefaultStructuralFPParams);
		ExplicitBitVect *rxnqFP = StructuralFingerprintChemReaction(*rxnq, DefaultStructuralFPParams);
		ExplicitBitVect *rxnq2FP = StructuralFingerprintChemReaction(*rxnq2, DefaultStructuralFPParams);

		TEST_ASSERT(!AllProbeBitsMatch(*rxnqFP,*rxnFP));
		TEST_ASSERT(AllProbeBitsMatch(*rxnq2FP,*rxnFP));

		delete rxn;
		delete rxnq;
		delete rxnq2;
		delete rxnFP;
		delete rxnqFP;
		delete rxnq2FP;
	}
}

void testDifferenceFingerprintsReaction(){

	BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
	BOOST_LOG(rdErrorLog) << "    Test Reaction DifferenceFingerprint" << std::endl;
	{
		std::string reaction1, reaction2;
		ChemicalReaction* rxn1, *rxn2;

		reaction1 = "C1CCCCC1>>C1CCCCC1";
		reaction2 = "C1CCCCC1>>C1CCNCC1";

		rxn1 = RxnSmartsToChemicalReaction(reaction1,0,true);
		rxn2 = RxnSmartsToChemicalReaction(reaction2,0,true);
		TEST_ASSERT(rxn1);
		TEST_ASSERT(rxn2);

		SparseIntVect<boost::uint32_t> *rxn1FP = DifferenceFingerprintChemReaction(*rxn1);
		SparseIntVect<boost::uint32_t> *rxn2FP = DifferenceFingerprintChemReaction(*rxn2);

		TEST_ASSERT(TanimotoSimilarity(*rxn1FP,*rxn2FP) == 0.0);;

		delete rxn1;
		delete rxn2;
		delete rxn1FP;
		delete rxn2FP;
	}
	{
		std::string reaction1, reaction2;
		ChemicalReaction* rxn1, *rxn2;

		reaction1 = "C1CCCCC1>>C1CCOCC1";
		reaction2 = "C1CCCCC1>>C1CCNCC1";

		rxn1 = RxnSmartsToChemicalReaction(reaction1,0,true);
		rxn2 = RxnSmartsToChemicalReaction(reaction2,0,true);
		TEST_ASSERT(rxn1);
		TEST_ASSERT(rxn2);

		SparseIntVect<boost::uint32_t> *rxn1FP = DifferenceFingerprintChemReaction(*rxn1);
		SparseIntVect<boost::uint32_t> *rxn2FP = DifferenceFingerprintChemReaction(*rxn2);

		TEST_ASSERT(TanimotoSimilarity(*rxn1FP,*rxn2FP) > 0.0);
	    TEST_ASSERT(TanimotoSimilarity(*rxn1FP,*rxn2FP) <= 1.0);

		delete rxn1;
		delete rxn2;
		delete rxn1FP;
		delete rxn2FP;
	}
	{
		std::string reaction1, reaction2;
		ChemicalReaction* rxn1, *rxn2;

		reaction1 = "c1ccccc1>>c1ccncn1";
		reaction2 = "c1ccccc1>>c1ccncc1";
		rxn1 = RxnSmartsToChemicalReaction(reaction1,0,true);
		rxn2 = RxnSmartsToChemicalReaction(reaction2,0,true);
		TEST_ASSERT(rxn1);
		TEST_ASSERT(rxn2);

		SparseIntVect<boost::uint32_t> *rxn1FP = DifferenceFingerprintChemReaction(*rxn1);
		SparseIntVect<boost::uint32_t> *rxn2FP = DifferenceFingerprintChemReaction(*rxn2);

		TEST_ASSERT(TanimotoSimilarity(*rxn1FP,*rxn2FP) > 0.0);
	    TEST_ASSERT(TanimotoSimilarity(*rxn1FP,*rxn2FP) <= 1.0);

		delete rxn1;
		delete rxn2;
		delete rxn1FP;
		delete rxn2FP;
	}
	{
		std::string reaction1, reaction2;
		ChemicalReaction* rxn1, *rxn2;

		reaction1 = "c1ccccc1>>c1ccncn1";
		reaction2 = "c1ccccc1>>c1ccncc1";
		rxn1 = RxnSmartsToChemicalReaction(reaction1,0,true);
		rxn2 = RxnSmartsToChemicalReaction(reaction2,0,true);
		TEST_ASSERT(rxn1);
		TEST_ASSERT(rxn2);

		ReactionFingerprintParams params;
		params.fpType = MorganFP;

		SparseIntVect<boost::uint32_t> *rxn1FP = DifferenceFingerprintChemReaction(*rxn1, params);
		SparseIntVect<boost::uint32_t> *rxn2FP = DifferenceFingerprintChemReaction(*rxn2, params);

		TEST_ASSERT(TanimotoSimilarity(*rxn1FP,*rxn2FP) > 0.0);
	    TEST_ASSERT(TanimotoSimilarity(*rxn1FP,*rxn2FP) <= 1.0);

		delete rxn1;
		delete rxn2;
		delete rxn1FP;
		delete rxn2FP;
	}
	{
		std::string reaction1, reaction2;
		ChemicalReaction* rxn1, *rxn2;

		reaction1 = "c1ccccc1>>c1ccncn1";
		reaction2 = "c1ccccc1>>c1ccncc1";
		rxn1 = RxnSmartsToChemicalReaction(reaction1,0,true);
		rxn2 = RxnSmartsToChemicalReaction(reaction2,0,true);
		TEST_ASSERT(rxn1);
		TEST_ASSERT(rxn2);

		ReactionFingerprintParams params;
		params.fpType = TopologicalTorsion;

		SparseIntVect<boost::uint32_t> *rxn1FP = DifferenceFingerprintChemReaction(*rxn1, params);
		SparseIntVect<boost::uint32_t> *rxn2FP = DifferenceFingerprintChemReaction(*rxn2, params);

		TEST_ASSERT(TanimotoSimilarity(*rxn1FP,*rxn2FP) > 0.0);
	    TEST_ASSERT(TanimotoSimilarity(*rxn1FP,*rxn2FP) <= 1.0);

		delete rxn1;
		delete rxn2;
		delete rxn1FP;
		delete rxn2FP;
	}
}


int main() { 
  RDLog::InitLogs();
    
  BOOST_LOG(rdInfoLog) << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Testing Chemical Reaction Fingerprints \n";

  testStructuralFingerprintsReaction();
  testDifferenceFingerprintsReaction();

  BOOST_LOG(rdInfoLog) << "*******************************************************\n";
  return(0);
}


