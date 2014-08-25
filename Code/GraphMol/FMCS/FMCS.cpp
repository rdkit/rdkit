//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <list>
#include <algorithm>
#include <math.h>
#include "SubstructMatchCustom.h"
#include "MaximumCommonSubgraph.h"

namespace RDKit {
    MCSResult findMCS (const std::vector<ROMOL_SPTR>& mols, const MCSParameters* params) {
        MCSParameters p;
        if(0 == params)
            params = &p;
        RDKit::FMCS::MaximumCommonSubgraph fmcs(params);
        return fmcs.find(mols);
    }

    bool MCSProgressCallbackTimeout(const MCSProgressData& stat, const MCSParameters &params, void* userData) {
        unsigned long long* t0 = (unsigned long long*)userData;
        unsigned long long  t  = nanoClock();
        return t - *t0 <= params.Timeout*1000000ULL;
    }

    // PREDEFINED FUNCTORS:

    //=== ATOM COMPARE ========================================================

    bool MCSAtomCompareAny      (const MCSAtomCompareParameters& p, const ROMol& mol1, unsigned int atom1, const ROMol& mol2, unsigned int atom2, void* ) {
        return true;
    }

    bool MCSAtomCompareElements (const MCSAtomCompareParameters& p, const ROMol& mol1, unsigned int atom1, const ROMol& mol2, unsigned int atom2, void* ) {
        const Atom& a1 = *mol1.getAtomWithIdx(atom1);
        const Atom& a2 = *mol2.getAtomWithIdx(atom2);
        if(a1.getAtomicNum() != a2.getAtomicNum())
          return false;
        if(p.MatchValences && a1.getTotalValence() != a2.getTotalValence())
            return false;
        return true;
    }

    bool MCSAtomCompareIsotopes (const MCSAtomCompareParameters& p, const ROMol& mol1, unsigned int atom1, const ROMol& mol2, unsigned int atom2, void* ud) {
        // ignore everything except isotope information:
        //if( ! MCSAtomCompareElements (p, mol1, atom1, mol2, atom2, ud))
        //    return false;
        const Atom& a1 = *mol1.getAtomWithIdx(atom1);
        const Atom& a2 = *mol2.getAtomWithIdx(atom2);
        return a1.getIsotope() == a2.getIsotope();
    }

    //=== BOND COMPARE ========================================================

    class BondMatchOrderMatrix {
        bool MatchMatrix [Bond::OTHER+1] [Bond::OTHER+1];
    public:
        BondMatchOrderMatrix(bool ignoreAromatization) {
            memset(MatchMatrix, 0, sizeof(MatchMatrix));
            for(size_t i=0; i <= Bond::OTHER; i++) { // fill cells of the same and unspecified type
                MatchMatrix[i][i] = true;
                MatchMatrix[Bond::UNSPECIFIED][i] = MatchMatrix[i][Bond::UNSPECIFIED] = true;
                MatchMatrix[Bond::OTHER][i] = MatchMatrix[i][Bond::OTHER] = true;
            }
            if(ignoreAromatization) {
                MatchMatrix[Bond::SINGLE][Bond::AROMATIC] = MatchMatrix[Bond::AROMATIC][Bond::SINGLE] = true;
                MatchMatrix[Bond::SINGLE][Bond::ONEANDAHALF] = MatchMatrix[Bond::ONEANDAHALF][Bond::SINGLE] = true;
                MatchMatrix[Bond::DOUBLE][Bond::TWOANDAHALF] = MatchMatrix[Bond::TWOANDAHALF][Bond::DOUBLE] = true;
                MatchMatrix[Bond::TRIPLE][Bond::THREEANDAHALF] = MatchMatrix[Bond::THREEANDAHALF][Bond::TRIPLE] = true;
                MatchMatrix[Bond::QUADRUPLE][Bond::FOURANDAHALF] = MatchMatrix[Bond::FOURANDAHALF][Bond::QUADRUPLE] = true;
                MatchMatrix[Bond::QUINTUPLE][Bond::FIVEANDAHALF] = MatchMatrix[Bond::FIVEANDAHALF][Bond::QUINTUPLE] = true;
            }
        }
        inline bool isEqual(unsigned i, unsigned j)const {
            return MatchMatrix [i] [j];
        }
    };


    bool checkRingMatch(const MCSBondCompareParameters& p, const ROMol& mol1, unsigned int bond1, const ROMol& mol2, unsigned int bond2, void* v_ringMatchMatrixSet) {
        if(!v_ringMatchMatrixSet)
            throw "v_ringMatchMatrixSet is NULL";   // never
        FMCS::RingMatchTableSet* ringMatchMatrixSet = static_cast<FMCS::RingMatchTableSet*>(v_ringMatchMatrixSet);

        const std::vector<size_t>& ringsIdx1 = ringMatchMatrixSet->getQueryBondRings (bond1); //indeces of rings
        const std::vector<size_t>& ringsIdx2 = ringMatchMatrixSet->getTargetBondRings(&mol2, bond2); //indeces of rings
        bool bond1inRing = !ringsIdx1.empty();
        bool bond2inRing = !ringsIdx2.empty();

        if(bond1inRing != bond2inRing)
            return false;
        if((! bond1inRing) )//the same: && (! bond2inRing))  // both bonds are NOT in rings
            return true;
        // both bonds are in rings
        if(p.CompleteRingsOnly) {
            const RingInfo::VECT_INT_VECT& r1 = mol1.getRingInfo()->bondRings();
            const RingInfo::VECT_INT_VECT& r2 = mol2.getRingInfo()->bondRings();
            // for each query ring contains bond1
            for(std::vector<size_t>::const_iterator r1i = ringsIdx1.begin(); r1i != ringsIdx1.end(); r1i++) {
                const INT_VECT& br1 = r1[*r1i]; // ring contains bond1
                // check all target rings contained bond2
                for(std::vector<size_t>::const_iterator r2i = ringsIdx2.begin(); r2i != ringsIdx2.end(); r2i++) {
                    const INT_VECT& br2 = r2[*r2i]; // ring contains bond2
                    if(br1.size() != br2.size())  // rings are different
                        continue;
                    // compare rings as substructures
                    if(ringMatchMatrixSet->isEqual(&br1, &br2, &mol2))  // EQUAL Rings found
                        return true;
                }
            }
            // all rings are different
            return false;
        } else
            return true; // bond1inRing && bond2inRing;   // both bonds are in rings
    }

    bool MCSBondCompareAny (const MCSBondCompareParameters& p, const ROMol& mol1, unsigned int bond1, const ROMol& mol2, unsigned int bond2, void* ud) {
        if(p.RingMatchesRingOnly)
            return checkRingMatch(p, mol1, bond1, mol2, bond2, ud);
        return true;
    }

    bool MCSBondCompareOrder (const MCSBondCompareParameters& p, const ROMol& mol1, unsigned int bond1, const ROMol& mol2, unsigned int bond2, void* ud) {
        static const BondMatchOrderMatrix match(true); // ignore Aromatization
        const Bond* b1 = mol1.getBondWithIdx(bond1);
        const Bond* b2 = mol2.getBondWithIdx(bond2);
        Bond::BondType t1 = b1->getBondType();
        Bond::BondType t2 = b2->getBondType();
        if(match.isEqual(t1,t2)) {
            if(p.RingMatchesRingOnly)
                return checkRingMatch(p, mol1, bond1, mol2, bond2, ud);
            else
                return true;
        }
        return false;
    }

    bool MCSBondCompareOrderExact (const MCSBondCompareParameters& p, const ROMol& mol1, unsigned int bond1, const ROMol& mol2, unsigned int bond2, void* ud) {
        static const BondMatchOrderMatrix match(false); // AROMATIC != SINGLE
        const Bond* b1 = mol1.getBondWithIdx(bond1);
        const Bond* b2 = mol2.getBondWithIdx(bond2);
        Bond::BondType t1 = b1->getBondType();
        Bond::BondType t2 = b2->getBondType();
        if(match.isEqual(t1,t2)) {
            if(p.RingMatchesRingOnly)
                return checkRingMatch(p, mol1, bond1, mol2, bond2, ud);
            else
                return true;
        }
        return false;
    }

}   // namespace RDKit
