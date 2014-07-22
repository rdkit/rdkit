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

namespace RDKit {
    namespace FMCS {
        class RingMatchTableSet {
            class RingMatchTable {
                FMCS::MatchTable                    MatchMatrix;
                std::map<const INT_VECT*, unsigned> RingIndex;
            public:
                inline void clear() {
                    MatchMatrix.clear();
                    RingIndex.clear();
                }
                inline void resize(unsigned s1, unsigned s2) {
                    MatchMatrix.resize(s1, s2);
                    for(size_t  i = 0; i < s1; i++)
                        for(size_t j = 0; j < s2; j++)
                            MatchMatrix.set(i, j, false);
                }
                inline void makeRingIndex(const ROMol* mol2) {
                    unsigned i=0;
                    // for each TARGET ring
                    const RingInfo::VECT_INT_VECT& rings2 = mol2->getRingInfo()->bondRings();
                    for(RingInfo::VECT_INT_VECT::const_iterator r2 = rings2.begin(); r2 != rings2.end(); r2++)
                        RingIndex[&*r2] = i++;
                }
                inline bool isEqual(unsigned i, const INT_VECT* r2)const {
                    return MatchMatrix.at(i, getRingIndex(r2));
                }
                inline void setMatch(unsigned i, const INT_VECT* r2) {
                    MatchMatrix.set(i, getRingIndex(r2), true);
                }
            private:
                inline unsigned getRingIndex(const INT_VECT* r2)const {
                    std::map<const INT_VECT*, unsigned>::const_iterator j = RingIndex.find(r2);
                    if(RingIndex.end() == j)
                        throw -1;
                    return j->second;
                }
            };

        private:
            std::vector<std::vector<size_t> >* QueryBondRingsIndeces;
            std::map<const ROMol*, std::vector<std::vector<size_t> > > TargetBondRingsIndecesSet; // by target molecules

            std::map<const ROMol*, RingMatchTable> MatchMatrixSet; // by target molecules
            std::map<const INT_VECT*, unsigned>    QueryRingIndex;

        public:
            RingMatchTableSet() : QueryBondRingsIndeces(0) {}

            inline void clear() {
                if(QueryBondRingsIndeces)
                    QueryBondRingsIndeces->clear();
                TargetBondRingsIndecesSet.clear();
                MatchMatrixSet.clear();
                QueryRingIndex.clear();

            }

            inline bool isQueryBondInRing(unsigned bi)const {
                return (*QueryBondRingsIndeces)[bi].empty();
            }
            inline const std::vector<size_t>& getQueryBondRings(unsigned bi)const {
                return (*QueryBondRingsIndeces)[bi];
            }

            inline bool isTargetBondInRing(const ROMol* target, unsigned bi)const {
                std::map<const ROMol*, std::vector<std::vector<size_t> > >::const_iterator
                i = TargetBondRingsIndecesSet.find(target);
                if(TargetBondRingsIndecesSet.end() == i)
                    throw -1; //never
                return i->second[bi].empty();
            }
            inline const std::vector<size_t>& getTargetBondRings(const ROMol* target, unsigned bi)const {
                std::map<const ROMol*, std::vector<std::vector<size_t> > >::const_iterator
                i = TargetBondRingsIndecesSet.find(target);
                if(TargetBondRingsIndecesSet.end() == i)
                    throw -1; //never
                return i->second[bi];
            }

            inline bool isEqual(const INT_VECT* r1, const INT_VECT* r2, const ROMol* mol2)const {
                const RingMatchTable& m = getTargetMatchMatrix(mol2);
                unsigned i = getQueryRingIndex(r1);
                return m.isEqual(i,r2);
            }

            void init(const ROMol* query) {
                MatchMatrixSet.clear();
                // fill out QueryRingIndex
                unsigned i=0;
                const RingInfo::VECT_INT_VECT& rings = query->getRingInfo()->bondRings();
                for(RingInfo::VECT_INT_VECT::const_iterator r = rings.begin(); r != rings.end(); r++)
                    QueryRingIndex[&*r] = i++;
                TargetBondRingsIndecesSet.clear();
                QueryBondRingsIndeces = &TargetBondRingsIndecesSet[query];
                QueryBondRingsIndeces->resize(query->getNumBonds());
                size_t ri = 0;
                for(RingInfo::VECT_INT_VECT::const_iterator r = rings.begin(); r != rings.end(); r++, ri++)
                    for(INT_VECT::const_iterator bi = r->begin(); bi != r->end(); bi++) // all bonds in the ring
                        (*QueryBondRingsIndeces)[*bi].push_back(ri);
            }
            inline void addTargetBondRingsIndeces(const ROMol* mol2) {
                std::vector<std::vector<size_t> >& m = TargetBondRingsIndecesSet[mol2];
                m.resize(mol2->getNumBonds());

                size_t ri = 0;
                const RingInfo::VECT_INT_VECT& rings = mol2->getRingInfo()->bondRings();
                for(RingInfo::VECT_INT_VECT::const_iterator r = rings.begin(); r != rings.end(); r++, ri++)
                    for(INT_VECT::const_iterator bi = r->begin(); bi != r->end(); bi++) // all bonds in the ring
                        m[*bi].push_back(ri);
            }

            void computeRingMatchTable(const ROMol* query, const ROMol* targetMolecule, const MCSParameters& parameters) { // call it for all targets
                const RingInfo::VECT_INT_VECT& rings1 = query->getRingInfo()->bondRings();
                const RingInfo::VECT_INT_VECT& rings2 = targetMolecule->getRingInfo()->bondRings();
                RingMatchTable& m = addTargetMatchMatrix(targetMolecule, rings1.size(), rings2.size());
                unsigned i=0;
                // for each query ring
                for(RingInfo::VECT_INT_VECT::const_iterator r1 = rings1.begin(); r1 != rings1.end(); r1++, i++) {
                    FMCS::Graph graph1;
                    makeRingGraph(graph1, *r1, query); // for each query ring bond ADD all atoms and bonds

                    // for each TARGET ring
                    for(RingInfo::VECT_INT_VECT::const_iterator r2 = rings2.begin(); r2 != rings2.end(); r2++) {
                        if(r1->size() != r2->size())  // rings are different
                            continue;
                        FMCS::Graph graph2;
                        makeRingGraph(graph2, *r2, targetMolecule); // for each TAG ring bond ADD all atoms and bonds

                        // check ring substruct match
                        MCSBondCompareParameters bp = parameters.BondCompareParameters;
                        bp.RingMatchesRingOnly = false;
                        bp.CompleteRingsOnly   = false;
                        bool match =
#ifdef NEVER_xxx_PRECOMPUTED_TABLES_MATCH // not computed yet, because MatchTable computation usees this ring info table
                            FMCS::SubstructMatchCustomTable(graph2, graph1, tag->AtomMatchTable, tag->BondMatchTable);
#else //noticable slowly:
                            FMCS::SubstructMatchCustom( graph2, *targetMolecule , graph1, *query
                                                        , parameters.AtomTyper, parameters.BondTyper, parameters.AtomCompareParameters, bp, NULL);
#endif
                        if(match)
                            m.setMatch(i, &*r2);
                    }
                }
            }

        private:
            void makeRingGraph(FMCS::Graph& g, const INT_VECT& ring, const ROMol* mol)const { // ADD all atoms and bonds
                std::map<const Atom*, unsigned> atomMap;

                for(size_t i=0; i < ring.size(); i++) {
                    const Bond* bond = mol->getBondWithIdx(ring[i]);
                    const Atom* atom1= bond->getBeginAtom();
                    const Atom* atom2= bond->getEndAtom();
                    unsigned j1 = -1;
                    unsigned j2 = -1;
                    std::map<const Atom*, unsigned>::const_iterator ai;
                    ai = atomMap.find(atom1);
                    if(atomMap.end() != ai)
                        j1 = ai->second;
                    ai = atomMap.find(atom2);
                    if(atomMap.end() != ai)
                        j2 = ai->second;
                    if(-1==j1) {
                        j1 = g.m_vertices.size();
                        atomMap[atom1] = j1;
                        g.addAtom(atom1->getIdx());
                    }
                    if(-1==j2) {
                        j2 = g.m_vertices.size();
                        atomMap[atom2] = j2;
                        g.addAtom(atom2->getIdx());
                    }
                    g.addBond(ring[i], j1, j2);
                }
            }

            inline unsigned getQueryRingIndex(const INT_VECT* r1)const {
                std::map<const INT_VECT*, unsigned>::const_iterator i = QueryRingIndex.find(r1);
                if(QueryRingIndex.end() == i)
                    throw -1;   //never
                return i->second;
            }
            inline const RingMatchTable& getTargetMatchMatrix(const ROMol* mol2)const {
                std::map<const ROMol*, RingMatchTable>::const_iterator mi = MatchMatrixSet.find(mol2);
                if(MatchMatrixSet.end() == mi)
                    throw -1;   //never
                return mi->second;
            }

            inline RingMatchTable& addTargetMatchMatrix(const ROMol* mol2, unsigned s1, unsigned s2) {
                RingMatchTable& m = MatchMatrixSet[mol2];
                m.resize(s1, s2);
                m.makeRingIndex(mol2);
                return m;
            }
        };
    }
}   // namespace RDKit
