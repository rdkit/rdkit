//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <map>
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include "SubstructMatchCustom.h"

#include "../Substruct/vf2.hpp"

namespace RDKit {
    namespace FMCS {
        class MolMatchFinalCheckFunctor {   // always TRUE. We DO NOT check Chirality
        public:
            MolMatchFinalCheckFunctor(const ROMol &mol1, const ROMol &mol2)
            {};
            bool operator()(const boost::detail::node_id[], const boost::detail::node_id[])const {
                return true;
            }
        };
        //=================================================================================================
        // PRECOMPUTED_TABLES_MATCH much faster in overall even with very simple compare functions

        class AtomTableCompareFunctor {
            const FMCS::Graph&      QueryTopology;
            const FMCS::Graph&      TargetTopology;
            const FMCS::MatchTable& MatchTable;
        public:
            AtomTableCompareFunctor(const FMCS::Graph& query, const FMCS::Graph& target, const FMCS::MatchTable& targetMatch)
                : QueryTopology(query), TargetTopology(target), MatchTable(targetMatch) {};
            bool operator()(unsigned int i,unsigned int j)const {
                return MatchTable.at(QueryTopology[i], TargetTopology[j]);
            }
        };

        class BondTableCompareFunctor {
            const FMCS::Graph& QueryTopology;
            const FMCS::Graph& TargetTopology;
            const FMCS::MatchTable& MatchTable;
        public:
            BondTableCompareFunctor(const FMCS::Graph& query, const FMCS::Graph& target, const FMCS::MatchTable& targetMatch)
                : QueryTopology(query), TargetTopology(target), MatchTable(targetMatch) {};
            bool operator()(FMCS::Graph::edge_descriptor i, FMCS::Graph::edge_descriptor j)const {
                return MatchTable.at(QueryTopology[i], TargetTopology[j]);
            }
        };

        bool SubstructMatchCustomTable(const FMCS::Graph& target, const FMCS::Graph& query, const MatchTable& atomMatchTable, const MatchTable& bondMatchTable, match_V_t* match) {
            if(query.m_vertices.size() > target.m_vertices.size()   // query > target
                    ||query.m_edges.size() > target.m_edges.size())
                return false;
            ROMol mo;
            MolMatchFinalCheckFunctor mc(mo, mo); // unused dummy item, just required by vf2() external implementation

            AtomTableCompareFunctor ac(query, target, atomMatchTable);
            BondTableCompareFunctor bc(query, target, bondMatchTable);

            match_V_t dummy_match;
            if(!match)
                match = &dummy_match;
            return boost::vf2(query, target, ac, bc, mc, *match);
        }

        //=========================================================================
        // slow implementation with absolutely the same functionality
        //=========================================================================

        class AtomLabelFunctor {
            const FMCS::Graph& QueryTopology;
            const FMCS::Graph& TargetTopology;
            const ROMol &d_query;
            const ROMol &d_mol;
            MCSAtomCompareFunction AtomCompare;
            const MCSAtomCompareParameters& Parameters;
            void* UserData;
        public:
            AtomLabelFunctor( const FMCS::Graph& query, const FMCS::Graph& target
                              , const ROMol &querySrc
                              , const ROMol &mol    // target
                              , MCSAtomCompareFunction atomCompare, const MCSAtomCompareParameters& p, void* ud
                            )
                : QueryTopology(query), TargetTopology(target)
                , d_query(querySrc)
                , d_mol(mol)
                , AtomCompare(atomCompare), Parameters(p), UserData(ud)
            {};
            bool operator()(unsigned int i,unsigned int j)const {
                return AtomCompare(Parameters, d_query, QueryTopology[i], d_mol, TargetTopology[j], UserData);
            }
        };

        class BondLabelFunctor {
            const FMCS::Graph& QueryTopology;
            const FMCS::Graph& TargetTopology;
            const ROMol &d_query;
            const ROMol &d_mol;
            MCSBondCompareFunction BondCompare;
            const MCSBondCompareParameters& Parameters;
            void* UserData;
        public:
            BondLabelFunctor(const FMCS::Graph& query, const FMCS::Graph& target
                             , const ROMol &querySrc
                             , const ROMol &mol
                             , MCSBondCompareFunction bondCompare
                             , const MCSBondCompareParameters& p, void* ud)
                : QueryTopology(query), TargetTopology(target)
                , d_query(querySrc)
                , d_mol(mol)
                , BondCompare(bondCompare), Parameters(p), UserData(ud)
            {};

            bool operator()(FMCS::Graph::edge_descriptor i, FMCS::Graph::edge_descriptor j)const {
                unsigned  ii =  QueryTopology[i]; //take actual Idx value for full source query molecule from index list
                unsigned  jj = TargetTopology[j]; // the same Idx
                return BondCompare(Parameters, d_query, ii, d_mol, jj, UserData);
            }
        };

        bool SubstructMatchCustom(const FMCS::Graph& target, const ROMol &mol
                                  , const FMCS::Graph&  query, const ROMol &querySrc // seed and full source query molecule
                                  , MCSAtomCompareFunction atomCompare, MCSBondCompareFunction bondCompare
                                  , const MCSAtomCompareParameters& acp
                                  , const MCSBondCompareParameters& bcp
                                  , void* ud
                                  , match_V_t* match
                                 ) {
            MolMatchFinalCheckFunctor matchChecker(querySrc, mol);
            AtomLabelFunctor atomLabeler(query, target, querySrc, mol, atomCompare, acp, ud);
            BondLabelFunctor bondLabeler(query, target, querySrc, mol, bondCompare, bcp, ud);

            match_V_t dummy_match;
            if(!match)
                match = &dummy_match;
            return boost::vf2(query, target, atomLabeler, bondLabeler, matchChecker, *match);
        }
    }
}
