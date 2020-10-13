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
#include <boost/graph/adjacency_list.hpp>
#include "vf2-orig.hpp"

namespace RDKit {
namespace FMCS {

class MolMatchFinalCheckFunctor {
  const FMCS::Graph& QueryTopology;
  const FMCS::Graph& TargetTopology;
  const ROMol& d_query;
  const ROMol& d_mol;
  const MCSParameters* Parameters;

 public:
  MolMatchFinalCheckFunctor(const FMCS::Graph& query, const FMCS::Graph& target,
                            const ROMol& querySrc, const ROMol& mol  // target
                            ,
                            const MCSParameters* parameters)
      : QueryTopology(query),
        TargetTopology(target),
        d_query(querySrc),
        d_mol(mol),
        Parameters(parameters){};

  bool operator()(const boost::detail::node_id c1[],
                  const boost::detail::node_id c2[]) const {
    if ((unsigned)c1[0] >= boost::num_vertices(QueryTopology)) {
      return false;  // invalid index - match failed, see v2f implementation
    }
    MCSFinalMatchCheckFunction compare =
        Parameters ? Parameters->FinalMatchChecker : nullptr;
    return compare ? compare(c1, c2, d_query, QueryTopology, d_mol,
                             TargetTopology, Parameters)
                   : true;
  }
};
//=================================================================================================
// PRECOMPUTED_TABLES_MATCH much faster in overall even with very simple compare
// functions

class AtomTableCompareFunctor {
  const FMCS::Graph& QueryTopology;
  const FMCS::Graph& TargetTopology;
  const FMCS::MatchTable& MatchTable;

 public:
  AtomTableCompareFunctor(const FMCS::Graph& query, const FMCS::Graph& target,
                          const FMCS::MatchTable& targetMatch)
      : QueryTopology(query), TargetTopology(target), MatchTable(targetMatch){};
  bool operator()(unsigned int i, unsigned int j) const {
    return MatchTable.at(QueryTopology[i], TargetTopology[j]);
  }
};

class BondTableCompareFunctor {
  const FMCS::Graph& QueryTopology;
  const FMCS::Graph& TargetTopology;
  const FMCS::MatchTable& MatchTable;

 public:
  BondTableCompareFunctor(const FMCS::Graph& query, const FMCS::Graph& target,
                          const FMCS::MatchTable& targetMatch)
      : QueryTopology(query), TargetTopology(target), MatchTable(targetMatch){};
  bool operator()(FMCS::Graph::edge_descriptor i,
                  FMCS::Graph::edge_descriptor j) const {
    return MatchTable.at(QueryTopology[i], TargetTopology[j]);
  }
};

bool SubstructMatchCustomTable(const FMCS::Graph& target, const ROMol& mol,
                               const FMCS::Graph& query, const ROMol& querySrc,
                               const MatchTable& atomMatchTable,
                               const MatchTable& bondMatchTable,
                               const MCSParameters* p, match_V_t* match) {
  if (query.m_vertices.size() > target.m_vertices.size()  // query > target
      || query.m_edges.size() > target.m_edges.size()) {
    return false;
  }

  MolMatchFinalCheckFunctor mc(query, target, querySrc, mol, p);

  AtomTableCompareFunctor ac(query, target, atomMatchTable);
  BondTableCompareFunctor bc(query, target, bondMatchTable);

  match_V_t dummy_match;
  if (!match) {
    match = &dummy_match;
  }
  return boost::vf2(query, target, ac, bc, mc, *match);
}

//=========================================================================
// slow implementation with absolutely the same functionality
//=========================================================================

class AtomLabelFunctor {
  const FMCS::Graph& QueryTopology;
  const FMCS::Graph& TargetTopology;
  const ROMol& d_query;
  const ROMol& d_mol;
  MCSAtomCompareFunction AtomCompare;
  const MCSAtomCompareParameters& Parameters;
  void* UserData;

 public:
  AtomLabelFunctor(const FMCS::Graph& query, const FMCS::Graph& target,
                   const ROMol& querySrc, const ROMol& mol  // target
                   ,
                   MCSAtomCompareFunction atomCompare,
                   const MCSAtomCompareParameters& p, void* ud)
      : QueryTopology(query),
        TargetTopology(target),
        d_query(querySrc),
        d_mol(mol),
        AtomCompare(atomCompare),
        Parameters(p),
        UserData(ud){};
  bool operator()(unsigned int i, unsigned int j) const {
    return AtomCompare(Parameters, d_query, QueryTopology[i], d_mol,
                       TargetTopology[j], UserData);
  }
};

class BondLabelFunctor {
  const FMCS::Graph& QueryTopology;
  const FMCS::Graph& TargetTopology;
  const ROMol& d_query;
  const ROMol& d_mol;
  MCSBondCompareFunction BondCompare;
  const MCSBondCompareParameters& Parameters;
  void* UserData;

 public:
  BondLabelFunctor(const FMCS::Graph& query, const FMCS::Graph& target,
                   const ROMol& querySrc, const ROMol& mol,
                   MCSBondCompareFunction bondCompare,
                   const MCSBondCompareParameters& p, void* ud)
      : QueryTopology(query),
        TargetTopology(target),
        d_query(querySrc),
        d_mol(mol),
        BondCompare(bondCompare),
        Parameters(p),
        UserData(ud){};

  bool operator()(FMCS::Graph::edge_descriptor i,
                  FMCS::Graph::edge_descriptor j) const {
    unsigned ii = QueryTopology[i];   // take actual Idx value for full source
                                      // query molecule from index list
    unsigned jj = TargetTopology[j];  // the same Idx
    return BondCompare(Parameters, d_query, ii, d_mol, jj, UserData);
  }
};

bool SubstructMatchCustom(
    const FMCS::Graph& target, const ROMol& mol, const FMCS::Graph& query,
    const ROMol& querySrc  // seed and full source query molecule
    ,
    MCSAtomCompareFunction atomCompare, MCSBondCompareFunction bondCompare,
    MCSFinalMatchCheckFunction finalCompare,
    const MCSAtomCompareParameters& acp, const MCSBondCompareParameters& bcp,
    void* ud, match_V_t* match) {
  RDUNUSED_PARAM(finalCompare);
  MolMatchFinalCheckFunctor matchChecker(query, target, querySrc, mol, nullptr);
  AtomLabelFunctor atomLabeler(query, target, querySrc, mol, atomCompare, acp,
                               ud);
  BondLabelFunctor bondLabeler(query, target, querySrc, mol, bondCompare, bcp,
                               ud);

  match_V_t dummy_match;
  if (!match) {
    match = &dummy_match;
  }
  return boost::vf2(query, target, atomLabeler, bondLabeler, matchChecker,
                    *match);
}
}
}
