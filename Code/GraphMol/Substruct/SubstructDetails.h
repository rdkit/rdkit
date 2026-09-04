//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_SUBSTRUCTDETAILS_H
#define RD_SUBSTRUCTDETAILS_H

// std bits
#include <vector>
#include <map>

#include <GraphMol/QueryAtom.h>

namespace RDKit {
class ROMol;
class RecursiveStructureQuery;
class SubstructMatchParameters;

namespace detail {
using SUBQUERY_MAP = std::map<unsigned int, QueryAtom::QUERYATOM_QUERY *>;
RDKIT_SUBSTRUCTMATCH_EXPORT void MatchSubqueries(
    const ROMol &mol, QueryAtom::QUERYATOM_QUERY *q,
    const SubstructMatchParameters &params, SUBQUERY_MAP &subqueryMap,
    std::vector<RecursiveStructureQuery *> &locked);
}  // namespace detail
}  // namespace RDKit

#endif
