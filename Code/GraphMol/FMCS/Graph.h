//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// graph topology in terms of indices in source molecule
#include <RDGeneral/export.h>
#pragma once
#include <RDGeneral/BoostStartInclude.h>
#include <boost/graph/adjacency_list.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {
namespace FMCS {
typedef unsigned int AtomIdx_t;
typedef unsigned int BondIdx_t;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              AtomIdx_t, BondIdx_t>
    Graph_t;

class RDKIT_FMCS_EXPORT Graph : public Graph_t {
 public:
  typedef edge_iterator EDGE_ITER;
  typedef std::pair<EDGE_ITER, EDGE_ITER> BOND_ITER_PAIR;

  void addAtom(unsigned int atom) {
    Graph::vertex_descriptor which = boost::add_vertex(*this);
    (*this)[which] = atom;
  }
  void addBond(unsigned int bond, unsigned int beginAtom,
               unsigned int endAtom) {
    bool res;
    Graph_t::edge_descriptor which;
    boost::tie(which, res) = boost::add_edge(beginAtom, endAtom, *this);
    (*this)[which] = bond;
  }
};
}  // namespace FMCS
}  // namespace RDKit
