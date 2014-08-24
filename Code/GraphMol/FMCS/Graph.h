//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// graph topology in terms of indeces in source molecule
#pragma once
#include <boost/graph/adjacency_list.hpp>

namespace RDKit {
    namespace FMCS {
        typedef unsigned AtomIdx_t;
        typedef unsigned BondIdx_t;
        typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, AtomIdx_t, BondIdx_t> Graph_t;

        class Graph : public Graph_t {
        public:
            void addAtom(unsigned atom) {
                Graph::vertex_descriptor which = boost::add_vertex(*this);
                (*this)[which] = atom;
            }
            void addBond(unsigned bond, unsigned beginAtom, unsigned endAtom) {
                bool res;
                Graph_t::edge_descriptor which;
                boost::tie(which, res) = boost::add_edge(beginAtom, endAtom, *this);
                (*this)[which] = bond;
            }
        };

    }
}
