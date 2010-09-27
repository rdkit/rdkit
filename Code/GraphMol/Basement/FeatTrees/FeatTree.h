//
//  Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef _RD_FEATTREE_H_
#define _RD_FEATTREE_H_

#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map.hpp>
#include <boost/shared_ptr.hpp>
#include <set>

namespace RDKit {
  class ROMol;

  namespace FeatTrees {
    typedef std::set<unsigned int> UINT_SET;

    // Each node of the feature tree topology contains:
    //   - a record of the atom indices that are lumped into
    //     that node
    struct FeatTreeNode_t {
      enum { num=1027 };
      typedef boost::vertex_property_tag kind;
    };
    typedef boost::property<FeatTreeNode_t,UINT_SET> FeatTreeNode;

    // Each edge of the feature tree topology contains:
    //   - an indicator of the number of rings at the ends
    //     (0, 1, or 2)
    struct FeatTreeEdge_t {
      enum { num=1028 };
      typedef boost::edge_property_tag kind;
    };
    typedef boost::property<FeatTreeEdge_t,unsigned int> FeatTreeEdge;

    typedef boost::adjacency_list < boost::vecS, boost::vecS,
				    boost::undirectedS,
				    FeatTreeNode,
				    FeatTreeEdge > FeatTreeGraph;
    typedef boost::shared_ptr<FeatTreeGraph> FeatTreeGraphSPtr;
    
    typedef boost::property_map<FeatTreeGraph,FeatTreeEdge_t>::type FeatTreeEdgePMap;
    typedef boost::property_map<FeatTreeGraph,FeatTreeNode_t>::type FeatTreeNodePMap;


    /*!

    */
    FeatTreeGraphSPtr molToBaseTree(const ROMol &mol);

    void baseTreeToFeatTree(FeatTreeGraph &baseTree);

    
  }

  
}
#endif
