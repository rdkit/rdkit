// $Id$
//
//  Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include "FeatTreeUtils.h"
#include "FeatTree.h"
#include <vector>
#include <algorithm>
#include <boost/graph/graph_utility.hpp>

namespace RDKit {
  namespace FeatTrees {
    typedef boost::property_map<FeatTreeGraph,FeatTreeEdge_t>::type FeatTreeEdgePMap;

    FeatTreeGraphSPtr molToBaseTree(const ROMol &mol){
      FeatTreeGraphSPtr resGraph(new FeatTreeGraph());


      // EFF: This can almost certainly be done with a single pass through
      // the atoms and rings, for now we're just trying to get it right

      // ------ ------ ------ ------
      // add a node to the initial feature tree for each of the SSSR rings:
      // ------ ------ ------ ------
      addRingsAndConnectors(mol,*resGraph);
      

      // ------ ------ ------ ------
      // Reduce any fused ring systems that form cycles in
      // the FeatTree to single points:
      // ------ ------ ------ ------
      replaceCycles(*resGraph);


      // ------ ------ ------ ------
      // We have not yet added bonds between ring systems, like
      // the N-N bond in:
      //    C1CN1-N1CC1
      // So loop over the rings that we currently have and add
      // those bonds as appropriate:
      // ------ ------ ------ ------
      addRingRingBonds(mol,*resGraph);
      
      // ------ ------ ------ ------
      // Now go through the molecule and add every non-ring atom that
      // has more than one neighbor; we include Hs in the neighbor count
      // ------ ------ ------ ------
      std::vector<unsigned int> atomIndices;
      atomIndices=addNonringAtoms(mol,*resGraph);

      // ------ ------ ------ ------
      // At this point all atom and ring nodes have been added to
      // the graph.
      // Add the atom-atom and atom-ring edges:
      // ------ ------ ------ ------
      addBondsFromNonringAtoms(mol,*resGraph,atomIndices);

      // ------ ------ ------ ------
      // Go through the graph and add zero nodes
      //  (points not associated with atoms that are used to
      //  break cycles caused by non-ring atoms bonding to rings)
      // ------ ------ ------ ------
      addZeroNodes(*resGraph);
      
      return resGraph;
    }

    void baseTreeToFeatTree(FeatTreeGraph &baseTree) {

    }

    
  } // end of namespace FeatTrees

}
