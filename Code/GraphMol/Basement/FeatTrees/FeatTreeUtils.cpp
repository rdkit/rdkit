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

#include "FeatTree.h"

#include <boost/graph/biconnected_components.hpp>
#include <boost/property_map.hpp>
#include <map>
#include <set>

namespace RDKit {
  namespace FeatTrees {
    typedef std::vector<boost::graph_traits<FeatTreeGraph>::vertex_descriptor> NodeVect;
    typedef std::set< boost::graph_traits<FeatTreeGraph>::vertex_descriptor > NodeSet;
    typedef std::set< boost::graph_traits<FeatTreeGraph>::edge_descriptor > EdgeSet;

    // local utility namespace
    namespace {
      /*!
	 This function replaces the elements of a system of
	 fused rings that form a cycle in the feature tree with
	 a single node. This is, in essence, a merge operation.

	 \param featGraph: the graph to be modified
	 \param featGraphCopy: a copy of the modified graph
	 \param components: a property map from edge->component id,
	        this property map refers to featGraphCopy
	 \param componentIdx: component id to be replaced

      */
      void mergeRingCycle(FeatTreeGraph &featGraph,FeatTreeGraph &featGraphCopy,
			  FeatTreeEdgePMap &components,
			  unsigned int componentIdx){
	UINT_SET atomIndices;
	FeatTreeNodePMap nodeMap = boost::get(FeatTreeNode_t(),featGraph);
	NodeVect nodeVect;
	
	boost::graph_traits<FeatTreeGraph>::edge_iterator edge,endEdges;
	for(boost::tie(edge,endEdges)=boost::edges(featGraphCopy);
	    edge!=endEdges;++edge){
	  if(components[*edge]==componentIdx){
	    boost::graph_traits<FeatTreeGraph>::vertex_descriptor node;

	    node = boost::source(*edge,featGraphCopy);
	    atomIndices.insert(nodeMap[node].begin(),nodeMap[node].end());
	    if(std::find(nodeVect.begin(),nodeVect.end(),node)==nodeVect.end()){
	      nodeVect.push_back(node);
	    }

	    node = boost::target(*edge,featGraphCopy);
	    atomIndices.insert(nodeMap[node].begin(),nodeMap[node].end());
	    if(std::find(nodeVect.begin(),nodeVect.end(),node)==nodeVect.end()){
	      nodeVect.push_back(node);
	    }
	  }
	}
	// ------ ------ ------ ------ ------ ------ ------
	// We now have a set of all the nodes involved in the cycle
	// as well as the full set of atom indices that are going
	// to be associated with the new node. 
	// We're going to:
	//   1) construct the new node and add it to the graph
	//   2) loop over each of the cycle nodes and:
	//       2.1) attach each of their outer edges to the new node
	//       2.2) remove the cycle node from the graph     
	boost::graph_traits<FeatTreeGraph>::vertex_descriptor newNode;
	newNode=boost::add_vertex(FeatTreeNode(atomIndices),featGraph);

	// we need to have the vertices sorted in inverse order (from
	// highest to lowest) before we start removing so that we don't
	// invalidate any of the "descriptors" that we have:
	std::sort(nodeVect.begin(),nodeVect.end());
	std::reverse(nodeVect.begin(),nodeVect.end());

	for(NodeVect::const_iterator nodeIt=nodeVect.begin();
	    nodeIt!=nodeVect.end();++nodeIt){
	  FeatTreeEdgePMap edgeMap = boost::get(FeatTreeEdge_t(),featGraph);
	  boost::graph_traits<FeatTreeGraph>::out_edge_iterator edge,endEdges;
	  for(boost::tie(edge,endEdges)=boost::out_edges(*nodeIt,featGraph);
	      edge!=endEdges;++edge){
	    // figure out which of the edge's nodes is the neighbor:
	    if(boost::source(*edge,featGraph)==*nodeIt){
	      // make sure the neighbor isn't in the nodeVect:
	      if(!std::binary_search(nodeVect.rbegin(),nodeVect.rend(),
				     boost::target(*edge,featGraph))){
		boost::add_edge(newNode,boost::target(*edge,featGraph),
				FeatTreeEdge(0),featGraph);
	      }
	    } else if(boost::target(*edge,featGraph)==*nodeIt) {
	      // make sure the neighbor isn't in the nodeSet:
	      if(!std::binary_search(nodeVect.rbegin(),nodeVect.rend(),
				     boost::source(*edge,featGraph))){
		boost::add_edge(boost::source(*edge,featGraph),newNode,
				FeatTreeEdge(0),featGraph);
	      }
	    } else {
	      CHECK_INVARIANT(false,
			      "inconsistent state encounted when replacing a cycle");
	    }
	  }
	  
	  // now remove the node from the graph:
	  boost::clear_vertex(*nodeIt,featGraph);
	  boost::remove_vertex(*nodeIt,featGraph);
	}
	
      }
      
    } // end of local utility namespace


    // -----------------------------------------------------------------------
    //
    // -----------------------------------------------------------------------
    void addRingsAndConnectors(const ROMol &mol,FeatTreeGraph &resGraph){
      RingInfo *ringInfo=mol.getRingInfo();
      unsigned int ringIdxI=0;
      for(VECT_INT_VECT::const_iterator ringItI=ringInfo->atomRings().begin();
	  ringItI != ringInfo->atomRings().end();++ringItI,++ringIdxI){
	UINT_SET s;
	s.insert(ringItI->begin(),ringItI->end());
	boost::add_vertex(FeatTreeNode(s),resGraph);

	// ------ ------ ------ ------
	// Add any relevant ring-ring connectors:
	// ------ ------ ------ ------
	// sort a copy of this ring's atoms so that it's easier to 
	// search for overlaps:
	INT_VECT ringI=*ringItI;
	std::sort(ringI.begin(),ringI.end());
	unsigned int ringIdxJ=0;
	for(VECT_INT_VECT::const_iterator ringItJ=ringInfo->atomRings().begin();
	    ringItJ != ringItI;
	    ++ringItJ,++ringIdxJ){
	  for(INT_VECT::const_iterator ringJElem=ringItJ->begin();
	      ringJElem != ringItJ->end();
	      ++ringJElem){
	    if(std::binary_search(ringI.begin(),ringI.end(),*ringJElem)){
	      // these two rings share a common atom, so set up an
	      // edge between them in the feature tree:
	      boost::add_edge(ringIdxI,ringIdxJ,
			      FeatTreeEdge(2),
			      resGraph);

	      // no point continuing the search for ringJ:
	      break;
	    }
	  }
	}
      }
    }

    // -----------------------------------------------------------------------
    //
    // -----------------------------------------------------------------------
    void addRingRingBonds(const ROMol &mol,FeatTreeGraph &resGraph){
      FeatTreeNodePMap nodeMap=boost::get(FeatTreeNode_t(),resGraph);

      boost::graph_traits<FeatTreeGraph>::vertex_iterator nodeI,endNodes;
      for(boost::tie(nodeI,endNodes)=boost::vertices(resGraph);
	  nodeI!=endNodes;++nodeI){
	boost::graph_traits<FeatTreeGraph>::vertex_iterator nodeJ;
	nodeJ=nodeI;
	while(++nodeJ!=endNodes){
	  boost::graph_traits<FeatTreeGraph>::edge_descriptor tmp;
	  bool found=false;
	  boost::tie(tmp,found) = boost::edge(*nodeI,*nodeJ,resGraph);
	  for(UINT_SET::const_iterator setI=nodeMap[*nodeI].begin();
	      setI!=nodeMap[*nodeI].end() && !found;
	      setI++){
	    for(UINT_SET::const_iterator setJ=nodeMap[*nodeJ].begin();
		setJ!=nodeMap[*nodeJ].end() && !found;
		setJ++){
	      if(mol.getBondBetweenAtoms(*setI,*setJ)){
		// set up the bond:
		boost::add_edge(*nodeI,*nodeJ,
				FeatTreeEdge(1),
				resGraph);
		found=true;
		break;
	      }
	    }
	  }
	}
      }
    }
    

    // -----------------------------------------------------------------------
    //
    // -----------------------------------------------------------------------
    std::vector<unsigned int> addNonringAtoms(const ROMol &mol,FeatTreeGraph &resGraph){
      RingInfo *ringInfo=mol.getRingInfo();
      unsigned int nAts = mol.getNumAtoms();
      std::vector<unsigned int> atomIndices(nAts,nAts+1);
      for(ROMol::ConstAtomIterator atomIt=mol.beginAtoms();
	  atomIt!=mol.endAtoms();++atomIt){
	const Atom *atom=*atomIt;
	if( (atom->getDegree()>1 || atom->getDegree()+atom->getTotalNumHs()>1) &&
	    !ringInfo->numAtomRings(atom->getIdx())){
	  UINT_SET s;
	  s.insert(atom->getIdx());
	  unsigned int idx = boost::add_vertex(FeatTreeNode(s),resGraph);
	  atomIndices[atom->getIdx()] = idx;
	}
      }
      return atomIndices;
    }

    // -----------------------------------------------------------------------
    //
    // -----------------------------------------------------------------------
    void addBondsFromNonringAtoms(const ROMol &mol,FeatTreeGraph &resGraph,
				  std::vector<unsigned int> &atomIndices){
      FeatTreeNodePMap nodeMap=boost::get(FeatTreeNode_t(),resGraph);
      unsigned int nAts=mol.getNumAtoms();
      for(ROMol::ConstAtomIterator atomIt=mol.beginAtoms();
	  atomIt!=mol.endAtoms();++atomIt){
	const Atom *atom=*atomIt;
	if(atomIndices[atom->getIdx()]<nAts){
	  // this atom has already been added as a free-standing atom:
	  ROMol::ADJ_ITER nbrIdx,endNbrs;
	  for(boost::tie(nbrIdx,endNbrs)=mol.getAtomNeighbors(atom);
	      nbrIdx!=endNbrs;++nbrIdx){
	    if(atomIndices[*nbrIdx]<nAts ) {
	      if(*nbrIdx>atom->getIdx()){
		// the neighbor has already been added, and it has a higher
		// index than ours (we make this check to avoid duplicated
		// bonds), so add a bond to it:
		boost::add_edge(atomIndices[atom->getIdx()],
				atomIndices[*nbrIdx],
				FeatTreeEdge(0),
				resGraph);
	      }
	    } else if(mol.getAtomWithIdx(*nbrIdx)->getDegree()>1) {
	      // the neighbor hasn't been added to the graph
	      // on its own and the degree is above 1, so it's
	      // got to be a ring atom, find the rings it's present in by
	      // looping over nodes already in the graph.  Note that we can't
	      // use the ringInfo structure anymore because we may have combined
	      // some rings in the call to replaceCycles() above:
	      boost::graph_traits<FeatTreeGraph>::vertex_iterator node,endNodes;
	      for(boost::tie(node,endNodes)=boost::vertices(resGraph);
		  node!=endNodes;++node){
		if(atomIndices[*node]!=*nbrIdx && nodeMap[*node].count(*nbrIdx)){
		  boost::add_edge(atomIndices[atom->getIdx()],
				  *node,
				  FeatTreeEdge(1),
				  resGraph);
		}
	      }
	    }
	  }
	}
      }      

    }

  
    // -----------------------------------------------------------------------
    //
    // -----------------------------------------------------------------------
    void replaceCycles(FeatTreeGraph &featGraph) {
      FeatTreeGraph featGraphCopy=featGraph;
      
      FeatTreeEdgePMap edgeMap = boost::get(FeatTreeEdge_t(),featGraph);
      FeatTreeEdgePMap componentMap=boost::get(FeatTreeEdge_t(),featGraphCopy);

      // ------ ------ ------ ------ ------ ------ ------ 
      // Start by finding the biconnected components:
      unsigned int numComponents=boost::biconnected_components(featGraphCopy,
							       componentMap);
      if(numComponents==boost::num_vertices(featGraph)){
	// no cycles, our work here is done, so go ahead and return
	return;
      }

      // ------ ------ ------ ------ ------ ------ ------ 
      // loop over the elements of the biconnected-components map and count:
      //    - how many times each index occurs
      std::vector<unsigned int> memberCount(numComponents,0);
      boost::graph_traits<FeatTreeGraph>::edge_iterator edge,endEdges;
      boost::graph_traits<FeatTreeGraph>::edge_iterator cpyEdge,endCpyEdges;
      for((boost::tie(edge,endEdges)=boost::edges(featGraph)),
	    (boost::tie(cpyEdge,endCpyEdges)=boost::edges(featGraphCopy));
	  edge!=endEdges;++edge,++cpyEdge){
	memberCount[componentMap[*cpyEdge]] += 1;
	CHECK_INVARIANT(edgeMap[*edge]==2,
			"replaceCycles() should be called with only ring-ring edges in the graph");
      }

      // ------ ------ ------ ------ ------ ------ ------ 
      // Replace any cycles that exist in the graph:
      for(unsigned int i=0; i<numComponents; ++i){
	if(memberCount[i]>1){
	  mergeRingCycle(featGraph,featGraphCopy,componentMap,i);


	}
      }
    }

    // -----------------------------------------------------------------------
    //
    // -----------------------------------------------------------------------
    void addZeroNodes(FeatTreeGraph &featGraph) {
      FeatTreeEdgePMap edgeMap = boost::get(FeatTreeEdge_t(),featGraph);
      FeatTreeGraph featGraphCopy=featGraph;
      FeatTreeEdgePMap componentMap=boost::get(FeatTreeEdge_t(),featGraphCopy);

      // ------ ------ ------ ------ ------ ------ ------ 
      // Start by finding the biconnected components:
      unsigned int numComponents=boost::biconnected_components(featGraphCopy,
							       componentMap);

      if(numComponents==boost::num_edges(featGraph)){
	// no cycles, our work here is done, so go ahead and return
	return;
      }

      // ------ ------ ------ ------ ------ ------ ------ 
      // loop over the elements of the biconnected-components map and:
      //    1) count how many times each index occurs
      //    2) determine if each component has a non-ring->ring bond
      std::vector<unsigned int> memberCount(numComponents,0);
      std::vector<bool> weight1EdgeObserved(numComponents,false);
      boost::graph_traits<FeatTreeGraph>::edge_iterator edge,endEdges;
      boost::graph_traits<FeatTreeGraph>::edge_iterator cpyEdge,endCpyEdges;
      for((boost::tie(edge,endEdges)=boost::edges(featGraph)),
	    (boost::tie(cpyEdge,endCpyEdges)=boost::edges(featGraphCopy));
	  edge!=endEdges;++edge,++cpyEdge){
	memberCount[componentMap[*cpyEdge]] += 1;
	if(edgeMap[*edge]==1){
	  weight1EdgeObserved[componentMap[*cpyEdge]]=true;
	}
      }

      for(unsigned int i=0;i<numComponents;i++){
	if(memberCount[i]>2){
	  CHECK_INVARIANT(weight1EdgeObserved[i],"internal inconsistency");
	  // add the zero node:
	  unsigned int newIdx;
	  UINT_SET emptySet;
	  newIdx=boost::add_vertex(FeatTreeNode(emptySet),featGraph);

	  // figure out who it needs to be connected to:
	  NodeSet nodesToConnect;
	  for((boost::tie(edge,endEdges)=boost::edges(featGraph)),
		(boost::tie(cpyEdge,endCpyEdges)=boost::edges(featGraphCopy));
	      cpyEdge!=endCpyEdges;
	      ++edge,++cpyEdge){
	    if(componentMap[*cpyEdge]==i){
	      nodesToConnect.insert(boost::source(*edge,featGraph));
	      nodesToConnect.insert(boost::target(*edge,featGraph));
	      // we can't actually remove edges at this stage, because we need
	      // the iterators for featGraph edges and featGraphCopy edges to
	      // stay in sync. So just mark this edge as one that should be
	      // removed and we'll take care of it later.
	      edgeMap[*edge]=23;
	    }
	  }
	  for(NodeSet::iterator nodeI=nodesToConnect.begin();
	      nodeI!=nodesToConnect.end();
	      nodeI++){
	    boost::add_edge(newIdx,*nodeI,4,featGraph);
	  }
	}
      }
      edgeMap = boost::get(FeatTreeEdge_t(),featGraph);
      boost::tie(edge,endEdges)=boost::edges(featGraph);
      while(edge!=endEdges){
	if(edgeMap[*edge]==23){
	  boost::graph_traits<FeatTreeGraph>::edge_iterator nextEdge=edge;
	  ++nextEdge;
	  boost::remove_edge(*edge,featGraph);
	  edge = nextEdge;
	} else {
	  ++edge;
	}
      }
      
    }      


    
    
  } //end of namespace FeatTrees
} // end of namespace RDKit

