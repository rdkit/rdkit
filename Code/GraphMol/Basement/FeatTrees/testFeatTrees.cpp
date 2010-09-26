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
#include <RDGeneral/RDLog.h>
#include <boost/log/functions.hpp>
#include <GraphMol/RDKitBase.h>
#include <string>
#include <iostream>
#include <GraphMol/FeatTrees/FeatTree.h>
#include <GraphMol/FeatTrees/FeatTreeUtils.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <boost/graph/graph_utility.hpp>

using namespace RDKit;


using namespace FeatTrees;
void showFeatTree(FeatTreeGraph &featGraph){
  typedef boost::property_map<FeatTreeGraph,FeatTreeNode_t>::type FeatTreeNodePMap;
  FeatTreeNodePMap nodes=boost::get(FeatTreeNode_t(),featGraph);
  boost::graph_traits<FeatTreeGraph>::vertex_iterator vtx,endVtx;
  unsigned int i;
  for(boost::tie(vtx,endVtx)=boost::vertices(featGraph),i=0;
      vtx!=endVtx;++vtx,++i){
    UINT_SET s = nodes[*vtx];
    std::cout << "\t" << i << ": [";
    std::copy(s.begin(),s.end(),std::ostream_iterator<int>(std::cout,","));
    std::cout << "]" << std::endl;
  }
}

void test1() {
  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test1: building FeatTreeGraphs \n\n";

  std::string smi;
  ROMol *m;
  FeatTreeGraphSPtr fGraph;
  std::vector<unsigned int> atomIndices;

  smi="CCC";
  std::cout << smi << std::endl;
  m=SmilesToMol(smi);
  TEST_ASSERT(m);
  fGraph.reset(new FeatTreeGraph());
  addRingsAndConnectors(*m,*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==0);
  TEST_ASSERT(boost::num_edges(*fGraph)==0);
  replaceCycles(*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==0);
  TEST_ASSERT(boost::num_edges(*fGraph)==0);
  addRingRingBonds(*m,*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==0);
  TEST_ASSERT(boost::num_edges(*fGraph)==0);
  atomIndices=addNonringAtoms(*m,*fGraph);
  TEST_ASSERT(atomIndices.size()==3);
  TEST_ASSERT(boost::num_vertices(*fGraph)==3);
  TEST_ASSERT(boost::num_edges(*fGraph)==0);
  addBondsFromNonringAtoms(*m,*fGraph,atomIndices);
  TEST_ASSERT(boost::num_vertices(*fGraph)==3);
  TEST_ASSERT(boost::num_edges(*fGraph)==2);
  addZeroNodes(*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==3);
  TEST_ASSERT(boost::num_edges(*fGraph)==2);
  fGraph = molToBaseTree(*m);
  TEST_ASSERT(boost::num_vertices(*fGraph)==3);
  TEST_ASSERT(boost::num_edges(*fGraph)==2);

  smi="C1CC1OOC1CCC1";
  std::cout << smi << std::endl;
  m=SmilesToMol(smi);
  TEST_ASSERT(m);
  fGraph.reset(new FeatTreeGraph());
  addRingsAndConnectors(*m,*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==2);
  TEST_ASSERT(boost::num_edges(*fGraph)==0);
  replaceCycles(*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==2);
  TEST_ASSERT(boost::num_edges(*fGraph)==0);
  addRingRingBonds(*m,*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==2);
  TEST_ASSERT(boost::num_edges(*fGraph)==0);
  atomIndices=addNonringAtoms(*m,*fGraph);
  TEST_ASSERT(atomIndices.size()==9);
  TEST_ASSERT(boost::num_vertices(*fGraph)==4);
  TEST_ASSERT(boost::num_edges(*fGraph)==0);
  addBondsFromNonringAtoms(*m,*fGraph,atomIndices);
  TEST_ASSERT(boost::num_vertices(*fGraph)==4);
  TEST_ASSERT(boost::num_edges(*fGraph)==3);
  addZeroNodes(*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==4);
  TEST_ASSERT(boost::num_edges(*fGraph)==3);
  fGraph = molToBaseTree(*m);
  TEST_ASSERT(boost::num_vertices(*fGraph)==4);
  TEST_ASSERT(boost::num_edges(*fGraph)==3);
  
  
  delete m;
  smi="C1CC2C1CCC2C(=O)";
  std::cout << smi << std::endl;
  m=SmilesToMol(smi);
  TEST_ASSERT(m);
  fGraph.reset(new FeatTreeGraph());
  addRingsAndConnectors(*m,*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==2);
  TEST_ASSERT(boost::num_edges(*fGraph)==1);
  replaceCycles(*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==2);
  TEST_ASSERT(boost::num_edges(*fGraph)==1);
  addRingRingBonds(*m,*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==2);
  TEST_ASSERT(boost::num_edges(*fGraph)==1);
  atomIndices=addNonringAtoms(*m,*fGraph);
  TEST_ASSERT(atomIndices.size()==9);
  TEST_ASSERT(boost::num_vertices(*fGraph)==3);
  TEST_ASSERT(boost::num_edges(*fGraph)==1);
  addBondsFromNonringAtoms(*m,*fGraph,atomIndices);
  TEST_ASSERT(boost::num_vertices(*fGraph)==3);
  TEST_ASSERT(boost::num_edges(*fGraph)==2);
  addZeroNodes(*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==3);
  TEST_ASSERT(boost::num_edges(*fGraph)==2);
  fGraph = molToBaseTree(*m);
  TEST_ASSERT(boost::num_vertices(*fGraph)==3);
  TEST_ASSERT(boost::num_edges(*fGraph)==2);

  delete m;
  smi="CC12CCCCC1CCCC2";
  std::cout << smi << std::endl;
  m=SmilesToMol(smi);
  TEST_ASSERT(m);
  fGraph.reset(new FeatTreeGraph());
  addRingsAndConnectors(*m,*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==2);
  TEST_ASSERT(boost::num_edges(*fGraph)==1);
  replaceCycles(*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==2);
  TEST_ASSERT(boost::num_edges(*fGraph)==1);
  addRingRingBonds(*m,*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==2);
  TEST_ASSERT(boost::num_edges(*fGraph)==1);
  atomIndices=addNonringAtoms(*m,*fGraph);
  TEST_ASSERT(atomIndices.size()==11);
  TEST_ASSERT(boost::num_vertices(*fGraph)==3);
  TEST_ASSERT(boost::num_edges(*fGraph)==1);
  addBondsFromNonringAtoms(*m,*fGraph,atomIndices);
  TEST_ASSERT(boost::num_vertices(*fGraph)==3);
  TEST_ASSERT(boost::num_edges(*fGraph)==3);
  addZeroNodes(*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==4);
  TEST_ASSERT(boost::num_edges(*fGraph)==3);
  fGraph = molToBaseTree(*m);
  TEST_ASSERT(boost::num_vertices(*fGraph)==4);
  TEST_ASSERT(boost::num_edges(*fGraph)==3);

  delete m;
  smi="C1CC2CCCCC2CC1C1CC(CC2)CCC21";
  std::cout << smi << std::endl;
  m=SmilesToMol(smi);
  TEST_ASSERT(m);
  fGraph.reset(new FeatTreeGraph());
  addRingsAndConnectors(*m,*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==5);
  TEST_ASSERT(boost::num_edges(*fGraph)==4);
  replaceCycles(*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==3);
  TEST_ASSERT(boost::num_edges(*fGraph)==1);
  addRingRingBonds(*m,*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==3);
  TEST_ASSERT(boost::num_edges(*fGraph)==2);
  atomIndices=addNonringAtoms(*m,*fGraph);
  TEST_ASSERT(atomIndices.size()==18);
  TEST_ASSERT(boost::num_vertices(*fGraph)==3);
  TEST_ASSERT(boost::num_edges(*fGraph)==2);
  addBondsFromNonringAtoms(*m,*fGraph,atomIndices);
  TEST_ASSERT(boost::num_vertices(*fGraph)==3);
  TEST_ASSERT(boost::num_edges(*fGraph)==2);
  addZeroNodes(*fGraph);
  TEST_ASSERT(boost::num_vertices(*fGraph)==3);
  TEST_ASSERT(boost::num_edges(*fGraph)==2);
  fGraph = molToBaseTree(*m);
  TEST_ASSERT(boost::num_vertices(*fGraph)==3);
  TEST_ASSERT(boost::num_edges(*fGraph)==2);



  BOOST_LOG(rdInfoLog) << "\tDone\n";
}



int main() { 
  RDLog::InitLogs();
  //boost::logging::disable_logs("rdApp.info");
  
  //BOOST_LOG(rdInfoLog) << "********************************************************\n";
  //BOOST_LOG(rdInfoLog) << "Testing FeatTrees\n";

  test1();

  //BOOST_LOG(rdInfoLog) << "*******************************************************\n";
  return(0);
}
