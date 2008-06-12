// $Id$
//
//  Copyright (C) 2001-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Atom.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>

#include <vector>
#include <algorithm> 

#include <boost/graph/connected_components.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/property_map.hpp>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>


const int ci_LOCAL_INF=static_cast<int>(1e8);

namespace RDKit{
  namespace MolOps {
    void cleanUp(RWMol &mol) {
      ROMol::AtomIterator ai; 
      int aid;
      bool aromHolder;
      for (ai = mol.beginAtoms(); ai != mol.endAtoms(); ai++) {
        switch( (*ai)->getAtomicNum() ){
        case 7:
          // convert neutral 5 coordinate Ns with double bonds to Os to the
          // zwitterionic form.  e.g.:
          //  CN(=O)=O -> C[N+](=O)[O-]
          // and:
          //  C1=CC=CN(=O)=C1 -> C1=CC=C[N+]([O-])=C1

          // we only want to do neutrals so that things like this don't get
          // munged: 
          //  O=[n+]1occcc1
          // this was sf.net issue 1811276
          if((*ai)->getFormalCharge()){
            continue;
          }

          // we need to play this little aromaticity game because the
          // explicit valence code modifies its results for aromatic
          // atoms.
          aromHolder = (*ai)->getIsAromatic();
          (*ai)->setIsAromatic(0);
          // NOTE that we are calling calcExplicitValence() here, we do
          // this because we cannot be sure that it has already been
          // called on the atom (cleanUp() gets called pretty early in
          // the sanitization process):
          if((*ai)->calcExplicitValence(false)==5 ) {
            aid = (*ai)->getIdx();
            RWMol::ADJ_ITER nid1,end1;
            boost::tie(nid1, end1) = mol.getAtomNeighbors(*ai);
            while (nid1 != end1) {
              if ((mol.getAtomWithIdx(*nid1)->getAtomicNum() == 8) &&
                  (mol.getBondBetweenAtoms(aid, *nid1)->getBondType() == Bond::DOUBLE)) {
                // here's the double bonded oxygen
                Bond *b = mol.getBondBetweenAtoms(aid, *nid1);
                b->setBondType(Bond::SINGLE);
                (*ai)->setFormalCharge(1);
                mol.getAtomWithIdx(*nid1)->setFormalCharge(-1);
                break;
              }
              nid1++;
            } // end of loop over the first neigh
          } // if this atom is 5 coordinate nitrogen
          // force a recalculation of the explicit valence here
          (*ai)->setIsAromatic(aromHolder);
          (*ai)->calcExplicitValence(false);
          break;
        case 17:
          // recognize perchlorate and convert it from:
          //    Cl(=O)(=O)(=O)[O-]
          // to:
          //    [Cl+3]([O-])([O-])([O-])[O-]
          if((*ai)->calcExplicitValence(false)==7 && (*ai)->getFormalCharge()==0){
            aid = (*ai)->getIdx();
            bool neighborsAllO=true;
            RWMol::ADJ_ITER nid1,end1;
            boost::tie(nid1, end1) = mol.getAtomNeighbors(*ai);
            while (nid1 != end1) {
              if(mol.getAtomWithIdx(*nid1)->getAtomicNum() != 8){
                neighborsAllO = false;
                break;
              }
              nid1++;
            }
            if(neighborsAllO){
              (*ai)->setFormalCharge(3);
              boost::tie(nid1, end1) = mol.getAtomNeighbors(*ai);
              while (nid1 != end1) {
                Bond *b = mol.getBondBetweenAtoms(aid, *nid1);
                if(b->getBondType()==Bond::DOUBLE){
                  b->setBondType(Bond::SINGLE);
                  Atom *otherAtom=mol.getAtomWithIdx(*nid1);
                  otherAtom->setFormalCharge(-1);
                  otherAtom->calcExplicitValence(false);
                }
                nid1++;
              }
              (*ai)->calcExplicitValence(false);
            }
          }
          break;
        }
      }
    }
                
    void adjustHs(RWMol &mol) {
      //
      //  Go through and adjust the number of implicit and explicit Hs
      //  on each atom in the molecule.
      //
      //  Atoms that do not *need* explicit Hs
      //
      //  Assumptions: this is called after the molecule has been
      //  sanitized, aromaticity has been perceived, and the implicit
      //  valence of everything has been calculated.
      //
      ROMol::AtomIterator ai; 
      for (ai = mol.beginAtoms(); ai != mol.endAtoms(); ai++) {
        int origImplicitV = (*ai)->getImplicitValence();
        (*ai)->calcExplicitValence();
        int origExplicitV = (*ai)->getNumExplicitHs();
        //(*ai)->setNumExplicitHs(0);
        
        int newImplicitV = (*ai)->calcImplicitValence();
        //
        //  Case 1: The disappearing Hydrogen
        //    Smiles:  O=C1NC=CC2=C1C=CC=C2
        //
        //    after perception is done, the N atom has two aromatic
        //    bonds to it and a single implict H.  When the Smiles is
        //    written, we get: n1ccc2ccccc2c1=O.  Here the nitrogen has
        //    no implicit Hs (because there are two aromatic bonds to
        //    it, giving it a valence of 3).  Also: this SMILES is bogus
        //    (un-kekulizable).  The correct SMILES would be:
        //    [nH]1ccc2ccccc2c1=O.  So we need to loop through the atoms
        //    and find those that have lost implicit H; we'll add those
        //    back as explict Hs.
        //
        //    <phew> that takes way longer to comment than it does to
        //    write:
        if(newImplicitV < origImplicitV){
          (*ai)->setNumExplicitHs(origExplicitV+(origImplicitV-newImplicitV));
          (*ai)->calcExplicitValence();
        }
      }
    }
                
                
    void sanitizeMol(RWMol &mol) {

      // clear out any cached properties
      mol.clearComputedProps();

      // clean up things like nitro groups
      cleanUp(mol);

      // update computed properties on atoms and bonds:
      mol.updatePropertyCache();

      // first do the kekulizations
      Kekulize(mol);
    
      // then do aromaticity perception
      setAromaticity(mol);
    
      // set conjugation
      setConjugation(mol);
    
      // set hybridization
      setHybridization(mol);

      // remove bogus chirality specs:
      cleanupChirality(mol);

      // adjust Hydrogen counts:
      adjustHs(mol);

    }

    unsigned int getMolFrags(const ROMol &mol, INT_VECT &mapping) {
      mapping.resize(mol.getNumAtoms());
      const MolGraph *G_p = mol.getTopology();
      return boost::connected_components(*G_p,&mapping[0]);
    };

    unsigned int getMolFrags(const ROMol &mol, VECT_INT_VECT &frags) {
      frags.clear();
      INT_VECT mapping;
      getMolFrags(mol, mapping);
  
      INT_INT_VECT_MAP comMap;
      for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
        int mi = mapping[i];
        if(comMap.find(mi)==comMap.end()){
          INT_VECT comp;
          comMap[mi] = comp;
        }
        comMap[mi].push_back(i);
      }

      for (INT_INT_VECT_MAP_CI mci = comMap.begin();
               mci != comMap.end();
               mci++) {
           frags.push_back((*mci).second);
      }
      return frags.size();
    }

    void findSpanningTree(const ROMol &mol,INT_VECT &mst){
      //
      //  The BGL provides Prim's and Kruskal's algorithms for finding
      //  the MST of a graph.  Prim's is O(n2) (n=# of atoms) while
      //  Kruskal's is O(e log e) (e=# of bonds).  For molecules, where
      //  e << n2, Kruskal's should be a win.
      //
      const MolGraph *mgraph = mol.getTopology();
      MolGraph *molGraph = const_cast<MolGraph *> (mgraph);
      ROMol::GRAPH_MOL_BOND_PMAP::const_type pMap = mol.getBondPMap();
    
      std::vector<MolGraph::edge_descriptor> treeEdges;
      treeEdges.reserve(boost::num_vertices(*molGraph));
    
      boost::property_map < MolGraph, edge_wght_t >::type w = boost::get(edge_wght_t(), *molGraph);
      boost::property_map < MolGraph, edge_bond_t>::type bps = boost::get(edge_bond_t(), *molGraph);
      boost::graph_traits < MolGraph >::edge_iterator e, e_end;
      Bond* bnd;
      for (boost::tie(e, e_end) = boost::edges(*molGraph); e != e_end; ++e) {
        bnd = bps[*e];
      
        if(!bnd->getIsAromatic()){
          w[*e] = (bnd->getBondTypeAsDouble());
        } else {
          w[*e] = 3.0/2.0;
        }
      }
    
      // FIX: this is a hack due to problems with MSVC++
#if 1
      typedef boost::graph_traits<MolGraph>::vertices_size_type size_type;
      typedef boost::graph_traits<MolGraph>::vertex_descriptor vertex_t;
      typedef boost::property_map<MolGraph,boost::vertex_index_t>::type index_map_t;
      boost::graph_traits<MolGraph>::vertices_size_type
        n = boost::num_vertices(*molGraph);
      std::vector<size_type> rank_map(n);
      std::vector<vertex_t> pred_map(n);

      boost::detail::kruskal_mst_impl
        (*molGraph, std::back_inserter(treeEdges),
         boost::make_iterator_property_map(rank_map.begin(),
                                           boost::get(boost::vertex_index, *molGraph),
                                           rank_map[0]),
         boost::make_iterator_property_map(pred_map.begin(),
                                           boost::get(boost::vertex_index, *molGraph),
                                           pred_map[0]),
         w);
  
#else  
      boost::kruskal_minimum_spanning_tree(*molGraph,std::back_inserter(treeEdges),
                                           w, *molGraph);
      //boost::weight_map(static_cast<boost::property_map<MolGraph,edge_wght_t>::const_type>(boost::get(edge_wght_t(),*molGraph))));
#endif
      mst.resize(0);
      for(std::vector<MolGraph::edge_descriptor>::iterator edgeIt=treeEdges.begin();
          edgeIt!=treeEdges.end();edgeIt++){
        mst.push_back(pMap[*edgeIt]->getIdx());
      }
    }

    int getFormalCharge(const ROMol &mol){
      int accum = 0;
      for(ROMol::ConstAtomIterator atomIt=mol.beginAtoms();
          atomIt!=mol.endAtoms();
          atomIt++){
        accum += (*atomIt)->getFormalCharge();
      }
      return accum;
    };
  }; // end of namespace MolOps
}; // end of namespace RDKit
