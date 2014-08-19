// $Id$
//
//  Copyright (C) 2001-2012 Greg Landrum and Rational Discovery LLC
//  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Atom.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/PeriodicTable.h>

#include <vector>
#include <algorithm> 

#include <boost/graph/connected_components.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 104000
#include <boost/property_map/property_map.hpp>
#else
#include <boost/property_map.hpp>
#endif

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <GraphMol/ROMol.h>

const int ci_LOCAL_INF=static_cast<int>(1e8);

namespace RDKit{
  namespace MolOps {
    namespace {
      void nitrogenCleanup(RWMol &mol,Atom *atom){
        // conversions here:
        // - neutral 5 coordinate Ns with double bonds to Os to the
        //   zwitterionic form.  e.g.:
        //   CN(=O)=O -> C[N+](=O)[O-]
        //   and:
        //   C1=CC=CN(=O)=C1 -> C1=CC=C[N+]([O-])=C1
        // - neutral 5 coordinate Ns with triple bonds to Ns to the
        //   zwitterionic form.  e.g.:
        //   C-N=N#N -> C-N=[N+]=[N-] 

        PRECONDITION(atom,"bad atom");
        bool aromHolder;

        // we only want to do neutrals so that things like this don't get
        // munged: 
        //  O=[n+]1occcc1
        // this was sf.net issue 1811276
        if(atom->getFormalCharge()) return;

        // we need to play this little aromaticity game because the
        // explicit valence code modifies its results for aromatic
        // atoms.
        aromHolder = atom->getIsAromatic();
        atom->setIsAromatic(0);
        // NOTE that we are calling calcExplicitValence() here, we do
        // this because we cannot be sure that it has already been
        // called on the atom (cleanUp() gets called pretty early in
        // the sanitization process):
        if(atom->calcExplicitValence(false)==5 ) {
          unsigned int aid = atom->getIdx();
          RWMol::ADJ_ITER nid1,end1;
          boost::tie(nid1, end1) = mol.getAtomNeighbors(atom);
          while (nid1 != end1) {
            if ((mol.getAtomWithIdx(*nid1)->getAtomicNum() == 8) &&
                (mol.getAtomWithIdx(*nid1)->getFormalCharge() == 0) &&
                (mol.getBondBetweenAtoms(aid, *nid1)->getBondType() == Bond::DOUBLE)) {
              // here's the double bonded oxygen
              Bond *b = mol.getBondBetweenAtoms(aid, *nid1);
              b->setBondType(Bond::SINGLE);
              atom->setFormalCharge(1);
              mol.getAtomWithIdx(*nid1)->setFormalCharge(-1);
              break;
            } else if ((mol.getAtomWithIdx(*nid1)->getAtomicNum() == 7) &&
                       (mol.getAtomWithIdx(*nid1)->getFormalCharge() == 0) &&
                       (mol.getBondBetweenAtoms(aid, *nid1)->getBondType() == Bond::TRIPLE)) {
              // here's the triple bonded nitrogen
              Bond *b = mol.getBondBetweenAtoms(aid, *nid1);
              b->setBondType(Bond::DOUBLE);
              atom->setFormalCharge(1);
              mol.getAtomWithIdx(*nid1)->setFormalCharge(-1);
              break;
            }
            ++nid1;
          } // end of loop over the first neigh
        } // if this atom is 5 coordinate nitrogen
          // force a recalculation of the explicit valence here
        atom->setIsAromatic(aromHolder);
        atom->calcExplicitValence(false);
      }
    }

    void cleanUp(RWMol &mol) {
      ROMol::AtomIterator ai; 
      for (ai = mol.beginAtoms(); ai != mol.endAtoms(); ++ai) {
        switch( (*ai)->getAtomicNum() ){
        case 7:
          nitrogenCleanup(mol,*ai);
          break;
        case 17:
          // recognize perchlorate and convert it from:
          //    Cl(=O)(=O)(=O)[O-]
          // to:
          //    [Cl+3]([O-])([O-])([O-])[O-]
          if((*ai)->calcExplicitValence(false)==7 && (*ai)->getFormalCharge()==0){
            unsigned int aid = (*ai)->getIdx();
            bool neighborsAllO=true;
            RWMol::ADJ_ITER nid1,end1;
            boost::tie(nid1, end1) = mol.getAtomNeighbors(*ai);
            while (nid1 != end1) {
              if(mol.getAtomWithIdx(*nid1)->getAtomicNum() != 8){
                neighborsAllO = false;
                break;
              }
              ++nid1;
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
                ++nid1;
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
      for (ROMol::AtomIterator ai = mol.beginAtoms();
	   ai != mol.endAtoms(); ++ai) {
        int origImplicitV = (*ai)->getImplicitValence();
        (*ai)->calcExplicitValence();
        int origExplicitV = (*ai)->getNumExplicitHs();
        
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

    void assignRadicals(RWMol &mol){
      for (ROMol::AtomIterator ai = mol.beginAtoms();
           ai != mol.endAtoms(); ++ai) {
        // we only put automatically assign radicals to things that
        // don't have them already and don't have implicit Hs:
        if( !(*ai)->getNoImplicit() || (*ai)->getNumRadicalElectrons() || !(*ai)->getAtomicNum()){
          continue;
        }
        double accum = 0.0;
        RWMol::OEDGE_ITER beg,end;
        boost::tie(beg,end) = mol.getAtomBonds(*ai);
        while(beg!=end){
          accum += mol[*beg]->getValenceContrib(*ai);
          ++beg;
        }
        accum += (*ai)->getNumExplicitHs();
        int totalValence = static_cast<int>(accum+0.1);
        int chg = (*ai)->getFormalCharge();
        int nOuter = PeriodicTable::getTable()->getNouterElecs((*ai)->getAtomicNum());
        int baseCount=8;
        if((*ai)->getAtomicNum()==1){
          baseCount=2;
        }

        // applies to later (more electronegative) elements:
        int numRadicals = std::max(baseCount - nOuter - totalValence + chg,0);
        // applies to earlier elements:
        int numRadicals2 = nOuter - totalValence - chg;
        if(numRadicals2>=0){
          numRadicals = std::min(numRadicals,numRadicals2);
        }
        (*ai)->setNumRadicalElectrons(numRadicals);
      }
    }
                
    void sanitizeMol(RWMol &mol){
      unsigned int failedOp=0;
      sanitizeMol(mol,failedOp,SANITIZE_ALL);
    }
    void sanitizeMol(RWMol &mol,
                     unsigned int &operationThatFailed,
                     unsigned int sanitizeOps){
      // clear out any cached properties
      mol.clearComputedProps();

      operationThatFailed=SANITIZE_CLEANUP;
      if(sanitizeOps & operationThatFailed){
        // clean up things like nitro groups
        cleanUp(mol);
      }

      // update computed properties on atoms and bonds:
      operationThatFailed = SANITIZE_PROPERTIES;
      if(sanitizeOps & operationThatFailed){
        mol.updatePropertyCache(true);
      } else {
        mol.updatePropertyCache(false);
      }
               
      operationThatFailed = SANITIZE_SYMMRINGS;
      if(sanitizeOps & operationThatFailed){
        VECT_INT_VECT arings;
        MolOps::symmetrizeSSSR(mol, arings);
      }

      // kekulizations
      operationThatFailed = SANITIZE_KEKULIZE;
      if(sanitizeOps & operationThatFailed){
        Kekulize(mol);
      }

      // look for radicals:
      // We do this now because we need to know
      // that the N in [N]1C=CC=C1 has a radical
      // before we move into setAromaticity().
      // It's important that this happen post-Kekulization
      // because there's no way of telling what to do 
      // with the same molecule if it's in the form
      // [n]1cccc1
      operationThatFailed = SANITIZE_FINDRADICALS;
      if(sanitizeOps & operationThatFailed){
        assignRadicals(mol);
      }
      
      // then do aromaticity perception
      operationThatFailed = SANITIZE_SETAROMATICITY;
      if(sanitizeOps & operationThatFailed){
        setAromaticity(mol);
      }
    
      // set conjugation
      operationThatFailed = SANITIZE_SETCONJUGATION;
      if(sanitizeOps & operationThatFailed){
        setConjugation(mol);
      }
    
      // set hybridization
      operationThatFailed = SANITIZE_SETHYBRIDIZATION;
      if(sanitizeOps & operationThatFailed){
        setHybridization(mol);
      }

      // remove bogus chirality specs:
      operationThatFailed = SANITIZE_CLEANUPCHIRALITY;
      if(sanitizeOps & operationThatFailed){
        cleanupChirality(mol);
      }

      // adjust Hydrogen counts:
      operationThatFailed = SANITIZE_ADJUSTHS;
      if(sanitizeOps & operationThatFailed){
        adjustHs(mol);
      }

      operationThatFailed = 0;
    }

    std::vector<ROMOL_SPTR> getMolFrags(const ROMol &mol,bool sanitizeFrags,
                                        INT_VECT *frags){
      bool ownIt=false;
      INT_VECT *mapping;
      if(frags){
        mapping=frags;
      } else {
        mapping = new INT_VECT;
        ownIt=true;
      }
      unsigned int nFrags=getMolFrags(mol,*mapping);
      std::vector<ROMOL_SPTR> res;
      if(nFrags==1){
        ROMol *tmp=new ROMol(mol);
        ROMOL_SPTR sptr(tmp);
        res.push_back(sptr);
      } else {
        res.reserve(nFrags);
        for(unsigned int fragIdx=0;fragIdx<nFrags;++fragIdx){
          // copy the molecule, then remove the atoms that
          // aren't in this fragment
          RWMol *tmp=new RWMol(mol);
          // loop over the atoms backwards so that indices
          // stay valid:
          for(int i=mol.getNumAtoms()-1;i>=0;--i){
            if((*mapping)[i]!=static_cast<int>(fragIdx)) tmp->removeAtom(i);
          }
          if(sanitizeFrags) sanitizeMol(*tmp);
          ROMOL_SPTR sptr(static_cast<ROMol *>(tmp));
          res.push_back(sptr);
        }
      }
      if(ownIt){
        delete mapping;
      }
      return res;
    }

    unsigned int getMolFrags(const ROMol &mol, INT_VECT &mapping) {
      mapping.resize(mol.getNumAtoms());
      return boost::connected_components(mol.getTopology(),&mapping[0]);
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


    template <typename T>
    std::map<T,boost::shared_ptr<ROMol> > getMolFragsWithQuery(const ROMol &mol,
                                                               T (*query)(const ROMol &,const Atom *),
                                                               bool sanitizeFrags,
                                                               const std::vector<T> *whiteList,
                                                               bool negateList){
      PRECONDITION(query,"no query");

      std::vector<T> assignments(mol.getNumAtoms());
      std::vector<int> ids(mol.getNumAtoms(),-1);
      std::map<T,boost::shared_ptr<ROMol> > res;
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        T where=query(mol,mol.getAtomWithIdx(i));
        if(whiteList){
          bool found=std::find(whiteList->begin(),whiteList->end(),where)!=whiteList->end();
          if(!found && !negateList) continue;
          else if (found && negateList) continue;
        }
        assignments[i]=where;
        if(res.find(where)==res.end()){
          res[where]=boost::shared_ptr<ROMol>(new ROMol());
        }
        RWMol *frag=static_cast<RWMol *>(res[where].get());
        ids[i]=frag->addAtom(mol.getAtomWithIdx(i)->copy(),false,true);
        // loop over neighbors and add bonds in the fragment to all atoms
        // that are already in the same fragment
        ROMol::ADJ_ITER nbrIdx,endNbrs;
        boost::tie(nbrIdx,endNbrs) = mol.getAtomNeighbors(mol.getAtomWithIdx(i));
        while(nbrIdx!=endNbrs){
          if(*nbrIdx<i && assignments[*nbrIdx]==where){
            Bond *nBond=mol.getBondBetweenAtoms(i,*nbrIdx)->copy();
            nBond->setOwningMol(static_cast<ROMol *>(frag));
            nBond->setBeginAtomIdx(ids[nBond->getBeginAtomIdx()]);
            nBond->setEndAtomIdx(ids[nBond->getEndAtomIdx()]);
            frag->addBond(nBond,true);
          }
          ++nbrIdx;
        }        
      }
      // update conformers
      for(ROMol::ConstConformerIterator cit=mol.beginConformers();
          cit!=mol.endConformers();++cit){
        for(typename std::map<T,boost::shared_ptr<ROMol> >::iterator iter=res.begin();
            iter!=res.end();++iter){
          ROMol *newM=iter->second.get();
          Conformer *conf=new Conformer(newM->getNumAtoms());
          conf->setId((*cit)->getId());
          conf->set3D((*cit)->is3D());
          newM->addConformer(conf);
        }
        for(unsigned int i=0;i<mol.getNumAtoms();++i){
          if(ids[i]<0) continue;
          res[assignments[i]]->getConformer((*cit)->getId()).setAtomPos(ids[i],(*cit)->getAtomPos(i));
        }
      }
      if(sanitizeFrags){
        for(typename std::map<T,boost::shared_ptr<ROMol> >::iterator iter=res.begin();
            iter!=res.end();++iter){
          sanitizeMol(*static_cast<RWMol *>(iter->second.get()));
        }
      }
      return res;
    }
    template std::map<std::string,boost::shared_ptr<ROMol> > getMolFragsWithQuery(const ROMol &mol,
                                                                                  std::string (*query)(const ROMol &,const Atom *),
                                                                                  bool sanitizeFrags,
                                                                                  const std::vector<std::string> *,
                                                                                  bool);
    template std::map<int,boost::shared_ptr<ROMol> > getMolFragsWithQuery(const ROMol &mol,
                                                                          int (*query)(const ROMol &,const Atom *),
                                                                          bool sanitizeFrags,
                                                                          const std::vector<int> *,
                                                                          bool);
    template std::map<unsigned int,boost::shared_ptr<ROMol> > getMolFragsWithQuery(const ROMol &mol,
                                                                                   unsigned int (*query)(const ROMol &,const Atom *),
                                                                                   bool sanitizeFrags,
                                                                                   const std::vector<unsigned int> *,
                                                                                   bool);

#if 0
    void findSpanningTree(const ROMol &mol,INT_VECT &mst){
      //
      //  The BGL provides Prim's and Kruskal's algorithms for finding
      //  the MST of a graph.  Prim's is O(n2) (n=# of atoms) while
      //  Kruskal's is O(e log e) (e=# of bonds).  For molecules, where
      //  e << n2, Kruskal's should be a win.
      //
      const MolGraph *mgraph = &mol.getTopology();
      MolGraph *molGraph = const_cast<MolGraph *> (mgraph);
    
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
        mst.push_back(mol[*edgeIt]->getIdx());
      }
    }
#endif
    int getFormalCharge(const ROMol &mol){
      int accum = 0;
      for(ROMol::ConstAtomIterator atomIt=mol.beginAtoms();
          atomIt!=mol.endAtoms();
          ++atomIt){
        accum += (*atomIt)->getFormalCharge();
      }
      return accum;
    };

    unsigned getNumAtomsWithDistinctProperty(const ROMol& mol, std::string prop)
    {
      unsigned numPropAtoms=0;
      for (ROMol::ConstAtomIterator ai = mol.beginAtoms();
      	  ai != mol.endAtoms(); ++ai) {
        if((*ai)->hasProp(prop)){
      	  ++numPropAtoms;
        }
      }
      return numPropAtoms;
    }
  }; // end of namespace MolOps
}; // end of namespace RDKit
