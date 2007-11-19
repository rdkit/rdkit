// $Id$
//
//  Copyright (C) 2003-2006  Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "RDKitBase.h"
#include <list>
#include "QueryAtom.h"
#include "QueryOps.h"
#include <Geometry/Transform3D.h>
#include <Geometry/point.h>

namespace RDKit{

  // Local utility functionality:
  namespace {
    Atom *getAtomNeighborNot(ROMol *mol,const Atom *atom,const Atom *other){
      PRECONDITION(mol,"bad molecule");
      PRECONDITION(atom,"bad atom");
      PRECONDITION(atom->getDegree()>1,"bad degree");
      PRECONDITION(other,"bad atom");
      Atom *res=0;

      ROMol::ADJ_ITER nbrIdx,endNbrs;
      boost::tie(nbrIdx,endNbrs) = mol->getAtomNeighbors(atom);
      while(nbrIdx!=endNbrs){
        if(*nbrIdx != other->getIdx()){
          res = mol->getAtomWithIdx(*nbrIdx);
          break;
        }
        nbrIdx++;
      }

      POSTCONDITION(res,"no neighbor found");
      return res;
    }

    void setHydrogenCoords(ROMol *mol,unsigned int hydIdx,unsigned int heavyIdx){
      // we will loop over all the coordinates 
      PRECONDITION(mol,"bad molecule");
      PRECONDITION(heavyIdx!=hydIdx,"degenerate atoms");
      Atom *hydAtom = mol->getAtomWithIdx(hydIdx);
      PRECONDITION(mol->getAtomDegree(hydAtom)==1,"bad atom degree");
      const Bond *bond=mol->getBondBetweenAtoms(heavyIdx,hydIdx);
      PRECONDITION(bond,"no bond between atoms");

      const Atom *heavyAtom = mol->getAtomWithIdx(heavyIdx);
      double bondLength = PeriodicTable::getTable()->getRb0(1) +
        PeriodicTable::getTable()->getRb0(heavyAtom->getAtomicNum());
    
      RDGeom::Point3D dirVect(0,0,0);

      RDGeom::Point3D perpVect,rotnAxis,nbrPerp;
      RDGeom::Point3D nbr1Vect,nbr2Vect,nbr3Vect;
      RDGeom::Transform3D tform;
      RDGeom::Point3D heavyPos, hydPos; 

      const Atom *nbr1=0,*nbr2=0,*nbr3=0;
      const Bond *nbrBond;
      ROMol::ADJ_ITER nbrIdx,endNbrs;

      switch(heavyAtom->getDegree()){
      case 1:
        // --------------------------------------------------------------------------
        //   No other atoms present:
        // --------------------------------------------------------------------------
        dirVect.z = 1;
        // loop over the conformations and set the coordinates
        for (ROMol::ConformerIterator cfi = mol->beginConformers();
             cfi != mol->endConformers(); cfi++) {
          heavyPos = (*cfi)->getAtomPos(heavyIdx);
          hydPos = heavyPos + dirVect*bondLength;
          (*cfi)->setAtomPos(hydIdx, hydPos);
        }
        break;

      case 2:
        // --------------------------------------------------------------------------
        //  One other neighbor:
        // --------------------------------------------------------------------------
        nbr1=getAtomNeighborNot(mol,heavyAtom,hydAtom);
        for (ROMol::ConformerIterator cfi = mol->beginConformers();
             cfi != mol->endConformers(); cfi++) {
          heavyPos = (*cfi)->getAtomPos(heavyIdx);
          RDGeom::Point3D nbr1Pos = (*cfi)->getAtomPos(nbr1->getIdx());
          // get a normalized vector pointing away from the neighbor:
          nbr1Vect = heavyPos.directionVector(nbr1Pos);
          nbr1Vect *= -1;

          // ok, nbr1Vect points away from the other atom, figure out where
          // this H goes:
          switch(heavyAtom->getHybridization()){
          case Atom::SP3:
            // get a perpendicular to nbr1Vect:
            perpVect=nbr1Vect.getPerpendicular();
            // and move off it:
            tform.SetRotation((180-109.471)*M_PI/180.,perpVect);
            dirVect = tform*nbr1Vect;
            hydPos = heavyPos + dirVect*bondLength;
            (*cfi)->setAtomPos(hydIdx, hydPos);
            break;
          case Atom::SP2:
            // default position is to just take an arbitrary perpendicular:
            perpVect = nbr1Vect.getPerpendicular();
          
            if(nbr1->getDegree()>1){
              // can we use the neighboring atom to establish a perpendicular?
              nbrBond=mol->getBondBetweenAtoms(heavyIdx,nbr1->getIdx());
              if(nbrBond->getIsAromatic() || nbrBond->getBondType()==Bond::DOUBLE){
                nbr2=getAtomNeighborNot(mol,nbr1,heavyAtom);
                nbr2Vect=nbr1Pos.directionVector((*cfi)->getAtomPos(nbr2->getIdx()));
                perpVect = nbr2Vect.crossProduct(nbr1Vect);
              }
            }
            perpVect.normalize();
            // rotate the nbr1Vect 60 degrees about perpVect and we're done:
            tform.SetRotation(60.*M_PI/180.,perpVect);
            dirVect = tform*nbr1Vect;
            hydPos = heavyPos + dirVect*bondLength;
            (*cfi)->setAtomPos(hydIdx, hydPos);
            break;
          case Atom::SP:
            // just lay the H along the vector:
            dirVect=nbr1Vect;
            hydPos = heavyPos + dirVect*bondLength;
            (*cfi)->setAtomPos(hydIdx, hydPos);
            break;
          default:
            // FIX: handle other hybridizations
            // for now, just lay the H along the vector:
            dirVect=nbr1Vect;
            hydPos = heavyPos + dirVect*bondLength;
            (*cfi)->setAtomPos(hydIdx, hydPos);
          }
        }
        break;
      case 3:
        // --------------------------------------------------------------------------
        // Two other neighbors:
        // --------------------------------------------------------------------------
        boost::tie(nbrIdx,endNbrs) = mol->getAtomNeighbors(heavyAtom);
        while(nbrIdx!=endNbrs){
          if(*nbrIdx != hydIdx){
            if(!nbr1) nbr1 = mol->getAtomWithIdx(*nbrIdx);
            else nbr2 = mol->getAtomWithIdx(*nbrIdx);
          }
          nbrIdx++;
        }
        TEST_ASSERT(nbr1);
        TEST_ASSERT(nbr2);
        for (ROMol::ConformerIterator cfi = mol->beginConformers();
             cfi != mol->endConformers(); cfi++) {
          // start along the average of the two vectors:
          heavyPos = (*cfi)->getAtomPos(heavyIdx);
          nbr1Vect = (*cfi)->getAtomPos(nbr1->getIdx()).directionVector(heavyPos);
          nbr2Vect = (*cfi)->getAtomPos(nbr2->getIdx()).directionVector(heavyPos);
          dirVect = nbr1Vect+nbr2Vect;
          dirVect.normalize();

          switch(heavyAtom->getHybridization()){
          case Atom::SP3:
            // get the perpendicular to the neighbors:
            nbrPerp = nbr1Vect.crossProduct(nbr2Vect);
            // and the perpendicular to that:
            rotnAxis = nbrPerp.crossProduct(dirVect);
            // and then rotate about that:
            rotnAxis.normalize();
            tform.SetRotation((109.471/2)*M_PI/180.,rotnAxis);
            dirVect = tform*dirVect;
            hydPos = heavyPos + dirVect*bondLength;
            (*cfi)->setAtomPos(hydIdx, hydPos);
            break;
          case Atom::SP2:
            // don't need to do anything here, the H atom goes right on the
            // direction vector
            hydPos = heavyPos + dirVect*bondLength;
            (*cfi)->setAtomPos(hydIdx, hydPos);
            break;
          default:
            // FIX: handle other hybridizations
            // for now, just lay the H along the neighbor vector;
            hydPos = heavyPos + dirVect*bondLength;
            (*cfi)->setAtomPos(hydIdx, hydPos);
            break;
          } 
        }
        break;
      case 4:
        // --------------------------------------------------------------------------
        // Three other neighbors:
        // --------------------------------------------------------------------------
        boost::tie(nbrIdx,endNbrs) = mol->getAtomNeighbors(heavyAtom);
        while(nbrIdx!=endNbrs){
          if(*nbrIdx != hydIdx){
            if(!nbr1) nbr1 = mol->getAtomWithIdx(*nbrIdx);
            else if(!nbr2) nbr2 = mol->getAtomWithIdx(*nbrIdx);
            else nbr3 = mol->getAtomWithIdx(*nbrIdx);
          }
          nbrIdx++;
        }
        TEST_ASSERT(nbr1);
        TEST_ASSERT(nbr2);
        TEST_ASSERT(nbr3);
        for (ROMol::ConformerIterator cfi = mol->beginConformers();
             cfi != mol->endConformers(); cfi++) {
          // use the average of the three vectors:
          heavyPos = (*cfi)->getAtomPos(heavyIdx);
          nbr1Vect = (*cfi)->getAtomPos(nbr1->getIdx()).directionVector(heavyPos);
          nbr2Vect = (*cfi)->getAtomPos(nbr2->getIdx()).directionVector(heavyPos);
          nbr3Vect = (*cfi)->getAtomPos(nbr3->getIdx()).directionVector(heavyPos);
          dirVect = nbr1Vect+nbr2Vect+nbr3Vect;
          dirVect.normalize();
          hydPos = heavyPos + dirVect*bondLength;
          (*cfi)->setAtomPos(hydIdx, hydPos);
        }
        break;
      default:
        // --------------------------------------------------------------------------
        // FIX: figure out what to do here
        // --------------------------------------------------------------------------
        hydPos = heavyPos + dirVect*bondLength;
        for (ROMol::ConformerIterator cfi = mol->beginConformers();
             cfi != mol->endConformers(); cfi++) {
          (*cfi)->setAtomPos(hydIdx, hydPos);
        }
        break;
      }
    }
  }  // end of unnamed namespace


  namespace MolOps {
    // NOTE that we do not need to check for chiral atoms when adding Hs
    // because the bond order goes from hydrogen implicitly first to
    // hydrogen explicitly last. This is a cyclic permutation, so it
    // doesn't affect the chirality.  (Of course there can only be one H
    // on the atom... otherwise it wouldn't be chiral!)
    ROMol *addHs(const ROMol &mol,bool explicitOnly,bool addCoords){
      RWMol *res = new RWMol(mol);

      // when we hit each atom, clear its computed properties
      // NOTE: it is essential that we not clear the ring info in the
      // molecule's computed properties.  We don't want to have to
      // regenerate that.  This caused Issue210 and Issue212:
      res->clearComputedProps(false);

      // precompute the number of hydrogens we are going to add so that we can
      // pre-allocate the necessary space on the conformations of the molecule
      // for their coordinates
      unsigned int numAddHyds = 0;
      for(ROMol::ConstAtomIterator at=mol.beginAtoms();at!=mol.endAtoms();++at){
        numAddHyds += (*at)->getNumExplicitHs();
        if (!explicitOnly) {
          numAddHyds += (*at)->getNumImplicitHs();
        }
      }
      unsigned int nSize = mol.getNumAtoms() + numAddHyds;
      // now if we want to add coordinates to the hydrogens;
      // loop over the conformations of the molecule and allocate new space
      // for their locations
      for (ROMol::ConformerIterator cfi = res->beginConformers();
           cfi != res->endConformers(); ++cfi) {
        (*cfi)->reserve(nSize);
      }

      for(ROMol::ConstAtomIterator at=mol.beginAtoms();at!=mol.endAtoms();++at){
        unsigned int oldIdx,newIdx;
        oldIdx = (*at)->getIdx();
        Atom *newAt = res->getAtomWithIdx(oldIdx);
        newAt->clearComputedProps();
        // always convert explicit Hs
        for(unsigned int i=0;i<(*at)->getNumExplicitHs();i++){
          newIdx=res->addAtom(new Atom(1),false,true);
          res->addBond(oldIdx,newIdx,Bond::SINGLE);
          res->getAtomWithIdx(newIdx)->updatePropertyCache();
          if(addCoords) setHydrogenCoords(res,newIdx,newAt->getIdx());
        }
        // clear the local property
        newAt->setNumExplicitHs(0);

        if(!explicitOnly){
          // take care of implicits
          for(unsigned int i=0;i<(*at)->getNumImplicitHs();i++){
            newIdx=res->addAtom(new Atom(1),false,true);
            res->addBond(oldIdx,newIdx,Bond::SINGLE);
            // set the isImplicit label so that we can strip these back
            // off later if need be.
            res->getAtomWithIdx(newIdx)->setProp("isImplicit",1, true);
            res->getAtomWithIdx(newIdx)->updatePropertyCache();
            if(addCoords) setHydrogenCoords(res,newIdx,newAt->getIdx());
          }
          // be very clear about implicits not being allowed in this representation
          newAt->setProp("origNoImplicit",newAt->getNoImplicit(), true);
          newAt->setNoImplicit(true);
        }
        // update the atom's derived properties (valence count, etc.)
        newAt->updatePropertyCache();
      }
      return static_cast<ROMol *>(res);
    };


    //
    //  This routine removes hydrogens (and bonds to them) from the molecular graph.
    //  Other Atom and bond indices may be affected by the removal.
    //
    //  NOTES:
    //   - Hydrogens which aren't connected to a heavy atom will not be
    //     removed.  This prevents molecules like "[H][H]" from having
    //     all atoms removed.
    //   - Labelled hydrogen (e.g. atoms with atomic number=1, but mass > 1),
    //     will not be removed.
    //
    ROMol *removeHs(const ROMol &mol,bool implicitOnly,bool updateExplicitCount,bool sanitize){
      Atom *atom;
      unsigned int currIdx=0,origIdx=0;
      RWMol *res = new RWMol(mol);
      while(currIdx < res->getNumAtoms()){
        atom = res->getAtomWithIdx(currIdx);
        atom->setProp("_origAtomIdx",origIdx);
        origIdx++;
        if(atom->getAtomicNum()==1){
          bool removeIt=false;

          if(atom->hasProp("isImplicit")){
            removeIt=true;
          } else if(!implicitOnly && atom->getMass()<1.1){
            ROMol::ADJ_ITER begin,end;
            boost::tie(begin,end) = res->getAtomNeighbors(atom);
            while(begin!=end){
              if(res->getAtomWithIdx(*begin)->getAtomicNum() != 1){
                removeIt=true;
                break;
              }
              begin++;
            }
          }

          if(removeIt){
            ROMol::OEDGE_ITER beg,end;
            boost::tie(beg,end) = res->getAtomBonds(atom);
            ROMol::GRAPH_MOL_BOND_PMAP::type pMap = res->getBondPMap();
            // note the assumption that the H only has one neighbor... I
            // feel no need to handle the case of hypervalent hydrogen!
            // :-) 
            Bond const *bond = pMap[*beg];
            Atom *heavyAtom =bond->getOtherAtom(atom);

            // we'll update the atom's explicit H count if we were told to
            // *or* if the atom is chiral, in which case the H is needed
            // in order to complete the coordination
            // *or* if the atom has the noImplicit flag set:
            if( updateExplicitCount || heavyAtom->getNoImplicit() || 
                heavyAtom->getChiralTag()!=Atom::CHI_UNSPECIFIED ){
              heavyAtom->setNumExplicitHs(heavyAtom->getNumExplicitHs()+1);
            } else {
              // this is a special case related to Issue 228 and the
              // "disappearing Hydrogen" problem discussed in MolOps::adjustHs
              //
              // If we remove a hydrogen from an aromatic N, we need to
              // be *sure* to increment the explicit count, even if the
              // H itself isn't marked as explicit
              if(heavyAtom->getAtomicNum()==7 && heavyAtom->getIsAromatic()
                 && heavyAtom->getFormalCharge()==0){
                heavyAtom->setNumExplicitHs(heavyAtom->getNumExplicitHs()+1);
              }
            }

            // One other consequence of removing the H from the graph is
            // that we may change the ordering of the bonds about a
            // chiral center.  This may change the chiral label at that
            // atom.  We deal with that by explicitly checking here:
            if(heavyAtom->getChiralTag()!=Atom::CHI_UNSPECIFIED){
              INT_LIST neighborIndices;
              unsigned int atomsBeforeHeavy=0;

              boost::tie(beg,end) = res->getAtomBonds(heavyAtom);
              while(beg!=end){
                if(pMap[*beg]->getIdx()!=bond->getIdx()){
                  neighborIndices.push_back(pMap[*beg]->getIdx());
                  if(pMap[*beg]->getOtherAtom(heavyAtom)->getIdx()<heavyAtom->getIdx()){
                    ++atomsBeforeHeavy;
                  }
                }
                ++beg;
              }
              if(atomsBeforeHeavy){
                neighborIndices.insert(++neighborIndices.begin(),bond->getIdx());
              } else {
                neighborIndices.push_front(bond->getIdx());
              }
              
              int nSwaps = heavyAtom->getPerturbationOrder(neighborIndices);
              //std::cerr << " swaps: " << nSwaps << " " << atomsBeforeHeavy << std::endl;
              if(nSwaps%2){
                if(heavyAtom->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW){
                  heavyAtom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
                } else if(heavyAtom->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW){
                  heavyAtom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
                }
              }
            }     
          
            res->removeAtom(atom);
          } else {
            // only increment the atom idx if we don't remove the atom
            currIdx++;
          }
        } else {
          // only increment the atom idx if we don't remove the atom
          currIdx++;
          if(atom->hasProp("origNoImplicit")){
            // we'll get in here if we haven't already processed the atom's implicit
            //  hydrogens. (this is protection for the case that removeHs() is called
            //  multiple times on a single molecule without intervening addHs() calls)
            bool tmpBool;
            atom->getProp("origNoImplicit",tmpBool);    
            atom->setNoImplicit(tmpBool);
            atom->clearProp("origNoImplicit");
          }
        }
      }
      //
      //  If we didn't only remove implicit Hs, which are guaranteed to
      //  be the highest numbered atoms, we may have altered atom indices.
      //  This can screw up derived properties (such as ring members), so
      //  we'll resanitize to ensure we stay in a consistent state.
      //
      if(!implicitOnly && sanitize){
        sanitizeMol(*res);
      }
      return static_cast<ROMol *>(res);
    };

    //
    //  This routine removes explicit hydrogens (and bonds to them) from
    //  the molecular graph and adds them as queries to the heavy atoms
    //  to which they are bound.  If the heavy atoms (or atom queries)
    //  already have hydrogen-count queries, they will be updated.
    //
    //  NOTE:
    //   - Hydrogens which aren't connected to a heavy atom will not be
    //     removed.  This prevents molecules like "[H][H]" from having
    //     all atoms removed.
    //
    ROMol *mergeQueryHs(const ROMol &mol){
      unsigned int currIdx=0;
      RWMol *res = new RWMol(mol);
      while(currIdx < res->getNumAtoms()){
        Atom *atom = res->getAtomWithIdx(currIdx);
        if(atom->getAtomicNum()==1){
          bool removeIt=false;
          ROMol::ADJ_ITER begin,end;
          boost::tie(begin,end) = res->getAtomNeighbors(atom);
          while(begin!=end){
            if(res->getAtomWithIdx(*begin)->getAtomicNum() > 1){
              removeIt=true;
              break;
            }
            begin++;
          }
          if(removeIt){
            ROMol::ADJ_ITER begin,end;
            boost::tie(begin,end) = res->getAtomNeighbors(atom);
            Atom *nbr = res->getAtomWithIdx(*begin);
            //
            //  We've found the neighbor.
            //   1) If the neighbor has no H query already:
            //        - add a generic H query
            //      else:
            //        - do nothing
            //   2) Remove the atom from the molecular graph
            //
            //  Examples:
            //    C[H] -> [C;!H0]
            //    [C;H1][H] -> [C;H1]
            //    [C;H2][H] -> [C;H2]
            //
            // FIX: this is going to behave oddly in the case of a contradictory
            //  SMARTS like: [C;H0][H], where it will give the equivalent of:
            //  [C;H0]  I think this is actually correct, but I can be persuaded
            //  otherwise.
            //
            //  First we'll search for an H query:
            bool hasHQuery=false;
            std::list<QueryAtom::QUERYATOM_QUERY::CHILD_TYPE> childStack;
            QueryAtom::QUERYATOM_QUERY::CHILD_VECT_CI child1;
            for(child1=nbr->getQuery()->beginChildren();
                child1!=nbr->getQuery()->endChildren();
                child1++){
              childStack.push_back(*child1);
            }
            while( !hasHQuery && childStack.size() ){
              QueryAtom::QUERYATOM_QUERY::CHILD_TYPE query = childStack.front();
              childStack.pop_front();
              if(query->getDescription()=="AtomHCount"){
                hasHQuery=true;
              } else {
                for(child1=query->beginChildren();
                    child1!=query->endChildren();
                    child1++){
                  childStack.push_back(*child1);
                }
              }
            }

            if(!hasHQuery){
              ATOM_EQUALS_QUERY *tmp=makeAtomHCountQuery(0);
              tmp->setNegation(true);
              nbr->expandQuery(tmp);
            }
            res->removeAtom(atom);
          } else {
            // only increment the atom idx if we don't remove the atom
            currIdx++;
          }
        } else {
          currIdx++;
        }
      }
      return static_cast<ROMol *>(res);
    };
  }; // end of namespace MolOps
}; // end of namespace RDKit
