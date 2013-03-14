// $Id$
//
//  Copyright (C) 2003-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RDDepictor.h"
#include <RDGeneral/types.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
#include <math.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Rings.h>
#include <Geometry/point.h>
#include <Geometry/Transform2D.h>
#include "EmbeddedFrag.h"
#include "DepictUtils.h"
#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>

namespace RDDepict {
  namespace DepictorLocal {
    // arings: indices of atoms in rings
    void embedFusedSystems(const RDKit::ROMol &mol,
                           const RDKit::VECT_INT_VECT &arings,
                           std::list<EmbeddedFrag> &efrags) {
      RDKit::INT_INT_VECT_MAP neighMap;
      RingUtils::makeRingNeighborMap(arings, neighMap);
      
      RDKit::INT_VECT fused;
      int cnrs = arings.size();
      boost::dynamic_bitset<> fusDone(cnrs);
      int curr = 0;
      
      while (curr < cnrs) {
        // embed all ring and fused ring systems
        fused.resize(0);
        RingUtils::pickFusedRings(curr, neighMap, fused, fusDone);
        RDKit::VECT_INT_VECT frings;
        frings.reserve(fused.size());
        for (RDKit::INT_VECT_CI rid = fused.begin();
             rid != fused.end(); ++rid) {
          frings.push_back(arings[*rid]);
        }
        EmbeddedFrag efrag(&mol, frings);
        efrag.setupNewNeighs(); 
        efrags.push_back(efrag);
        int rix;
        for (rix = 0; rix < cnrs; ++rix) {
          if (!fusDone[rix]) {
            curr = rix;
            break;
          }
        }
        if (rix == cnrs) {
          break;
        }
      }
      
    }
    
    void embedCisTransSystems(const RDKit::ROMol &mol,
                              std::list<EmbeddedFrag> &efrags) {
      for (RDKit::ROMol::ConstBondIterator cbi = mol.beginBonds();
           cbi != mol.endBonds(); ++cbi) {
        // check if this bond is in a cis/trans double bond
        // and it is not a ring bond
        if (((*cbi)->getBondType() == RDKit::Bond::DOUBLE) // this is a double bond
            && ((*cbi)->getStereo() > RDKit::Bond::STEREOANY) // and has stereo chemistry specified
            && (!(*cbi)->getOwningMol().getRingInfo()->numBondRings((*cbi)->getIdx())) ){ // not in a ring
          if((*cbi)->getStereoAtoms().size() != 2 ){
            BOOST_LOG(rdWarningLog)<<"WARNING: bond found with stereo spec but no stereo atoms"<<std::endl;
            continue;
          }
          EmbeddedFrag efrag(*cbi);
          efrag.setupNewNeighs(); 
          efrags.push_back(efrag);
        }
      }
    }
    
    RDKit::INT_LIST getNonEmbeddedAtoms(const RDKit::ROMol &mol,
                                        const std::list<EmbeddedFrag> &efrags) {
      RDKit::INT_LIST res;
      boost::dynamic_bitset<> done(mol.getNumAtoms());
      for (std::list<EmbeddedFrag>::const_iterator efri = efrags.begin();
           efri != efrags.end(); efri++) {
        const INT_EATOM_MAP &oatoms = efri->GetEmbeddedAtoms();
        for (INT_EATOM_MAP_CI ori = oatoms.begin(); ori != oatoms.end(); ori++) {
          done[ori->first]=1;
        }
      }
      for (RDKit::ROMol::ConstAtomIterator ai = mol.beginAtoms();
           ai != mol.endAtoms(); ai++) {
        int aid = (*ai)->getIdx();
        if (!done[aid]){
          res.push_back(aid);
        }
      }
      return res;
    }
    
    // find the largest fragments that is not done yet (
    //  i.e. merged with the master fragments)
    // if do not find anything we return efrags.end()
    std::list<EmbeddedFrag>::iterator
    _findLargestFrag(std::list<EmbeddedFrag> &efrags) {
      std::list<EmbeddedFrag>::iterator mfri;
      int msiz = 0;
      for (std::list<EmbeddedFrag>::iterator efri = efrags.begin();
           efri != efrags.end(); efri++) {
        if ((!efri->isDone()) && (efri->Size() > msiz) ) {
          msiz = efri->Size();
          mfri = efri;
        }
      }
      if (msiz == 0) {
        mfri = efrags.end();
      }
      return mfri;
    }
    
    void _shiftCoords(std::list<EmbeddedFrag> &efrags) {
      // shift the coordinates if there are multiple fragments
      // so that the fargments do not overlap each other
      
      for(std::list<EmbeddedFrag>::iterator efi=efrags.begin();
          efi!=efrags.end();efi++){
        efi->computeBox();
      }
      std::list<EmbeddedFrag>::iterator eri = efrags.begin();
      double xmax = eri->getBoxPx();
      double xmin = eri->getBoxNx();
      double ymax = eri->getBoxPy();
      double ymin = eri->getBoxNy();
      
      ++eri;
      while (eri != efrags.end()) {
        bool xshift = true;
        if (xmax+xmin > ymax+ymin) {
          xshift = false;
        }
        double xn = eri->getBoxNx();
        double xp = eri->getBoxPx();
        double yn = eri->getBoxNy();
        double yp = eri->getBoxPy();
        RDGeom::Point2D shift(0.0,0.0);
        if (xshift) {
          shift.x = xmax + xn + 1.0;
          shift.y = 0.0;
          xmax += xp + xn + 1.0;
        } else {
          shift.x = 0.0;
          shift.y = ymax + yn + 1.0;
          ymax += yp + yn + 1.0; 
        }
        eri->Translate(shift);
        
        ++eri;
      }
    }
  }

  void computeInitialCoords(RDKit::ROMol &mol,
                            const RDGeom::INT_POINT2D_MAP *coordMap, 
                            std::list<EmbeddedFrag> &efrags) {
    RDKit::INT_VECT atomRanks;
    atomRanks.resize(mol.getNumAtoms());

    RDKit::MolOps::assignStereochemistry(mol, false);

    efrags.clear();
    RDKit::VECT_INT_VECT arings;

    // first find all the rings
    RDKit::MolOps::symmetrizeSSSR(mol, arings);

    // user specfied coordinates exist
    bool preSpec = false;
    // first embed any atoms for which the coordinates have been specified.
    if ((coordMap) && (coordMap->size() > 1) ) {
      EmbeddedFrag efrag(&mol, *coordMap);
      // add this to the list of embedded fragments 
      efrags.push_back(efrag);
      preSpec = true;
    }
    
    if (arings.size() > 0) {
      // first deal with the fused rings
      DepictorLocal::embedFusedSystems(mol, arings, efrags);
    }
    // deal with any cis/trans systems
    DepictorLocal::embedCisTransSystems(mol, efrags);
    // now get the atoms that are not yet embedded in either a cis/trans system
    // or a ring system (or simply the first atom)
    RDKit::INT_LIST nratms = DepictorLocal::getNonEmbeddedAtoms(mol, efrags);
    std::list<EmbeddedFrag>::iterator mri; 
    if (preSpec) {
      // if the user specified coordinates on some of the atoms use that as 
      // as the starting fragment and it should be at the beginning of the vector
      mri = efrags.begin();
    } else {
      // otherwise - find the largest fragment that was embedded
      mri = DepictorLocal::_findLargestFrag(efrags);
    }

    while ((mri != efrags.end()) || (nratms.size() > 0)) {
      if (mri == efrags.end()) {
        // we are out of embedded fragments, if there are any 
        // non embedded atoms use them to start a fragment
        int mrank, rank;
        mrank = static_cast<int>(RDKit::MAX_INT);
        RDKit::INT_LIST_I nri, mnri;
        for (nri = nratms.begin(); nri != nratms.end(); nri++) {
          rank=atomRanks[*nri];
          rank *= 1000;
          // use the atom index as well so that we at least
          // get reproduceable depictions in cases where things
          // have identical ranks.
          rank += *nri;
          if (rank < mrank) {
            mrank = rank;
            mnri = nri;
          }
        }
        EmbeddedFrag efrag((*mnri), &mol);
        nratms.erase(mnri);
        efrags.push_back(efrag);
        mri = efrags.end();
        mri--;
      }
      mri->markDone();
      mri->expandEfrag(nratms, efrags); 
      mri = DepictorLocal::_findLargestFrag(efrags);
    }
    // at this point any remaining efrags should belong individual fragments in the molecule
  }

  unsigned int copyCoordinate(RDKit::ROMol &mol, std::list<EmbeddedFrag> &efrags, bool clearConfs) {
    // create a conformation to store the coordinates and add it to the molecule 
    RDKit::Conformer *conf = new RDKit::Conformer(mol.getNumAtoms());
    conf->set3D(false);
    std::list<EmbeddedFrag>::iterator eri;
    for (eri = efrags.begin(); eri != efrags.end(); eri++) {
      const INT_EATOM_MAP &eatoms = eri->GetEmbeddedAtoms();
      INT_EATOM_MAP_CI eai;
      for (eai = eatoms.begin(); eai != eatoms.end(); eai++) {
        int aid = eai->first;
        RDGeom::Point2D cr = eai->second.loc;
        RDGeom::Point3D fcr(cr.x, cr.y, 0.0);
        conf->setAtomPos(aid, fcr);
      }
    }
    unsigned int confId = 0;
    if (clearConfs) {
      // clear all the conformation on the molecules and assign conf ID 0 to this 
      // conformation
      mol.clearConformers();
      conf->setId(confId);
      // conf ID has already been set in this case to 0 - not other
      // confs on the molecule at this point
      mol.addConformer(conf);
    } else {
      // let add conf assign a conformation ID for the conformation
      confId = mol.addConformer(conf, true);
    }
    return confId;
  }
  //
  //
  // 50,000 foot algorithm:
  //   1) Find rings
  //   2) Find fused systems
  //   3) embed largest fused system
  //   4) foreach unfinished atom:
  //      1) find neighbors
  //      2) if neighbor is non-ring atom, embed it; otherwise merge the
  //         ring system
  //      3) add all atoms just merged/embedded to unfinished atom list
  //      
  //
  unsigned int compute2DCoords(RDKit::ROMol &mol,
                               const RDGeom::INT_POINT2D_MAP *coordMap,
                               bool canonOrient, bool clearConfs,
                               unsigned int nFlipsPerSample,
                               unsigned int nSamples,
                               int sampleSeed, bool permuteDeg4Nodes) {
    
    // storage for pieces of a molecule/s that are embedded in 2D
    std::list<EmbeddedFrag> efrags;
    computeInitialCoords(mol, coordMap, efrags);

    std::list<EmbeddedFrag>::iterator eri;
    // perform random sampling here to improve the density
    for (eri = efrags.begin(); eri != efrags.end(); eri++) {
      // either sample the 2D space by randomly flipping rotatable
      // bonds in the structure or flip onyl bonds along the shortest
      // path between colliding atoms - don't do both
      if ((nSamples > 0) && (nFlipsPerSample > 0)) {
        eri->randomSampleFlipsAndPermutations(nFlipsPerSample, nSamples,
                                              sampleSeed, 0, 0.0,
                                              permuteDeg4Nodes);
      } else {
        eri->removeCollisionsBondFlip();
      }
    }
    for (eri = efrags.begin(); eri != efrags.end(); eri++) {
      // if there are any remaining collisions
      eri->removeCollisionsOpenAngles();
      eri->removeCollisionsShortenBonds();
    }
    if (!coordMap || !coordMap->size() ) {
      if (canonOrient && efrags.size() ) {
        // if we do not have any prespecified coordinates - canonicalize
        // the orientation of the fragment so that the longest axes fall
        // along the x-axis etc.
        for (eri = efrags.begin(); eri != efrags.end(); eri++) {
          eri->canonicalizeOrientation();
        }
      }
    }
    DepictorLocal::_shiftCoords(efrags);
    // create a confomation on the moelcule and copy the coodinates
    unsigned int cid = copyCoordinate(mol, efrags, clearConfs);

    // special case for a single-atom coordMap template 
    if ((coordMap) && (coordMap->size() == 1) ) {
      RDKit::Conformer &conf=mol.getConformer(cid);
      RDGeom::INT_POINT2D_MAP::const_iterator cRef=coordMap->begin();
      RDGeom::Point3D confPos=conf.getAtomPos(cRef->first);
      RDGeom::Point2D refPos=cRef->second;
      refPos.x -= confPos.x;
      refPos.y -= confPos.y;
      for(unsigned int i=0;i<conf.getNumAtoms();++i){
        confPos = conf.getAtomPos(i);
        confPos.x += refPos.x;
        confPos.y += refPos.y;
        conf.setAtomPos(i,confPos);
      }
    }


    return cid;
  }

  //! \brief Compute the 2D coordinates such the interatom distances
  //!        mimic those in a distance matrix
  /*!  
    This function generates 2D coordinates such that the inter atom
    distance mimic those specified via dmat. This is done by randomly
    sampling(flipping) the rotatable bonds in the molecule and
    evaluating a cost function which contains two components. The
    first component is the sum of inverse of the squared inter-atom
    distances, this helps in spreading the atoms far from each
    other. The second component is the sum of squares of the
    difference in distance between those in dmat and the generated
    structure.  The user can adjust the relative importance of the two
    components via a adjustable paramter (see below)

    ARGUMENTS:
    \param mol - molecule involved in the frgament

    \param dmat - the distance matrix we want to mimic, this is
                  symmteric N by N matrix when N is the number of
                  atoms in mol. All ngative entries in dmat are
                  ignored.

    \param canonOrient - canonicalze the orientation after the 2D
                         embedding is done

    \param clearConfs - clear any previously existing conformations on
                        mol before adding a conformation

    \param weightDistMat - A value between 0.0 and 1.0, this
                           determines the importance of mimicing the
                           the inter atoms distances in dmat. (1.0 -
                           weightDistMat) is the weight associated to
                           spreading out the structure (density) in
                           the cost function

    \param nFlipsPerSample - the number of rotatable bonds that are
                             randomly flipped for each sample

    \param nSample - the number of samples

    \param sampleSeed - seed for the random sampling process
  */
  unsigned int compute2DCoordsMimicDistMat(RDKit::ROMol &mol,
                                           const DOUBLE_SMART_PTR *dmat,
                                           bool canonOrient,
                                           bool clearConfs, double weightDistMat,
                                           unsigned int nFlipsPerSample,
                                           unsigned int nSamples,
                                           int sampleSeed,
                                           bool permuteDeg4Nodes){
    // storage for pieces of a molecule/s that are embedded in 2D
    std::list<EmbeddedFrag> efrags;
    computeInitialCoords(mol, 0, efrags);
    
    // now perform random flips of rotatable bonds so taht we can sample the space
    // and try to mimic the distances in dmat
    std::list<EmbeddedFrag>::iterator eri;
    for (eri = efrags.begin(); eri != efrags.end(); eri++) {
      eri->randomSampleFlipsAndPermutations(nFlipsPerSample, nSamples, sampleSeed, dmat, 
                                            weightDistMat, permuteDeg4Nodes);
    }
    if (canonOrient && efrags.size()) {
      // canonicalize the orientation of the fragment so that the
      // longest axes fall along the x-axis etc.
      for (eri = efrags.begin(); eri != efrags.end(); eri++) {
        eri->canonicalizeOrientation();
      }
    }

    DepictorLocal::_shiftCoords(efrags);
    // create a confomation on the moelcule and copy the coodinates
    unsigned int cid = copyCoordinate(mol, efrags, clearConfs);
    return cid;
  }
    
    
    
}
          





