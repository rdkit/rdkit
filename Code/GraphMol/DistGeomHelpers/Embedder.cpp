// $Id$
//
//  Copyright (C) 2004-2012 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "Embedder.h"
#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/DistGeomUtils.h>
#include <DistGeom/TriangleSmooth.h>
#include "BoundsMatrixBuilder.h"
#include <ForceField/ForceField.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/AtomIterators.h>

#include <GraphMol/Conformer.h>
#include <RDGeneral/types.h>
#include <RDGeneral/RDLog.h>
#include <RDBoost/Exceptions.h>

#include <Geometry/Transform3D.h>
#include <Numerics/Alignment/AlignPoints.h>
#include <DistGeom/ChiralSet.h>
#include <GraphMol/MolOps.h>
#include <boost/dynamic_bitset.hpp>
#include <iomanip>

#define ERROR_TOL 0.00001

namespace RDKit {
  namespace DGeomHelpers {
    typedef std::pair<int,int> INT_PAIR;
    typedef std::vector<INT_PAIR> INT_PAIR_VECT;
    
    bool _embedPoints(RDGeom::PointPtrVect &positions, 
                      const DistGeom::BoundsMatPtr mmat,
                      bool useRandomCoords,double boxSizeMult,
                      bool randNegEig, 
                      unsigned int numZeroFail, double optimizerForceTol,
                      double basinThresh, int seed, unsigned int maxIterations,
                      const DistGeom::VECT_CHIRALSET &chiralCenters){
      unsigned int nat = positions.size();
      if(maxIterations==0){
        maxIterations=10*nat;
      }
      RDNumeric::DoubleSymmMatrix distMat(nat, 0.0);

      // The basin threshold just gets us into trouble when we're using
      // random coordinates since it ends up ignoring 1-4 (and higher)
      // interactions. This causes us to get folded-up (and self-penetrating)
      // conformations for large flexible molecules
      if(useRandomCoords) basinThresh=1e8;

      RDKit::double_source_type *rng=0;
      RDKit::rng_type *generator;
      RDKit::uniform_double *distrib;
      if(seed>0){
        generator=new RDKit::rng_type(42u);
        generator->seed(seed);
        distrib=new RDKit::uniform_double(0.0,1.0);
        rng = new RDKit::double_source_type(*generator,*distrib);
      } else {
        rng = &RDKit::getDoubleRandomSource();
      }
      
      bool gotCoords = false;
      unsigned int iter = 0;
      double largestDistance=-1.0;
      while ((gotCoords == false) && (iter < maxIterations)) {
        ++iter;
        if(!useRandomCoords){
          largestDistance=DistGeom::pickRandomDistMat(*mmat, distMat, *rng);
          gotCoords = DistGeom::computeInitialCoords(distMat, positions,*rng,
                                                     randNegEig, numZeroFail);
        } else {
          double boxSize;
          if(boxSizeMult>0){
            boxSize=5.*boxSizeMult;
          } else {
            boxSize=-1*boxSizeMult;
          }
          gotCoords = DistGeom::computeRandomCoords(positions,boxSize,*rng);
        }
      }
      if(seed>0 && rng){
        delete rng;
        delete generator;
        delete distrib;
      }
      if (gotCoords) {
        ForceFields::ForceField *field = DistGeom::constructForceField(*mmat, positions,
                                                                       chiralCenters,
                                                                       1.0, 0.1, 
                                                                       0,basinThresh);
        unsigned int nPasses=0;
        field->initialize();
        //std::cerr<<"FIELD E: "<<field->calcEnergy()<<std::endl;
        if(field->calcEnergy() > ERROR_TOL){
          int needMore = 1;
          while(needMore){
            needMore = field->minimize(200,optimizerForceTol);
            //needMore = field->minimize(200,1e-6);
            ++nPasses;
          }
        }
        //std::cerr<<"   "<<field->calcEnergy()<<" after npasses: "<<nPasses<<std::endl;
        // now redo the minimization if we have a chiral center, this
        // time removing the chiral constraints and
        // increasing the weight on the fourth dimension
        if (chiralCenters.size()>0 || useRandomCoords) {
          ForceFields::ForceField *field2 = DistGeom::constructForceField(*mmat, positions,
                                                                          chiralCenters,
                                                                          0.1, 1.0, 0,
                                                                          basinThresh);
          field2->initialize();
          //std::cerr<<"FIELD2 E: "<<field2->calcEnergy()<<std::endl;
          if(field2->calcEnergy() > ERROR_TOL){
            int needMore = 1;
            int nPasses2=0;
            while(needMore){
              needMore = field2->minimize(200,optimizerForceTol);
              ++nPasses2;
            }
            //std::cerr<<"   "<<field2->calcEnergy()<<" after npasses: "<<nPasses2<<std::endl;
          }
          delete field2;
        }
        delete field;
      }
      return gotCoords;
    }
    
    void _findChiralSets(const ROMol &mol, DistGeom::VECT_CHIRALSET &chiralCenters) {
      ROMol::ConstAtomIterator ati;
      INT_PAIR_VECT nbrs;
      ROMol::OEDGE_ITER beg,end;
      Atom *oatom;
      for (ati = mol.beginAtoms(); ati != mol.endAtoms(); ati++) {
        if ((*ati)->getAtomicNum() != 1) { //skip hydrogens
          if ((*ati)->hasProp("_CIPCode")) { 
            // make a chiral set from the neighbors
            nbrs.clear();
            nbrs.reserve(4);
            // find the neighbors of this atom and enter them into the
            // nbr list along with their CIPRanks
            boost::tie(beg,end) = mol.getAtomBonds(*ati);
            while (beg != end) {
              oatom = mol[*beg]->getOtherAtom(*ati);
              int rank;
              oatom->getProp("_CIPRank", rank);
              INT_PAIR rAid(rank, oatom->getIdx());
              nbrs.push_back(rAid);
              ++beg;
            }
            // if we have less than 4 heavy atoms as neighbors,
            // we need to include the chiral center into the mix
            // we should at least have 3 though
            bool includeSelf = false;
            CHECK_INVARIANT(nbrs.size() >= 3, "Cannot be a chiral center");

            std::sort(nbrs.begin(), nbrs.end());
            if (nbrs.size() < 4) {
              int rank;
              (*ati)->getProp("_CIPRank", rank);
              INT_PAIR rAid(rank, (*ati)->getIdx());
              nbrs.insert(nbrs.begin(), rAid); 
              includeSelf = true;
            }
                            
            // now create a chiral set and set the upper and lower bound on the volume
            std::string cipCode;
            (*ati)->getProp("_CIPCode", cipCode);
            
            if (cipCode == "S") { 
              // postive chiral volume
              DistGeom::ChiralSet *cset = new DistGeom::ChiralSet(nbrs[0].second,
                                                                  nbrs[1].second,
                                                                  nbrs[2].second,
                                                                  nbrs[3].second,
                                                                  5.0, 100.0);
              DistGeom::ChiralSetPtr cptr(cset);
              chiralCenters.push_back(cptr);
            } else {
              DistGeom::ChiralSet *cset = new DistGeom::ChiralSet(nbrs[0].second,
                                                                  nbrs[1].second,
                                                                  nbrs[2].second,
                                                                  nbrs[3].second,
                                                                  -100.0, -5.0);
              DistGeom::ChiralSetPtr cptr(cset);
              chiralCenters.push_back(cptr);
            }
          } // if block -chirality check
        } // if block - heavy atom check
      } // for loop over atoms
    } // end of _findChiralSets

    void _fillAtomPositions(RDGeom::Point3DConstPtrVect &pts, const Conformer &conf) {
      unsigned int na = conf.getNumAtoms();
      pts.clear();
      unsigned int ai;
      pts.reserve(na);
      for (ai = 0; ai < na; ++ai) {
        pts.push_back(&conf.getAtomPos(ai));
      }
    }
        
    bool _isConfFarFromRest(const ROMol &mol, const Conformer &conf,
                            double threshold) {
      // NOTE: it is tempting to use some triangle inequality to prune
      // conformations here but some basic testing has shown very
      // little advantage and given that the time for pruning fades in
      // comparison to embedding - we will use a simple for loop below
      // over all conformation until we find a match
      ROMol::ConstConformerIterator confi;

      RDGeom::Point3DConstPtrVect refPoints, prbPoints;
      _fillAtomPositions(refPoints, conf);

      bool res = true;
      unsigned int na = conf.getNumAtoms();
      double ssrThres = na*threshold*threshold;

      RDGeom::Transform3D trans;
      double ssr;
      for (confi = mol.beginConformers(); confi != mol.endConformers(); confi++) {
        _fillAtomPositions(prbPoints, *(*confi));
        ssr = RDNumeric::Alignments::AlignPoints(refPoints, prbPoints, trans);
        if (ssr < ssrThres) {
          res = false;
          break;
        }
      }
      return res;
    }

    int EmbedMolecule(ROMol &mol, unsigned int maxIterations, int seed,
                      bool clearConfs,
                      bool useRandomCoords,double boxSizeMult,
                      bool randNegEig, unsigned int numZeroFail,
                      const std::map<int,RDGeom::Point3D> *coordMap,
                      double optimizerForceTol,
                      bool ignoreSmoothingFailures,
                      double basinThresh){

      INT_VECT confIds;
      confIds=EmbedMultipleConfs(mol,1,maxIterations,seed,clearConfs,
                                 useRandomCoords,boxSizeMult,randNegEig,
                                 numZeroFail,-1.0,coordMap,optimizerForceTol,
                                 ignoreSmoothingFailures,basinThresh);

      int res;
      if(confIds.size()){
        res=confIds[0];
      } else {
        res=-1;
      }
      return res;
    }

    void adjustBoundsMatFromCoordMap(DistGeom::BoundsMatPtr mmat,unsigned int nAtoms,
                                     const std::map<int,RDGeom::Point3D> *coordMap){
      // std::cerr<<std::endl;
      // for(unsigned int i=0;i<nAtoms;++i){
      //   for(unsigned int j=0;j<nAtoms;++j){
      //     std::cerr<<"  "<<std::setprecision(3)<<mmat->getVal(i,j);
      //   }
      //   std::cerr<<std::endl;
      // }
      // std::cerr<<std::endl;
      for(std::map<int,RDGeom::Point3D>::const_iterator iIt=coordMap->begin();
          iIt!=coordMap->end();++iIt){
        int iIdx=iIt->first;
        const RDGeom::Point3D &iPoint=iIt->second;
        
        std::map<int,RDGeom::Point3D>::const_iterator jIt=iIt;
        while(++jIt != coordMap->end()){
          int jIdx=jIt->first;
          const RDGeom::Point3D &jPoint=jIt->second;
          double dist=(iPoint-jPoint).length();
          mmat->setUpperBound(iIdx,jIdx,dist);
          mmat->setLowerBound(iIdx,jIdx,dist);
        }
      }
      // std::cerr<<std::endl;
      // for(unsigned int i=0;i<nAtoms;++i){
      //   for(unsigned int j=0;j<nAtoms;++j){
      //     std::cerr<<"  "<<std::setprecision(3)<<mmat->getVal(i,j);
      //   }
      //   std::cerr<<std::endl;
      // }
      // std::cerr<<std::endl;
    }

    
    INT_VECT EmbedMultipleConfs(ROMol &mol, unsigned int numConfs,
                                unsigned int maxIterations, 
                                int seed, bool clearConfs, 
                                bool useRandomCoords,double boxSizeMult,
                                bool randNegEig, unsigned int numZeroFail,
                                double pruneRmsThresh,
                                const std::map<int,RDGeom::Point3D>  *coordMap,
                                double optimizerForceTol,
                                bool ignoreSmoothingFailures,
                                double basinThresh){
      if(!mol.getNumAtoms()){
        throw ValueErrorException("molecule has no atoms");
      }

      INT_VECT fragMapping;
      std::vector<ROMOL_SPTR> molFrags=MolOps::getMolFrags(mol,true,&fragMapping);
      if(molFrags.size()>1 && coordMap){
        BOOST_LOG(rdWarningLog)<<"Constrained conformer generation (via the coordMap argument) does not work with molecules that have multiple fragments."<<std::endl;
        coordMap=0;
      }
      std::vector< Conformer * > confs;
      confs.reserve(numConfs);
      for(unsigned int i=0;i<numConfs;++i){
        confs.push_back(new Conformer(mol.getNumAtoms()));
      }
      boost::dynamic_bitset<> confsOk(numConfs);
      confsOk.set();

      if (clearConfs) {
        mol.clearConformers();
      }
      INT_VECT res;
      
      for(unsigned int fragIdx=0;fragIdx<molFrags.size();++fragIdx){
        ROMOL_SPTR piece=molFrags[fragIdx];
        unsigned int nAtoms = piece->getNumAtoms();
        DistGeom::BoundsMatrix *mat = new DistGeom::BoundsMatrix(nAtoms);
        DistGeom::BoundsMatPtr mmat(mat);
        initBoundsMat(mmat);
      
        double tol=0.0;
        setTopolBounds(*piece, mmat, true, false);
        if(coordMap){
          adjustBoundsMatFromCoordMap(mmat,nAtoms,coordMap);
          tol=0.05;
        }
        if (!DistGeom::triangleSmoothBounds(mmat,tol)) {
          // ok this bound matrix failed to triangle smooth - re-compute the bounds matrix 
          // without 15 bounds and with VDW scaling
          initBoundsMat(mmat);
          setTopolBounds(*piece, mmat, false, true);

          if(coordMap){
            adjustBoundsMatFromCoordMap(mmat,nAtoms,coordMap);
          }

          // try triangle smoothing again 
          if (!DistGeom::triangleSmoothBounds(mmat,tol)) {
            // ok, we're not going to be able to smooth this,
            if(ignoreSmoothingFailures){
              // proceed anyway with the more relaxed bounds matrix
              initBoundsMat(mmat);
              setTopolBounds(*piece, mmat, false, true);

              if(coordMap){
                adjustBoundsMatFromCoordMap(mmat,nAtoms,coordMap);
              }
            } else {
              BOOST_LOG(rdWarningLog)<<"Could not triangle bounds smooth molecule."<<std::endl;
              return res;
            }
          }
        }
#if 0
        for(unsigned int li=0;li<piece->getNumAtoms();++li){
          for(unsigned int lj=li+1;lj<piece->getNumAtoms();++lj){
            std::cerr<<" ("<<li<<","<<lj<<"): "<<mat->getLowerBound(li,lj)<<" -> "<<mat->getUpperBound(li,lj)<<std::endl;
          }
        }
#endif
        // find all the chiral centers in the molecule
        DistGeom::VECT_CHIRALSET chiralCenters;
        MolOps::assignStereochemistry(*piece);
        _findChiralSets(*piece, chiralCenters);

        // if we have any chiral centers or are using random coordinates, we will 
        // first embed the molecule in four dimensions, otherwise we will use 3D 
        bool fourD = false;
        if (useRandomCoords || chiralCenters.size() > 0) {
          fourD = true;
        }
        RDGeom::PointPtrVect positions;
        for (unsigned int i = 0; i < nAtoms; ++i) {
          if(fourD){
            positions.push_back(new RDGeom::PointND(4));
          } else {
            positions.push_back(new RDGeom::Point3D());
          }
        }
        for (unsigned int ci=0; ci<numConfs; ci++) {
          if(!confsOk[ci]){
            // if one of the fragments here has already failed, there's no
            // sense in embedding this one
            continue;
          }
          bool gotCoords = _embedPoints(positions, mmat,
                                        useRandomCoords,boxSizeMult,
                                        randNegEig, numZeroFail,
                                        optimizerForceTol,
                                        basinThresh, (ci+1)*seed,
                                        maxIterations, chiralCenters);
          if (gotCoords) {
            Conformer *conf = confs[ci];
            unsigned int fragAtomIdx=0;
            for (unsigned int i = 0; i < mol.getNumAtoms();++i){
              if(fragMapping[i]==static_cast<int>(fragIdx) ){
                conf->setAtomPos(i, RDGeom::Point3D((*positions[fragAtomIdx])[0],
                                                    (*positions[fragAtomIdx])[1],
                                                    (*positions[fragAtomIdx])[2]));
                ++fragAtomIdx;
              }
            }
          } else {
            confsOk[ci]=0;
          }
        }
        for (unsigned int i = 0; i < nAtoms; ++i) {
          delete positions[i];
        }
      }
      for(unsigned int ci=0;ci<confs.size();++ci){
        Conformer *conf = confs[ci];
        if(confsOk[ci]){
          // check if we are pruning away conformations and 
          // a closeby conformation has already been chosen :
          if (pruneRmsThresh > 0.0 && 
              !_isConfFarFromRest(mol, *conf, pruneRmsThresh)) { 
            delete conf;
          } else {
            int confId = (int)mol.addConformer(conf, true);
            res.push_back(confId);
          }
        } else {
          delete conf;
        }
      }
      return res;
    }
  }
}
    
