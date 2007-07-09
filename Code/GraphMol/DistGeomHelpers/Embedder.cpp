//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
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

#include <Geometry/Transform3D.h>
#include <Numerics/Alignment/AlignPoints.h>
#include <DistGeom/EmbedObject.h>
#include <GraphMol/MolOps.h>

#define ERROR_TOL 0.00001

namespace RDKit {
  namespace DGeomHelpers {
    typedef std::pair<int,int> INT_PAIR;
    typedef std::vector<INT_PAIR> INT_PAIR_VECT;

    
    bool _embedPoints(RDGeom::PointPtrVect &positions, 
		      const DistGeom::BoundsMatPtr mmat, bool randNegEig, 
		      unsigned int numZeroFail, double optimizerForceTol,
		      double basinThresh, int seed, unsigned int maxIterations,
		      const DistGeom::VECT_CHIRALSET &chiralCenters) {
      
      
      unsigned int nat = positions.size();
      RDNumeric::DoubleSymmMatrix distMat(nat, 0.0);
      
      bool gotCoords = false;
      unsigned int iter = 0;
      
      while ((gotCoords == false) && (iter < maxIterations)) {
	iter += 1;
	if (seed > 0) {
	  pickRandomDistMat(*mmat, distMat, iter*seed);
	} else {
	  pickRandomDistMat(*mmat, distMat);
	}
	gotCoords = DistGeom::computeInitialCoords(distMat, positions,
						   randNegEig, numZeroFail);
      }
      if (gotCoords) {
	ForceFields::ForceField *field = DistGeom::constructForceField(*mmat, positions, chiralCenters,
								       1.0, 0.1, 
								       0,basinThresh);
	if (field) {
          field->initialize();
          if(field->calcEnergy() > ERROR_TOL){
            int needMore = 1;
            while(needMore){
              needMore = field->minimize(200,optimizerForceTol);
            }
          }
          delete field;
        }
	// now if we have a chiral center redo the minimization but this time removing the chiral constraints and
	// increasing the weight on the fourth dimension
	if (chiralCenters.size() > 0) {
	  ForceFields::ForceField *field = DistGeom::constructForceField(*mmat, positions, chiralCenters,
									 0.0, 1.0, 0,basinThresh);
	  if (field) {
	    field->initialize();
	    if(field->calcEnergy() > ERROR_TOL){
	      int needMore = 1;
	      while(needMore){
		needMore = field->minimize(200,optimizerForceTol);
	      }
	    }
	    delete field;
	  }
	}
	  
      }
      return gotCoords;
    }
    
    
    void findChiralSets(const ROMol &mol, DistGeom::VECT_CHIRALSET &chiralCenters) {

      
      ROMol::ConstAtomIterator ati;
      INT_PAIR_VECT nbrs;
      ROMol::OEDGE_ITER beg,end;
      ROMol::GRAPH_MOL_BOND_PMAP::const_type pMap = mol.getBondPMap();
      Atom *oatom;
      for (ati = mol.beginAtoms(); ati != mol.endAtoms(); ati++) {
	if ((*ati)->getAtomicNum() != 1) { //skip hydrogens
	  if ((*ati)->hasProp("_CIPCode")) { 
	    // make a chiral set from the neighbors
	    nbrs.clear();
	    nbrs.reserve(4);
	    // find the neighbors of this atom and enter them into the nbr list along with their CIPRanks
	    boost::tie(beg,end) = mol.getAtomBonds(*ati);
	    while (beg != end) {
	      oatom = pMap[*beg]->getOtherAtom(*ati);
	      //if (oatom->getAtomicNum() != 1) { // skip hydrogens
	      int rank;
	      oatom->getProp("_CIPRank", rank);
	      INT_PAIR rAid(rank, oatom->getIdx());
	      nbrs.push_back(rAid);
	      //}
	      beg++;
	    }
	    // if we have less than 4 heavy atoms as neighbors, we need to include the chiral center into the mix
	    // we should atleast have 3 though
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
	      DistGeom::ChiralSet *cset = new DistGeom::ChiralSet(nbrs[0].second, nbrs[1].second, nbrs[2].second,
								  nbrs[3].second, 3.0, 100.0);
	      DistGeom::ChiralSetPtr cptr(cset);
	      chiralCenters.push_back(cptr);
	    } else {
	      DistGeom::ChiralSet *cset = new DistGeom::ChiralSet(nbrs[0].second, nbrs[1].second, nbrs[2].second,
								  nbrs[3].second, -100.0, -3.0);
	      DistGeom::ChiralSetPtr cptr(cset);
	      chiralCenters.push_back(cptr);
	    }
	  } // if block -chirality check
	} // if block - heavy atom check
      } // for loop over atoms
      
    }
	      
    int EmbedMolecule(ROMol &mol, unsigned int maxIterations, int seed, bool clearConfs,
                      bool randNegEig, unsigned int numZeroFail,
		      double optimizerForceTol,double basinThresh) {
      unsigned int nat = mol.getNumAtoms();
      DistGeom::BoundsMatrix *mat = new DistGeom::BoundsMatrix(nat);
      DistGeom::BoundsMatPtr mmat(mat);
      initBoundsMat(mmat);

      // set the bounds using topology - including 15 bounds
      setTopolBounds(mol, mmat, true, false);
      
      if (!DistGeom::triangleSmoothBounds(mmat)) {
        // ok this bound matrix failed to triangle smooth - re-compute the bounds matrix 
        // with out 15 bound and with VDW scaling
        initBoundsMat(mmat);
        setTopolBounds(mol, mmat, false, true);
        // try triangle smoothing again - give up if we fail
        if (!DistGeom::triangleSmoothBounds(mmat)) {
          //BOOST_LOG(rdDebugLog) << "failed 14 and vdw scaling\n";
          return -1;
        }
      }
      
      // find all the chiral centers in the molecule
      DistGeom::VECT_CHIRALSET chiralCenters;
      MolOps::assignAtomChiralCodes(mol);
      findChiralSets(mol, chiralCenters);
      // if we have any chiral centers we will first embed the molecule in four dimensions
      // other we will use 3D 
      RDGeom::PointPtrVect positions;
      bool fourD = false;
      unsigned int i;
      if (chiralCenters.size() > 0) {
	fourD = true;
	for (i = 0; i < nat; ++i) {
	  RDGeom::PointND *pt = new RDGeom::PointND(4);
	  positions.push_back(pt);
	}
      } else {
	for (i = 0; i < nat; ++i) {
	  RDGeom::Point3D *pt = new RDGeom::Point3D();
	  positions.push_back(pt);
	}
      }
      bool gotCoords = _embedPoints(positions, mmat, randNegEig, numZeroFail, optimizerForceTol,
				    basinThresh, seed, maxIterations, chiralCenters);
      
      int confId = -1;
      if (gotCoords) {
	Conformer *conf = new Conformer(nat);
	
	for (i = 0; i < nat; ++i) {
	  conf->setAtomPos(i, RDGeom::Point3D((*positions[i])[0], (*positions[i])[1], (*positions[i])[2]));
	  delete positions[i];
	}
	
        if (clearConfs) {
	  mol.clearConformers();
	  conf->setId(0);
          confId = (int)mol.addConformer(conf);
        } else {
          confId = (int)mol.addConformer(conf, true);
        }
      } 
      /*
      else {
        delete conf;
	}*/
      return confId;
    }

    
    void _fillAtomPositions(RDGeom::Point3DConstPtrVect &pts, const Conformer &conf) {
      unsigned int na = conf.getNumAtoms();
      pts.clear();
      unsigned int ai;
      pts.reserve(na);
      for (ai = 0; ai < na; ++ai) {
	pts.push_back(&conf.getAtomPos(ai));
      }
    }
        
    bool _isConfFarFromRest(const ROMol &mol, const Conformer &conf, double threshold) {
      // NOTE: it is tempting to use some triangle inequality to prune conformations here 
      // but some basic testing has shown very little advantage and given that the time for 
      // pruning fades in comparison to embedding - we will use a simple for loop below over all
      // conformation untill we find a match
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

    INT_VECT EmbedMultipleConfs(ROMol &mol, unsigned int numConfs, unsigned int maxIterations, 
                                int seed, bool clearConfs, 
                                bool randNegEig, unsigned int numZeroFail,
				double optimizerForceTol,double basinThresh, double pruneRmsThresh) {

      unsigned int nat = mol.getNumAtoms();
      DistGeom::BoundsMatrix *mat = new DistGeom::BoundsMatrix(nat);
      DistGeom::BoundsMatPtr mmat(mat);
      initBoundsMat(mmat);
      
      INT_VECT res;
      setTopolBounds(mol, mmat, true, false);
      
      if (!DistGeom::triangleSmoothBounds(mmat)) {
        // ok this bound matrix failed to triangle smooth - re-compute the bounds matrix 
        // with out 15 bound and with VDW scaling
        initBoundsMat(mmat);
        setTopolBounds(mol, mmat, false, true);
        // try triangle smoothing again - give up if we fail
        if (!DistGeom::triangleSmoothBounds(mmat)) {
          return res;
        }
      }
      
      if (clearConfs) {
        mol.clearConformers();
      }

      // find all the chiral centers in the molecule
      DistGeom::VECT_CHIRALSET chiralCenters;
      MolOps::assignAtomChiralCodes(mol);
      findChiralSets(mol, chiralCenters);

      // if we have any chiral centers we will first embed the molecule in four dimensions
      // otherwise we will use 3D 
      RDGeom::PointPtrVect positions;
      bool fourD = false;
      unsigned int i;
      if (chiralCenters.size() > 0) {
	fourD = true;
	for (i = 0; i < nat; ++i) {
	  RDGeom::PointND *pt = new RDGeom::PointND(4);
	  positions.push_back(pt);
	}
      } else {
	for (i = 0; i < nat; ++i) {
	  RDGeom::Point3D *pt = new RDGeom::Point3D();
	  positions.push_back(pt);
	}
      }

      unsigned int ci;
      //RDNumeric::DoubleSymmMatrix distMat(nat, 0.0);
      bool gotCoords;
      
      for (ci = 0; ci < numConfs; ci++) {
	
	gotCoords = _embedPoints(positions, mmat, randNegEig, numZeroFail, optimizerForceTol,
				 basinThresh, (ci+1)*seed, maxIterations, chiralCenters);
		
	
        if (gotCoords) {
	  Conformer *conf = new Conformer(nat);
	
	  for (i = 0; i < nat; ++i) {
	    conf->setAtomPos(i, RDGeom::Point3D((*positions[i])[0], (*positions[i])[1], (*positions[i])[2]));
	  }

	  bool addConf = true; // add the conformation to the molecule by default
	  if (pruneRmsThresh > 0.0) { // check if we are pruning away conformations
	    if (!_isConfFarFromRest(mol, *conf, pruneRmsThresh)) { // check if a closeby conformation has already been chosen
	      addConf = false;
	      delete conf;
	    }
	  }
	  if (addConf) {
	    int confId = (int)mol.addConformer(conf, true);
	    res.push_back(confId);
	  }
        } 
      }
      for (i = 0; i < nat; ++i) {
	delete positions[i];
      }
      return res;
    } 
  }
}
    
