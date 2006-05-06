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
#include <GraphMol/Conformer.h>
#include <RDGeneral/types.h>
#include <RDGeneral/RDLog.h>

#define ERROR_TOL 0.00001

namespace RDKit {
  namespace DGeomHelpers {
    
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
        
      bool gotCoords = false;

      // get pointers to atom position into an array so that they can be passed around
      Conformer *conf = new Conformer(nat);
      DistGeom::PointPtrVect positions;
      unsigned int i;
      for (i = 0; i < nat; i++) {
	positions.push_back(&conf->getAtomPos(i));
      }

      int confId = -1;
      if (clearConfs) {
        mol.clearConformers();
        conf->setId(0);
      }
      
      RDNumeric::DoubleSymmMatrix distMat(nat, 0.0);
      unsigned int iter = 0;
      while ((gotCoords == false) && (iter < maxIterations)) {
        if (seed > 0) {
          // we will seed the distance matrix picker by the the iteration ID
          // so that we will always get the same coordinate for a given topology
          pickRandomDistMat(*mmat, distMat, (iter+1)*seed); //(mol.getNumConformers()+1));
        } else {
          pickRandomDistMat(*mmat, distMat);
        }
	
        // now compute the coordinates
        gotCoords = DistGeom::computeInitialCoords(distMat, positions,
						   randNegEig, numZeroFail);
        iter += 1;
      }
      //BOOST_LOG(rdDebugLog) << iter << "\n";
      if (gotCoords) { 
        // if we managed to get an initial embedding, 
        // let us minimize the distance violation error function
        ForceFields::ForceField *field = DistGeom::constructForceField(*mmat, positions,
								       0,basinThresh);
        
        if (field) {
          field->initialize();
          if(field->calcEnergy() > ERROR_TOL){
            int needMore = 1;
            while(needMore){
              needMore = field->minimize(2,optimizerForceTol);
            }
          }
          delete field;
        }
      }
      
      if (gotCoords) {
        if (clearConfs) {
          confId = (int)mol.addConformer(conf);
        } else {
          confId = (int)mol.addConformer(conf, true);
        }
      } else {
        delete conf;
      }
      return confId;
    }


    INT_VECT EmbedMultipleConfs(ROMol &mol, unsigned int numConfs, unsigned int maxIterations, 
                                int seed, bool clearConfs, 
                                bool randNegEig, unsigned int numZeroFail,
				double optimizerForceTol,double basinThresh) {

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

      unsigned int ci;
      RDNumeric::DoubleSymmMatrix distMat(nat, 0.0);
      bool gotCoords;
      unsigned int numDPicks = 0;
      for (ci = 0; ci < numConfs; ci++) {
        Conformer *conf = new Conformer(nat);
        DistGeom::PointPtrVect positions;
        unsigned int i;
        for (i = 0; i < nat; i++) {
          positions.push_back(&conf->getAtomPos(i));
        }
        gotCoords = false;
        unsigned int iter = 0;

        while ((gotCoords == false) && (iter < maxIterations)) {
	  numDPicks++;
	  // update the distance matrix picker seed: by the the iteration ID
	  if(seed>0){
	    pickRandomDistMat(*mmat, distMat, seed*numDPicks);
	  } else {
	    pickRandomDistMat(*mmat, distMat);
	  }

          // now compute the coordinates
          gotCoords = DistGeom::computeInitialCoords(distMat, positions, randNegEig, numZeroFail);
          iter += 1;
        }

        if (gotCoords) {
          ForceFields::ForceField *field = DistGeom::constructForceField(*mmat, positions,0,basinThresh);
        
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
        if (gotCoords) {
          int confId = (int)mol.addConformer(conf, true);
          res.push_back(confId);
        } else {
          delete conf;
        }
      }
      return res;
    } 
  }
}
    
