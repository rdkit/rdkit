/* 
* $Id$
*
*  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
*  All rights reserved.
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met: 
*
*     * Redistributions of source code must retain the above copyright 
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above
*       copyright notice, this list of conditions and the following 
*       disclaimer in the documentation and/or other materials provided 
*       with the distribution.
*     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
*       nor the names of its contributors may be used to endorse or promote 
*       products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
%{
#include <DistGeom/ChiralSet.h>
#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/DistGeomUtils.h>
#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <DistGeom/TriangleSmooth.h>
#include <GraphMol/ROMol.h>
%}

// This prevents duplicate definitions in Java code (due to 2 C++ functions resolving to the same Java function)
%ignore RDKit::DGeomHelpers::initBoundsMat(DistGeom::BoundsMatPtr, double, double);
%ignore RDKit::DGeomHelpers::initBoundsMat(DistGeom::BoundsMatPtr,double);
%ignore RDKit::DGeomHelpers::initBoundsMat(DistGeom::BoundsMatPtr);

// This conflicts with definitions in UFF::, so ignore and make it a method on the bounds matrix
%ignore DistGeom::constructForceField;

%include <DistGeom/ChiralSet.h>
%include <DistGeom/BoundsMatrix.h>
%include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>
%include <DistGeom/DistGeomUtils.h>
%include <DistGeom/TriangleSmooth.h>

// A class to hang special distance geometry methods on.
%inline {
  namespace DistGeom {
    class DistanceGeom {
      public:

      /* From Embedder.h */
      static int EmbedMolecule(RDKit::ROMol &mol,
        unsigned int maxIterations=0, int seed=-1,
        bool clearConfs=true,
        bool useRandomCoords=false,double boxSizeMult=2.0,
        bool randNegEig=true,
        unsigned int numZeroFail=1,
        const std::map<int,RDGeom::Point3D> *coordMap=0,
        double optimizerForceTol=1e-3,
        double basinThresh=5.0) {
        
        return RDKit::DGeomHelpers::EmbedMolecule(mol, 
          maxIterations, seed, 
          clearConfs, 
          useRandomCoords, boxSizeMult,
          randNegEig, 
          numZeroFail, 
          coordMap, 
          optimizerForceTol, 
          basinThresh);
      }

      static RDKit::INT_VECT EmbedMultipleConfs(RDKit::ROMol &mol,
        unsigned int numConfs=10,
        unsigned int maxIterations=30, 
        int seed=-1, bool clearConfs=true, 
      	bool useRandomCoords=false,double boxSizeMult=2.0,
        bool randNegEig=true, unsigned int numZeroFail=1,
        double pruneRmsThresh=-1.0,
        const std::map<int,RDGeom::Point3D> *coordMap=0,
        double optimizerForceTol=1e-3,double basinThresh=5.0) {

        return RDKit::DGeomHelpers::EmbedMultipleConfs(mol,
          numConfs,
          maxIterations, 
          seed, clearConfs, 
          useRandomCoords,boxSizeMult,
          randNegEig, numZeroFail,
          pruneRmsThresh, 
          coordMap,
          optimizerForceTol, basinThresh);
      }    

      static void SetTopolBounds(RDKit::ROMol &mol,
        DistGeom::BoundsMatrix* mmat,
        bool set15bounds=true, bool scaleVDW=false) {

        DistGeom::BoundsMatPtr* mmatAsSmartPtr = new DistGeom::BoundsMatPtr(mmat);
        RDKit::DGeomHelpers::setTopolBounds(mol, *mmatAsSmartPtr, set15bounds, scaleVDW);
      }

      static ForceFields::ForceField *ConstructForceField(DistGeom::BoundsMatrix& mmat,
         RDGeom::PointPtrVect &positions, const VECT_CHIRALSET &csets,
         double weightChiral=1.0,
         double weightFourthDim=0.1,
         std::map< std::pair<int,int>,double> *extraWeights=0,
         double basinSizeTol=5.0) {

        return DistGeom::constructForceField(mmat, positions, csets,
          weightChiral,
          weightFourthDim,
          extraWeights,
          basinSizeTol);
      } 

      // This is hard to override in Java, so do it here.
      static bool ComputeInitialCoords(const RDNumeric::SymmMatrix<double> &distmat,  
        std::vector<RDGeom::Point3D *> &positions, bool randNegEig=false, 
        unsigned int numZeroFail=2) {
        
        return DistGeom::computeInitialCoords(distmat, (RDGeom::PointPtrVect &) positions, randNegEig, numZeroFail);
      }
    };
  }
}


