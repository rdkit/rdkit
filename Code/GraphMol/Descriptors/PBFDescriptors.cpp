//
//  Copyright (c) 2012, Institue of Cancer Research.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Institue of Cancer Research.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// For more information on the Plane of Best Fit please see http://pubs.acs.org/doi/abs/10.1021/ci300293f
//
//  If this code has been useful to you, please include the reference
//  in any work which has made use of it:

//  Plane of Best Fit: A Novel Method to Characterize the Three-Dimensionality of Molecules, Nicholas C. Firth, Nathan Brown, and Julian Blagg, Journal of Chemical Information and Modeling 2012 52 (10), 2516-2525

//
//
// Created by Nicholas Firth, November 2011
// Modified by Greg Landrum for inclusion in the RDKit distribution November 2012
// Modified by David Hahn for inclusion in the RDKit distribution August 2014
//

#include "PBFDescriptors.h"
#include <vector>

#ifdef RDK_USE_EIGEN3
#include <Eigen/Dense>
#else
#include <Numerics/EigenSolvers/JacobiEigenSolver.h>
#endif

using namespace RDKit;
using namespace RDNumeric;

double _distanceFromAPlane(const RDGeom::Point3D &pt,const std::vector<double> &plane, double denom){
  double numer=0.0;
  numer = std::abs(pt.x*plane[0]+pt.y*plane[1]+pt.z*plane[2]+plane[3]);

  return numer/denom;
}

void _getBestFitPlane(const std::vector<RDGeom::Point3D> &points,
                     std::vector<double> &plane,
                     const std::vector<double> *weights) {
  PRECONDITION((!weights || weights->size()>=points.size()),"bad weights vector");
  RDGeom::Point3D origin(0,0,0);
  double wSum=0.0;

  for(unsigned int i=0;i<points.size();++i){
    if(weights){
      double w=(*weights)[i];
      wSum+=w;
      origin+=points[i]*w;
    } else {
      wSum+=1;
      origin+=points[i];
    }
  }
  origin /= wSum;

  double sumXX=0,sumXY=0,sumXZ=0,sumYY=0,sumYZ=0,sumZZ=0;
  for(unsigned int i=0;i<points.size();++i){
    RDGeom::Point3D delta=points[i]-origin;
    if(weights){
      double w=(*weights)[i];
      delta *= w;
    }
    sumXX += delta.x*delta.x;
    sumXY += delta.x*delta.y;
    sumXZ += delta.x*delta.z;
    sumYY += delta.y*delta.y;
    sumYZ += delta.y*delta.z;
    sumZZ += delta.z*delta.z;
  }
  sumXX/=wSum;
  sumXY/=wSum;
  sumXZ/=wSum;
  sumYY/=wSum;
  sumYZ/=wSum;
  sumZZ/=wSum;

  RDGeom::Point3D normal;

#ifdef RDK_USE_EIGEN3
  Eigen::Matrix3d mat;
  mat << sumXX, sumXY, sumXZ,
         sumXY, sumYY, sumYZ,
         sumXZ, sumYZ, sumZZ;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(mat);
  if(eigensolver.info() != Eigen::Success){
    throw ValueErrorException("Eigen3: eigenvalue calculation did not converge");
  }

  normal.x=eigensolver.eigenvectors()(0,0);
  normal.y=eigensolver.eigenvectors()(1,0);
  normal.z=eigensolver.eigenvectors()(2,0);
#else
  RDNumeric::DoubleSymmMatrix mat(3);
  mat.setVal(0,0,sumXX);
  mat.setVal(1,0,sumXY);
  mat.setVal(2,0,sumXZ);
  mat.setVal(1,1,sumYY);
  mat.setVal(1,2,sumYZ);
  mat.setVal(2,2,sumZZ);

  DoubleVector eigenValues(3);
  DoubleSquareMatrix eigenVectors(3);

  unsigned int iter = EigenSolvers::jacobiEigenSolver(mat,eigenValues,eigenVectors, 100);
  if ( iter == 100 ) {
    throw ValueErrorException("JacobiEigenSolver: eigenvalue calculation did not converge");
  }

  normal.x=eigenVectors.getVal(0,0);
  normal.y=eigenVectors.getVal(0,1);
  normal.z=eigenVectors.getVal(0,2);
#endif
  plane[0] = normal.x;
  plane[1] = normal.y;
  plane[2] = normal.z;
  plane[3] = -1*normal.dotProduct(origin);
}

double Descriptors::calcPBF(const ROMol &mol,int confId){
  PRECONDITION(mol.getNumConformers()>=1,"molecule has no conformers")
  int numAtoms = mol.getNumAtoms();
  if(numAtoms<4) return 0;

  const Conformer &conf = mol.getConformer(confId);
  if(!conf.is3D()) return 0 ;

  std::vector<RDGeom::Point3D> points;
  points.reserve(numAtoms);
  for(unsigned int i=0; i<numAtoms; ++i){
    points.push_back(conf.getAtomPos(i));
  } 
    
  std::vector<double> plane(4);
  _getBestFitPlane(points,plane,0);

  double denom=0.0;
  for(unsigned int i=0; i<3; ++i){
    denom += plane[i]*plane[i];
  }
  denom = pow(denom,0.5);
    
  double res=0.0;
  for(unsigned int i=0; i<numAtoms; ++i){
    res += _distanceFromAPlane(points[i], plane, denom);
  }
  res /= numAtoms;

  return res;
}

