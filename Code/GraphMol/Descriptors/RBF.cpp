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
// Further modified by Greg Landrum for inclusion in the RDKit core September 2016
// Adding RBF descriptors to 3D descriptors by Guillaume Godin

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

#include "RBF.h"
#include <Numerics/Matrix.h>
#include <Numerics/SquareMatrix.h>
#include <Numerics/SymmMatrix.h>
#include <boost/foreach.hpp>
#include <math.h>
#include <Eigen/Dense>

namespace RDKit {
namespace Descriptors{
namespace {


std::vector<double> getG(int n){
  std::vector<double> res = new std::vector<double>();
  for (int i=0;i<n;i++) {
    res[i] = 1+i*n/2;
  }  
  return res;
}



double getAtomDistance(const RDGeom::Point3D x1, const RDGeom::Point3D x2){
  double res=0;
  for (int i=0;i<3;i++) {
    res+=pow(x1[i]-x2[i],2);
  }
  return sqrt(res);
}



std::vector<std::vector<double>> GetGeometricalDistanceMatrix(const std::vector<RDGeom::Point3D> &points){
    int numAtoms= points.size();

    std::vector<std::vector<double>> res(numAtoms,std::vector<double>(numAtoms,0));
    for(unsigned int i=0; i<numAtoms; ++i){
        for(unsigned int j=i+1; j<numAtoms; ++j){
            res[i][j]=getAtomDistance(points[i], points[j]);
            res[j][i]=res[i][j];
    }

    return res;

}



std::vector<double> CalculateUnweightRDF(const Conformer &conf,const std::vector<RDGeom::Point3D> &points){
   std::vector<double>  R = getG(30);
   std::vector<double>  RDFres(std::vector<double>(numAtoms,0));

   int numAtoms = conf.getNumAtoms();

   std::vector<RDGeom::Point3D> points;
     points.reserve(numAtoms);
     for(unsigned int i=0; i<numAtoms; ++i){
          points.push_back(conf.getAtomPos(i));
     }

  std::vector<std::vector<double>> DM = GetGeometricalDistanceMatrix(&points);

  for (int i=0;i<30;i++) {
      double res=0;
      for (int j=0;j<numAtoms-1;j++)  {
        for (int k=j+1;k<numAtoms;k++)  {
          res+=exp(-100*pow(R[i]-DM[j][k],2));
        }
      }

      RDFres.push_back(res);
  }

  return RDFres;
}


} //end of anonymous namespace
} // end of Descriptors namespace
} // end of RDKit namespace
