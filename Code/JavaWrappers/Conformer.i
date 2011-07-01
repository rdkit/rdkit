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
#include <GraphMol/Conformer.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <Geometry/point.h>
%}

%typemap(javacode) RDKit::Conformer %{
  public void setSwigCMemOwn(boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
  }
%}


%include <GraphMol/Conformer.h>

%extend RDKit::Conformer {
   /* From MolTransforms.h */

  RDGeom::Point3D computeCentroid(bool ignoreHs=true) {
	return MolTransforms::computeCentroid(*($self), ignoreHs);
  }


  RDNumeric::DoubleSymmMatrix *computeCovarianceMatrix(const RDGeom::Point3D &center,
                                                       bool normalize=false, 
                                                       bool ignoreHs=true) {
	return MolTransforms::computeCovarianceMatrix(*($self), center, normalize, ignoreHs);
  }

  RDGeom::Transform3D *computeCanonicalTransform(const RDGeom::Point3D *center=0,
                                                 bool normalizeCovar=false,
                                                 bool ignoreHs=true) {
	return MolTransforms::computeCanonicalTransform(*(self), center, normalizeCovar, ignoreHs);
  }

  void transformConformer(const RDGeom::Transform3D &trans) {
	MolTransforms::transformConformer(*($self), trans);
  }

  void canonicalizeConformer(const RDGeom::Point3D *center=0,
                             bool normalizeCovar=false, bool ignoreHs=true) {
	MolTransforms::canonicalizeConformer(*($self), center, normalizeCovar, ignoreHs);
  }

}
