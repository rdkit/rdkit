//
//  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
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

/*! \file USRDescriptor.h

  \brief Contains the USR descriptor. Use MolDescriptors.h in client code.

*/
#include <RDGeneral/export.h>
#ifndef __RD_USR_H__
#define __RD_USR_H__

#include <Geometry/point.h>
#include <Numerics/Vector.h>

namespace RDKit {
class ROMol;
class Conformer;
namespace Descriptors {
/*!
  Calculates the ultra-fast shape recognition (USR) descriptor

  Reference: P. J. Ballester, W. G. Richards, JCC (2007), 28, 1711 - 1723.

  Derived from RDKit Python implementation of Jan Domanski
  who derived his code from Adrian Schreyer's code:
  http://hg.adrianschreyer.eu/usrcat/src/70e075d93cd2?at=default

  \param mol          the molecule of interest
  \param descriptor   storage for the computed USR descriptor
  \param confId       the conformer Id

*/
RDKIT_DESCRIPTORS_EXPORT void USR(const ROMol &mol,
                                  std::vector<double> &descriptor,
                                  int confId = -1);

/*!
  Calculates the ultra-fast shape recognition with CREDO atom types (USRCAT)
  descriptor

  Reference: A. M. Schreyer, T. Blundell, J. Cheminf. (2012), 4, 27.

  Derived from Python implementation Adrian Schreyer:
  http://hg.adrianschreyer.eu/usrcat/src/70e075d93cd2?at=default

  \param mol          the molecule of interest
  \param descriptor   storage for the computed USR descriptor
  \param confId       the conformer Id

*/
RDKIT_DESCRIPTORS_EXPORT void USRCAT(
    const ROMol &mol, std::vector<double> &descriptor,
    std::vector<std::vector<unsigned int>> &atomIds, int confId = -1);
/*!
  Calculates the four distance distributions for the USR descriptor

  \param coords      the atom coordinates
  \param dist        storage for the four distance distributions
  \param points      storage for the four points

*/
RDKIT_DESCRIPTORS_EXPORT void calcUSRDistributions(
    const RDGeom::Point3DConstPtrVect &coords,
    std::vector<std::vector<double>> &dist,
    std::vector<RDGeom::Point3D> &points);

/*!
  Calculates the four distance distributions for the USR descriptor

  \param coords      the atom coordinates
  \param points      vector with the points
  \param dist   storage for the distance distributions

*/
RDKIT_DESCRIPTORS_EXPORT void calcUSRDistributionsFromPoints(
    const RDGeom::Point3DConstPtrVect &coords,
    const std::vector<RDGeom::Point3D> &points,
    std::vector<std::vector<double>> &dist);

/*!
  Calculates the USR descriptor from the four distance distributions

  \param dist        vector with the four distance distributions
  \param descriptor  storage for the computed USR descriptor

*/
RDKIT_DESCRIPTORS_EXPORT void calcUSRFromDistributions(
    const std::vector<std::vector<double>> &dist,
    std::vector<double> &descriptor);

/*!
  Calculates the score between two USRCAT descriptors with weights

  \param d1       descriptor 1
  \param d2       descriptor 2
  \param weights  the weights for each subset of moments

  \return the score
*/
RDKIT_DESCRIPTORS_EXPORT double calcUSRScore(
    const std::vector<double> &d1, const std::vector<double> &d2,
    const std::vector<double> &weights);

}  // end of namespace Descriptors
}  // end of namespace RDKit

#endif
