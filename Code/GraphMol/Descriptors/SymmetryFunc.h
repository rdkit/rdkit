//
//  Copyright (c) 2020, Manan Goel
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

#ifndef SymmetryFuncRDKIT_H_JUNE2020
#define SymmetryFuncRDKIT_H_JUNE2020
#ifdef RDK_HAS_EIGEN3

#ifdef RDK_BUILD_DESCRIPTORS3D
#include <Eigen/Dense>
namespace RDKit {
class ROMol;
namespace Descriptors {
const std::string CoulombMatVersion = "1.0.0";
RDKIT_DESCRIPTORS_EXPORT Eigen::ArrayXXd cosine_cutoff(Eigen::ArrayXXd distances, double cutoff);
RDKIT_DESCRIPTORS_EXPORT Eigen::VectorXd generateSpeciesVector(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT Eigen::ArrayXd neighbor_pairs(Eigen::ArrayXXd coordinates, Eigen::VectorXd species, double cutoff);
RDKIT_DESCRIPTORS_EXPORT Eigen::ArrayXXd radial_terms(double cutoff, Eigen::ArrayXXd distances);
RDKIT_DESCRIPTORS_EXPORT Eigen::ArrayXXd angular_terms(double cutoff, Eigen::ArrayXXd vectors12);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> cumsum_from_zero(std::vector<int> count);
RDKIT_DESCRIPTORS_EXPORT std::pair<std::vector<int>, Eigen::ArrayXXd> triple_by_molecules(Eigen::ArrayXXd atom_index12_angular);
RDKIT_DESCRIPTORS_EXPORT Eigen::ArrayXXd triu_index(int num_species);
RDKIT_DESCRIPTORS_EXPORT Eigen::ArrayXXd SymmetryFunc(const ROMol &mol, int confId=-1);
RDKIT_DESCRIPTORS_EXPORT Eigen::ArrayXXd index_select(Eigen::ArrayXXd vector1, Eigen::ArrayXXd vector2, Eigen::ArrayXd index, int dim);
RDKIT_DESCRIPTORS_EXPORT Eigen::ArrayXXd index_add(Eigen::ArrayXXd vector1, Eigen::ArrayXXd vector2, Eigen::ArrayXXd index, int multi, int numAtoms);
}
}
#endif
#endif
#endif