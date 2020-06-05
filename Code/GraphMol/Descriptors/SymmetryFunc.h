
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
}
}
#endif
#endif
#endif