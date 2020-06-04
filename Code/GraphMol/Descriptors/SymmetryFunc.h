
#ifndef SymmetryFuncRDKIT_H_JUNE2020
#define SymmetryFuncRDKIT_H_JUNE2020
#ifdef RDK_HAS_EIGEN3

#ifdef RDK_BUILD_DESCRIPTORS3D
#include <eigen3/Eigen/Dense>
namespace RDKit {
class ROMol;
namespace Descriptors {
const std::string CoulombMatVersion = "1.0.0";
RDKIT_DESCRIPTORS_EXPORT Eigen::ArrayXXd SymmetryFunc(const ROMol &mol, std::vector<std::vector<double>> &res, 
					 int confId=-1);
}
}
#endif
#endif
#endif