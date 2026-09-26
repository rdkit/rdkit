#ifndef RDKIT_OSMORDRED_H
#define RDKIT_OSMORDRED_H

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>

#include <Eigen/Dense>  // we should try to remove those...

#if defined(_MSC_VER) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#endif

#include <lapacke.h>

namespace RDKit {
namespace Descriptors {
namespace Osmordred {
template <class T>
double InfoEntropy(const std::vector<T> &data) {
  T nInstances = 0;
  double accum = 0.0, d;

  for (const auto &val : data) {
    nInstances += val;
  }

  if (nInstances != 0) {
    for (const auto &val : data) {
      d = static_cast<double>(val) / nInstances;
      if (d != 0) {
        accum += -d * std::log(d);
      }
    }
  }
  return accum / std::log(2.0);
}

template <class T>
double WeightedInfoEntropy(const std::vector<T> &data,
                           const std::vector<double> &w) {
  T nInstances = 0;
  double accum = 0.0, d;

  for (const auto &val : data) {
    nInstances += val;
  }

  if (nInstances != 0) {
    for (size_t i = 0; i < data.size(); ++i) {
      const auto &val = data[i];
      const auto &wi = w[i];
      d = static_cast<double>(val) / nInstances;
      if (d != 0) {
        accum += -d * std::log(d) * wi;
      }
    }
  }
  return accum / std::log(2.0);
}

template <class T>
double WeightedCrossInfoEntropy(const std::vector<T> &data,
                                const std::vector<T> &w) {
  T nInstances = 0;
  double accum = 0.0, d;

  for (const auto &val : data) {
    nInstances += val;
  }

  if (nInstances != 0) {
    for (size_t i = 0; i < data.size(); ++i) {
      const auto &val = data[i];
      const auto &wi = w[i] * data[i];
      d = static_cast<double>(val) / nInstances;
      if (d != 0) {
        accum += -d * std::log(d) * wi;
      }
    }
  }
  return accum / std::log(2.0);
}

Eigen::MatrixXd calculateAdjacencyMatrix(const RDKit::ROMol &mol);
Eigen::MatrixXd calculateDistanceMatrix(const ROMol &mol);
std::vector<std::vector<double>> calculateDistanceMatrixL(
    const RDKit::ROMol &mol);
Eigen::MatrixXd calculateChargeTermMatrix(const Eigen::MatrixXd &A,
                                          const Eigen::MatrixXd &D);
Eigen::VectorXd calculateEccentricity(const ROMol &mol);

std::vector<double> calcValence(const RDKit::ROMol &mol);

std::vector<double> calcEStateIndices(const RDKit::ROMol &mol);
std::vector<double> calcIStateIndices(const RDKit::ROMol &mol);
std::vector<double> CalcHEStateIndices(const RDKit::ROMol &mol);

double getValenceElectrons(const Atom &atom);
double getSigmaElectrons(const Atom &atom);
double getIntrinsicState(const Atom &atom);
int GetPrincipalQuantumNumber(int atomicNum);

const std::map<int, double> &McGowanVolumAtomicMap();
const std::map<int, double> &Polarizability78AtomicMap();
const std::map<int, double> &Polarizability94AtomicMap();
const std::map<int, double> &VdWAtomicMap();
const std::map<int, double> &SandersonENAtomicMap();
const std::map<int, double> &PaulingENAtomicMap();
const std::map<int, double> &Allred_rocow_ENAtomicMap();
const std::map<int, double> &ionizationEnergyAtomicMap();

inline double vdw_volume(double r) {
  return (4.0 / 3.0) * M_PI * std::pow(r, 3);
}

// Enum for ChiType
enum class ChiType {
  Path = 1,
  Cluster,
  PathCluster,
  Chain
};

// Function to convert ChiType to string
inline std::string toString(ChiType type) {
  switch (type) {
    case ChiType::Path:
      return "Path";
    case ChiType::Cluster:
      return "Cluster";
    case ChiType::PathCluster:
      return "PathCluster";
    case ChiType::Chain:
      return "Chain";
    default:
      return "Unknown";
  }
}

ChiType classifySubgraph(const std::set<int> &degrees, bool isChain);
ChiType classifySubgraph(const RDKit::ROMol &mol,
                         const std::vector<int> &bondPath);

std::vector<std::tuple<std::vector<int>, std::set<int>, ChiType>>
extractAndClassifyPaths(const RDKit::ROMol &mol, unsigned int targetLength,
                        bool useHs);

void solveLinearSystem(const ROMol &mol, std::vector<double> &A,
                       std::vector<double> &B, int n, int nrhs, bool &success);

void compute_eigenvalues_and_eigenvectors(const Eigen::MatrixXd &matrix,
                                          Eigen::VectorXd &eigenvalues,
                                          Eigen::MatrixXd &eigenvectors);

double spAbs(const Eigen::VectorXd &eigenvalues);
double spMax(const Eigen::VectorXd &eigenvalues);
double spDiam(const Eigen::VectorXd &eigenvalues);
double spMean(const Eigen::VectorXd &eigenvalues);
double spAD(const Eigen::VectorXd &eigenvalues, double mean);
double logEE(const Eigen::VectorXd &eigenvalues);
double SM1(const Eigen::MatrixXd &matrix);

// Coefficient Sum of the Last Eigenvector (VE1)
double VE1(const Eigen::MatrixXd &matrix, Eigen::VectorXd &eigenvalues,
           Eigen::MatrixXd &eigenvectors);  // Average Coefficient of the Last
                                            // Eigenvector (VE2)
double VE2(const Eigen::MatrixXd &matrix, int numAtoms,
           Eigen::VectorXd &eigenvalues, Eigen::MatrixXd &eigenvectors);
// Logarithmic Coefficient Sum of the Last Eigenvector (VE3)
double VE3(const Eigen::MatrixXd &matrix, int numAtoms,
           Eigen::VectorXd &eigenvalues, Eigen::MatrixXd &eigenvectors);

// Randic-like Eigenvector-Based Index (VR1)
double VR1(const Eigen::MatrixXd &matrix,
           const std::vector<std::pair<int, int>> &bonds,
           Eigen::VectorXd &eigenvalues, Eigen::MatrixXd &eigenvectors);
// Normalized Randic-like Eigenvector-Based Index (VR2)
double VR2(const Eigen::MatrixXd &matrix,
           const std::vector<std::pair<int, int>> &bonds, int numAtoms,
           Eigen::VectorXd &eigenvalues, Eigen::MatrixXd &eigenvectors);
// Logarithmic Randic-like Eigenvector-Based Index (VR3)
double VR3(const Eigen::MatrixXd &matrix,
           const std::vector<std::pair<int, int>> &bonds, int numAtoms,
           Eigen::VectorXd &eigenvalues, Eigen::MatrixXd &eigenvectors);

void compute_eigenvalues_and_eigenvectorsL(
    std::vector<std::vector<double>> &matrix, std::vector<double> &eigenvalues,
    std::vector<std::vector<double>> &eigenvectors);
// Spectral Absolute Sum
double spAbsL(const std::vector<double> &eigenvalues);
// Leading Eigenvalue
double spMaxL(const std::vector<double> &eigenvalues);
// Spectral Diameter
double spDiamL(const std::vector<double> &eigenvalues);
// Mean of Eigenvalues
double spMeanL(const std::vector<double> &eigenvalues);
// Spectral Absolute Deviation
double spADL(const std::vector<double> &eigenvalues, double mean);
double logEEL(const std::vector<double> &eigenvalues);
double logEE_stable(const std::vector<double> &eigenvalues,
                    double threshold = 1e-10);
// Trace of Matrix
double SM1L(const std::vector<std::vector<double>> &matrix);
// Coefficient Sum of the Last Eigenvector
double VE1L(const std::vector<std::vector<double>> &eigenvectors);
// Average Coefficient of the Last Eigenvector
double VE2L(double ve1, int numAtoms);
// Logarithmic Coefficient Sum of the Last Eigenvector
double VE3L(double ve1, int numAtoms);
// Randic-like Eigenvector-Based Index
double VR1L(const std::vector<std::vector<double>> &eigenvectors,
            const std::vector<std::pair<int, int>> &bonds);
// Normalized Randic-like Eigenvector-Based Index
double VR2L(double vr1, int numAtoms);
// Logarithmic Randic-like Eigenvector-Based Index
double VR3L(double vr1, int numAtoms);

// Floyd Warshall shortest paths algorithms
Eigen::MatrixXd floydWarshall(Eigen::MatrixXd &A);
std::vector<std::vector<double>> floydWarshallL(
    std::vector<std::vector<double>> &matrix);

template<class MOL>
const RingInfo & getRings(const MOL &mol) {
  if (!mol.getRingInfo() || !mol.getRingInfo()->isSssrOrBetter()) {
    RDKit::MolOps::findSSSR(mol);
  }
  return *mol.getRingInfo();
}

}  // namespace Osmordred
}  // namespace Descriptors
}  // namespace RDKit

#endif
