#ifndef RDKIT_OSMORDRED_H
#define RDKIT_OSMORDRED_H

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>

#include <Eigen/Dense>  // we should try to remove those...

#include <GraphMol/Substruct/SubstructMatch.h>

#include <memory>
#include <optional>

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

//! Per-molecule intermediates shared by the descriptor blocks of one
//! calcOsmordred call ("compute once, use everywhere"). Not exported: the
//! public calcXxx(mol) functions build a local context, calcOsmordred builds
//! one for all blocks. Every intermediate is built lazily, on first use, from
//! the molecule as it is at that point, exactly as the blocks built it before.
class OsmordredContext {
 public:
  explicit OsmordredContext(const ROMol &mol) : d_mol(mol) {}
  OsmordredContext(const OsmordredContext &) = delete;
  OsmordredContext &operator=(const OsmordredContext &) = delete;

  const ROMol &mol() const { return d_mol; }

  //! MolOps::addHs(mol). The copy is shared only when the input already has
  //! SSSR-or-better ring information: blocks never modify that, so the H
  //! molecule (which copies the ring information) is the same whenever it is
  //! built. Otherwise a fresh copy is made on every call, as before. Blocks
  //! may set computed properties (Gasteiger charges, cached distance matrix)
  //! on it but must not modify its structure.
  const ROMol &molWithHs();

  //! unique matches of every extended EState atom-type query on mol (index =
  //! query index); filled by the EState code on first use
  std::unique_ptr<std::vector<std::vector<MatchVectType>>> estateExtMatches;

  //! Descriptors::calcNumHBA / calcNumHBD of mol, filled on first use
  std::optional<unsigned int> numHBA;
  std::optional<unsigned int> numHBD;

  //! calcPol(mol), also the basis of calcMR; filled on first use
  std::optional<double> pol;

 private:
  const ROMol &d_mol;
  std::unique_ptr<ROMol> d_molWithHs;
  bool d_molWithHsShared = false;
};

// Internal context overloads of the public descriptor blocks (same results
// as the public functions, which call them with a local context).
std::vector<int> calcAtomCounts(OsmordredContext &ctx);
std::vector<int> calcBondCounts(OsmordredContext &ctx);
std::vector<double> calcWeight(OsmordredContext &ctx);
double calcVdwVolumeABC(OsmordredContext &ctx);
std::vector<int> calcLipinskiGhose(OsmordredContext &ctx);
double calcMcGowanVolume(OsmordredContext &ctx);
std::vector<double> calcPolarizability(OsmordredContext &ctx);
std::vector<double> calcConstitutional(OsmordredContext &ctx);
std::vector<double> calcRNCG_RPCG(OsmordredContext &ctx);
std::vector<double> calcAutoCorrelation(OsmordredContext &ctx);
double calcFramework(OsmordredContext &ctx);
std::vector<double> calcEStateDescs(OsmordredContext &ctx, bool extended);
std::vector<double> calcHydrogenBond(OsmordredContext &ctx);
double calcPol(OsmordredContext &ctx);
double calcMR(OsmordredContext &ctx);
std::vector<double> calcBEStateDescs(OsmordredContext &ctx);
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

//! value stored for \c atomicNum in one of the atomic property maps above,
//! or 0 when the element is missing (what std::map::operator[] used to insert)
inline double atomicMapValue(const std::map<int, double> &atomicMap,
                             int atomicNum) {
  const auto it = atomicMap.find(atomicNum);
  return it == atomicMap.end() ? 0.0 : it->second;
}

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

//! (begin atom, end atom) of every bond, indexed by bond index
std::vector<std::pair<int, int>> getBondAtoms(const RDKit::ROMol &mol);

//! Classifies a connected bond subgraph (as returned by
//! findAllSubgraphsOfLengthN) as Chain/Path/PathCluster/Cluster.
//! \c bondAtoms comes from getBondAtoms(). On return \c atoms holds the
//! subgraph's atom indices in ascending order. \c degreeScratch must have at
//! least getNumAtoms() zero entries; it is left zeroed so it can be reused
//! across calls without reallocation.
ChiType classifyBondSubgraph(const std::vector<std::pair<int, int>> &bondAtoms,
                             const std::vector<int> &bondPath,
                             std::vector<int> &degreeScratch,
                             std::vector<int> &atoms);

std::vector<std::tuple<std::vector<int>, std::set<int>, ChiType>>
extractAndClassifyPaths(const RDKit::ROMol &mol, unsigned int targetLength,
                        bool useHs);

//! Conservative substructure screen: returns false only if some atom of
//! queryMol is compatible (Query::Match, the test the matcher applies) with
//! no atom of mol, in which case SubstructMatch(mol, queryMol) finds nothing.
//! Recursive SMARTS are screened through their own query molecules.
bool queryMolMayMatch(const ROMol &mol, const ROMol &queryMol);

//! Number of matches SubstructMatch(mol, queryMol, matches, true) returns
//! (unique matches, default parameters), skipping the matcher when
//! queryMolMayMatch() rules a match out and counting single-atom queries
//! without recursive SMARTS directly.
unsigned int countUniqueMatches(const ROMol &mol, const ROMol &queryMol);

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
//! in-place Floyd-Warshall on a dense row-major distance matrix
void floydWarshallL(std::vector<std::vector<double>> &matrix);

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
