#include "OsmordredHelpers.h"
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

// Function to get the Adjacency Matrix
Eigen::MatrixXd calculateAdjacencyMatrix(const RDKit::ROMol &mol) {
  unsigned int nAtoms = mol.getNumAtoms();
  Eigen::MatrixXd adjMatrix(nAtoms, nAtoms);
  adjMatrix.setZero();

  // Populate the adjacency matrix using RDKit's getAdjacencyMatrix is this
  // faster than calling
  for (unsigned int i = 0; i < nAtoms; ++i) {
    for (unsigned int j = 0; j < nAtoms; ++j) {
      const RDKit::Bond *bond = mol.getBondBetweenAtoms(i, j);
      adjMatrix(i, j) = (bond != nullptr) ? 1.0 : 0.0;
    }
  }

  return adjMatrix;
}

std::vector<std::vector<double>> calculateDistanceMatrixL(
    const RDKit::ROMol &mol) {
  unsigned int nAtoms = mol.getNumAtoms();

  // Get the distance matrix using RDKit's MolOps::getDistanceMat
  double *distanceMat = MolOps::getDistanceMat(
      mol, false, false, false);  // No bond order, no weights, no hydrogens

  // Convert the raw pointer to a 2D vector
  std::vector<std::vector<double>> distMatrix(nAtoms,
                                              std::vector<double>(nAtoms, 0.0));
  for (unsigned int i = 0; i < nAtoms; ++i) {
    for (unsigned int j = 0; j < nAtoms; ++j) {
      distMatrix[i][j] = distanceMat[i * nAtoms + j];
    }
  }

  return distMatrix;
}

Eigen::MatrixXd calculateDistanceMatrix(const ROMol &mol) {
  unsigned int nAtoms = mol.getNumAtoms();

  // Get the distance matrix using RDKit's MolOps::getDistanceMat
  double *distanceMat = MolOps::getDistanceMat(
      mol, false, false, false);  // no bond order, no weights, no hydrogens

  // Convert the raw pointer to an Eigen MatrixXd (distance matrix)
  Eigen::MatrixXd distMatrix(nAtoms, nAtoms);

  for (unsigned int i = 0; i < nAtoms; ++i) {
    for (unsigned int j = i; j < nAtoms; ++j) {
      distMatrix(i, j) = distanceMat[i * nAtoms + j];
      if (j > i) {
        distMatrix(j, i) =
            distMatrix(i, j);  // diagonal is already set only the symetrical
                               // part need to be clone
      }
    }
  }

  return distMatrix;
}

Eigen::MatrixXd calculateChargeTermMatrix(const Eigen::MatrixXd &A,
                                          const Eigen::MatrixXd &D) {
  // Step 1: Invert non-zero elements of D^2 and zero out the diagonal
  Eigen::MatrixXd D2 = (D.array() != 0).select(D.array().pow(-2), 0.0);
  D2.diagonal().setZero();

  // Step 2: Perform matrix multiplication and antisymmetric operation
  Eigen::MatrixXd M = A * D2;
  return M - M.transpose();
}

// Function to calculate eccentricity (max distance for each atom)
Eigen::VectorXd calculateEccentricity(const ROMol &mol) {
  unsigned int nAtoms = mol.getNumAtoms();

  // Calculate the distance matrix using Eigen
  Eigen::MatrixXd distanceMatrix = calculateDistanceMatrix(mol);

  // For each atom, calculate the maximum distance to any other atom
  // (eccentricity)
  Eigen::VectorXd eccentricity(nAtoms);

  // Eigen's rowwise max operation
  eccentricity = distanceMatrix.colwise().maxCoeff();

  return eccentricity;
}

// ZagrebIndex this is correct! we can use the Eigen too...
std::vector<double> calcValence(const RDKit::ROMol &mol) {
  unsigned int nAtoms = mol.getNumAtoms();
  double *adjMat =
      RDKit::MolOps::getAdjacencyMatrix(mol, false, false, false, "noBO");

  // Sum along each column to get the valence of each atom
  std::vector<double> valences(nAtoms, 0);
  for (unsigned int i = 0; i < nAtoms; ++i) {
    valences[i] = static_cast<double>(
        std::accumulate(&adjMat[i * nAtoms], &adjMat[(i + 1) * nAtoms], 0.0));
  }

  return valences;
}

// Calculate EState indices for a molecule
std::vector<double> calcEStateIndices(const RDKit::ROMol &mol) {
  std::vector<double> Is = calcIStateIndices(mol);
  int numAtoms = mol.getNumAtoms();

  // Compute distance matrix
  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"

  std::vector<double> accum(numAtoms, 0.0);

  // Compute accumulative EState contributions
  for (int i = 0; i < numAtoms; ++i) {
    for (int j = i + 1; j < numAtoms; ++j) {
      double p =
          distances[i * numAtoms + j] + 1.;  // Fix : use "p = distance + 1."
      if (p < 1e6) {                         // Valid distance
        double tmp = (Is[i] - Is[j]) / (p * p);
        accum[i] += tmp;
        accum[j] -= tmp;
      }
    }
  }

  // Combine initial EState values and accumulative contributions
  std::vector<double> res(numAtoms, 0.0);
  for (int i = 0; i < numAtoms; ++i) {
    res[i] = accum[i] + Is[i];
    mol.getAtomWithIdx(i)->setProp("EState", res[i]);
  }
  return res;
}

std::vector<double> calcIStateIndices(const RDKit::ROMol &mol) {
  const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
  int numAtoms = mol.getNumAtoms();

  std::vector<double> Is(numAtoms, 0.0);

  // Compute initial EState values
  for (int i = 0; i < numAtoms; ++i) {
    const auto *atom = mol.getAtomWithIdx(i);
    int degree = atom->getDegree();
    if (degree > 0) {
      int atomicNum = atom->getAtomicNum();
      int dv = tbl->getNouterElecs(atomicNum) - atom->getTotalNumHs();
      int N = GetPrincipalQuantumNumber(atomicNum);
      Is[i] = (4.0 / (N * N) * dv + 1.0) / degree;
    }
    mol.getAtomWithIdx(i)->setProp("IState", Is[i]);
  }
  return Is;
}

namespace {
double getRKHE(const Atom &atom) {
  const PeriodicTable *tbl = PeriodicTable::getTable();

  int hi = atom.getTotalNumHs();

  int Z = atom.getAtomicNum();

  int Sigma = atom.getDegree();

  int SigmaV = tbl->getNouterElecs(Z) - hi;

  int N = GetPrincipalQuantumNumber(Z);

  return (SigmaV - Sigma) / static_cast<double>(N * N);
}

}  // namespace

std::vector<double> CalcHEStateIndices(const RDKit::ROMol &mol) {
  size_t nAtoms = mol.getNumAtoms();

  int numAtoms = mol.getNumAtoms();

  std::vector<double> RKHE(numAtoms, 0.0);
  std::vector<double> hasHs(numAtoms, 0.0);

  for (int i = 0; i < numAtoms; ++i) {
    const auto *atom = mol.getAtomWithIdx(i);
    int hi = atom->getTotalNumHs();
    if (hi > 0) {
      hasHs[i] = 1.;
      RKHE[i] = getRKHE(*atom);
    }
  }

  // Get distance matrix
  double *distMat = MolOps::getDistanceMat(
      mol, false, false, false);  // no bond order, no weights, no hydrogens
  ;
  std::vector<std::vector<double>> invDists(nAtoms,
                                            std::vector<double>(nAtoms, 0.0));

  for (size_t i = 0; i < nAtoms; ++i) {
    for (size_t j = 0; j < nAtoms; ++j) {
      double dist = distMat[i * nAtoms + j] + 1.;
      invDists[i][j] = dist > 0 ? 1.0 / (dist * dist) : 0.0;
    }
  }
  // Compute adjusted HEState
  std::vector<double> heStateIndices(nAtoms, 0.0);
  for (size_t i = 0; i < nAtoms; ++i) {
    double accum = 0.0;
    for (size_t j = 0; j < nAtoms; ++j) {
      if (i != j) {
        accum += RKHE[j] * invDists[i][j];
      }
    }
    heStateIndices[i] = hasHs[i] * (accum + 2.0 * RKHE[i] + 0.2);
  }

  return heStateIndices;
}

// caution use the Valence Method upper faster and cheaper!
double getValenceElectrons(const Atom &atom) {
  unsigned int N = atom.getAtomicNum();
  if (N == 1) {
    return 0.0;  // Hydrogen has 0 valence electrons in this context
  }

  const PeriodicTable *tbl = PeriodicTable::getTable();
  double Zv = tbl->getNouterElecs(N) - atom.getFormalCharge();
  double Z = atom.getAtomicNum() - atom.getFormalCharge();
  unsigned int hi = atom.getTotalNumHs();

  unsigned int he = std::count_if(
      atom.getOwningMol().atomNeighbors(&atom).begin(),
      atom.getOwningMol().atomNeighbors(&atom).end(),
      [](const Atom *neighbor) { return neighbor->getAtomicNum() == 1; });
  unsigned int h = hi + he;

  return (Zv - h) / (Z - Zv - 1.);
}

double getSigmaElectrons(const Atom &atom) {
  // Retrieve the molecule owning the atom
  const ROMol &mol = atom.getOwningMol();

  // Get the neighbors of the atom
  auto neighbors = mol.atomNeighbors(&atom);

  // Count the number of neighbors that are not hydrogen
  unsigned int sigmaElectrons = std::count_if(
      neighbors.begin(), neighbors.end(), [](const Atom *neighbor) {
        return neighbor->getAtomicNum() != 1;  // Exclude hydrogens
      });

  return static_cast<double>(sigmaElectrons);
}

// Function to get intrinsic state (equivalent to get_intrinsic_state)
double getIntrinsicState(const Atom &atom) {
  unsigned int i = atom.getAtomicNum();
  double d = getSigmaElectrons(atom);
  double dv = getValenceElectrons(atom);

  if (d == 0) {
    return 0.;
  }

  // Use GetPrincipalQuantumNumber to get the principal quantum number
  int n = GetPrincipalQuantumNumber(i);

  return ((2.0 / n) * (2.0 / n) * dv + 1) / d;
}

// Get principal quantum number for an atomic number
int GetPrincipalQuantumNumber(int atomicNum) {
  if (atomicNum <= 2) return 1;
  if (atomicNum <= 10) return 2;
  if (atomicNum <= 18) return 3;
  if (atomicNum <= 36) return 4;
  if (atomicNum <= 54) return 5;
  if (atomicNum <= 86) return 6;
  return 7;
}

const std::map<int, double> &McGowanVolumAtomicMap() {
  static const std::map<int, double> atomicProperties = {
      {1, 8.71},    {2, 6.75},    {3, 22.23},  {4, 20.27},  {5, 18.31},
      {6, 16.35},   {7, 14.39},   {8, 12.43},  {9, 10.47},  {10, 8.51},
      {11, 32.71},  {12, 30.75},  {13, 28.79}, {14, 26.83}, {15, 24.87},
      {16, 22.91},  {17, 20.95},  {18, 18.99}, {19, 51.89}, {20, 50.28},
      {21, 48.68},  {22, 47.07},  {23, 45.47}, {24, 43.86}, {25, 42.26},
      {26, 40.65},  {27, 39.05},  {28, 37.44}, {29, 35.84}, {30, 34.23},
      {31, 32.63},  {32, 31.02},  {33, 29.42}, {34, 27.81}, {35, 26.21},
      {36, 24.60},  {37, 60.22},  {38, 58.61}, {39, 57.01}, {40, 55.40},
      {41, 53.80},  {42, 52.19},  {43, 50.59}, {44, 48.98}, {45, 47.38},
      {46, 45.77},  {47, 44.17},  {48, 42.56}, {49, 40.96}, {50, 39.35},
      {51, 37.75},  {52, 36.14},  {53, 34.54}, {54, 32.93}, {55, 77.25},
      {56, 76.00},  {57, 74.75},  {58, 73.49}, {59, 72.24}, {60, 70.99},
      {61, 69.74},  {62, 68.49},  {63, 67.23}, {64, 65.98}, {65, 64.73},
      {66, 63.48},  {67, 62.23},  {68, 60.97}, {69, 59.72}, {70, 58.47},
      {71, 57.22},  {72, 55.97},  {73, 54.71}, {74, 53.46}, {75, 52.21},
      {76, 50.96},  {77, 49.71},  {78, 48.45}, {79, 47.20}, {80, 45.95},
      {81, 44.70},  {82, 43.45},  {83, 42.19}, {84, 40.94}, {85, 39.69},
      {86, 38.44},  {87, 75.59},  {88, 74.34}, {89, 73.09}, {90, 71.83},
      {91, 70.58},  {92, 69.33},  {93, 68.08}, {94, 66.83}, {95, 65.57},
      {96, 64.32},  {97, 63.07},  {98, 61.82}, {99, 60.57}, {100, 59.31},
      {101, 58.06}, {102, 56.81}, {103, 55.56}};

  return atomicProperties;
}

const std::map<int, double> &Polarizability78AtomicMap() {
  static std::map<int, double> atomicProperties = {
      {1, 0.666793}, {2, 0.204956}, {3, 24.3},    {4, 5.6},    {5, 3.03},
      {6, 1.76},     {7, 1.1},      {8, 0.802},   {9, 0.557},  {10, 0.3956},
      {11, 23.6},    {12, 10.6},    {13, 6.8},    {14, 5.38},  {15, 3.63},
      {16, 2.9},     {17, 2.18},    {18, 1.6411}, {19, 43.4},  {20, 22.8},
      {21, 17.8},    {22, 14.6},    {23, 12.4},   {24, 11.6},  {25, 9.4},
      {26, 8.4},     {27, 7.5},     {28, 6.8},    {29, 6.1},   {30, 7.1},
      {31, 8.12},    {32, 6.07},    {33, 4.31},   {34, 3.77},  {35, 3.05},
      {36, 2.4844},  {37, 47.3},    {38, 27.6},   {39, 22.7},  {40, 17.9},
      {41, 15.7},    {42, 12.8},    {43, 11.4},   {44, 9.6},   {45, 8.6},
      {46, 4.8},     {47, 7.2},     {48, 7.2},    {49, 10.2},  {50, 7.7},
      {51, 6.6},     {52, 5.5},     {53, 5.35},   {54, 4.044}, {55, 59.6},
      {56, 39.7},    {57, 31.1},    {58, 29.6},   {59, 28.2},  {60, 31.4},
      {61, 30.1},    {62, 28.8},    {63, 27.7},   {64, 23.5},  {65, 25.5},
      {66, 24.5},    {67, 23.6},    {68, 22.7},   {69, 21.8},  {70, 21.0},
      {71, 21.9},    {72, 16.2},    {73, 13.1},   {74, 11.1},  {75, 9.7},
      {76, 8.5},     {77, 7.6},     {78, 6.5},    {79, 5.8},   {80, 5.7},
      {81, 7.6},     {82, 6.8},     {83, 7.4},    {84, 6.8},   {85, 6.0},
      {86, 5.3},     {87, 48.7},    {88, 38.3},   {89, 32.1},  {90, 32.1},
      {91, 25.4},    {92, 27.4},    {93, 24.8},   {94, 24.5},  {95, 23.3},
      {96, 23.0},    {97, 22.7},    {98, 20.5},   {99, 19.7},  {100, 23.8},
      {101, 18.2},   {102, 17.5}};

  return atomicProperties;
}

const std::map<int, double> &Polarizability94AtomicMap() {
  static std::map<int, double> atomicProperties = {
      {1, 0.666793}, {2, 0.2050522}, {3, 24.33},   {4, 5.60},   {5, 3.03},
      {6, 1.67},     {7, 1.10},      {8, 0.802},   {9, 0.557},  {10, 0.39432},
      {11, 24.11},   {12, 10.6},     {13, 6.8},    {14, 5.53},  {15, 3.63},
      {16, 2.90},    {17, 2.18},     {18, 1.6411}, {19, 43.06}, {20, 22.8},
      {21, 17.8},    {22, 14.6},     {23, 12.4},   {24, 11.6},  {25, 9.4},
      {26, 8.4},     {27, 7.5},      {28, 6.8},    {29, 6.2},   {30, 5.75},
      {31, 8.12},    {32, 5.84},     {33, 4.31},   {34, 3.77},  {35, 3.05},
      {36, 2.4844},  {37, 47.24},    {38, 23.5},   {39, 22.7},  {40, 17.9},
      {41, 15.7},    {42, 12.8},     {43, 11.4},   {44, 9.6},   {45, 8.6},
      {46, 4.8},     {47, 6.78},     {48, 7.36},   {49, 10.2},  {50, 7.84},
      {51, 6.6},     {52, 5.5},      {53, 5.35},   {54, 4.044}, {55, 59.42},
      {56, 39.7},    {57, 31.1},     {58, 29.6},   {59, 28.2},  {60, 31.4},
      {61, 30.1},    {62, 28.8},     {63, 27.7},   {64, 23.5},  {65, 25.5},
      {66, 24.5},    {67, 23.6},     {68, 22.7},   {69, 21.8},  {70, 20.9},
      {71, 21.9},    {72, 16.2},     {73, 13.1},   {74, 11.1},  {75, 9.7},
      {76, 8.5},     {77, 7.6},      {78, 6.5},    {79, 5.8},   {80, 5.02},
      {81, 7.6},     {82, 7.01},     {83, 7.4},    {84, 6.8},   {85, 6.0},
      {86, 5.3},     {87, 48.6},     {88, 38.3},   {89, 32.1},  {90, 32.1},
      {91, 25.4},    {92, 24.9},     {93, 24.8},   {94, 24.5},  {95, 23.3},
      {96, 23.0},    {97, 22.7},     {98, 20.5},   {99, 19.7},  {100, 23.8},
      {101, 18.2},   {102, 16.4},    {112, 4.06},  {114, 4.59}, {119, 24.26}};

  return atomicProperties;
}
const std::map<int, double> &VdWAtomicMap() {
  static std::map<int, double> VdW_atomicRadii = {
      {1, 1.10},  {2, 1.40},  {3, 1.82},  {4, 1.53},   {5, 1.92},   {6, 1.70},
      {7, 1.55},  {8, 1.52},  {9, 1.47},  {10, 1.54},  {11, 2.27},  {12, 1.73},
      {13, 1.84}, {14, 2.10}, {15, 1.80}, {16, 1.80},  {17, 1.75},  {18, 1.88},
      {19, 2.75}, {20, 2.31}, {21, 2.15}, {22, 2.11},  {23, 2.07},  {24, 2.06},
      {25, 2.05}, {26, 2.04}, {27, 2.00}, {28, 1.97},  {29, 1.96},  {30, 2.01},
      {31, 1.87}, {32, 2.11}, {33, 1.85}, {34, 1.90},  {35, 1.85},  {36, 2.02},
      {37, 3.03}, {38, 2.49}, {39, 2.32}, {40, 2.23},  {41, 2.18},  {42, 2.17},
      {43, 2.16}, {44, 2.13}, {45, 2.10}, {46, 2.10},  {47, 2.11},  {48, 2.18},
      {49, 1.93}, {50, 2.17}, {51, 2.06}, {52, 2.06},  {53, 1.98},  {54, 2.16},
      {55, 3.43}, {56, 2.68}, {57, 2.43}, {58, 2.42},  {59, 2.40},  {60, 2.39},
      {61, 2.38}, {62, 2.36}, {63, 2.35}, {64, 2.34},  {65, 2.33},  {66, 2.31},
      {67, 2.30}, {68, 2.29}, {69, 2.27}, {70, 2.26},  {71, 2.24},  {72, 2.23},
      {73, 2.22}, {74, 2.18}, {75, 2.16}, {76, 2.16},  {77, 2.13},  {78, 2.13},
      {79, 2.14}, {80, 2.23}, {81, 1.96}, {82, 2.02},  {83, 2.07},  {84, 1.97},
      {85, 2.02}, {86, 2.20}, {87, 3.48}, {88, 2.83},  {89, 2.47},  {90, 2.45},
      {91, 2.43}, {92, 2.41}, {93, 2.39}, {94, 2.43},  {95, 2.44},  {96, 2.45},
      {97, 2.44}, {98, 2.45}, {99, 2.45}, {100, 2.45}, {101, 2.46}, {102, 2.46},
      {103, 2.46}};

  return VdW_atomicRadii;
}

const std::map<int, double> &SandersonENAtomicMap() {
  static std::map<int, double> SandersonElectronnegativityAtomicMap = {
      {1, 2.592},  {3, 0.670},  {4, 1.810},  {5, 2.275},  {6, 2.746},
      {7, 3.194},  {8, 3.654},  {9, 4.000},  {10, 4.5},   {11, 0.560},
      {12, 1.318}, {13, 1.714}, {14, 2.138}, {15, 2.515}, {16, 2.957},
      {17, 3.475}, {18, 3.31},  {19, 0.445}, {20, 0.946}, {21, 1.02},
      {22, 1.09},  {23, 1.39},  {24, 1.66},  {25, 2.2},   {26, 2.2},
      {27, 2.56},  {28, 1.94},  {29, 2.033}, {30, 2.223}, {31, 2.419},
      {32, 2.618}, {33, 2.816}, {34, 3.014}, {35, 3.219}, {36, 2.91},
      {37, 0.312}, {38, 0.721}, {39, 0.65},  {40, 0.9},   {41, 1.42},
      {42, 1.15},  {47, 1.826}, {48, 1.978}, {49, 2.138}, {50, 2.298},
      {51, 2.458}, {52, 2.618}, {53, 2.778}, {54, 2.34},  {55, 0.220},
      {56, 0.651}, {74, 0.98},  {80, 2.195}, {81, 2.246}, {82, 2.291},
      {83, 2.342}};

  return SandersonElectronnegativityAtomicMap;
}

const std::map<int, double> &PaulingENAtomicMap() {
  static std::map<int, double> PaulingelectronegativityAtomicMap = {
      {1, 2.2},   {3, 0.98},  {4, 1.57},  {5, 2.04},  {6, 2.55},  {7, 3.04},
      {8, 3.44},  {9, 3.98},  {11, 0.93}, {12, 1.31}, {13, 1.61}, {14, 1.9},
      {15, 2.19}, {16, 2.58}, {17, 3.16}, {19, 0.82}, {20, 1.0},  {21, 1.36},
      {22, 1.54}, {23, 1.63}, {24, 1.66}, {25, 1.55}, {26, 1.83}, {27, 1.88},
      {28, 1.91}, {29, 1.9},  {30, 1.65}, {31, 1.81}, {32, 2.01}, {33, 2.18},
      {34, 2.55}, {35, 2.96}, {36, 3.0},  {37, 0.82}, {38, 0.95}, {39, 1.22},
      {40, 1.33}, {41, 1.6},  {42, 2.16}, {43, 1.9},  {44, 2.2},  {45, 2.28},
      {46, 2.2},  {47, 1.93}, {48, 1.69}, {49, 1.78}, {50, 1.96}, {51, 2.05},
      {52, 2.1},  {53, 2.66}, {54, 2.6},  {55, 0.79}, {56, 0.89}, {57, 1.1},
      {58, 1.12}, {59, 1.13}, {60, 1.14}, {62, 1.17}, {64, 1.2},  {66, 1.22},
      {67, 1.23}, {68, 1.24}, {69, 1.25}, {71, 1.27}, {72, 1.3},  {73, 1.5},
      {74, 2.36}, {75, 1.9},  {76, 2.2},  {77, 2.2},  {78, 2.28}, {79, 2.54},
      {80, 2.0},  {81, 1.62}, {82, 2.33}, {83, 2.02}, {84, 2.0},  {85, 2.2},
      {87, 0.7},  {88, 0.9},  {89, 1.1},  {90, 1.3},  {91, 1.5},  {92, 1.38},
      {93, 1.36}, {94, 1.28}, {95, 1.3},  {96, 1.3},  {97, 1.3},  {98, 1.3},
      {99, 1.3},  {100, 1.3}, {101, 1.3}, {102, 1.3}};

  return PaulingelectronegativityAtomicMap;
}

const std::map<int, double> &Allred_rocow_ENAtomicMap() {
  static std::map<int, double> allred_rocow_electron_negativityAtomicMap = {
      {1, 2.20},  {3, 0.97},  {4, 1.47},  {5, 2.01},  {6, 2.50},  {7, 3.07},
      {8, 3.50},  {9, 4.10},  {11, 1.01}, {12, 1.23}, {13, 1.47}, {14, 1.74},
      {15, 2.06}, {16, 2.44}, {17, 2.83}, {19, 0.91}, {20, 1.04}, {21, 1.20},
      {22, 1.32}, {23, 1.45}, {24, 1.56}, {25, 1.60}, {26, 1.64}, {27, 1.70},
      {28, 1.75}, {29, 1.75}, {30, 1.66}, {31, 1.82}, {32, 2.02}, {33, 2.20},
      {34, 2.48}, {35, 2.74}, {37, 0.89}, {38, 0.99}, {39, 1.11}, {40, 1.22},
      {41, 1.23}, {42, 1.30}, {43, 1.36}, {44, 1.42}, {45, 1.45}, {46, 1.35},
      {47, 1.42}, {48, 1.46}, {49, 1.49}, {50, 1.72}, {51, 1.82}, {52, 2.01},
      {53, 2.21}, {55, 0.86}, {56, 0.97}, {57, 1.08}, {72, 1.23}, {73, 1.33},
      {74, 1.40}, {75, 1.46}, {76, 1.52}, {77, 1.55}, {78, 1.44}, {79, 1.42},
      {80, 1.44}, {81, 1.44}, {82, 1.55}, {83, 1.67}, {84, 1.76}, {85, 1.90}};

  return allred_rocow_electron_negativityAtomicMap;
}

const std::map<int, double> &ionizationEnergyAtomicMap() {
  static std::map<int, double> ionizationEnergyAtomicMap = {
      {1, 13.598443}, {2, 24.587387},  {3, 5.391719},   {4, 9.32270},
      {5, 8.29802},   {6, 11.26030},   {7, 14.5341},    {8, 13.61805},
      {9, 17.4228},   {10, 21.56454},  {11, 5.139076},  {12, 7.646235},
      {13, 5.985768}, {14, 8.15168},   {15, 10.48669},  {16, 10.36001},
      {17, 12.96763}, {18, 15.759610}, {19, 4.3406633}, {20, 6.11316},
      {21, 6.56149},  {22, 6.82812},   {23, 6.74619},   {24, 6.76651},
      {25, 7.43402},  {26, 7.9024},    {27, 7.88101},   {28, 7.6398},
      {29, 7.72638},  {30, 9.394199},  {31, 5.999301},  {32, 7.89943},
      {33, 9.7886},   {34, 9.75239},   {35, 11.8138},   {36, 13.99961},
      {37, 4.177128}, {38, 5.69485},   {39, 6.2173},    {40, 6.63390},
      {41, 6.75885},  {42, 7.09243},   {43, 7.28},      {44, 7.36050},
      {45, 7.45890},  {46, 8.3369},    {47, 7.57623},   {48, 8.99382},
      {49, 5.78636},  {50, 7.34392},   {51, 8.60839},   {52, 9.0096},
      {53, 10.45126}, {54, 12.12984},  {55, 3.893905},  {56, 5.211664},
      {57, 5.5769},   {58, 5.5387},    {59, 5.473},     {60, 5.5250},
      {61, 5.582},    {62, 5.6437},    {63, 5.67038},   {64, 6.14980},
      {65, 5.8638},   {66, 5.9389},    {67, 6.0215},    {68, 6.1077},
      {69, 6.18431},  {70, 6.25416},   {71, 5.42586},   {72, 6.82507},
      {73, 7.54957},  {74, 7.86403},   {75, 7.83352},   {76, 8.43823},
      {77, 8.96702},  {78, 8.9588},    {79, 9.22553},   {80, 10.4375},
      {81, 6.108194}, {82, 7.41663},   {83, 7.2855},    {84, 8.414},
      {86, 10.7485},  {87, 4.072741},  {88, 5.278423},  {89, 5.17},
      {90, 6.3067},   {91, 5.89},      {92, 6.1941},    {93, 6.2657},
      {94, 6.0260},   {95, 5.9738},    {96, 5.9914},    {97, 6.1979},
      {98, 6.2817},   {99, 6.42},      {100, 6.50},     {101, 6.58},
      {102, 6.65},    {103, 4.9},      {104, 6.0}};

  return ionizationEnergyAtomicMap;
}

namespace {
void performDFS(const RDKit::ROMol &mol, int startAtomIdx,
                const std::vector<int> &path, std::set<int> &visitedNodes,
                std::set<std::pair<int, int>> &visitedEdges,
                std::set<int> &degrees, bool &isChain) {
  std::set<int> pathBonds(path.begin(),
                          path.end());  // Bonds in the path for quick lookup
  std::unordered_map<int, std::set<int>>
      neighbors;  // Neighbors in the subgraph

  // Populate neighbors for the subgraph
  for (int bondIdx : path) {
    const auto *bond = mol.getBondWithIdx(bondIdx);
    int begin = bond->getBeginAtomIdx();
    int end = bond->getEndAtomIdx();
    neighbors[begin].insert(end);
    neighbors[end].insert(begin);
  }

  // Perform DFS
  std::stack<int> stack;
  std::unordered_map<int, int> parent;  // To track parent nodes in DFS
  stack.push(startAtomIdx);
  parent[startAtomIdx] = -1;  // Root node has no parent

  while (!stack.empty()) {
    int node = stack.top();
    stack.pop();

    if (visitedNodes.count(node)) {
      continue;
    }

    visitedNodes.insert(node);

    // Calculate degree for this node based on subgraph neighbors
    int degree = neighbors[node].size();
    degrees.insert(degree);  // Add degree to the set

    // Traverse neighbors in the subgraph
    for (int neighbor : neighbors[node]) {
      std::pair<int, int> edge = std::minmax(node, neighbor);

      if (!visitedNodes.count(neighbor)) {
        stack.push(neighbor);
        parent[neighbor] = node;  // Set parent for the neighbor
        visitedEdges.insert(edge);
      } else if (parent[node] != neighbor) {  // Detect back edge
        isChain = true;                       // Cycle detected
      }
    }
  }
}

bool allDegreesAreOneOrTwo(const std::set<int> &degrees) {
  return std::all_of(degrees.begin(), degrees.end(),
                     [](int d) { return d == 1 || d == 2; });
}

}  // namespace
ChiType classifySubgraph(const std::set<int> &degrees, bool isChain) {
  // Classification logic based on degrees and cycle detection
  if (isChain) {
    return ChiType::Chain;
  } else if (allDegreesAreOneOrTwo(degrees)) {
    return ChiType::Path;
  } else if (degrees.count(2)) {
    return ChiType::PathCluster;
  } else {
    return ChiType::Cluster;
  }
}

bool detectCycle(const std::vector<int> &atomPath) {
  std::unordered_set<int> visited;
  for (int atomIdx : atomPath) {
    if (visited.count(atomIdx)) {
      return true;  // Cycle detected
    }
    visited.insert(atomIdx);
  }
  return false;  // No cycle
}

ChiType classifySubgraph(const RDKit::ROMol &mol,
                         const std::vector<int> &bondPath) {
  // Map to store atom index mapping between subgraph and parent molecule
  INT_MAP_INT atomIdxMap;

  // Create a submolecule for the given bond path
  std::unique_ptr<RDKit::ROMol> subMol(
      Subgraphs::pathToSubmol(mol, bondPath, false, atomIdxMap));
  // Collect degrees in the subgraph
  std::set<int> degrees;
  for (const auto &atom : subMol->atoms()) {
    degrees.insert(atom->getDegree());
  }

  // Determine if the subgraph is a chain
  bool isChain = detectCycle(bondPath);

  // Apply classification logic
  if (isChain) {
    return ChiType::Chain;
  } else if (allDegreesAreOneOrTwo(degrees)) {
    return ChiType::Path;
  } else if (degrees.count(2)) {
    return ChiType::PathCluster;
  } else {
    return ChiType::Cluster;
  }
}

// Main function to extract and classify subgraphs
std::vector<std::tuple<std::vector<int>, std::set<int>, ChiType>>
extractAndClassifyPaths(const RDKit::ROMol &mol, unsigned int targetLength,
                        bool useHs) {
  std::vector<std::tuple<std::vector<int>, std::set<int>, ChiType>> results;

  // Get all subgraphs of the given length
  auto paths = findAllSubgraphsOfLengthN(
      mol, targetLength, useHs,
      -1);  // return both path and complex subgraphs AllPath return only true
            // Path ...! maybe we can leverage that except if it is too
            // expensive...

  for (const auto &path : paths) {
    // Prepare sets for DFS traversal
    std::set<int> visitedNodes;
    std::set<std::pair<int, int>> visitedEdges;
    std::set<int> degrees;
    bool isChain = false;

    // Start DFS from the first bond in the path
    if (!path.empty()) {
      int startAtomIdx = mol.getBondWithIdx(path.front())->getBeginAtomIdx();
      performDFS(mol, startAtomIdx, path, visitedNodes, visitedEdges, degrees,
                 isChain);
    }

    // If a cycle is detected, it's a Chain Path by definitin of the isChain
    // bool flag from DFS code this is a decision tree: Chain first than only 1
    // and 2 => Path than has 2 => Path Cluster else Cluster!
    ChiType type;
    if (isChain) {
      type = ChiType::Chain;
    } else if (allDegreesAreOneOrTwo(degrees)) {
      type = ChiType::Path;
    } else if (degrees.count(2)) {
      type = ChiType::PathCluster;
    } else {
      type = ChiType::Cluster;
    }
    results.emplace_back(path, visitedNodes, type);
  }
  return results;
}

void solveLinearSystem(const ROMol &mol, std::vector<double>& A, std::vector<double>& B,
		       int n, int nrhs, bool& success) {
    int lda = n; // Leading dimension of A
    int ldb = n; // Leading dimension of B
    int info;

    success = false; // Initialize success flag

    // CRITICAL FIX v2.0: Save original RHS before any LAPACK calls modify B
    // LAPACK routines modify B in-place, even when they fail!
    std::vector<double> B_original = B;

    // First, try dposv (Cholesky factorization for positive definite matrices)
    std::vector<double> A_copy = A; // Copy A because LAPACK modifies it
    info = LAPACKE_dposv(LAPACK_COL_MAJOR, 'U', n, nrhs, A_copy.data(), lda, B.data(), ldb);

    if (info == 0) {
        success = true;
        return;
    } else {
        // dposv failed; fall back to dgesv (LU factorization)
        // CRITICAL FIX v2.0: Restore original RHS before calling dgesv
        // dposv modified B even though it failed!
        B = B_original;
        
        std::vector<int> ipiv(n); // Pivot array for dgesv
        A_copy = A; // Reset A because it was modified by dposv
        info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, A_copy.data(), lda, ipiv.data(), B.data(), ldb);

        if (info == 0) {
            success = true;
            return;
        } else {
            // dgesv failed (singular matrix); fall back to dgelss (pseudo-inverse via SVD)
            // CRITICAL FIX v2.0: Added dgelss fallback for singular matrices
            // This provides a minimum-norm least-squares solution when exact solution doesn't exist
            B = B_original; // Restore original RHS values
            
            std::vector<double> A_copy2 = A; // Fresh copy for dgelss
            std::vector<double> B_copy = B; // Copy B because dgelss modifies it
            
            // Allocate workspace for dgelss
            std::vector<double> s(n); // Singular values
            int rank; // Rank of matrix
            double rcond = 1e-15; // Condition number threshold
            
            // dgelss computes least-squares solution: min ||Ax - b||_2
            info = LAPACKE_dgelss(LAPACK_COL_MAJOR, n, n, nrhs,
                                   A_copy2.data(), lda, B_copy.data(), ldb,
                                   s.data(), rcond, &rank);
            
            if (info == 0) {
                // dgelss succeeded - copy solution back to B
                B = B_copy;
                success = true;
                return;
            } else {
                // All solvers failed - this is a true error
                std::string outputSmiles = RDKit::MolToSmiles(mol);
                std::cerr << "ERROR: All LAPACK solvers failed (dposv, dgesv, dgelss): info=" 
                          << info << ", Smiles:" << outputSmiles << "\n";
            }
        }
    }
}

////// Barysz Matrixes Eigen style

// Function to compute eigenvalues and eigenvectors
void compute_eigenvalues_and_eigenvectors(const Eigen::MatrixXd &matrix,
                                          Eigen::VectorXd &eigenvalues,
                                          Eigen::MatrixXd &eigenvectors) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(
      matrix);  // hermitian (because we have symetrical matrix! critial)
  eigenvalues = solver.eigenvalues().real();  // Take real parts if complex
  eigenvectors = solver.eigenvectors().real();

  // Canonicalize eigenvector signs: force first non-zero entry positive per
  // column
  const double canonical_eps = 1e-12;
  for (int j = 0; j < eigenvectors.cols(); ++j) {
    for (int i = 0; i < eigenvectors.rows(); ++i) {
      double v = eigenvectors(i, j);
      if (std::abs(v) > canonical_eps) {
        if (v < 0.0) {
          eigenvectors.col(j) *= -1.0;
        }
        break;
      }
    }
  }
}

// Spectral Absolute Sum: sum of absolute eigenvalues
double spAbs(const Eigen::VectorXd &eigenvalues) {
  return eigenvalues.array().abs().sum();
}

// Leading Eigenvalue: largest eigenvalue
double spMax(const Eigen::VectorXd &eigenvalues) {
  return eigenvalues.maxCoeff();
}

// Spectral Diameter: difference between largest and smallest eigenvalues
double spDiam(const Eigen::VectorXd &eigenvalues) {
  return eigenvalues.maxCoeff() - eigenvalues.minCoeff();
}

// Mean of Eigenvalues: average eigenvalue
double spMean(const Eigen::VectorXd &eigenvalues) { return eigenvalues.mean(); }

// Spectral Absolute Deviation: sum of absolute deviations from the mean
double spAD(const Eigen::VectorXd &eigenvalues, double mean) {
  return (eigenvalues.array() - mean).abs().sum();
}

// Estrada-like index: log sum of exponentials of eigenvalues
double logEE(const Eigen::VectorXd &eigenvalues) {
  double a = std::max(eigenvalues.maxCoeff(), 0.0);
  double sx = (eigenvalues.array() - a).exp().sum() + std::exp(-a);
  return a + std::log(sx);
}

double SM1(const Eigen::MatrixXd &matrix) { return matrix.trace(); }

// Coefficient Sum of the Last Eigenvector (VE1)
double VE1(const Eigen::MatrixXd &matrix, Eigen::VectorXd &eigenvalues,
           Eigen::MatrixXd &eigenvectors) {
  compute_eigenvalues_and_eigenvectors(
      matrix, eigenvalues, eigenvectors);  // Reuse the eigenvalue computation
  Eigen::VectorXd eigenvector =
      eigenvectors.col(eigenvalues.size() - 1);  // Last eigenvector
  return eigenvector.array().abs().sum();        // Sum of the absolute values
}

// Average Coefficient of the Last Eigenvector (VE2)
double VE2(const Eigen::MatrixXd &matrix, int numAtoms,
           Eigen::VectorXd &eigenvalues, Eigen::MatrixXd &eigenvectors) {
  double VE1_value = VE1(matrix, eigenvalues, eigenvectors);
  return VE1_value / numAtoms;
}

// Logarithmic Coefficient Sum of the Last Eigenvector (VE3)
double VE3(const Eigen::MatrixXd &matrix, int numAtoms,
           Eigen::VectorXd &eigenvalues, Eigen::MatrixXd &eigenvectors) {
  double VE1_value = VE1(matrix, eigenvalues, eigenvectors);
  return std::log(0.1 * numAtoms * VE1_value);
}

// Randic-like Eigenvector-Based Index (VR1)
double VR1(const Eigen::MatrixXd &matrix,
           const std::vector<std::pair<int, int>> &bonds,
           Eigen::VectorXd &eigenvalues, Eigen::MatrixXd &eigenvectors) {
  compute_eigenvalues_and_eigenvectors(
      matrix, eigenvalues, eigenvectors);  // Reuse the eigenvalue computation
  Eigen::VectorXd eigenvector =
      eigenvectors.col(eigenvalues.size() - 1);  // Last eigenvector
  double result = 0.0;

  for (const auto &bond : bonds) {
    int i = bond.first;
    int j = bond.second;
    result +=
        std::pow(std::abs(eigenvector(i)) * std::abs(eigenvector(j)), -0.5);
  }

  return result;
}

// Normalized Randic-like Eigenvector-Based Index (VR2)
double VR2(const Eigen::MatrixXd &matrix,
           const std::vector<std::pair<int, int>> &bonds, int numAtoms,
           Eigen::VectorXd &eigenvalues, Eigen::MatrixXd &eigenvectors) {
  double VR1_value = VR1(matrix, bonds, eigenvalues, eigenvectors);
  return VR1_value / numAtoms;
}

// Logarithmic Randic-like Eigenvector-Based Index (VR3)
double VR3(const Eigen::MatrixXd &matrix,
           const std::vector<std::pair<int, int>> &bonds, int numAtoms,
           Eigen::VectorXd &eigenvalues, Eigen::MatrixXd &eigenvectors) {
  double VR1_value = VR1(matrix, bonds, eigenvalues, eigenvectors);
  return std::log(0.1 * numAtoms * VR1_value);
}

void compute_eigenvalues_and_eigenvectorsL(
    std::vector<std::vector<double>> &matrix, std::vector<double> &eigenvalues,
    std::vector<std::vector<double>> &eigenvectors) {
  int n = matrix.size();
  eigenvalues.resize(n);
  eigenvectors = matrix;  // Copy matrix to preserve the original

  // Convert the 2D vector to a 1D array in column-major order for LAPACK
  std::vector<double> flatMatrix(n * n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) flatMatrix[j * n + i] = eigenvectors[i][j];

  // Call LAPACKE_dsyev to compute eigenvalues and eigenvectors
  int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', n, flatMatrix.data(), n,
                           eigenvalues.data());
  if (info != 0) {
    throw std::runtime_error("Error in LAPACKE_dsyev: " + std::to_string(info));
  }

  // Reshape the flatMatrix back into eigenvectors
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) eigenvectors[i][j] = flatMatrix[j * n + i];

  // Canonicalize eigenvector signs: force first non-zero entry positive per
  // column
  const double canonical_eps = 1e-12;
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      double v = eigenvectors[i][j];
      if (std::abs(v) > canonical_eps) {
        if (v < 0.0) {
          for (int k = 0; k < n; ++k) {
            eigenvectors[k][j] = -eigenvectors[k][j];
          }
        }
        break;
      }
    }
  }
}

// Spectral Absolute Sum
double spAbsL(const std::vector<double> &eigenvalues) {
  return std::accumulate(
      eigenvalues.begin(), eigenvalues.end(), 0.0,
      [](double acc, double val) { return acc + std::abs(val); });
}

// Leading Eigenvalue
double spMaxL(const std::vector<double> &eigenvalues) {
  return *std::max_element(eigenvalues.begin(), eigenvalues.end());
}

// Spectral Diameter
double spDiamL(const std::vector<double> &eigenvalues) {
  auto [minIt, maxIt] =
      std::minmax_element(eigenvalues.begin(), eigenvalues.end());
  return *maxIt - *minIt;
}

// Mean of Eigenvalues
double spMeanL(const std::vector<double> &eigenvalues) {
  return std::accumulate(eigenvalues.begin(), eigenvalues.end(), 0.0) /
         eigenvalues.size();
}

// Spectral Absolute Deviation
double spADL(const std::vector<double> &eigenvalues, double mean) {
  return std::accumulate(
      eigenvalues.begin(), eigenvalues.end(), 0.0,
      [mean](double acc, double val) { return acc + std::abs(val - mean); });
}

double logEEL(const std::vector<double> &eigenvalues) {
  if (eigenvalues.empty()) {
    throw std::runtime_error("Eigenvalues vector is empty");
  }

  double maxVal = *std::max_element(eigenvalues.begin(), eigenvalues.end());
  double sumExp = 0.0;

  for (double val : eigenvalues) {
    sumExp += std::exp(val - maxVal);
  }

  // Avoid issues with log(0) if sumExp is very small
  if (sumExp < 1e-10) {
    return maxVal;  // If sumExp is tiny, Log_EE reduces to maxVal
  }

  return maxVal + std::log(sumExp);
}

double logEE_stable(const std::vector<double> &eigenvalues, double threshold) {
  // Handle small eigenvalues close to zero
  std::vector<double> adjustedEigenvalues = eigenvalues;
  for (auto &val : adjustedEigenvalues) {
    if (std::abs(val) < threshold) {
      val = 0.0;
    }
  }

  double maxVal =
      *std::max_element(adjustedEigenvalues.begin(), adjustedEigenvalues.end());
  double sumExp =
      std::accumulate(adjustedEigenvalues.begin(), adjustedEigenvalues.end(),
                      0.0, [maxVal](double acc, double val) {
                        return acc + std::exp(val - maxVal);
                      });
  return maxVal + std::log(sumExp);
}

// Trace of Matrix
double SM1L(const std::vector<std::vector<double>> &matrix) {
  double trace = 0.0;
  for (size_t i = 0; i < matrix.size(); ++i) {
    trace += matrix[i][i];
  }
  return trace;
}

// Coefficient Sum of the Last Eigenvector
double VE1L(const std::vector<std::vector<double>> &eigenvectors) {
  size_t n = eigenvectors.size();
  // const auto& lastEigenvector = eigenvectors[n - 1];

  std::vector<double> lastEigenvector(n);
  for (unsigned int i = 0; i < n; ++i) {
    lastEigenvector[i] = eigenvectors[i][n - 1];
  }

  return std::accumulate(
      lastEigenvector.begin(), lastEigenvector.end(), 0.0,
      [](double acc, double val) { return acc + std::abs(val); });
}

// Average Coefficient of the Last Eigenvector
double VE2L(double ve1, int numAtoms) { return ve1 / numAtoms; }

// Logarithmic Coefficient Sum of the Last Eigenvector
double VE3L(double ve1, int numAtoms) { return std::log(0.1 * numAtoms * ve1); }

// Randic-like Eigenvector-Based Index
double VR1L(const std::vector<std::vector<double>> &eigenvectors,
            const std::vector<std::pair<int, int>> &bonds) {
  size_t n = eigenvectors.size();

  std::vector<double> lastEigenvector(n);
  for (size_t i = 0; i < n; ++i) {
    lastEigenvector[i] = eigenvectors[i][n - 1];
  }

  // const auto& lastEigenvector = eigenvectors[n - 1];
  double result = 0.0;

  for (const auto &bond : bonds) {
    int i = bond.first;
    int j = bond.second;
    double product =
        std::abs(lastEigenvector[i]) * std::abs(lastEigenvector[j]);
    result += 1.0 / std::sqrt(product);
  }

  return result;
}

// Normalized Randic-like Eigenvector-Based Index
double VR2L(double vr1, int numAtoms) { return vr1 / numAtoms; }

// Logarithmic Randic-like Eigenvector-Based Index
double VR3L(double vr1, int numAtoms) { return std::log(0.1 * numAtoms * vr1); }

// BaryszMatrix deps

Eigen::MatrixXd floydWarshall(Eigen::MatrixXd &A) {
  int n = A.rows();

  // Set diagonal elements to 0
  A.diagonal().setZero();

  // Perform the relaxation to get shortest path BFS equivalent faster for
  // complex case than BFS
  for (int i = 0; i < n; ++i) {
    // Broadcast the row i and column i sums and apply the min operation
    Eigen::RowVectorXd rowSum = A.row(i);  // i-th row
    Eigen::VectorXd colSum = A.col(i);     // i-th column

    // Ensure correct dimensions for broadcasting
    Eigen::MatrixXd rowpart = rowSum.replicate(n, 1);
    Eigen::MatrixXd colpart = colSum.replicate(1, n);

    Eigen::MatrixXd tempMatrix = rowpart + colpart;

    A = A.cwiseMin(tempMatrix);
  }

  return A;
}

std::vector<std::vector<double>> floydWarshallL(
    std::vector<std::vector<double>> &matrix) {
  int n = matrix.size();

  for (int k = 0; k < n; ++k) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (matrix[i][k] < std::numeric_limits<double>::infinity() &&
            matrix[k][j] < std::numeric_limits<double>::infinity()) {
          matrix[i][j] = std::min(matrix[i][j], matrix[i][k] + matrix[k][j]);
        }
      }
    }
  }

  return matrix;
}

}  // namespace Osmordred
}  // namespace Descriptors
}  // namespace RDKit
