//
//  Copyright (C) 2022 Sreya Gogineni and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "DetermineBonds.h"
#include <GraphMol/RDKitBase.h>
#ifdef RDK_BUILD_YAEHMOP_SUPPORT
#include <YAeHMOP/EHTTools.h>
#endif
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <unordered_map>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <GraphMol/FileParsers/ProximityBonds.h>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>
    Graph;
using boost::multiprecision::uint1024_t;

namespace {

// see http://phrogz.net/lazy-cartesian-product
template <typename T>
struct LazyCartesianProduct {
  std::vector<std::vector<T>> d_listOfLists;
  std::vector<uint1024_t> d_divs;
  std::vector<uint1024_t> d_mods;
  uint1024_t d_maxSize;
  uint1024_t d_currentPos;

  explicit LazyCartesianProduct(const std::vector<std::vector<T>> &input)
      : d_listOfLists(input), d_currentPos(0) {
    auto size = d_listOfLists.size();
    d_divs.resize(size);
    d_mods.resize(size);
    d_maxSize = 1;

    for (int i = size - 1; i >= 0; --i) {
      uint1024_t items(d_listOfLists[i].size());
      d_divs[i] = d_maxSize;
      d_mods[i] = items;
      d_maxSize *= items;
    }
  }

  std::vector<T> entryAt(uint1024_t pos) const;
  std::vector<T> next() { return entryAt(d_currentPos++); }
  bool atEnd() const { return d_currentPos >= d_maxSize; }
};

template <typename T>
std::vector<T> LazyCartesianProduct<T>::entryAt(uint1024_t pos) const {
  auto length = d_listOfLists.size();
  std::vector<T> res(length);
  for (auto i = 0u; i < length; ++i) {
    res[i] = d_listOfLists[i][static_cast<size_t>(
        static_cast<uint1024_t>(pos / d_divs[i]) % d_mods[i])];
  }
  return res;
}

std::vector<unsigned int> possibleValences(
    const RDKit::Atom *atom,
    const std::unordered_map<int, std::vector<unsigned int>> &atomicValence) {
  auto atomNum = atom->getAtomicNum() - atom->getFormalCharge();
  auto numBonds = atom->getDegree();

  std::vector<unsigned int> avalences;
  auto valences = atomicValence.find(atomNum);
  if (valences != atomicValence.end()) {
    avalences = valences->second;
  } else {
    for (auto v : RDKit::PeriodicTable::getTable()->getValenceList(atomNum)) {
      if (v >= 0) {
        avalences.push_back(v);
      }
    }
  }
  std::vector<unsigned int> possible;
  for (const auto &valence : avalences) {
    if (numBonds <= valence) {
      possible.push_back(valence);
    }
  }
  return possible;
}

LazyCartesianProduct<unsigned int> getValenceCombinations(
    const RDKit::RWMol &mol) {
  auto numAtoms = mol.getNumAtoms();
  const std::unordered_map<int, std::vector<unsigned int>> atomicValence = {
      {1, {1}},  {5, {3, 4}}, {6, {4}},     {7, {3, 4}},     {8, {2, 1, 3}},
      {9, {1}},  {14, {4}},   {15, {5, 3}}, {16, {6, 3, 2}}, {17, {1}},
      {32, {4}}, {35, {1}},   {53, {1}}};
  std::vector<std::vector<unsigned int>> possible(numAtoms);
  for (unsigned int i = 0; i < numAtoms; i++) {
    possible[i] = possibleValences(mol.getAtomWithIdx(i), atomicValence);
    if (possible[i].empty()) {
      const auto atom = mol.getAtomWithIdx(i);
      std::vector<unsigned int> valences =
          atomicValence.at(atom->getAtomicNum());
      std::stringstream ss;
      ss << "Valence of atom " << i << " is " << atom->getDegree()
         << ", which is larger than the allowed maximum, "
         << valences[valences.size() - 1];
      throw ValueErrorException(ss.str());
    }
  }

  return LazyCartesianProduct<unsigned int>(possible);
}

}  // namespace

namespace RDKit {

#ifdef RDK_BUILD_YAEHMOP_SUPPORT
void connectivityHueckel(RWMol &mol, int charge) {
  auto numAtoms = mol.getNumAtoms();
  mol.getAtomWithIdx(0)->setFormalCharge(charge);
  EHTTools::EHTResults res;
  bool success = runMol(mol, res);
  RDUNUSED_PARAM(success);
  // as of this writing runMol() always returns true, so we ignore the return
  // value.
  double *mat = res.reducedOverlapPopulationMatrix.get();
  int matInd = 0;
  for (unsigned int i = 0; i < numAtoms; i++) {
    for (unsigned int j = 0; j < i + 1; j++) {
      if (i != j && mat[matInd] >= 0.15) {
        mol.addBond(i, j, Bond::BondType::SINGLE);
      }
      matInd++;
    }
  }
}  // connectivityHueckel()
#else
void connectivityHueckel(RWMol &, int) {
  CHECK_INVARIANT(0, "YAeHMOP support not available");
}
#endif

void connectivityVdW(RWMol &mol, double covFactor) {
  auto numAtoms = mol.getNumAtoms();
  double *distMat = MolOps::get3DDistanceMat(mol);

  std::vector<double> rcov(numAtoms);
  for (unsigned int i = 0; i < numAtoms; i++) {
    rcov[i] = covFactor * PeriodicTable::getTable()->getRcovalent(
                              mol.getAtomWithIdx(i)->getAtomicNum());
  }
  for (unsigned int i = 0; i < numAtoms; i++) {
    for (unsigned int j = i + 1; j < numAtoms; j++) {
      if (distMat[i * numAtoms + j] <= (rcov[i] + rcov[j])) {
        mol.addBond(i, j, Bond::BondType::SINGLE);
      }
    }
  }
}  // connectivityVdW()

void determineConnectivity(RWMol &mol, bool useHueckel, int charge,
                           double covFactor, bool useVdw) {
#ifndef RDK_BUILD_YAEHMOP_SUPPORT
  if (useHueckel) {
    throw ValueErrorException(
        "The RDKit was not compiled with YAeHMOP support");
  }
#endif
  auto numAtoms = mol.getNumAtoms();
  for (unsigned int i = 0; i < numAtoms; i++) {
    for (unsigned int j = i + 1; j < numAtoms; j++) {
      mol.removeBond(i, j);
      mol.getAtomWithIdx(i)->setNoImplicit(true);
      mol.getAtomWithIdx(j)->setNoImplicit(true);
    }
  }
  if (useHueckel) {
    connectivityHueckel(mol, charge);
  } else if (useVdw) {
    connectivityVdW(mol, covFactor);
  } else {
    ConnectTheDots(&mol, ctdIGNORE_H_H_CONTACTS);
  }
}  // determineConnectivity()

void getUnsaturated(const std::vector<unsigned int> &order,
                    const std::vector<unsigned int> &valency,
                    std::vector<unsigned int> &unsat) {
  for (unsigned int i = 0; i < order.size(); i++) {
    if (order[i] > valency[i]) {
      unsat.push_back(i);
    }
  }
}

void getUnsaturatedPairs(
    const std::vector<std::vector<unsigned int>> &ordMat,
    const std::vector<unsigned int> &unsat,
    std::vector<std::pair<unsigned int, unsigned int>> &unsatPairs) {
  for (unsigned int i = 0; i < unsat.size(); i++) {
    for (unsigned int j = i + 1; j < unsat.size(); j++) {
      if (ordMat[unsat[i]][unsat[j]]) {
        unsatPairs.push_back(std::make_pair(unsat[i], unsat[j]));
      }
    }
  }
}

bool checkValency(const std::vector<unsigned int> &order,
                  const std::vector<unsigned int> &valency) {
  for (unsigned int i = 0; i < valency.size(); i++) {
    if (valency[i] > order[i]) {
      return false;
    }
  }
  return true;
}

int getAtomicCharge(int atom, unsigned int valence) {
  if (atom == 1) {
    return 1 - valence;
  } else if (atom == 5) {
    return 3 - valence;
  } else if (atom == 15 && valence == 5) {
    return 0;
  } else if (atom == 16 && valence == 6) {
    return 0;
  } else {
    return PeriodicTable::getTable()->getNouterElecs(atom) - 8 + valence;
  }
}

bool checkCharge(RWMol &mol, const std::vector<unsigned int> &valency,
                 int charge) {
  int molCharge = 0;
  for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    const auto atom = mol.getAtomWithIdx(i);
    int atomCharge = getAtomicCharge(atom->getAtomicNum(), valency[i]);
    molCharge += atomCharge;
    if (atom->getAtomicNum() == 6) {
      if (atom->getDegree() == 2 && valency[i] == 2) {
        molCharge += 1;
        atomCharge = 2;
      } else if (atom->getDegree() == 3 && (molCharge + 1 < charge)) {
        molCharge += 2;
        atomCharge = 1;
      }
    }
  }
  return molCharge == charge;
}

void setAtomicCharges(RWMol &mol, const std::vector<unsigned int> &valency,
                      int charge) {
  int molCharge = 0;
  for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    auto atom = mol.getAtomWithIdx(i);
    int atomCharge = getAtomicCharge(
        atom->getAtomicNum() - atom->getFormalCharge(), valency[i]);
    molCharge += atomCharge;
    if (atom->getAtomicNum() == 6) {
      if (atom->getDegree() == 2 && valency[i] == 2) {
        molCharge += 1;
        atomCharge = 0;
      } else if (atom->getDegree() == 3 && (molCharge + 1 < charge)) {
        molCharge += 2;
        atomCharge = 1;
      }
    }
    if (atomCharge != 0) {
      atom->setFormalCharge(atomCharge);
    }
  }
}

void setAtomicRadicals(RWMol &mol, const std::vector<unsigned int> &valency,
                       int charge) {
  for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    auto atom = mol.getAtomWithIdx(i);
    int atomCharge = getAtomicCharge(atom->getAtomicNum(), valency[i]);
    if (atomCharge != 0) {
      atom->setNumRadicalElectrons(std::abs(charge));
    }
  }
}

bool checkSaturation(const std::vector<unsigned int> &order,
                     const std::vector<unsigned int> &valency) {
  std::vector<unsigned int> unsat;
  getUnsaturated(order, valency, unsat);
  return unsat.empty();
}

void setAtomMap(RWMol &mol) {
  for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    auto atom = mol.getAtomWithIdx(i);
    atom->setAtomMapNum(i + 1);
  }
}

void setChirality(RWMol &mol) {
  MolOps::sanitizeMol(mol);
  MolOps::setDoubleBondNeighborDirections(mol, &mol.getConformer());
  MolOps::assignStereochemistryFrom3D(mol);
}

void addBondOrdering(RWMol &mol,
                     const std::vector<std::vector<unsigned int>> &ordMat,
                     const std::vector<unsigned int> &valency,
                     bool allowChargedFragments, bool embedChiral,
                     bool useAtomMap, int charge) {
  auto numAtoms = mol.getNumAtoms();

  for (unsigned int i = 0; i < numAtoms; i++) {
    for (unsigned int j = i + 1; j < numAtoms; j++) {
      if (ordMat[i][j] == 0) {
        continue;
      } else if (ordMat[i][j] == 1) {
        mol.getBondBetweenAtoms(i, j)->setBondType(Bond::BondType::SINGLE);
      } else if (ordMat[i][j] == 2) {
        mol.getBondBetweenAtoms(i, j)->setBondType(Bond::BondType::DOUBLE);
      } else if (ordMat[i][j] == 3) {
        mol.getBondBetweenAtoms(i, j)->setBondType(Bond::BondType::TRIPLE);
      } else {
        mol.getBondBetweenAtoms(i, j)->setBondType(Bond::BondType::SINGLE);
      }
    }
  }

  if (allowChargedFragments) {
    setAtomicCharges(mol, valency, charge);
  } else {
    setAtomicRadicals(mol, valency, charge);
  }

  if (MolOps::getFormalCharge(mol) != charge) {
    mol.debugMol(std::cerr);
    std::stringstream ss;
    ss << "Final molecular charge (" << charge << ") does not match input ("
       << MolOps::getFormalCharge(mol)
       << "); could not find valid bond ordering";
    throw ValueErrorException(ss.str());
  }

  if (useAtomMap) {
    setAtomMap(mol);
  }

  if (embedChiral) {
    setChirality(mol);
  }
}

void determineBondOrders(RWMol &mol, int charge, bool allowChargedFragments,
                         bool embedChiral, bool useAtomMap) {
  auto numAtoms = mol.getNumAtoms();

  std::vector<std::vector<unsigned int>> conMat(
      numAtoms, std::vector<unsigned int>(numAtoms, 0));
  std::vector<unsigned int> origValency(numAtoms, 0);
  for (unsigned int i = 0; i < numAtoms; i++) {
    for (unsigned int j = i + 1; j < numAtoms; j++) {
      if (mol.getBondBetweenAtoms(i, j)) {
        conMat[i][j]++;
        origValency[i]++;
        origValency[j]++;
      }
    }
  }

  std::vector<std::vector<unsigned int>> best(conMat);
  std::vector<unsigned int> bestValency(origValency);
  int bestSum = std::accumulate(origValency.begin(), origValency.end(), 0);

  auto valenceCombos = getValenceCombinations(mol);

  bool valencyValid = false;
  bool chargeValid = false;
  bool saturationValid = false;

  while (!valenceCombos.atEnd()) {
    auto order = valenceCombos.next();
    std::vector<unsigned int> unsat;
    getUnsaturated(order, origValency, unsat);
    // checks whether the atomic connectivity is valid for the current set of
    // atomic valences
    if (unsat.empty()) {
      valencyValid = checkValency(order, origValency);
      chargeValid = checkCharge(mol, origValency, charge);
      saturationValid = checkSaturation(order, origValency);

      if (valencyValid && chargeValid && saturationValid) {
        addBondOrdering(mol, conMat, origValency, allowChargedFragments,
                        embedChiral, useAtomMap, charge);
        return;
      } else {
        continue;
      }
    }

    std::vector<std::vector<unsigned int>> ordMat(conMat);
    std::vector<unsigned int> valency(origValency);
    bool newBonds = false;
    do {
      newBonds = false;
      std::vector<unsigned int> unsat;
      getUnsaturated(order, valency, unsat);
      std::vector<std::pair<unsigned int, unsigned int>> unsatPairs;
      getUnsaturatedPairs(conMat, unsat, unsatPairs);

      if (!unsatPairs.empty()) {
        Graph graph(unsatPairs.begin(), unsatPairs.end(), numAtoms);
        std::vector<boost::graph_traits<Graph>::vertex_descriptor> mate(
            numAtoms);
        edmonds_maximum_cardinality_matching(graph, &mate[0]);

        boost::graph_traits<Graph>::vertex_iterator vi, viEnd;
        for (boost::tie(vi, viEnd) = vertices(graph); vi != viEnd; ++vi) {
          if (mate[*vi] != boost::graph_traits<Graph>::null_vertex() &&
              *vi < mate[*vi]) {
            newBonds = true;
            ordMat[*vi][mate[*vi]]++;
            valency[*vi]++;
            valency[mate[*vi]]++;
          }
        }
      }

    } while (newBonds == true);

    valencyValid = checkValency(order, valency);
    chargeValid = checkCharge(mol, valency, charge);
    saturationValid = checkSaturation(order, valency);

    if (valencyValid && chargeValid) {
      if (saturationValid) {
        addBondOrdering(mol, ordMat, valency, allowChargedFragments,
                        embedChiral, useAtomMap, charge);
        return;
      } else {
        int sum = std::accumulate(valency.begin(), valency.end(), 0);
        ;
        if (sum > bestSum) {
          best = ordMat;
          bestSum = sum;
          bestValency = valency;
        }
      }
    }
  }

  addBondOrdering(mol, best, bestValency, allowChargedFragments, embedChiral,
                  useAtomMap, charge);
  return;
}  // determineBondOrdering()

void determineBonds(RWMol &mol, bool useHueckel, int charge, double covFactor,
                    bool allowChargedFragments, bool embedChiral,
                    bool useAtomMap, bool useVdw) {
  if (mol.getNumAtoms() <= 1) {
    return;
  }
  determineConnectivity(mol, useHueckel, charge, covFactor, useVdw);
  determineBondOrders(mol, charge, allowChargedFragments, embedChiral,
                      useAtomMap);
}  // determineBonds()

}  // namespace RDKit
