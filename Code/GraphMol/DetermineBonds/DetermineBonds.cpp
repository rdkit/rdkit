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
#include <YAeHMOP/EHTTools.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <unordered_map>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

namespace RDKit {

void connectivityHueckel(RWMol &mol, int charge) {
    auto numAtoms = mol.getNumAtoms();
    mol.getAtomWithIdx(0)->setFormalCharge(charge);
    EHTTools::EHTResults res;
    bool success = runMol(mol, res);
    // as of this writing runMol() always returns true, so we ignore the return value.
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
}

void connectivityVdW(RWMol &mol, double covFactor) {
    auto numAtoms = mol.getNumAtoms();
    double *distMat = MolOps::get3DDistanceMat(mol);
    
    double rcov[numAtoms];
    for (unsigned int i = 0; i < numAtoms; i++) {
        rcov[i] = covFactor * PeriodicTable::getTable()->getRcovalent(mol.getAtomWithIdx(i)->getAtomicNum());
    }
    for (unsigned int i = 0; i < numAtoms; i++) {
        for (unsigned int j = i + 1; j < numAtoms; j++) {
            if (distMat[i * numAtoms + j] <= (rcov[i] + rcov[j])) {
                mol.addBond(i, j, Bond::BondType::SINGLE);
            }
        }
    }
}

void determineConnectivity(RWMol &mol, bool useHueckel, int charge, double covFactor) {
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
    } else {
        connectivityVdW(mol, covFactor);
    }
}

std::vector<int> possibleValences(Atom *atom, std::unordered_map<int,
                                  std::vector<int>> &atomicValence) {
    int atomNum = atom->getAtomicNum();
    int numBonds = atom->getDegree();
    std::vector<int> valences = atomicValence[atomNum];
    std::vector<int> possible;
    for (auto &valence : valences) {
        if (numBonds <= valence) {
            possible.push_back(valence);
        }
    }
    return possible;
}

void valenceCombinations(RWMol &mol, std::vector<std::vector<int>> &possible,
                         std::vector<std::vector<int>> &combos) {
    int numAtoms = mol.getNumAtoms();
    std::unordered_map<int, std::vector<int>> atomicValence =
    {
        {1, {1}},
        {5, {3,4}},
        {6, {4}},
        {7, {3,4}},
        {8, {2,1,3}},
        {9, {1}},
        {14, {4}},
        {15, {5,3}},
        {16, {6,3,2}},
        {17, {1}},
        {32, {4}},
        {35, {1}},
        {53, {1}}
    };

    possible.resize(numAtoms);
    int numCombos = 1;
    for (int i = 0; i < numAtoms; i++) {
        possible[i] = possibleValences(mol.getAtomWithIdx(i), atomicValence);
        if (possible[i].size() != 0) {
            numCombos *= possible[i].size();
        } else {
            auto atom = mol.getAtomWithIdx(i);
            std::vector<int> valences = atomicValence[atom->getAtomicNum()];
            std::stringstream ss;
            ss << "Valence of atom " << i << " is " << atom->getDegree()
            << ", which is larger than the allowed maximum, "
            << valences[valences.size() - 1];
            throw std::runtime_error(ss.str());
        }
    }

    int orig = numCombos;
    combos.resize(numCombos, std::vector<int>(numAtoms));
    for (int i = 0; i < numAtoms; i++) {
         int seg = numCombos / possible[i].size();
         for (int p = 0; p < orig / numCombos; p++) {
           for (int j = 0; j < possible[i].size(); j++) {
             for (int k = 0; k < seg; k++) {
                 combos[p*possible[i].size()*seg + j*seg + k][i] = possible[i][j];
             }
           }
         }
         numCombos /= possible[i].size();
     }
}

void getUnsaturated(std::vector<int> &order, std::vector<int> &valency,
                    std::vector<int> &unsat) {
    for (int i = 0; i < order.size(); i++) {
        int leftover = order[i] - valency[i];
        if (leftover > 0) {
            unsat.push_back(i);
        }
    }
}

void getUnsaturatedPairs(std::vector<std::vector<int>> &ordMat, std::vector<int> &unsat,
                         std::vector<std::pair<int, int>> &unsatPairs) {
    for (int i = 0; i < unsat.size(); i++) {
        for (int j = i + 1; j < unsat.size(); j++) {
            if (ordMat[unsat[i]][unsat[j]]) {
                unsatPairs.push_back(std::make_pair(unsat[i], unsat[j]));
            }
        }
    }
}

bool checkValency(std::vector<int> &order, std::vector<int> &valency) {
    for (int i = 0; i < valency.size(); i++) {
        if (valency[i] > order[i]) {
            return false;
        }
    }
    return true;
}

int getAtomicCharge(int atom, int valence) {
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

bool checkCharge(RWMol &mol, std::vector<int> &valency, int charge) {
    int molCharge = 0;
    for (int i = 0; i < mol.getNumAtoms(); i++) {
        auto atom = mol.getAtomWithIdx(i);
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

void setAtomicCharges(RWMol &mol, std::vector<int> &valency, int charge) {
    int molCharge = 0;
    for (int i = 0; i < mol.getNumAtoms(); i++) {
        auto atom = mol.getAtomWithIdx(i);
        int atomCharge = getAtomicCharge(atom->getAtomicNum(), valency[i]);
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

void setAtomicRadicals(RWMol &mol, std::vector<int> &valency, int charge) {
    for (int i = 0; i < mol.getNumAtoms(); i++) {
        auto atom = mol.getAtomWithIdx(i);
        int atomCharge = getAtomicCharge(atom->getAtomicNum(), valency[i]);
        if (atomCharge != 0) {
            atom->setNumRadicalElectrons(std::abs(charge));
        }
    }
}

bool checkSaturation(std::vector<int> &order, std::vector<int> &valency) {
    std::vector<int> unsat;
    getUnsaturated(order, valency, unsat);
    if (unsat.size() == 0) {
        return true;
    } else {
        return false;
    }
}

void setAtomMap(RWMol &mol) {
    for (int i = 0; i < mol.getNumAtoms(); i++) {
        auto atom = mol.getAtomWithIdx(i);
        atom->setAtomMapNum(i + 1);
    }
}

void setChirality(RWMol &mol) {
    MolOps::sanitizeMol(mol);
    MolOps::setDoubleBondNeighborDirections(mol, &mol.getConformer());
    MolOps::assignStereochemistryFrom3D(mol);
    MolOps::assignChiralTypesFrom3D(mol);
}

void addBondOrdering(RWMol &mol, std::vector<std::vector<int>> &ordMat,
                     std::vector<int> &valency, bool allowChargedFragments,
                     bool embedChiral, bool useAtomMap, int charge) {
    for (int i = 0; i < mol.getNumAtoms(); i++) {
        for (int j = i + 1; j < mol.getNumAtoms(); j++) {
            if (ordMat[i][j] == 0 || ordMat[i][j] == 1) {
                continue;
            } else if (ordMat[i][j] == 2) {
                mol.removeBond(i, j);
                mol.addBond(i, j, Bond::BondType::DOUBLE);
            } else if (ordMat[i][j] == 3) {
                mol.removeBond(i, j);
                mol.addBond(i, j, Bond::BondType::TRIPLE);
            } else {
                mol.removeBond(i, j);
                mol.addBond(i, j, Bond::BondType::SINGLE);
            }
        }
    }
    
    if (allowChargedFragments) {
        setAtomicCharges(mol, valency, charge);
    } else {
        setAtomicRadicals(mol, valency, charge);
    }
    
    if (useAtomMap) {
        setAtomMap(mol);
    }
    
    if (embedChiral) {
        setChirality(mol);
    }
    
    int successKek = MolOps::KekulizeIfPossible(mol);
    // TODO:: what to do if this fails?
    int successArom = MolOps::setAromaticity(mol);
    // TODO:: what to do if this fails?
}

void determineBondOrder(RWMol &mol, int charge, bool allowChargedFragments,
                        bool embedChiral, bool useAtomMap) {
    int numAtoms = mol.getNumAtoms();
    
    std::vector<std::vector<int>> conMat(numAtoms, std::vector<int>(numAtoms, 0));
    std::vector<int> origValency(numAtoms, 0);
    for (int i = 0; i < numAtoms; i++) {
        for (int j = i + 1; j < numAtoms; j++) {
            if (mol.getBondBetweenAtoms(i, j)) {
                conMat[i][j]++;
                origValency[i]++;
                origValency[j]++;
            }
        }
    }
    
    std::vector<std::vector<int>> best(conMat);
    std::vector<int> bestValency(origValency);
    int bestSum = std::accumulate(origValency.begin(), origValency.end(), 0);

    std::vector<std::vector<int>> possible;
    std::vector<std::vector<int>> orders;
    valenceCombinations(mol, possible, orders);

    bool valencyValid = false;
    bool chargeValid = false;
    bool saturationValid = false;
  
    for (auto &order : orders) {
        std::vector<int> unsat;
        getUnsaturated(order, origValency, unsat);
        // checks whether the atomic connectivity is valid for the current set of atomic valences
        if (unsat.size() == 0) {
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
        
        std::vector<std::vector<int>> ordMat(conMat);
        std::vector<int> valency(origValency);
        bool newBonds = false;
        do {
            newBonds = false;
            std::vector<int> unsat;
            getUnsaturated(order, valency, unsat);
            std::vector<std::pair<int, int>> unsatPairs;
            getUnsaturatedPairs(conMat, unsat, unsatPairs);
            
            if (unsatPairs.size() > 0) {
                Graph graph(unsatPairs.begin(), unsatPairs.end(), numAtoms);
                std::vector<boost::graph_traits<Graph>::vertex_descriptor> mate(numAtoms);
                edmonds_maximum_cardinality_matching(graph, &mate[0]);
                
                boost::graph_traits<Graph>::vertex_iterator vi, viEnd;
                for (boost::tie(vi, viEnd) = vertices(graph); vi != viEnd; ++vi) {
                    if (mate[*vi] != boost::graph_traits<Graph>::null_vertex()
                        && *vi < mate[*vi]) {
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
                int sum = std::accumulate(valency.begin(), valency.end(), 0);;
                if (sum > bestSum) {
                    best = ordMat;
                    bestSum = sum;
                    bestValency = valency;
                }
            }
                
        }
        
    }
    
    addBondOrdering(mol, best, bestValency, allowChargedFragments,
                    embedChiral, useAtomMap, charge);
    return;
} // determineBondOrdering()

} // namespace RDKit


