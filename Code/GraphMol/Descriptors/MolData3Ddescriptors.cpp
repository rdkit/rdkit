#include <cstdlib>
#include <iostream>
#include "MolData3Ddescriptors.h"
#include <GraphMol/RDKitBase.h>

#include "GraphMol/PartialCharges/GasteigerCharges.h"
#include "GraphMol/PartialCharges/GasteigerParams.h"

using namespace std;

MolData3Ddescriptors::MolData3Ddescriptors() {}

std::vector<double> MolData3Ddescriptors::GetUn(int numAtoms) {
  std::vector<double> u(numAtoms, 1.0);

  return u;
}

std::vector<double> MolData3Ddescriptors::GetRelativeMW(
    const RDKit::ROMol &mol) {
  double *relativeMw = data3D.getMW();
  int numAtoms = mol.getNumAtoms();

  std::vector<double> pol(numAtoms, 0.0);
  for (int i = 0; i < numAtoms; ++i) {
    pol[i] = relativeMw[mol.getAtomWithIdx(i)->getAtomicNum() - 1];
  }
  return pol;
}

std::vector<double> MolData3Ddescriptors::GetRelativePol(
    const RDKit::ROMol &mol) {
  int numAtoms = mol.getNumAtoms();
  double *relativePol = data3D.getPOL();

  std::vector<double> pol(numAtoms, 0.0);
  for (int i = 0; i < numAtoms; ++i) {
    pol[i] = relativePol[mol.getAtomWithIdx(i)->getAtomicNum() - 1];
  }
  return pol;
}

std::vector<double> MolData3Ddescriptors::GetRelativeVdW(
    const RDKit::ROMol &mol) {
  int numAtoms = mol.getNumAtoms();
  double *relativeVdW = data3D.getVDW();

  std::vector<double> vdw(numAtoms, 0.0);
  for (int i = 0; i < numAtoms; ++i) {
    vdw[i] = relativeVdW[mol.getAtomWithIdx(i)->getAtomicNum() - 1];
  }
  return vdw;
}

std::vector<double> MolData3Ddescriptors::GetRelativeRcov(
    const RDKit::ROMol &mol) {
  int numAtoms = mol.getNumAtoms();
  double *rcov = data3D.getRCOV();

  std::vector<double> wroc(numAtoms, 0.0);
  for (int i = 0; i < numAtoms; ++i) {
    wroc[i] = rcov[mol.getAtomWithIdx(i)->getAtomicNum() - 1] / rcov[5];
  }
  return wroc;
}

std::vector<double> MolData3Ddescriptors::GetRelativeENeg(
    const RDKit::ROMol &mol) {
  int numAtoms = mol.getNumAtoms();
  double *relativeNeg = data3D.getNEG();

  std::vector<double> neg(numAtoms, 0.0);
  for (int i = 0; i < numAtoms; ++i) {
    neg[i] = relativeNeg[mol.getAtomWithIdx(i)->getAtomicNum() - 1];
  }
  return neg;
}

std::vector<double> MolData3Ddescriptors::GetRelativeIonPol(
    const RDKit::ROMol &mol) {
  int numAtoms = mol.getNumAtoms();
  double *absionpol = data3D.getIonPOL();

  std::vector<double> ionpols(numAtoms, 0.0);
  for (int i = 0; i < numAtoms; ++i) {
    ionpols[i] = absionpol[mol.getAtomWithIdx(i)->getAtomicNum() - 1];
  }
  return ionpols;
}

std::vector<double> MolData3Ddescriptors::GetCustomAtomProp(
    const RDKit::ROMol &mol, const std::string &customAtomPropName) {
  int numAtoms = mol.getNumAtoms();

  std::vector<double> customAtomArray(numAtoms, 1.0);
  for (auto &atom : mol.atoms()) {
    atom->getPropIfPresent(customAtomPropName, customAtomArray[atom->getIdx()]);
  }
  return customAtomArray;
}

std::vector<double> MolData3Ddescriptors::GetCharges(const RDKit::ROMol &mol) {
  std::vector<double> charges(mol.getNumAtoms(), 0);
  // use 12 iterations... can be more
  RDKit::computeGasteigerCharges(mol, charges, 12, true);
  return charges;
}

int MolData3Ddescriptors::GetPrincipalQuantumNumber(int AtomicNum) {
  if (AtomicNum <= 2) {
    return 1;
  } else if (AtomicNum <= 10) {
    return 2;
  } else if (AtomicNum <= 18) {
    return 3;
  } else if (AtomicNum <= 36) {
    return 4;
  } else if (AtomicNum <= 54) {
    return 5;
  } else if (AtomicNum <= 86) {
    return 6;
  } else {
    return 7;
  }
}

std::vector<double> MolData3Ddescriptors::GetIState(const RDKit::ROMol &mol) {
  int numAtoms = mol.getNumAtoms();
  std::vector<double> Is(numAtoms, 1.0);  // values set to 1 for Hs

  for (int i = 0; i < numAtoms; ++i) {
    const RDKit::Atom *atom = mol.getAtomWithIdx(i);
    int atNum = atom->getAtomicNum();
    int degree = atom->getDegree();  // number of substituants (heavy of not?)
    if (degree > 0 && atNum > 1) {
      int h = atom->getTotalNumHs(
          true);  // caution getTotalNumHs(true) to count h !!!!
      int dv = RDKit::PeriodicTable::getTable()->getNouterElecs(atNum) -
               h;  // number of valence (explicit with Hs)
      int N = GetPrincipalQuantumNumber(atNum);  // principal quantum number
      double d = (double)degree - h;             // degree-h
      if (d > 0) {
        Is[i] = std::round(1000 * (4.0 / (N * N) * dv + 1.0) / d) / 1000;
      }
    }
  }

  return Is;
}

std::vector<double> MolData3Ddescriptors::GetIStateDrag(
    const RDKit::ROMol &mol) {
  int numAtoms = mol.getNumAtoms();
  std::vector<double> Is(numAtoms, 1.0);

  for (int i = 0; i < numAtoms; ++i) {
    const RDKit::Atom *atom = mol.getAtomWithIdx(i);
    int atNum = atom->getAtomicNum();
    int degree = atom->getDegree();  // number of substituants
    if (degree > 0 && atNum > 1) {
      int h = atom->getTotalNumHs(true);
      int Zv = RDKit::PeriodicTable::getTable()->getNouterElecs(
          atNum);                  // number of valence (explicit with Hs)
      double dv = (double)Zv - h;  // number of valence electron without Hs
      int N = GetPrincipalQuantumNumber(atNum);  // principal quantum number
      double d = (double)degree - h;             // degree-h
      if (d > 0) {
        Is[i] = std::round(1000 * (4.0 / (N * N) * dv + 1.0) / d) / 1000;
      }
    }
  }

  return Is;
}

// adaptation from EState.py
// we need the Is value only there
std::vector<double> MolData3Ddescriptors::GetEState(const RDKit::ROMol &mol) {
  int numAtoms = mol.getNumAtoms();

  std::vector<double> Is = GetIState(mol);

  double tmp, p;
  double *dist = RDKit::MolOps::getDistanceMat(mol, false, false);
  std::vector<double> accum(numAtoms, 0.0);

  for (int i = 0; i < numAtoms; i++) {
    for (int j = i + 1; j < numAtoms; j++) {
      p = dist[i * numAtoms + j] + 1;
      if (p < 1e6) {
        tmp = (Is[i] - Is[j]) / (p * p);
        accum[i] += tmp;
        accum[j] -= tmp;
      }
    }
  }

  for (int i = 0; i < numAtoms; i++) {
    Is[i] += accum[i];
  }

  return Is;
}

// modification of previous code to follow documentation from Padel code
std::vector<double> MolData3Ddescriptors::GetEState2(const RDKit::ROMol &mol) {
  int numAtoms = mol.getNumAtoms();

  std::vector<double> Si = GetIState(mol);

  // in WHIM definition it's write:
  double tmp, p, d;
  double *dist = RDKit::MolOps::getDistanceMat(mol, false, false);
  std::vector<double> accum(numAtoms, 0.0);

  for (int i = 0; i < numAtoms; i++) {
    for (int j = i + 1; j < numAtoms; j++) {
      d = dist[i * numAtoms + j];
      p = dist[i * numAtoms + j] + 1;
      if (d == 1) {
        tmp = (Si[i] - Si[j]) / (p * p);
        accum[i] += tmp;
        accum[j] -= tmp;
      }
    }
  }

  // add the Accum to the Si
  // WHIM Si values
  // electrotopological indices are scaled thus: Si'=Si + 7 => Si' > 0
  // In this case, only the nonhydrogen atoms are considered,
  // and the atomic electrotopological charge of each atom depends on its atom
  // neighbor.
  // So we should not use all the terms in the sum but only Adj matrix cases!

  // Correct the Si adding the rescaling parameter for WHIM only
  for (int i = 0; i < numAtoms; i++) {
    Si[i] += accum[i] + 7.0;
  }

  return Si;
}
