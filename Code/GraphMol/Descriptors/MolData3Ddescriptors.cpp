#include <cstdlib>
#include <iostream>
#include "MolData3Ddescriptors.h"
#include <GraphMol/RDKitBase.h>


using namespace std;

MolData3Ddescriptors::MolData3Ddescriptors(){}


std::vector<double> MolData3Ddescriptors::GetRelativeMW(const RDKit::ROMol& mol){
  double* relativeMw=data3D.getMW();
  int numAtoms= mol.getNumAtoms();

  std::vector<double> pol(numAtoms, 0.0);
  for( int i=0; i<numAtoms; ++i){

    pol[i]=relativeMw[mol.getAtomWithIdx(i)->getAtomicNum()-1];
  }
  return pol;
}

std::vector<double> MolData3Ddescriptors::GetRelativePol(const RDKit::ROMol& mol){
   int numAtoms= mol.getNumAtoms();
  double* relativePol=data3D.getPOL();

  std::vector<double> pol(numAtoms, 0.0);
  for( int i=0; i<numAtoms; ++i){

    pol[i]=relativePol[mol.getAtomWithIdx(i)->getAtomicNum()-1];

  }
  return pol;
}

std::vector<double> MolData3Ddescriptors::GetRelativeVdW(const RDKit::ROMol& mol){
   int numAtoms= mol.getNumAtoms();
   double* relativeVdW=data3D.getVDW();

  std::vector<double> vdw(numAtoms, 0.0);
  for( int i=0; i<numAtoms; ++i){

    vdw[i]=relativeVdW[mol.getAtomWithIdx(i)->getAtomicNum()-1];

  }
  return vdw;
}


std::vector<double> MolData3Ddescriptors::GetRelativeRcov(const RDKit::ROMol& mol){
  int numAtoms= mol.getNumAtoms();
  double* rcov=data3D.getRCOV();

  std::vector<double> wroc(numAtoms, 0.0);
  for( int i=0; i<numAtoms; ++i){

    wroc[i]=rcov[mol.getAtomWithIdx(i)->getAtomicNum()-1]/rcov[5];
  }
  return wroc;
}


std::vector<double> MolData3Ddescriptors::GetRelativeENeg(const RDKit::ROMol& mol){
   int numAtoms= mol.getNumAtoms();
   double* relativeNeg=data3D.getNEG();

  std::vector<double> neg(numAtoms, 0.0);
  for( int i=0; i<numAtoms; ++i){

    neg[i]=relativeNeg[mol.getAtomWithIdx(i)->getAtomicNum()-1];
  }
  return neg;
}


std::vector<double> MolData3Ddescriptors::GetRelativeIonPol(const RDKit::ROMol& mol){
   int numAtoms= mol.getNumAtoms();
	double* absionpol=data3D.getIonPOL();

  std::vector<double> ionpols(numAtoms, 0.0);
  for( int i=0; i<numAtoms; ++i){

    ionpols[i]=absionpol[mol.getAtomWithIdx(i)->getAtomicNum()-1]/absionpol[5];
  }
  return ionpols;
}