//
//  Copyright (c) 2012, Institue of Cancer Research.
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
//     * Neither the name of Institue of Cancer Research.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
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
//
// For more information on the Plane of Best Fit please see
// http://pubs.acs.org/doi/abs/10.1021/ci300293f
//
//  If this code has been useful to you, please include the reference
//  in any work which has made use of it:

//  Plane of Best Fit: A Novel Method to Characterize the Three-Dimensionality
//  of Molecules, Nicholas C. Firth, Nathan Brown, and Julian Blagg, Journal of
//  Chemical Information and Modeling 2012 52 (10), 2516-2525

//
//
// Created by Nicholas Firth, November 2011
// Modified by Greg Landrum for inclusion in the RDKit distribution November
// 2012
// Further modified by Greg Landrum for inclusion in the RDKit core September
// 2016
// Adding MORSE descriptors to 3D descriptors by Guillaume Godin

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include "MORSE.h"
#include "GraphMol/PartialCharges/GasteigerCharges.h"
#include "GraphMol/PartialCharges/GasteigerParams.h"
#include <Numerics/Matrix.h>
#include <Numerics/SquareMatrix.h>
#include <Numerics/SymmMatrix.h>
#include <boost/foreach.hpp>
#include <math.h>
#include <Eigen/Dense>

// data checked using book Todeschini R., Consonni V. - Molecular Descriptors for Chemoinformatics 2009 atomic properties page 21/22

double Pol1[] = {
    0.67,  0,    24.3, 5.60,  3.03, 1.76,  1.10, 0.80, 0.56, 0,    23.6, 10.6,
    6.80,  5.38, 3.63, 2.90,  2.18, 0,     43.4, 22.8, 0,    0,    0,    11.60,
    9.40,  8.40, 7.50, 6.80,  6.10, 7.10,  8.12, 6.07, 4.31, 3.77, 3.05, 0,
    47.3,  27.6, 0,    0,     0,    12.80, 0,    0,    0,    0,    7.20, 7.20,
    10.20, 7.70, 6.60, 5.50,  5.35, 0,     0,    0,    0,    0,    0,    0,
    0,     0,    0,    23.50, 0,    0,     0,    0,    0,    0,    0,    0,
    0,     0,    0,    0,     0,    6.50,  5.80, 5.70, 7.60, 6.80, 7.40};
double ElectroNeg1[] = {
    2.59, 0,    0.89, 1.81, 2.28, 2.75, 3.19, 3.65, 4.0,  0,    0.56, 1.32,
    1.71, 2.14, 2.52, 2.96, 3.48, 0,    0.45, 0.95, 0,    0,    0,    1.66,
    2.2,  2.2,  2.56, 1.94, 1.98, 2.23, 2.42, 2.62, 2.82, 3.01, 3.22, 0,
    0.31, 0.72, 0,    0,    0,    1.15, 0,    0,    0,    0,    1.83, 1.98,
    2.14, 2.3,  2.46, 2.62, 2.78, 0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    2.0,  0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    2.28, 2.54, 2.2,  2.25, 2.29, 2.34};

double VdW1[] = {
    6.71,  0,     25.25, 0.0,   17.88, 22.45, 15.6,  11.49, 9.20,   0,     49.0,
    21.69, 36.51, 31.98, 26.52, 24.43, 22.45, 0,     87.11, 0.0,   0,     0,
    0,     44.6,  43.4,  41.05, 35.04, 17.16, 11.49, 11.25, 27.39, 28.73, 26.52,
    28.73, 31.06, 0,     0.0,   0.0,   0,     0,     0,     33.51, 0,     0,
    0,     0,     21.31, 16.52, 30.11, 45.83, 38.79, 36.62, 32.52, 0,     0,
    0,     0,     0,     0,     0,     0,     0,     0,     72.78, 0,     0,
    0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    22.45, 19.16, 15.6,  31.54, 34.53, 38.79};

double IonPolarizability1[] = {
  13.598, 0, 5.392, 9.323, 8.298, 11.26, 14.534, 13.618, 17.423, 0, 5.139, 7.646,
  5.986, 8.152, 10.487, 10.360, 12.968, 0, 4.341, 6.113, 0, 0, 0, 6.767, 7.434,
  7.902, 7.881, 7.640, 7.723, 9.394, 5.999, 7.900, 9.815, 9.752, 11.814, 0, 4.177,
  5.695, 0, 0, 0, 7.092, 0, 0, 0, 0, 7.576, 8.994, 5.786, 7.344, 8.640, 9.010,
  10.451, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6.15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 9.000, 9.226, 10.438, 6.108, 7.417, 7.289};



namespace RDKit {
namespace Descriptors {

namespace {

std::vector<double> getG(int n) {
  std::vector<double> res(n);
  for (int i = 0; i < n; i++) {
    res[i] = i;
  }
  return res;
}

std::vector<double> GetCharges(const ROMol &mol) {
  std::vector<double> charges(mol.getNumAtoms(), 0);
  // use 12 iterations... can be more
  computeGasteigerCharges(mol, charges, 12, true);
  return charges;
}

std::vector<double> GetRelativePol(const ROMol &mol) {
  int numAtoms = mol.getNumAtoms();

  std::vector<double> pol(numAtoms, 0);
  for (int i = 0; i < numAtoms; ++i) {
    pol[i] = Pol1[mol.getAtomWithIdx(i)->getAtomicNum()-1] / Pol1[5];
  }

  return pol;
}

std::vector<double> GetRelativeElectroNeg(const ROMol &mol) {
  int numAtoms = mol.getNumAtoms();

  std::vector<double> REN(numAtoms, 0);
  for (int i = 0; i < numAtoms; ++i) {
    REN[i] =
        ElectroNeg1[mol.getAtomWithIdx(i)->getAtomicNum()-1] / ElectroNeg1[5];
  }

  return REN;
}

std::vector<double> GetRelativeIonPol(const ROMol &mol) {
  int numAtoms = mol.getNumAtoms();

  std::vector<double> ionpol(numAtoms, 0);
  for (int i = 0; i < numAtoms; ++i) {
    ionpol[i] = IonPolarizability1[mol.getAtomWithIdx(i)->getAtomicNum()-1] / IonPolarizability1[5];
  }

  return ionpol;
}


std::vector<double> GetRelativeVdW(const ROMol &mol) {
  int numAtoms = mol.getNumAtoms();

  std::vector<double> vdw(numAtoms, 0);
  for (int i = 0; i < numAtoms; ++i) {
    vdw[i] = VdW1[mol.getAtomWithIdx(i)->getAtomicNum()-1] / VdW1[5];
  }

  return vdw;
}

std::vector<double> GetAbsPol(const ROMol &mol) {
  int numAtoms = mol.getNumAtoms();
  std::vector<double> pol(numAtoms, 0);
  for (int i = 0; i < numAtoms; ++i) {
    pol[i] = Pol1[mol.getAtomWithIdx(i)->getAtomicNum()-1];
  }

  return pol;
}



std::vector<double> CalcUnweightedMORSE(const ROMol &mol,
    const Conformer &conf) {
  int numAtoms = conf.getNumAtoms();
  int confId = conf.getId();

  std::vector<double> R = getG(32);
  std::vector<double> RDFres;
  double *DM = MolOps::get3DDistanceMat(mol, confId);

  for (int i = 0; i < 32; i++) {
    double res = 0;
    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {

        if (i==0) { res += 1;}

        else

          res += sin(R[i] * DM[j * numAtoms + k]) / (R[i] * DM[j * numAtoms + k]);
      }
    }

    RDFres.push_back(round( 1000 * res) / 1000);
  }

  return RDFres;
}

std::vector<double> CalcChargeMORSE(
    const ROMol &mol, const Conformer &conf) {
  int numAtoms = conf.getNumAtoms();
  int confId = conf.getId();

  std::vector<double> R = getG(32);
  std::vector<double> RDFres;
  double *DM = MolOps::get3DDistanceMat(mol, confId);
  std::vector<double> charges = GetCharges(mol);

  for (int i = 0; i < 32; i++) {
    double res = 0;
    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {
       
        if (i==0) { res += charges[j] * charges[k];}

        else

        res += charges[j] * charges[k] * sin(R[i] * DM[j * numAtoms + k]) / (R[i] * DM[j * numAtoms + k]);
      }
    }

    RDFres.push_back(round( 1000 * res) / 1000);
  }

  return RDFres;
}

std::vector<double> CalcMassMORSE(const ROMol &mol, const Conformer &conf) {
  int numAtoms = conf.getNumAtoms();
  int confId = conf.getId();

  std::vector<double> R = getG(32);
  std::vector<double> RDFres;
  double *DM = MolOps::get3DDistanceMat(mol, confId);
  std::vector<double> mass;
  for (int p = 0; p < numAtoms; p++) {
    mass.push_back(mol.getAtomWithIdx(p)->getMass());

  }

  for (int i = 0; i < 32; i++) {
    double res = 0;
    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {

       
        if (i==0) { res += mass[j] * mass[k];}

        else

        res += mass[j] * mass[k] * sin(R[i] * DM[j * numAtoms + k]) / (R[i] * DM[j * numAtoms + k]);
      }
    }


    // 144 = 12*12 mass of the Carbon in rdkit MC=12.011 caution

    RDFres.push_back(round(1000 * res / 12.011 / 12.011) / 1000);
  }



  return RDFres;
}

std::vector<double> CalcAtomNumMORSE(
    const ROMol &mol, const Conformer &conf) {
  int numAtoms = conf.getNumAtoms();
  int confId = conf.getId();

  std::vector<double> R = getG(32);
  std::vector<double> RDFres;

  double *DM = MolOps::get3DDistanceMat(mol, confId);

  std::vector<int> AN;
  for (int p = 0; p < numAtoms; p++) {
    AN.push_back(mol.getAtomWithIdx(p)->getAtomicNum());
  }
  /* debug
    for (int p = 0; p < numAtoms; p++) {
        std::cout << AN[p] << ",";
    }
    std::cout << "\n";
  */

  for (int i = 0; i < 32; i++) {
    double res = 0;
    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {
       
        if (i==0) { res += AN[j] * AN[k];}

        else

        res += AN[j] * AN[k] * sin(R[i] * DM[j * numAtoms + k]) / (R[i] * DM[j * numAtoms + k]);
      }
    }
    // 36 = 6*6 Atomic number of the Carbon
    
    RDFres.push_back(round( 1000 * res / 36) / 1000);
  }

  return RDFres;
}

std::vector<double> CalcIonPolMORSE(const ROMol &mol, const Conformer &conf) {
  int numAtoms = conf.getNumAtoms();
  int confId = conf.getId();

  std::vector<double> R = getG(32);
  std::vector<double> RDFres;
  double *DM = MolOps::get3DDistanceMat(mol, confId);

  std::vector<double> RelativeIonPol = GetRelativeIonPol(mol);

  for (int i = 0; i < 32; i++) {
    double res = 0;
    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {
     
        if (i==0) { res += RelativeIonPol[j] * RelativeIonPol[k];}

        else
        res += RelativeIonPol[j] * RelativeIonPol[k] * sin(R[i] * DM[j * numAtoms + k]) /
               (R[i] * DM[j * numAtoms + k]);
      }
    }

    RDFres.push_back(round( 1000 * res) / 1000);
  }

  return RDFres;
}

std::vector<double> CalcPolMORSE(const ROMol &mol, const Conformer &conf) {
  int numAtoms = conf.getNumAtoms();
  int confId = conf.getId();

  std::vector<double> R = getG(32);
  std::vector<double> RDFres;
  double *DM = MolOps::get3DDistanceMat(mol, confId);

  std::vector<double> RelativePol = GetRelativePol(mol);

  for (int i = 0; i < 32; i++) {
    double res = 0;
    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {
     
        if (i==0) { res += RelativePol[j] * RelativePol[k];}

        else
        res += RelativePol[j] * RelativePol[k] * sin(R[i] * DM[j * numAtoms + k]) /
               (R[i] * DM[j * numAtoms + k]);
      }
    }

    RDFres.push_back(round( 1000 * res) / 1000);
  }

  return RDFres;
}


std::vector<double> CalcElectroNegMORSE(
    const ROMol &mol, const Conformer &conf) {
  int numAtoms = conf.getNumAtoms();
  int confId = conf.getId();

  std::vector<double> R = getG(32);
  std::vector<double> RDFres;
  double *DM = MolOps::get3DDistanceMat(mol, confId);

  std::vector<double> RelativeElectroNeg = GetRelativeElectroNeg(mol);

  for (int i = 0; i < 32; i++) {
    double res = 0;
    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {
   
        if (i==0) { res += RelativeElectroNeg[j] * RelativeElectroNeg[k];}

        else
        res += RelativeElectroNeg[j] * RelativeElectroNeg[k] *
               sin(R[i] * DM[j * numAtoms + k]) / (R[i] * DM[j * numAtoms + k]);
      }
    }

    RDFres.push_back(round( 1000 * res) / 1000);
  }

  return RDFres;
}

std::vector<double> CalcVdWvolMORSE(
    const ROMol &mol, const Conformer &conf) {
  int numAtoms = conf.getNumAtoms();
  int confId = conf.getId();

  std::vector<double> R = getG(32);
  std::vector<double> RDFres;
  double *DM = MolOps::get3DDistanceMat(mol, confId);

  std::vector<double> RelativeVdW = GetRelativeVdW(mol);

  for (int i = 0; i < 32; i++) {
    double res = 0;
    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {

        if (i==0) { res += RelativeVdW[j] * RelativeVdW[k];}

        else

        res += RelativeVdW[j] * RelativeVdW[k] * sin(R[i] * DM[j * numAtoms + k]) /(R[i] * DM[j * numAtoms + k]);
      }
    }

    RDFres.push_back(round( 1000 * res) / 1000);
  }

  return RDFres;
}

}  // end of anonymous namespace

std::vector<double> MORSE(const ROMol &mol, int confId) {
  PRECONDITION(mol.getNumConformers() >= 1, "molecule has no conformers")

  const Conformer &conf = mol.getConformer(confId);

  std::vector<double> res1 = CalcUnweightedMORSE(mol,conf);

  std::vector<double> res2 = CalcMassMORSE(mol, conf);
  res1.insert(res1.end(),res2.begin(), res2.end());

  std::vector<double> res6 = CalcVdWvolMORSE(mol, conf);
  res1.insert(res1.end(),res6.begin(), res6.end());

  std::vector<double> res5 = CalcElectroNegMORSE(mol, conf);
  res1.insert(res1.end(),res5.begin(), res5.end());

  std::vector<double> res4 = CalcPolMORSE(mol, conf);
  res1.insert(res1.end(),res4.begin(), res4.end());

  std::vector<double> res7 = CalcIonPolMORSE(mol, conf);
  res1.insert(res1.end(),res7.begin(), res7.end());

 // std::vector<double> res3 = CalcChargeMORSE(mol, conf);
 // res1.insert(res1.end(),res3.begin(), res3.end());

 // std::vector<double> res7 = CalcAtomNumMORSE(mol, conf);
 // res1.insert(res1.end(),res7.begin(), res7.end());

   
  /* debug  */
    for (int p = 0; p < res1.size(); p++) {
      std::cout << res1[p] << ",";
      if ((p+1) % 32 == 0)     std::cout << "\n"; 
      
    }
      std::cout << "\n"; 
  /**/
  return res1;
}

}  // end of Descriptors namespace
}  // end of RDKit namespace
