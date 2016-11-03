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
#include "MolData3Ddescriptors.h"

#include <Numerics/Matrix.h>
#include <Numerics/SquareMatrix.h>
#include <Numerics/SymmMatrix.h>
#include <boost/foreach.hpp>
#include <math.h>
#include <Eigen/Dense>

// data checked using book Todeschini R., Consonni V. - Molecular Descriptors for Chemoinformatics 2009 atomic properties page 21/22

namespace RDKit {
namespace Descriptors {

namespace {


  MolData3Ddescriptors moldata3D;

std::vector<double> getG(int n) {
  std::vector<double> res(n);
  for (int i = 0; i < n; i++) {
    res[i] = i;
  }
  return res;
}

/*
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
*/

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


std::vector<double> CalcMassMORSE(const ROMol &mol, const Conformer &conf) {
  int numAtoms = conf.getNumAtoms();
  int confId = conf.getId();

  std::vector<double> R = getG(32);
  std::vector<double> RDFres;
  double *DM = MolOps::get3DDistanceMat(mol, confId);

  
  std::vector<double> mass = moldata3D.GetRelativeMW(mol);

  for (int i = 0; i < 32; i++) {
    double res = 0;
    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {

       
        if (i==0) { res += mass[j] * mass[k];}

        else

        res += mass[j] * mass[k] * sin(R[i] * DM[j * numAtoms + k]) / (R[i] * DM[j * numAtoms + k]);
      }
    }

    RDFres.push_back(round(1000 * res ) / 1000);
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

  std::vector<double> RelativeIonPol = moldata3D.GetRelativeIonPol(mol);

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

  std::vector<double> RelativePol = moldata3D.GetRelativePol(mol);

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

  std::vector<double> RelativeElectroNeg = moldata3D.GetRelativeENeg(mol);

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

  std::vector<double> RelativeVdW = moldata3D.GetRelativeVdW(mol);


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

   
  return res1;
}

}  // end of Descriptors namespace
}  // end of RDKit namespace
