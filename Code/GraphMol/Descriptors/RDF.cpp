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
// Adding RBF descriptors to 3D descriptors by Guillaume Godin

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>

#include <GraphMol/MolTransforms/MolTransforms.h>

#include "RDF.h"
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
    res[i] = 1 + i * 0.5;
  }
  return res;
}

/*
std::vector<double> CalcChargeRDF(const ROMol &mol, const Conformer &conf) {
  int numAtoms = conf.getNumAtoms();
  int confId = conf.getId();

  std::vector<double> R = getG(30);
  std::vector<double> RDFres;
  //std::vector<std::vector<double> > DM = GetGeometricalDistanceMatrix(points);
  double *DM = MolOps::get3DDistanceMat(mol,confId);

  std::vector<double> charges = GetCharges(mol);

  for (int i = 0; i < 30; i++) {
    double res = 0;
    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {
        res += charges[j] * charges[k] * exp(-100 * pow(R[i] - DM[j * numAtoms + k], 2));
      }
    }

      RDFres.push_back(round( 1000 * res) / 1000);
  }

  return RDFres;

*/


std::vector<double> CalcUnweightedRDF(const ROMol &mol, const Conformer &conf) {
  int numAtoms = conf.getNumAtoms();
  int confId = conf.getId();

  std::vector<double> R = getG(30);

  std::vector<double> RDFres;

  double *DM = MolOps::get3DDistanceMat(mol,confId);

  for (int i = 0; i < 30; i++) {
    double res = 0;
    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {
        res += exp(-100 * pow(R[i] - DM[j * numAtoms + k], 2));
          
      }
    }

    RDFres.push_back(round( 1000 * res) / 1000);
  }

  return RDFres;
}



std::vector<double> CalcMassRDF(const ROMol &mol, const Conformer &conf) {
  int numAtoms = conf.getNumAtoms();
  int confId = conf.getId();

  std::vector<double> R = getG(30);
  std::vector<double> RDFres;

  double *DM = MolOps::get3DDistanceMat(mol,confId);


  std::vector<double> Mass = moldata3D.GetRelativeMW(mol);



  for (int i = 0; i < 30; i++) {
    double res = 0;
    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {
        res += Mass[j] * Mass[k] * exp(-100 * pow(R[i] - DM[j * numAtoms + k], 2));
      }
    }
      RDFres.push_back(round( 1000 * res) / 1000);
  }

  return RDFres;
}

std::vector<double> CalcPolRDF(const ROMol &mol, const Conformer &conf) {
  int numAtoms = conf.getNumAtoms();
  int confId = conf.getId();

  std::vector<double> R = getG(30);
  std::vector<double> RDFres;

  double *DM = MolOps::get3DDistanceMat(mol,confId);

  std::vector<double> RelativePol = moldata3D.GetRelativePol(mol);

  for (int i = 0; i < 30; i++) {
    double res = 0;
    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {
        res += RelativePol[j] * RelativePol[k] *
               exp(-100 * pow(R[i] - DM[j * numAtoms + k], 2));
      }
    }

      RDFres.push_back(round( 1000 * res) / 1000);
  }

  return RDFres;
}

std::vector<double> CalcIonPolRDF(const ROMol &mol, const Conformer &conf) {
  int numAtoms = conf.getNumAtoms();
  int confId = conf.getId();

  std::vector<double> R = getG(30);
  std::vector<double> RDFres;

  double *DM = MolOps::get3DDistanceMat(mol,confId);

  std::vector<double> IonPol = moldata3D.GetRelativeIonPol(mol);

  for (int i = 0; i < 30; i++) {
    double res = 0;
    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {
        res += IonPol[j] * IonPol[k] *
               exp(-100 * pow(R[i] - DM[j * numAtoms + k], 2));
      }
    }

      RDFres.push_back(round( 1000 * res) / 1000);
  }

  return RDFres;
}


std::vector<double> CalcIStateRDF(const ROMol &mol, const Conformer &conf) {
  int numAtoms = conf.getNumAtoms();
  int confId = conf.getId();

  std::vector<double> R = getG(30);
  std::vector<double> RDFres;
  //std::vector<std::vector<double> > DM = GetGeometricalDistanceMatrix(points);

  MolOps::removeHs(mol, false, false);

  double *DM = MolOps::get3DDistanceMat(mol,confId);

  std::vector<double> IState = moldata3D.GetEState2(mol);

  for (int i = 0; i < 30; i++) {
    double res = 0;
    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {
        res += IState[j] * IState[k] *
               exp(-100 * pow(R[i] - DM[j * numAtoms + k], 2));
      }
    }
      RDFres.push_back(round( 1000 * res) / 1000);
  }

  return RDFres;
}


std::vector<double> CalcElectroNegRDF(const ROMol &mol, const Conformer &conf) {
  int numAtoms = conf.getNumAtoms();
  int confId = conf.getId();

  std::vector<double> R = getG(30);
  std::vector<double> RDFres;

  double *DM = MolOps::get3DDistanceMat(mol,confId);

  std::vector<double> RelativeElectroNeg = moldata3D.GetRelativeENeg(mol);

  for (int i = 0; i < 30; i++) {
    double res = 0;
    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {
        res += RelativeElectroNeg[j] * RelativeElectroNeg[k] *
               exp(-100 * pow(R[i] - DM[j * numAtoms + k], 2));
      }
    }

      RDFres.push_back(round( 1000 * res) / 1000);
  }

  return RDFres;
}

std::vector<double> CalcVdWvolRDF(const ROMol &mol, const Conformer &conf) {
  int numAtoms = conf.getNumAtoms();
  int confId = conf.getId();

  std::vector<double> R = getG(30);
  std::vector<double> RDFres;
  double *DM = MolOps::get3DDistanceMat(mol, confId);

  std::vector<double> RelativeVdW = moldata3D.GetRelativeVdW(mol);

  for (int i = 0; i < 30; i++) {
    double res = 0;
    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {
        res += RelativeVdW[j] * RelativeVdW[k] *
               exp(-100 * pow(R[i] - DM[j * numAtoms + k], 2));
      }
    }

      RDFres.push_back(round( 1000 * res) / 1000);
  }

  return RDFres;
}

}  // end of anonymous namespace

std::vector<double> RDF(const ROMol &mol, int confId) {
  std::vector<double> reserror(std::vector<double>(30, 0));

// RDF010u  RDF015u RDF020u RDF025u RDF030u RDF035u RDF040u RDF045u RDF050u RDF055u RDF060u RDF065u RDF070u RDF075u RDF080u RDF085u RDF090u RDF095u RDF100u RDF105u RDF110u RDF115u RDF120u RDF125u RDF130u RDF135u RDF140u RDF145u RDF150u RDF155u 
// RDF010m RDF015m RDF020m RDF025m RDF030m RDF035m RDF040m RDF045m RDF050m RDF055m RDF060m RDF065m RDF070m RDF075m RDF080m RDF085m RDF090m RDF095m RDF100m RDF105m RDF110m RDF115m RDF120m RDF125m RDF130m RDF135m RDF140m RDF145m RDF150m RDF155m 
// RDF010v RDF015v RDF020v RDF025v RDF030v RDF035v RDF040v RDF045v RDF050v RDF055v RDF060v RDF065v RDF070v RDF075v RDF080v RDF085v RDF090v RDF095v RDF100v RDF105v RDF110v RDF115v RDF120v RDF125v RDF130v RDF135v RDF140v RDF145v RDF150v RDF155v 
// RDF010e RDF015e RDF020e RDF025e RDF030e RDF035e RDF040e RDF045e RDF050e RDF055e RDF060e RDF065e RDF070e RDF075e RDF080e RDF085e RDF090e RDF095e RDF100e RDF105e RDF110e RDF115e RDF120e RDF125e RDF130e RDF135e RDF140e RDF145e RDF150e RDF155e 
// RDF010p RDF015p RDF020p RDF025p RDF030p RDF035p RDF040p RDF045p RDF050p RDF055p RDF060p RDF065p RDF070p RDF075p RDF080p RDF085p RDF090p RDF095p RDF100p RDF105p RDF110p RDF115p RDF120p RDF125p RDF130p RDF135p RDF140p RDF145p RDF150p RDF155p 
// RDF010i RDF015i RDF020i RDF025i RDF030i RDF035i RDF040i RDF045i RDF050i RDF055i RDF060i RDF065i RDF070i RDF075i RDF080i RDF085i RDF090i RDF095i RDF100i RDF105i RDF110i RDF115i RDF120i RDF125i RDF130i RDF135i RDF140i RDF145i RDF150i RDF155i 
// RDF010s RDF015s RDF020s RDF025s RDF030s RDF035s RDF040s RDF045s RDF050s RDF055s RDF060s RDF065s RDF070s RDF075s RDF080s RDF085s RDF090s RDF095s RDF100s RDF105s RDF110s RDF115s RDF120s RDF125s RDF130s RDF135s RDF140s RDF145s RDF150s RDF155s


  PRECONDITION(mol.getNumConformers() >= 1, "molecule has no conformers")
  int numAtoms = mol.getNumAtoms();
  if (numAtoms < 4) return reserror;

  const Conformer &conf = mol.getConformer(confId);
  if (!conf.is3D()) return reserror;

  std::vector<double> res1 = CalcUnweightedRDF(mol,conf);

  std::vector<double> res2=CalcMassRDF(mol,conf);
  res1.insert(res1.end(),res2.begin(), res2.end());

  std::vector<double> res6=CalcVdWvolRDF(mol,conf);
  res1.insert(res1.end(),res6.begin(), res6.end());

  std::vector<double> res5=CalcElectroNegRDF(mol,conf);
  res1.insert(res1.end(),res5.begin(), res5.end());

  std::vector<double> res4=CalcPolRDF(mol,conf);
  res1.insert(res1.end(),res4.begin(), res4.end());


  std::vector<double> res7=CalcIonPolRDF(mol,conf);
  res1.insert(res1.end(),res7.begin(), res7.end());

  std::vector<double> res8=CalcIStateRDF(mol,conf);
  res1.insert(res1.end(),res8.begin(), res8.end());


 // std::vector<double> res3=CalcChargeRDF(mol,conf);
 // res1.insert(res1.end(),res3.begin(), res3.end());


  return res1;
}

}  // end of Descriptors namespace
}  // end of RDKit namespace
