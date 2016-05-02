//
//  Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
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
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
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
#include "Property.h"
#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/Descriptors/MolSurf.h>

namespace RDKit {
namespace Descriptors {

const char * Property::getNameForProp(PropType prop) {
  switch (prop) {
    case Property::MW: return "Molecular Weight";
    case Property::TPSA: return "Polar Surface Area";
    case Property::ALOGP: return "A logP";
    case Property::NumRotors: return "Number of Rotatable bonds";
    case Property::MissingStereo: return "Number of undefined stereo centers";
    case Property::LipinskiHBA: return "Number of Lipinski Hydrogen Bond Acceptors";
    case Property::LipinskiHBD: return "Number of Lipinski Hydrogen Bond Donors";
    case Property::HBA: return "Number of Hydrogen Bond Acceptors";
    case Property::HBD: return "Number of Hydrogen Bond Donors";
    case Property::NumHeteroAtoms: return "Number of Hetero Atoms";
    case Property::NumAmideBonds: return "Number of Amide Bonds";
    case Property::FractionCSP3: return "Fraction of SP3 carbons";
    case Property::NumRings: return "Number of Rings";
    case Property::NumAromaticRings: return "Number of Aromatic Rings";
    case Property::NumAliphaticRings: return "Number of Aliphatic Rings";
    case Property::NumSaturatedRings: return "Number of Saturated Rings";
    case Property::NumHeteroCycles: return "Number of Heterocycles";
    case Property::NumSaturatedHeteroCycles: return "Number of Saturated Heterocycles";
    case Property::NumAromaticHeteroCycles: return "Number of Aromatic HeteroCycles";
    case Property::NumAliphaticCarboCycles: return "Number of Aliphatic CarboCycles";
    case Property::NumAromaticCarboCycles: return "Number of Aromatic CarboCycles";
    case Property::NumSaturatedCarboCycles: return "Number of Saturated CarboCycles";
    case Property::NumAlpiphaticHeteroCycles: return "Number of Aliphatic Hetero Cycles";
    case Property::NumSpiroAtoms: return "Number of Spiro Atoms";
    case Property::NumBridgeheadAtoms: return "Number of Bridgehead Atoms";
    case Property::AMW: return "Average Molecular WeightW";
    case Property::LabuteASA: return "Labute Approximate Surface Area";
      
    default:
        return "User defined property";
  }
}

std::string Property::getVersionForProp(PropType prop) {
  switch (prop) {
    case Property::MW: return  exactmwVersion;
    case Property::TPSA: return tpsaVersion;
    case Property::ALOGP: return crippenVersion;
    case Property::NumRotors: return NumRotatableBondsVersion;
    case Property::MissingStereo: return "1.0";
    case Property::LipinskiHBA: return "1.0";
    case Property::LipinskiHBD: return "1.0";
    case Property::HBA: return NumHBAVersion;
    case Property::HBD: return NumHBDVersion;
    case Property::NumHeteroAtoms: return NumHeteroatomsVersion;
    case Property::NumAmideBonds: return NumAmideBondsVersion;
    case Property::FractionCSP3: return FractionCSP3Version;
    case Property::NumRings: return NumRingsVersion;
    case Property::NumAromaticRings: return NumAromaticRingsVersion;
    case Property::NumAliphaticRings: return NumAliphaticRingsVersion;
    case Property::NumSaturatedRings: return NumSaturatedRingsVersion;
    case Property::NumHeteroCycles: return NumHeterocyclesVersion;
    case Property::NumSaturatedHeteroCycles: return NumSaturatedHeterocyclesVersion;
    case Property::NumAromaticHeteroCycles: return NumAromaticHeterocyclesVersion;
    case Property::NumAromaticCarboCycles: return NumAromaticCarbocyclesVersion;
    case Property::NumSaturatedCarboCycles: return NumSaturatedCarbocyclesVersion;
    case Property::NumAlpiphaticHeteroCycles: return NumAliphaticCarbocyclesVersion;
    case Property::NumSpiroAtoms: return NumSpiroAtomsVersion;
    case Property::NumBridgeheadAtoms: return NumBridgeheadAtomsVersion;
    case Property::AMW: return amwVersion;
    case Property::LabuteASA: return labuteASAVersion;
      
    default:
        return "User defined property";
  }
}

bool Property::isAdditive() const {
  switch(m_proptype) {
    case Property::AMW: 
    case Property::LabuteASA:
      return false;
    case Property::User:
      PRECONDITION(m_propfxn.get(), "Null Property fxn");
      return m_propfxn->isAdditive();
    default:
      return true;
  }
}


namespace {
double numUnspecifiedStereoAtoms(const ROMol &mol) {
  double res=0.;
  for (ROMol::ConstAtomIterator atom = mol.beginAtoms(); atom != mol.endAtoms();
       ++atom) {
    if ((*atom)->hasProp(common_properties::_ChiralityPossible) &&
        (*atom)->getChiralTag() == Atom::CHI_UNSPECIFIED)
      res++;
  }
  return res;
}
}

double Property::computeProperty(const ROMol &mol) const {
  const bool ignoreCachedValues = true;
  const bool includeHs = true; // used for default values...
  double logp;
  double mr;
  
  switch (m_proptype) {
    case Property::MW: return calcExactMW(mol);
    case Property::TPSA: return calcTPSA(mol, ignoreCachedValues);
    case Property::ALOGP:
      calcCrippenDescriptors(mol, logp, mr, includeHs, ignoreCachedValues);
      return logp;
    case Property::NumRotors: return calcNumRotatableBonds(mol);
    case Property::MissingStereo: return numUnspecifiedStereoAtoms(mol);
    case Property::LipinskiHBA: return calcLipinskiHBA(mol);
    case Property::LipinskiHBD: return calcLipinskiHBD(mol);
    case Property::HBA: return calcNumHBA(mol);
    case Property::HBD: return calcNumHBD(mol);
    case Property::NumHeteroAtoms: return calcNumHeteroatoms(mol); 
    case Property::NumAmideBonds: return calcNumAmideBonds(mol);
    case Property::FractionCSP3: return calcFractionCSP3(mol);
    case Property::NumRings: return calcNumRings(mol);
    case Property::NumAromaticRings: return calcNumAromaticRings(mol);
    case Property::NumAliphaticRings: return calcNumAliphaticRings(mol);
    case Property::NumSaturatedRings: return calcNumSaturatedRings(mol);
    case Property::NumHeteroCycles: return calcNumHeterocycles(mol);
    case Property::NumSaturatedHeteroCycles: return calcNumSaturatedHeterocycles(mol);
    case Property::NumAromaticHeteroCycles: return calcNumAromaticHeterocycles(mol);
    case Property::NumAromaticCarboCycles:return calcNumAromaticCarbocycles(mol);
    case Property::NumAliphaticCarboCycles: return calcNumAliphaticCarbocycles(mol);
    case Property::NumSaturatedCarboCycles: return calcNumSaturatedCarbocycles(mol);
    case Property::NumAlpiphaticHeteroCycles: return calcNumAliphaticHeterocycles(mol);
    case Property::NumSpiroAtoms: return calcNumSpiroAtoms(mol);
    case Property::NumBridgeheadAtoms: return calcNumBridgeheadAtoms(mol);
    case Property::AMW: return calcAMW(mol);
    case Property::LabuteASA: return calcLabuteASA(mol);
      
    default:
      return m_propfxn.get()->compute(mol);
  }
  
}

}
}
