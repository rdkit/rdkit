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

#ifndef RDKIT_FREESASA_H
#define RDKIT_FREESASA_H

#include <GraphMol/RDKitBase.h>

namespace RDKit {
  namespace common_properties {
    namespace Atom {
      extern const std::string SASA; // Solvent Accessible Surface Area for atom- double
      extern const std::string SASAClass;     // Class type, 0,1,2... etc
      extern const std::string SASAClassName; // Class name, Polar, APolar etc...
    }
    namespace Molecule {
      extern const std::string SASA; // Total Solvent Accessible Surface area for molecule;
    }
  }
}

namespace FreeSASA {
struct SASAOpts {
  enum Algorithm { LeeRichards = 0, ShrakeRupley = 1 };
  enum Classifier { Protor = 0, NACCESS = 1, OONS = 2 };
  enum Classes { Unclassified = 0, APolar = 1, Polar = 2 };

  Algorithm d_alg;
  Classifier d_classifier;
  Classes d_classes;
  SASAOpts() : d_alg(LeeRichards), d_classifier(Protor) {}
};

// Classify atoms using standard freesaa classifiers
//  Classes end up in atom.getProp(common_properties::Atom::FreeSASAClass)
//  returns false if no atoms could be classified
bool classifyAtoms(RDKit::ROMol &mol, std::vector<double> &radii,
                   const SASAOpts &opts = SASAOpts());

// calculate all atom contributions
//  data is stored in atom.getProp(common_properites::Atom::SASA);
double calcSASA(const RDKit::ROMol &mol, const std::vector<double> &radii,
                const SASAOpts &opts = SASAOpts());

// calculate all atom contributions
//  data is stored in atom.getProp(common_properites::Atom::SASA);
double calcSASA(const RDKit::ROMol &mol, const std::vector<double> &radii,
                int confIdx, const SASAOpts &opts = SASAOpts());

// calculate all atom contributions from the specified query
//  data is stored in atom.getProp(common_properites::Atom::SASA);

double calcSASA(const RDKit::ROMol &mol, const std::vector<double> &radii,
                const RDKit::QueryAtom *query,
                const SASAOpts &opts = SASAOpts());

double calcSASA(const RDKit::ROMol &mol, const std::vector<double> &radii,
                const RDKit::QueryAtom *query,
                int confIdx,
                const SASAOpts &opts = SASAOpts());

const RDKit::QueryAtom * makeAPolarAtomQuery();
const RDKit::QueryAtom * makePolarAtomQuery();

}

#endif
