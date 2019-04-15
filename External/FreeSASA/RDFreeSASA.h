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

#include <RDGeneral/export.h>
#include <GraphMol/RDKitBase.h>

namespace RDKit {
namespace common_properties {
namespace Atom {
RDKIT_FREESASALIB_EXPORT extern const std::string
    SASA;  // Solvent Accessible Surface Area for atom- double
RDKIT_FREESASALIB_EXPORT extern const std::string
    SASAClass;  // Class type, 0,1,2... etc
RDKIT_FREESASALIB_EXPORT extern const std::string
    SASAClassName;  // Class name, Polar, APolar etc...
}  // namespace Atom
namespace Molecule {
RDKIT_FREESASALIB_EXPORT extern const std::string
    SASA;  // Total Solvent Accessible Surface area for molecule;
}
}  // namespace common_properties
}  // namespace RDKit

namespace FreeSASA {
struct RDKIT_FREESASALIB_EXPORT SASAOpts {
  enum Algorithm { LeeRichards = 0, ShrakeRupley = 1 };
  enum Classifier { Protor = 0, NACCESS = 1, OONS = 2 };
  enum Classes { Unclassified = 0, APolar = 1, Polar = 2 };

  Algorithm algorithm;
  Classifier classifier;
  double probeRadius;
  SASAOpts();
  SASAOpts(Algorithm alg, Classifier cls);
  SASAOpts(Algorithm alg, Classifier cls, double pr);
};

//! Classify atoms using standard freesaa classifiers
/*!
  Note:

    FreeSASA identified Classes end up in
  atom.getProp<int>(common_properties::Atom::SASAClassName)
    FreeSASA Class names end up in
  atom.getProp<string>(common_properties::Atom::SASAClassName)

    \param mol:    Molecule to analyze
    \param radii   output vector of radii where radii[idx] is the radius for
  atom with index idx
    \return false if no atoms could be classified
*/
RDKIT_FREESASALIB_EXPORT bool classifyAtoms(
    RDKit::ROMol &mol, std::vector<double> &radii,
    const FreeSASA::SASAOpts &opts = SASAOpts());

//! calculate the Solvent Accessible Surface Area using the FreeSASA library.
/*!
  SASA atom contribution data is stored in
  atom.getProp(common_properites::Atom::SASA);

  \param mol:    Molecule to analyze
  \param radii   vector of radii where radii[idx] is the radius for atom with
  index idx
                 These can be passed in or calculated with classifyAtoms for
  some proteins.
  \param confIdx specify the conformation [default -1]
  \param query    query atom to limit the number of atoms to the ones matching
  the query
                  precanned query atoms can be made with
  makeFreeSasaPolarAtomQuery and
                  makeFreeSasaAPolarAtomQuery for classified polar and apolar
  atoms respectively.

  \param opts     SASAOpts class specifying options.
  \return the requested solvent accessible surface area
*/
RDKIT_FREESASALIB_EXPORT double calcSASA(const RDKit::ROMol &mol,
                                         const std::vector<double> &radii,
                                         int confIdx = -1,
                                         const RDKit::QueryAtom *query = NULL,
                                         const SASAOpts &opts = SASAOpts());

//! Make a query atom returning the FreeSASA supplied apolar atom classification
/*!
    These are atoms that have the "SASAClassName" property set to "Apolar"
    after calling classifyAtoms.

    \return QueryAtom pointer
*/
RDKIT_FREESASALIB_EXPORT const RDKit::QueryAtom *makeFreeSasaAPolarAtomQuery();
//! Make a query atom returning the FreeSASA supplied polar atom classification
/*!
    These are atoms that have the "SASAClassName" property set to "Polar"
    after calling classifyAtoms.
    \return QueryAtom pointer
*/
RDKIT_FREESASALIB_EXPORT const RDKit::QueryAtom *makeFreeSasaPolarAtomQuery();
}  // namespace FreeSASA

#endif
