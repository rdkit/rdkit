//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file MolStandardize.h

	\brief Defines the CleanupParameters and some convenience functions.

*/
#ifndef __RD_MOLSTANDARDIZE_H__
#define __RD_MOLSTANDARDIZE_H__

#include <string>
#include <GraphMol/RDKitBase.h>

namespace RDKit {
class RWMol;
class ROMol;

namespace MolStandardize {

//! The CleanupParameters structure defines the default parameters for the
// cleanup process and also allows the user to customize the process by changing
// the parameters.
/*!

  <b>Notes:</b>
    - To customize the parameters, the stucture must be initialized first. 
		  (Another on the TODO list)
		- For this project, not all the parameters have been revealed. (TODO)

*/
struct CleanupParameters {
  // TODO reveal all parameters
  std::string rdbase = std::getenv("RDBASE");
  std::string normalizations;
  std::string acidbaseFile;
  std::string fragmentFile;
  // std::vector<std::string> chargeCorrections;
  std::string tautomerTransforms;
  // std::vector<std::string> TautomerScores;
  int maxRestarts;     // The maximum number of times to attempt to apply the
                       // series of normalizations (default 200).
  int maxTautomers;    // The maximum number of tautomers to enumerate (default
                       // 1000).
  bool preferOrganic;  // Whether to prioritize organic fragments when choosing
                       // fragment parent (default False).

  CleanupParameters()
      :  // TODO
         //			normalizations(""),//this->DEFAULT_TRANSFORMS),
        normalizations(rdbase + "/Data/MolStandardize/normalizations.txt"),
        acidbaseFile(rdbase + "/Data/MolStandardize/acid_base_pairs.txt"),
        fragmentFile(rdbase + "/Data/MolStandardize/fragmentPatterns.txt"),
        // chargeCorrections()
        tautomerTransforms(rdbase +
                           "/Data/MolStandardize/tautomerTransforms.in"),
        // TautomerScores()
        maxRestarts(200),
        maxTautomers(1000),
        preferOrganic(false) {}
};

extern const CleanupParameters defaultCleanupParameters;

//! The cleanup function is equivalent to the
// molvs.Standardizer().standardize(mol) function. It calls the same steps,
// namely: RemoveHs, RDKit SanitizeMol, MetalDisconnector, Normalizer, Reionizer, 
// RDKit AssignStereochemistry.
RWMol *cleanup(const RWMol &mol,
               const CleanupParameters &params = defaultCleanupParameters);

//! TODO not yet finished!
void tautomerParent(RWMol &mol,
                    const CleanupParameters &params = defaultCleanupParameters);

//! Returns the fragment parent of a given molecule. The fragment parent is the
// largest organic covalent unit in the molecule.
RWMol *fragmentParent(
    const RWMol &mol,
    const CleanupParameters &params = defaultCleanupParameters,
    bool skip_standardize = false);

// TODO
void stereoParent(RWMol &mol,
                  const CleanupParameters &params = defaultCleanupParameters);

// TODO
void isotopeParent(RWMol &mol,
                   const CleanupParameters &params = defaultCleanupParameters);

//! Returns the charge parent of a given molecule. The charge parent is the uncharged 
// version of the fragment parent.
RWMol *chargeParent(const RWMol &mol,
                    const CleanupParameters &params = defaultCleanupParameters,
                    bool skip_standardize = false);

// TODO Need to do tautomers first
void superParent(RWMol &mol,
                 const CleanupParameters &params = defaultCleanupParameters);

//! Works the same as Normalizer().normalize(mol)
RWMol *normalize(const RWMol *mol,
                 const CleanupParameters &params = defaultCleanupParameters);

//! Works the same as Reionizer().reionize(mol)
RWMol *reionize(const RWMol *mol,
                const CleanupParameters &params = defaultCleanupParameters);

//! Convenience function for quickly standardizing a single SMILES string.
// Returns a standardized canonical SMILES string given a SMILES string.
std::string standardizeSmiles(const std::string &smiles);

//! TODO
std::vector<std::string> enumerateTautomerSmiles(
    const std::string &smiles,
    const CleanupParameters &params = defaultCleanupParameters);
};  // namespace MolStandardize
}  // namespace RDKit
#endif
