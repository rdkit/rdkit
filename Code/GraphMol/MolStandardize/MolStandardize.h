//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_MOLSTANDARDIZE_H__
#define __RD_MOLSTANDARDIZE_H__

#include <string>
#include <GraphMol/RDKitBase.h>

namespace RDKit {
class RWMol;
class ROMol;

namespace MolStandardize {

struct CleanupParameters {
  // TODO
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

RWMol *cleanup(const RWMol &mol,
               const CleanupParameters &params = defaultCleanupParameters);

void tautomerParent(RWMol &mol,
                    const CleanupParameters &params = defaultCleanupParameters);

RWMol *fragmentParent(
    const RWMol &mol,
    const CleanupParameters &params = defaultCleanupParameters,
    bool skip_standardize = false);

void stereoParent(RWMol &mol,
                  const CleanupParameters &params = defaultCleanupParameters);

void isotopeParent(RWMol &mol,
                   const CleanupParameters &params = defaultCleanupParameters);

RWMol *chargeParent(const RWMol &mol,
                    const CleanupParameters &params = defaultCleanupParameters,
                    bool skip_standardize = false);

void superParent(RWMol &mol,
                 const CleanupParameters &params = defaultCleanupParameters);

RWMol *normalize(const RWMol *mol,
                 const CleanupParameters &params = defaultCleanupParameters);

RWMol *reionize(const RWMol *mol,
                const CleanupParameters &params = defaultCleanupParameters);

std::string standardizeSmiles(const std::string &smiles);

std::vector<std::string> enumerateTautomerSmiles(
    const std::string &smiles,
    const CleanupParameters &params = defaultCleanupParameters);
};  // namespace MolStandardize
}  // namespace RDKit
#endif
