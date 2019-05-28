//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file Normalize.h

        \brief Defines the Normalizer class.

*/
#include <RDGeneral/export.h>
#ifndef __RD_NORMALIZE_H__
#define __RD_NORMALIZE_H__

#include <Catalogs/Catalog.h>
#include <GraphMol/MolStandardize/TransformCatalog/TransformCatalogEntry.h>
#include <GraphMol/MolStandardize/TransformCatalog/TransformCatalogParams.h>
#include <GraphMol/MolStandardize/MolStandardize.h>

namespace RDKit {
class RWMol;
class ROMol;

namespace MolStandardize {
RDKIT_MOLSTANDARDIZE_EXPORT extern const CleanupParameters
    defaultCleanupParameters;

typedef RDCatalog::HierarchCatalog<TransformCatalogEntry,
                                   TransformCatalogParams, int>
    TransformCatalog;

//! The Normalizer class for applying Normalization transforms.
/*!

  <b>Notes:</b>
    - This class is typically used to apply a series of Normalization transforms
  to correct functional groups and recombine charges.
                - Each transform is repeatedly applied until no further changes
  occur.
*/

class RDKIT_MOLSTANDARDIZE_EXPORT Normalizer {
 public:
  Normalizer();
  //! Construct a Normalizer with a particular normalizeFile and maxRestarts
  Normalizer(const std::string normalizeFile, const unsigned int maxRestarts);
  //! Construct a Normalizer with a particular stream (with parameters) and
  //! maxRestarts
  Normalizer(std::istream &normalizeStream, const unsigned int maxRestarts);
  //! making Normalizer objects non-copyable
  Normalizer(const Normalizer &other) = delete;
  Normalizer &operator=(Normalizer const &) = delete;
  ~Normalizer();

  //! Apply a series of Normalization transforms to correct functional groups
  //! and recombine charges.
  /*!
    <b>Notes:</b>
      - A series of transforms are applied to the molecule. For each
    Normalization, the transform is applied repeatedly until no further changes
    occur.
      - If any changes occurred, we go back and start from the first
    Normalization again, in case the changes mean an earlier transform is now
    applicable.
                        - The molecule is returned once the entire series of
    Normalizations cause no further changes or if max_restarts (default 200) is
    reached.
  */
  ROMol *normalize(const ROMol &mol);
  struct Product {
    std::string Smiles;
    boost::shared_ptr<ROMol> Mol;
    Product(std::string smiles, boost::shared_ptr<ROMol> &mol)
        : Smiles(smiles), Mol(mol) {}

    // sorting products alphabetically by SMILES
    bool operator<(const Product &pdt) const { return (Smiles < pdt.Smiles); }
  };

 private:
  TransformCatalog *d_tcat;
  unsigned int MAX_RESTARTS;

  ROMol *normalizeFragment(
      const ROMol &mol,
      const std::vector<std::shared_ptr<ChemicalReaction>> &transforms);
  boost::shared_ptr<ROMol> applyTransform(const ROMol &mol,
                                          ChemicalReaction &rule);

};  // Normalizer class
}  // namespace MolStandardize
}  // namespace RDKit

#endif
