//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef __RD_FRAGMENT_REMOVER_H__
#define __RD_FRAGMENT_REMOVER_H__

#include <Catalogs/Catalog.h>
#include <GraphMol/MolStandardize/FragmentCatalog/FragmentCatalogEntry.h>
#include <GraphMol/MolStandardize/FragmentCatalog/FragmentCatalogParams.h>
#include <GraphMol/MolStandardize/MolStandardize.h>

namespace RDKit {
class ROMol;

namespace MolStandardize {

RDKIT_MOLSTANDARDIZE_EXPORT extern const CleanupParameters
    defaultCleanupParameters;

typedef RDCatalog::HierarchCatalog<FragmentCatalogEntry, FragmentCatalogParams,
                                   int>
    FragmentCatalog;

class RDKIT_MOLSTANDARDIZE_EXPORT FragmentRemover {
 public:
  FragmentRemover();
  FragmentRemover(const std::string fragmentFile, const bool leave_last);
  //  FragmentRemover(bool leave_last) : LEAVE_LAST(leave_last){};
  //! making FragmentRemover objects non-copyable
  FragmentRemover(const FragmentRemover &other) = delete;
  FragmentRemover &operator=(FragmentRemover const &) = delete;
  ~FragmentRemover();

  ROMol *remove(const ROMol &mol);

 private:
  // Setting leave_last to True will ensure at least one fragment
  //  is left in the molecule, even if it is matched by a
  //  FragmentPattern
  bool LEAVE_LAST;
  FragmentCatalog *d_fcat;

};  // class FragmentRemover

class RDKIT_MOLSTANDARDIZE_EXPORT LargestFragmentChooser {
 public:
  //  LargestFragmentChooser(){};
  LargestFragmentChooser(bool prefer_organic = false)
      : PREFER_ORGANIC(prefer_organic){};
  LargestFragmentChooser(const LargestFragmentChooser &other);
  ~LargestFragmentChooser(){};

  ROMol *choose(const ROMol &mol);
  struct Largest {
    Largest();
    Largest(std::string &smiles, const boost::shared_ptr<ROMol> &fragment,
            unsigned int &numatoms, double &weight, bool &organic);
    std::string Smiles;
    boost::shared_ptr<ROMol> Fragment;
    unsigned int NumAtoms;
    double Weight;
    bool Organic;
  };

 private:
  bool PREFER_ORGANIC;
};  // class LargestFragmentChooser
}  // namespace MolStandardize
}  // namespace RDKit

#endif
