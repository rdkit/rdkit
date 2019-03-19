
#include <RDGeneral/export.h>

#ifndef RDKIT_EDITABLESUBSTRUCTLIBRARY_TSP_H
#define RDKIT_EDITABLESUBSTRUCTLIBRARY_TSP_H

#include "EditableSubstructLibrary.h"

namespace RDKit {

class RDKIT_SUBSTRUCTLIBRARY_EXPORT EditableSubstructLibraryTrustedSmilesWithPattern
    : public EditableSubstructLibrary<std::string, CachedTrustedSmilesMolHolder,
                                      PatternHolder> {
 public:
  EditableSubstructLibraryTrustedSmilesWithPattern()
      : EditableSubstructLibrary() {}

  int addSmiles(const std::string &smi, const std::string &fp);

  int addSmiles(const std::vector<std::string> smilesVector,
                const std::vector<std::string> fpVector);

  ExplicitBitVect * makeFingerprint(const std::string &smiles) const;

  std::string makeStringFingerprint(const std::string &smiles) const;  

  };

}  // namespace RDKit

#endif
