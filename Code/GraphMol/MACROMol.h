//
//  Copyright (C) 2024 Tad Hurst and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_MACROMOL_H
#define RD_MACROMOL_H

#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include "FileParsers/FileParsers.h"
#include "FileParsers/FileParserUtils.h"

namespace RDKit {

class RDKIT_GRAPHMOL_EXPORT MACROMolTemplate {
 private:
  void init(std::unique_ptr<RDKit::ROMol> &mol, std::string className,
            std::vector<std::string> templateNames,
            std::vector<std::string> templateAttrs);

 public:
  MACROMolTemplate(std::unique_ptr<ROMol> &mol, std::string className,
                   std::vector<std::string> templateNames,
                   std::vector<std::string> templateAttrs);
  MACROMolTemplate(std::unique_ptr<ROMol> &mol, std::string className,
                   std::string templateName,
                   std::vector<std::string> templateAttrs);
  MACROMolTemplate() = delete;
  MACROMolTemplate(const MACROMolTemplate &other) = delete;
  MACROMolTemplate(MACROMolTemplate &&other) noexcept = delete;
  MACROMolTemplate &operator=(MACROMolTemplate &&other) noexcept = delete;
  MACROMolTemplate &operator=(const MACROMolTemplate &) =
      delete;  // disable assignment
  ~MACROMolTemplate() {}

  RDKit::ROMol *getMol() { return p_mol.get(); }

 private:
  std::unique_ptr<RDKit::ROMol> p_mol;
};

class RDKIT_GRAPHMOL_EXPORT MACROMol : public RWMol {
 private:
  // for libraries owned by this macromol
  std::vector<std::unique_ptr<MACROMolTemplate>> p_templateLibrary;

  // used for library access.  Could be a pointer to the internal library or
  // an external one that is shared across multiple macromols

  std::vector<std::unique_ptr<MACROMolTemplate>> *p_templateLibraryPtr;
  bool p_atomIdxToTemplateIdxIsStale;
  std::map<unsigned int, unsigned int> p_atomIdxToTemplateIdx;

 public:
  MACROMol() : p_atomIdxToTemplateIdxIsStale(true) {};
  MACROMol(const MACROMol &other) = delete;
  MACROMol(MACROMol &&other) noexcept = delete;
  MACROMol &operator=(MACROMol &&other) noexcept = delete;

  MACROMol &operator=(const MACROMol &) = delete;  // disable assignment
  ~MACROMol() {}

  void setTemplateLibrary(std::vector<std::unique_ptr<MACROMolTemplate>> &lib,
                          bool takeOwnership = true) {
    if (takeOwnership) {
      for (auto &libTemplate : lib) {
        p_templateLibrary.push_back(std::move(libTemplate));
      }
      p_templateLibraryPtr = &p_templateLibrary;
    }

    else {
      p_templateLibraryPtr = &lib;
    }
  }

  void addTemplate(std::unique_ptr<MACROMolTemplate> templateMol) {
    PRECONDITION(templateMol, "bad template molecule");
    p_templateLibrary.push_back(std::move(templateMol));
    p_templateLibraryPtr = &p_templateLibrary;
  }

  unsigned int getTemplateCount() const { return p_templateLibraryPtr->size(); }

  ROMol *getTemplate(unsigned int index) {
    return (*p_templateLibraryPtr)[index]->getMol();
  };

  unsigned int addMacroAtom(std::string className, std::string templateName);

  void addMacroBond(unsigned int fromAtomIdx, unsigned int toAtomIdx,
                    Bond::BondType bondType, std::string fromConnectionPoint,
                    std::string toConnectionPoint);

  unsigned int atomIdxToTemplateIdx(unsigned int atomIdx);
  ROMol *atomIdxToTemplateMol(unsigned int atomIdx);
};
typedef boost::shared_ptr<MACROMol> MACROMol_SPTR;
typedef boost::shared_ptr<MACROMolTemplate> MACROMolTemplate_SPTR;
}  // namespace RDKit

#endif
