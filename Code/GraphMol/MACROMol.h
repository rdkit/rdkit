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
#include <mutex>
#include <unordered_map>


namespace RDKit {

class RDKIT_GRAPHMOL_EXPORT MACROMolTemplate : public RDKit::RWMol {
 private:
  void init(std::string className,
            std::vector<std::string> templateNames,
            std::vector<std::pair<std::string,std::string>> templateAttrs);

  void findMainSgroupForTemplate(std::string className,
                                 std::string templateName) const;

 public:
  MACROMolTemplate(std::unique_ptr<RWMol> &mol, std::string className,
                   std::vector<std::string> templateNames,
                   std::vector<std::pair<std::string,std::string>> templateAttrs);

  MACROMolTemplate(std::unique_ptr<RWMol> &mol, std::string className,
                   std::string templateName,
                   std::vector<std::pair<std::string,std::string>> templateAttrs);

  MACROMolTemplate() = delete;
  MACROMolTemplate(const MACROMolTemplate &other);
  MACROMolTemplate(MACROMolTemplate &&other) noexcept = delete;
  MACROMolTemplate &operator=(MACROMolTemplate &&other) noexcept = delete;
  MACROMolTemplate &operator=(const MACROMolTemplate &) =
      delete;  // disable assignment
  ~MACROMolTemplate() {}

  RDKit::SubstanceGroup *getMainSgroup();

 private:
  mutable unsigned int p_mainSgroupIdx;
  mutable std::once_flag p_mainSgroupIdxOnceFlag;
};


   //! Key is (symbol, monomer_class)
class MACROMolTemplateKey: public std::pair<std::string, std::string>{};


class RDKIT_GRAPHMOL_EXPORT MACROMolTemplateLib : public std::vector<MACROMolTemplate *> {
  
    private:

    //! Hash function for MonomerKey
    struct MACROMolTemplateKeyHash {
        std::size_t operator()(const MACROMolTemplateKey& key) const {
            std::size_t h1 = std::hash<std::string>{}(key.first);
            std::size_t h2 = std::hash<std::string>{}(key.second);
            return h1 ^ (h2 << 1);
        }
    };

    std::unordered_map<MACROMolTemplateKey, unsigned int, MACROMolTemplateKeyHash> d_keyToIndex;

    // all templates USED are references.  Some may be owned by this lib, and some may be owned by
    // an external lib.  The ones owned by this lib are destroyed when this template is destroyed.
    // this allows one or more global libraries to be added, but the member template are owned by
    // externally. 
    //
    // Also, if the templates are owned externally, one or more local templates can be added.
    // THis allows, an an example, for "immediate" templates such the the MonomerMol smiles-type, which 
    // are defined as a smiles string in the template instantiation.

    std::vector<std::unique_ptr<MACROMolTemplate>> ownedTemplates;   

    public:
    void addTemplate(MACROMolTemplate *templateMol);    // does NOT take ownership
    void addTemplate(std::unique_ptr<MACROMolTemplate> &templateMol, bool takeOwnership=true);
    void addTemplateLib(MACROMolTemplateLib &libToAdd, bool takeOwnership=false);


    unsigned int getTemplateCount() const { return this->size(); }

    RDKit::MACROMolTemplate *getTemplate(unsigned int index) {
      return this->at(index);
    }
    
    unsigned int getMACROMolTemplateIndex(std::string templateClass, std::string templateName) {
      MACROMolTemplateKey key(std::pair(templateClass, templateName));
      if ( ! d_keyToIndex.contains(key)) {
         throw FileParseException("Template for macro atom not found");
      }
      return d_keyToIndex[key];
    }
  };


class RDKIT_GRAPHMOL_EXPORT MACROMol : public RWMol {
 private:
  // for libraries owned by this macromol
  MACROMolTemplateLib p_templateLibrary;

  // used for library access.  Could be a pointer to the internal library or
  // an external one that is shared across multiple macromols

  MACROMolTemplateLib *p_templateLibraryPtr;
  bool p_atomIdxToTemplateIdxIsStale;
  std::map<unsigned int, unsigned int> p_atomIdxToTemplateIdx;

 public:
  MACROMol() : p_atomIdxToTemplateIdxIsStale(true) {};
  MACROMol(const MACROMol &other) = delete;
  MACROMol(MACROMol &&other) noexcept = delete;
  MACROMol &operator=(MACROMol &&other) noexcept = delete;

  MACROMol &operator=(const MACROMol &) = delete;  // disable assignment

  MACROMol(std::unique_ptr<RWMol> &rwMol)
      : RWMol(std::move(*(rwMol.get()))), p_atomIdxToTemplateIdxIsStale(true) {
    rwMol.release();
  }

  ~MACROMol() {}

  MACROMolTemplateLib *getTemplateLibrary() {
    return p_templateLibraryPtr;
  }

  void setTemplateLibrary(MACROMolTemplateLib &lib,
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

  void addTemplate(std::unique_ptr<MACROMolTemplate> &templateMol) {
    PRECONDITION(templateMol, "bad template molecule");
    p_templateLibrary.addTemplate(templateMol);
    p_templateLibraryPtr = &p_templateLibrary;
  }

  unsigned int getTemplateCount() const { return p_templateLibraryPtr->getTemplateCount(); }

  RDKit::MACROMolTemplate *getTemplate(unsigned int index) {
    return p_templateLibraryPtr->getTemplate(index);
  };

  unsigned int addMacroAtom(std::string className, std::string templateName);

  void addMacroBond(unsigned int fromAtomIdx, unsigned int toAtomIdx,
                    Bond::BondType bondType, std::string fromConnectionPoint,
                    std::string toConnectionPoint);

  unsigned int atomIdxToTemplateIdx(unsigned int atomIdx);
  MACROMolTemplate *atomIdxToMACROMolTemplate(unsigned int atomIdx);
};
typedef boost::shared_ptr<MACROMol> MACROMol_SPTR;
typedef boost::shared_ptr<MACROMolTemplate> MACROMolTemplate_SPTR;
}  // namespace RDKit

#endif
