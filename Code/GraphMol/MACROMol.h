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
    void clearTemplateLib(){
        this->clear();
        ownedTemplates.clear();
        d_keyToIndex.clear();
    }
    void setTemplateLib(MACROMolTemplateLib &libToSet, bool takeOwnership=false) {
      clearTemplateLib();
      addTemplateLib(libToSet, takeOwnership);
    }


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
  // elements (MACROMolTemplate items) of the library are either owned by the library or owned externally.
  // Externally owned items support global libraies

  MACROMolTemplateLib p_templateLibrary;

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
    return &p_templateLibrary;
  }

  void clearTemplateLibrary() {
    p_templateLibrary.clear();
  }

  void addTemplateLibrary(MACROMolTemplateLib &lib,
                          bool takeOwnership = true) {
    p_templateLibrary.addTemplateLib(lib, takeOwnership);
  }


  void setTemplateLibrary(MACROMolTemplateLib &lib,
                          bool takeOwnership = true) {
    p_templateLibrary.setTemplateLib(lib, takeOwnership);
  }

  void addTemplate(std::unique_ptr<MACROMolTemplate> &templateMol, bool takeOwnership=true) {
    PRECONDITION(templateMol, "bad template molecule");
    p_templateLibrary.addTemplate(templateMol, takeOwnership);
  }

  unsigned int getTemplateCount() const { return p_templateLibrary.getTemplateCount(); }

  RDKit::MACROMolTemplate *getTemplate(unsigned int index) {
    return p_templateLibrary.getTemplate(index);
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
