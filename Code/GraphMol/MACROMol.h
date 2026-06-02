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
  MACROMolTemplate &operator=(const MACROMolTemplate &) = delete; 
  ~MACROMolTemplate() {}

  RDKit::SubstanceGroup *getMainSgroup();

 private:
  mutable unsigned int d_mainSgroupIdx;
  mutable std::once_flag d_mainSgroupIdxOnceFlag;
};


   //! Key is (symbol, monomer_class)
class MACROMolTemplateKey: public std::pair<std::string, std::string>{};


class RDKIT_GRAPHMOL_EXPORT MACROMolTemplateLib {
  
  
    private:
      std::vector<std::unique_ptr<MACROMolTemplate>> d_templates;

    //! Hash function for MonomerKey
    struct MACROMolTemplateKeyHash {
        std::size_t operator()(const MACROMolTemplateKey& key) const {
            std::size_t h1 = std::hash<std::string>{}(key.first);
            std::size_t h2 = std::hash<std::string>{}(key.second);
            return h1 ^ (h2 << 1);
        }
    };

    std::unordered_map<MACROMolTemplateKey, unsigned int, MACROMolTemplateKeyHash> d_keyToIndex;

    // all templates in the library are owned by the library

public:
    MACROMolTemplateLib() = default;
    MACROMolTemplateLib(const MACROMolTemplateLib &other) = delete;
    MACROMolTemplateLib(MACROMolTemplateLib &&other) noexcept = delete;
    MACROMolTemplateLib &operator=(MACROMolTemplateLib &&other) noexcept = delete;
    MACROMolTemplateLib &operator=(const MACROMolTemplateLib &) = delete;  
    ~MACROMolTemplateLib() {}

    std::vector<std::unique_ptr<MACROMolTemplate>>::iterator begin() {
      return d_templates.begin();
    }
    std::vector<std::unique_ptr<MACROMolTemplate>>::iterator end() {
      return d_templates.end();
    }std::vector<std::unique_ptr<MACROMolTemplate>>::const_iterator begin() const{
      return d_templates.begin();
    }
    std::vector<std::unique_ptr<MACROMolTemplate>>::const_iterator end() const{
      return d_templates.end();
    }

    void addTemplate(std::unique_ptr<MACROMolTemplate> &templateMol);
    void addTemplateLib(MACROMolTemplateLib &libToAdd);
    void copyTemplateLib(const MACROMolTemplateLib &libToCopy);
    void clearTemplateLib(){
        d_templates.clear();
        d_keyToIndex.clear();
    }
    void setTemplateLib(MACROMolTemplateLib &libToSet) {
      clearTemplateLib();
      addTemplateLib(libToSet);
    }

    unsigned int getNumTemplates(){
      return d_templates.size();
    }

    unsigned int getNumTemplates() const { 
      return d_templates.size(); 
    }

    RDKit::MACROMolTemplate *getTemplate(unsigned int index) const {
      return d_templates.at(index).get();
    }
    
    unsigned int getMACROMolTemplateIndex(std::string templateClass, std::string templateName) const {
      MACROMolTemplateKey key(std::pair(templateClass, templateName));
      if ( ! d_keyToIndex.contains(key)) {
         return UINT_MAX;
      }
      return d_keyToIndex.at(key);
    }

    RDKit::MACROMolTemplate *getTemplate(std::string templateClass, std::string templateName) const {
      auto index = getMACROMolTemplateIndex(templateClass, templateName);
      if (index == UINT_MAX) {
        return nullptr;
      }
      return d_templates.at(index).get();
    }

    bool libContains(std::string templateClass, std::string templateName) const {
      MACROMolTemplateKey key(std::pair(templateClass, templateName));
      return d_keyToIndex.contains(key);
    }

    bool doesLibHaveCoords() {
      for (auto const &macroTemplate : d_templates) {
        if (macroTemplate->getNumConformers() == 0) {
          return false;
        }
      }
      return true;
    }
  };


class RDKIT_GRAPHMOL_EXPORT MACROMol : public RWMol {
 private:
  // elements (MACROMolTemplate items) of the library are either owned by the library or owned externally.
  // Externally owned items support global libraies

  MACROMolTemplateLib d_templateLibrary;

  std::vector<const MACROMolTemplateLib *> d_externalTemplateLibs;

  mutable bool d_atomIdxToTemplatePtrIsStale;
  mutable std::map<unsigned int, MACROMolTemplate *> d_atomIdxToTemplatePtr;

 public:
  MACROMol() : d_atomIdxToTemplatePtrIsStale(true) {};
  MACROMol(const MACROMol &other) : RWMol((RWMol) other), d_atomIdxToTemplatePtrIsStale(true){
    for (const auto templateLib : other.d_externalTemplateLibs) {
      this->d_externalTemplateLibs.push_back(templateLib);
    }
    d_templateLibrary.copyTemplateLib(*other.getTemplateLibrary());
  }
  MACROMol(MACROMol &&other) noexcept = delete;
  MACROMol &operator=(MACROMol &&other) noexcept = delete;

  MACROMol &operator=(const MACROMol &) = delete;  // disable assignment

  MACROMol(std::unique_ptr<RWMol> &rwMol)
      : RWMol(std::move(*(rwMol))), d_atomIdxToTemplatePtrIsStale(true) {
    rwMol = nullptr;
  }

  ~MACROMol() {}

  MACROMolTemplate *atomIdxToTemplatePtr(unsigned int atomIdx) const;

  MACROMolTemplateLib *getTemplateLibrary() {
    return &d_templateLibrary;
  }
    const MACROMolTemplateLib *getTemplateLibrary() const {
    return &d_templateLibrary;
  }

  void clearTemplateLibrary() {
    d_templateLibrary.clearTemplateLib();
  }

  void clearExternalTemplateLibraries() {
    d_externalTemplateLibs.clear();
  }

  void addTemplateLibrary(const MACROMolTemplateLib *lib) {
   d_externalTemplateLibs.push_back(lib);
  }


   // the following adds a template to the internal libraty for this MACROMol
  void addTemplate(std::unique_ptr<MACROMolTemplate> &templateMol) {
    PRECONDITION(templateMol, "bad template molecule");
    d_templateLibrary.addTemplate(templateMol);
  }

  unsigned int getNumTemplates() const { return d_templateLibrary.getNumTemplates(); }

  MACROMolTemplate *getTemplate(unsigned int idx) { return d_templateLibrary.getTemplate(idx); }
  const MACROMolTemplate *getTemplate(unsigned int idx) const { return d_templateLibrary.getTemplate(idx); }

  unsigned int getNumExternalTemplateLibs() const { return d_externalTemplateLibs.size();}

  unsigned int addMacroAtom(std::string className, std::string templateName);

  void addMacroBond(unsigned int fromAtomIdx, unsigned int toAtomIdx,
                    Bond::BondType bondType, std::string fromConnectionPoint,
                    std::string toConnectionPoint);

};
typedef boost::shared_ptr<MACROMol> MACROMol_SPTR;
typedef boost::shared_ptr<MACROMolTemplate> MACROMolTemplate_SPTR;
}  // namespace RDKit

#endif
