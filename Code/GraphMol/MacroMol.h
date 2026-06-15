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

class RDKIT_GRAPHMOL_EXPORT MacroMolTemplate : public RDKit::RWMol {

 public:
  MacroMolTemplate(std::unique_ptr<RWMol> &mol, std::string className,
                   std::vector<std::string> templateNames,
                   std::vector<std::pair<std::string,std::string>> templateAttrs);

  MacroMolTemplate(std::unique_ptr<RWMol> &mol, std::string className,
                   std::string templateName,
                   std::vector<std::pair<std::string,std::string>> templateAttrs);

  MacroMolTemplate() = delete;
  MacroMolTemplate(const MacroMolTemplate &other);
  MacroMolTemplate(MacroMolTemplate &&other) noexcept = delete;
  MacroMolTemplate &operator=(MacroMolTemplate &&other) noexcept = delete;
  MacroMolTemplate &operator=(const MacroMolTemplate &) = delete; 
  ~MacroMolTemplate() {}

  RDKit::SubstanceGroup *getMainSgroup();
  const RDKit::SubstanceGroup *getMainSgroup() const;

 private:
  void init(std::string className,
            std::vector<std::string> templateNames,
            std::vector<std::pair<std::string,std::string>> templateAttrs);

  void findMainSgroupForTemplate(std::string className,
                                 std::string templateName) const;
  void initMainSgroupIdx() const;
  mutable unsigned int d_mainSgroupIdx;
  mutable std::once_flag d_mainSgroupIdxOnceFlag;
};


class RDKIT_GRAPHMOL_EXPORT MacroMolTemplateLib {

  private:
    friend class MacroMol;

    // All templates in the library are owned by the library. Consumers get
    // const access; the library mutates them only during construction.
    std::vector<std::unique_ptr<MacroMolTemplate>> d_templates;

    //! Key is (monomer_class, symbol)
    using MacroMolTemplateKey = std::pair<std::string, std::string>;

    //! Hash function for MacroMolTemplateKey
    struct MacroMolTemplateKeyHash {
        std::size_t operator()(const MacroMolTemplateKey& key) const {
            std::size_t h1 = std::hash<std::string>{}(key.first);
            std::size_t h2 = std::hash<std::string>{}(key.second);
            return h1 ^ (h2 << 1);
        }
    };

    std::unordered_map<MacroMolTemplateKey, unsigned int, MacroMolTemplateKeyHash> d_keyToIndex;

  protected:
  
    RDKit::MacroMolTemplate* findMutable(const std::string& templateClass, const std::string& templateName) {
      auto iter = d_keyToIndex.find({templateClass, templateName}); 
      if (iter == d_keyToIndex.end()) {
        return nullptr;
      }
      return d_templates.at(iter->second).get();
    }

  public:
    MacroMolTemplateLib() = default;
    MacroMolTemplateLib(const MacroMolTemplateLib &other) = delete;
    MacroMolTemplateLib(MacroMolTemplateLib &&other) noexcept = delete;
    MacroMolTemplateLib &operator=(MacroMolTemplateLib &&other) noexcept = delete;
    MacroMolTemplateLib &operator=(const MacroMolTemplateLib &) = delete;  
    ~MacroMolTemplateLib() {}

    std::vector<std::unique_ptr<MacroMolTemplate>>::const_iterator begin() const {
      return d_templates.begin();
    }
    std::vector<std::unique_ptr<MacroMolTemplate>>::const_iterator end() const {
      return d_templates.end();
    }

    void addTemplate(std::unique_ptr<MacroMolTemplate> &templateMol);
    void copyTemplateLib(const MacroMolTemplateLib &libToCopy);

    unsigned int size() const {
      return d_templates.size(); 
    }

    const RDKit::MacroMolTemplate* find(const std::string& templateClass, const std::string& templateName) const {
      auto iter = d_keyToIndex.find({templateClass, templateName}); 
      if (iter == d_keyToIndex.end()) {
        return nullptr;
      }
      return d_templates.at(iter->second).get();
    }

    bool contains(const std::string& templateClass, const std::string& templateName) const {
      return d_keyToIndex.contains({templateClass, templateName});
    }

    bool doesLibHaveCoords() const {
      for (auto const &macroTemplate : d_templates) {
        if (macroTemplate->getNumConformers() == 0) {
          return false;
        }
      }
      return true;
    }
  };


class RDKIT_GRAPHMOL_EXPORT MacroMol : public RWMol {
 private:
   friend class SCSRUtils;
  // elements (MacroMolTemplate items) of the library are owned by the library
  MacroMolTemplateLib d_templateLibrary;

  protected:
    // the following are used in SCSR parsing, where the templates are still being contructed, and should not be used by other callers
    MacroMolTemplate *getMutableTemplate(unsigned int atomIdx);

 public:
  MacroMol() = default;
  MacroMol(const MacroMol &other) : RWMol((RWMol) other) {
    d_templateLibrary.copyTemplateLib(*other.getTemplateLibrary());
  }
  MacroMol(MacroMol &&other) noexcept = delete;
  MacroMol &operator=(MacroMol &&other) noexcept = delete;

  MacroMol &operator=(const MacroMol &) = delete;  // disable assignment

  MacroMol(std::unique_ptr<RWMol> &rwMol)
      : RWMol(std::move(*(rwMol))) {
    rwMol = nullptr;
  }

  ~MacroMol() {}

  const MacroMolTemplate *getTemplate(unsigned int atomIdx) const;

  const MacroMolTemplateLib *getTemplateLibrary() const {
    return &d_templateLibrary;
  }

   // the following adds a template to the internal libraty for this MacroMol
  void addTemplate(std::unique_ptr<MacroMolTemplate> &templateMol) {
    PRECONDITION(templateMol, "bad template molecule");
    d_templateLibrary.addTemplate(templateMol);
  }

  unsigned int size() const { return d_templateLibrary.size(); }

  unsigned int addMacroAtom(std::string className, std::string templateName);

  void addMacroBond(unsigned int fromAtomIdx, unsigned int toAtomIdx,
                    Bond::BondType bondType, std::string fromConnectionPoint,
                    std::string toConnectionPoint);

};
typedef boost::shared_ptr<MacroMol> MacroMol_SPTR;
typedef boost::shared_ptr<MacroMolTemplate> MacroMolTemplate_SPTR;
}  // namespace RDKit

#endif
