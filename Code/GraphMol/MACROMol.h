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


class RDKIT_GRAPHMOL_EXPORT MACROMolTemplateLib {

  private:
    // All templates in the library are owned by the library. Consumers get
    // const access; the library mutates them only during construction.
    std::vector<std::unique_ptr<MACROMolTemplate>> d_templates;

    //! Key is (monomer_class, symbol)
    using MACROMolTemplateKey = std::pair<std::string, std::string>;

    //! Hash function for MACROMolTemplateKey
    struct MACROMolTemplateKeyHash {
        std::size_t operator()(const MACROMolTemplateKey& key) const {
            std::size_t h1 = std::hash<std::string>{}(key.first);
            std::size_t h2 = std::hash<std::string>{}(key.second);
            return h1 ^ (h2 << 1);
        }
    };

    std::unordered_map<MACROMolTemplateKey, unsigned int, MACROMolTemplateKeyHash> d_keyToIndex;

  public:
    MACROMolTemplateLib() = default;
    MACROMolTemplateLib(const MACROMolTemplateLib &other) = delete;
    MACROMolTemplateLib(MACROMolTemplateLib &&other) noexcept = delete;
    MACROMolTemplateLib &operator=(MACROMolTemplateLib &&other) noexcept = delete;
    MACROMolTemplateLib &operator=(const MACROMolTemplateLib &) = delete;  
    ~MACROMolTemplateLib() {}

    std::vector<std::unique_ptr<MACROMolTemplate>>::const_iterator begin() const {
      return d_templates.begin();
    }
    std::vector<std::unique_ptr<MACROMolTemplate>>::const_iterator end() const {
      return d_templates.end();
    }

    void addTemplate(std::unique_ptr<MACROMolTemplate> &templateMol);
    void copyTemplateLib(const MACROMolTemplateLib &libToCopy);

    unsigned int size() const {
      return d_templates.size(); 
    }

    const RDKit::MACROMolTemplate* find(const std::string& templateClass, const std::string& templateName) const {
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


class RDKIT_GRAPHMOL_EXPORT MACROMol : public RWMol {
 private:
  // elements (MACROMolTemplate items) of the library are owned by the library
  MACROMolTemplateLib d_templateLibrary;

 public:
  MACROMol() = default;
  MACROMol(const MACROMol &other) : RWMol((RWMol) other) {
    d_templateLibrary.copyTemplateLib(*other.getTemplateLibrary());
  }
  MACROMol(MACROMol &&other) noexcept = delete;
  MACROMol &operator=(MACROMol &&other) noexcept = delete;

  MACROMol &operator=(const MACROMol &) = delete;  // disable assignment

  MACROMol(std::unique_ptr<RWMol> &rwMol)
      : RWMol(std::move(*(rwMol))) {
    rwMol = nullptr;
  }

  ~MACROMol() {}

  const MACROMolTemplate *getTemplate(unsigned int atomIdx) const;

  const MACROMolTemplateLib *getTemplateLibrary() const {
    return &d_templateLibrary;
  }

   // the following adds a template to the internal libraty for this MACROMol
  void addTemplate(std::unique_ptr<MACROMolTemplate> &templateMol) {
    PRECONDITION(templateMol, "bad template molecule");
    d_templateLibrary.addTemplate(templateMol);
  }

  unsigned int size() const { return d_templateLibrary.size(); }

  unsigned int addMacroAtom(std::string className, std::string templateName);

  void addMacroBond(unsigned int fromAtomIdx, unsigned int toAtomIdx,
                    Bond::BondType bondType, std::string fromConnectionPoint,
                    std::string toConnectionPoint);

};
typedef boost::shared_ptr<MACROMol> MACROMol_SPTR;
typedef boost::shared_ptr<MACROMolTemplate> MACROMolTemplate_SPTR;
}  // namespace RDKit

#endif
