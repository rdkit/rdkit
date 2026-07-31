//
//  Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_MACROMOLTEMPLATE_H
#define RD_MACROMOLTEMPLATE_H

#include <RDGeneral/export.h>
#include <GraphMol/MacroAtomInfo.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SubstanceGroup.h>

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace RDKit {

//! Describes the atoms removed at one MacroMol template attachment point.
struct RDKIT_GRAPHMOL_EXPORT MacroMolLeavingGroup {
  std::vector<unsigned int> atomIdxs;
  unsigned int attachAtomIdx = 0;
  unsigned int leavingAtomIdx = 0;
  int attachPoint = 0;
};

class MacroMolTemplateBuilder;

//! An immutable, annotated molecule template for a macromolecule monomer.
/*!
  The typed main- and leaving-group definitions are authoritative. The owned
  molecule mirrors those definitions as SUP SGroups for RDKit I/O.
*/
class RDKIT_GRAPHMOL_EXPORT MacroMolTemplate final {
 public:
  MacroMolTemplate(const MacroMolTemplate &) = default;
  MacroMolTemplate(MacroMolTemplate &&) noexcept = default;
  MacroMolTemplate &operator=(const MacroMolTemplate &) = delete;
  MacroMolTemplate &operator=(MacroMolTemplate &&) = delete;

  //! Returns the immutable underlying template molecule.
  const ROMol &getMol() const { return d_mol; }
  //! Returns the monomer class used for lookup and the main SGroup CLASS.
  MonomerClass getMonomerClass() const { return d_monomerClass; }
  //! Returns the template name, e.g. "ALA".
  const std::string &getTemplateName() const { return d_templateName; }
  //! Returns the template symbol, e.g. "A" for alanine.
  const std::string &getSymbol() const { return d_symbol; }
  //! Returns the original template definition (SMILES, SDF, etc.).
  const std::string &getOriginalData() const { return d_originalData; }
  //! Returns the atom indices belonging to the retained monomer group.
  const std::vector<unsigned int> &getMainAtomIdxs() const {
    return d_mainAtomIdxs;
  }
  //! Returns the typed leaving-group definitions.
  const std::vector<MacroMolLeavingGroup> &getLeavingGroups() const {
    return d_leavingGroups;
  }
  //! Returns the main SUP SGroup mirrored onto the owned molecule.
  const SubstanceGroup &getMainSgroup() const;

 private:
  friend class MacroMolTemplateBuilder;

  MacroMolTemplate(RWMol mol, MonomerClass monomerClass,
                   std::string templateName, std::string symbol,
                   std::string originalData,
                   std::vector<unsigned int> mainAtomIdxs,
                   std::vector<MacroMolLeavingGroup> leavingGroups,
                   unsigned int mainSgroupIdx)
      : d_mol(std::move(mol)),
        d_monomerClass(monomerClass),
        d_templateName(std::move(templateName)),
        d_symbol(std::move(symbol)),
        d_originalData(std::move(originalData)),
        d_mainAtomIdxs(std::move(mainAtomIdxs)),
        d_leavingGroups(std::move(leavingGroups)),
        d_mainSgroupIdx(mainSgroupIdx) {}

  RWMol d_mol;
  MonomerClass d_monomerClass;
  std::string d_templateName;
  std::string d_symbol;
  std::string d_originalData;
  std::vector<unsigned int> d_mainAtomIdxs;
  std::vector<MacroMolLeavingGroup> d_leavingGroups;
  unsigned int d_mainSgroupIdx;
};

//! Builds and validates a completed MacroMolTemplate.
class RDKIT_GRAPHMOL_EXPORT MacroMolTemplateBuilder {
 public:
  MacroMolTemplateBuilder(const ROMol &mol, MonomerClass monomerClass,
                          std::string templateName, std::string symbol,
                          std::string originalData)
      : d_mol(mol),
        d_monomerClass(monomerClass),
        d_templateName(std::move(templateName)),
        d_symbol(std::move(symbol)),
        d_originalData(std::move(originalData)) {}

  //! Defines the atoms retained as the main monomer group.
  MacroMolTemplateBuilder &setMainGroup(
      std::vector<unsigned int> atomIdxs);
  //! Adds a leaving group and its attachment-point definition.
  MacroMolTemplateBuilder &addLeavingGroup(MacroMolLeavingGroup leavingGroup);
  //! Validates the complete definition and returns an immutable template.
  std::unique_ptr<MacroMolTemplate> build() &&;

 private:
  RWMol d_mol;
  MonomerClass d_monomerClass;
  std::string d_templateName;
  std::string d_symbol;
  std::string d_originalData;
  std::vector<unsigned int> d_mainAtomIdxs;
  std::vector<MacroMolLeavingGroup> d_leavingGroups;
  bool d_mainGroupSet = false;
};

class RDKIT_GRAPHMOL_EXPORT MacroMolTemplateLibrary {
 public:
  //! Adds a completed template and takes ownership of it.
  void addTemplate(std::unique_ptr<MacroMolTemplate> macroMolTemplate);

  //! Returns templates ordered by descending main-group size.
  const std::vector<const MacroMolTemplate *> &entries() const;

  //! Returns a matching template, or nullptr if none has been added.
  const MacroMolTemplate *getByName(
      MonomerClass monomerClass, const std::string &templateName) const;
  //! Returns a matching template, or nullptr if none has been added.
  const MacroMolTemplate *getBySymbol(
      MonomerClass monomerClass, const std::string &symbol) const;

 private:
  using MacroMolTemplateKey = std::pair<MonomerClass, std::string>;

  std::map<MacroMolTemplateKey, const MacroMolTemplate *> byTemplateName;
  std::map<MacroMolTemplateKey, const MacroMolTemplate *> bySymbol;
  std::vector<std::unique_ptr<const MacroMolTemplate>> ownedTemplates;
  std::vector<const MacroMolTemplate *> orderedEntries;
};

}  // namespace RDKit

#endif
