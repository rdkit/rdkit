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
#include <GraphMol/RWMol.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/MacroAtomInfo.h>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace RDKit {
class RDKIT_GRAPHMOL_EXPORT MacroMolTemplate : public RWMol {
 public:
  MacroMolTemplate() = default;
  explicit MacroMolTemplate(const RWMol &other) : RWMol(other) {}

  unsigned int setMainGroup(const std::vector<unsigned int> &atomIdxs,
                            const std::string &className);

  unsigned int addLeavingGroup(const std::vector<unsigned int> &atomIdxs,
                               unsigned int attachAtomIdx,
                               unsigned int leavingAtomIdx, const int attachPt);

  const SubstanceGroup *getMainSgroup() const;

  std::vector<const SubstanceGroup *> getLeavingGroups() const;

 private:
  SubstanceGroup *getMainSgroup();
};

struct RDKIT_GRAPHMOL_EXPORT MacroMolEntry {
  MonomerClass monomerClass = MonomerClass::OTHER;  // e.g., AA, NA, CHEM
  std::string templateName;  // Name of the template, e.g., "ALA"
  std::string symbol;        // e.g., "A" for Alanine
  std::string original_data;  // Original definition (SMILES, SDF, etc.)
  // Parsed, annotated template molecule
  std::shared_ptr<MacroMolTemplate> molTemplate;

  // Library-managed: set by MacroMolTemplateLibrary::addEntry from the
  // template's main group and used to order entries largest-first. Do not set
  // manually.
  size_t mainGroupSize = 0;
};

class RDKIT_GRAPHMOL_EXPORT MacroMolTemplateLibrary {
 public:
  void addEntry(const std::shared_ptr<MacroMolEntry> &macroMolEntry);
  const std::vector<std::shared_ptr<MacroMolEntry>> &entries() const;
  const std::shared_ptr<MacroMolEntry> &getByTemplateName(
      MonomerClass monomerClass, const std::string &templateName) const;
  const std::shared_ptr<MacroMolEntry> &getBySymbol(
      MonomerClass monomerClass, const std::string &symbol) const;

 private:
  using MacroMolTemplateKey = std::pair<MonomerClass, std::string>;

  std::map<MacroMolTemplateKey, std::shared_ptr<MacroMolEntry>> byTemplateName;
  std::map<MacroMolTemplateKey, std::shared_ptr<MacroMolEntry>> bySymbol;
  std::vector<std::shared_ptr<MacroMolEntry>> orderedEntries;
};

}  // namespace RDKit

#endif
