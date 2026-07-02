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
  std::string monomerClass;  // MonomerClass as string, e.g., "AA", "NA", "CHEM"
  std::string templateName;  // Name of the template, e.g., "ALA"
  std::string symbol;        // e.g., "A" for Alanine
  std::string original_data;  // Original definition (SMILES, SDF, etc.)
  // Parsed, annotated template molecule
  std::shared_ptr<MacroMolTemplate> molTemplate;
};

class RDKIT_GRAPHMOL_EXPORT MacroMolTemplateLibrary {
 public:
  void addEntry(const std::shared_ptr<MacroMolEntry> &macroMolEntry);
  const std::shared_ptr<MacroMolEntry> &getByTemplateName(
      const std::string &monomerClass, const std::string &templateName) const;
  const std::shared_ptr<MacroMolEntry> &getBySymbol(
      const std::string &monomerClass, const std::string &symbol) const;

 private:
  using MacroMolTemplateKey = std::pair<std::string, std::string>;

  std::map<MacroMolTemplateKey, std::shared_ptr<MacroMolEntry>> byTemplateName;
  std::map<MacroMolTemplateKey, std::shared_ptr<MacroMolEntry>> bySymbol;
};

}  // namespace RDKit

#endif
