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
#include <map>
#include <memory>
#include <string>
#include <utility>

namespace RDKit {

struct RDKIT_GRAPHMOL_EXPORT MacroMolTemplate {
  std::string monomerClass;  // MonomerClass as string, e.g., "AA", "NA", "CHEM"
  std::string templateName;  // Name of the template, e.g., "ALA"
  std::string symbol;        // e.g., "A" for Alanine
  std::string original_data;   // Original definition (SMILES, SDF, etc.)
  std::shared_ptr<RWMol> mol;  // Parsed molecule
};

class RDKIT_GRAPHMOL_EXPORT MacroMolTemplateLibrary {
 public:
  void addTemplate(const std::shared_ptr<MacroMolTemplate> &macroMolTemplate);
  const std::shared_ptr<MacroMolTemplate> &getByTemplateName(
      const std::string &monomerClass, const std::string &templateName) const;
  const std::shared_ptr<MacroMolTemplate> &getBySymbol(
      const std::string &monomerClass, const std::string &symbol) const;

 private:
  using MacroMolTemplateKey = std::pair<std::string, std::string>;

  std::map<MacroMolTemplateKey, std::shared_ptr<MacroMolTemplate>>
      byTemplateName;
  std::map<MacroMolTemplateKey, std::shared_ptr<MacroMolTemplate>> bySymbol;
};

}  // namespace RDKit

#endif
