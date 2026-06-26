//
//  Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "MacroMolTemplate.h"

namespace RDKit {

void MacroMolTemplateLibrary::addTemplate(
    const std::shared_ptr<MacroMolTemplate> &tmpl) {
  byTemplateName[{tmpl->monomerClass, tmpl->templateName}] = tmpl;
  bySymbol[{tmpl->monomerClass, tmpl->symbol}] = tmpl;
}

const std::shared_ptr<MacroMolTemplate> &
MacroMolTemplateLibrary::getByTemplateName(
    const std::string &monomerClass, const std::string &templateName) const {
  static std::shared_ptr<MacroMolTemplate> empty;
  auto it = byTemplateName.find({monomerClass, templateName});
  if (it != byTemplateName.end()) return it->second;
  return empty;
}

const std::shared_ptr<MacroMolTemplate> &MacroMolTemplateLibrary::getBySymbol(
    const std::string &monomerClass, const std::string &symbol) const {
  static std::shared_ptr<MacroMolTemplate> empty;
  auto it = bySymbol.find({monomerClass, symbol});
  if (it != bySymbol.end()) return it->second;
  return empty;
}

}  // namespace RDKit
