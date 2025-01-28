//
//  Copyright (C) 2023 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "Templates.h"

#include "RDDepictor.h"

namespace RDDepict {
void CoordinateTemplates::assertValidTemplate(RDKit::ROMol &mol,
                                              const std::string &smiles) {
  // template must have 2D coordinates
  if (mol.getNumConformers() == 0) {
    std::string msg = "Template missing coordinates: " + smiles;
    throw RDDepict::DepictException(msg);
  }
  if (mol.getConformer().is3D()) {
    std::string msg =
        "Template has 3D coordinates, 2D coordinates required: " + smiles;
    throw RDDepict::DepictException(msg);
  }

  // Make sure this template is a single ring system (spiro'd ring systems are
  // OK). We can check this by ensuring that every bond is in a ring and that
  // there is only one connected component in the molecular graph.
  if (RDKit::MolOps::getMolFrags(mol).size() != 1) {
    std::string msg =
        "Template consists of multiple fragments, single fragment required: " +
        smiles;
    throw RDDepict::DepictException(msg);
  }

  if (mol.getNumAtoms() == 1) {
    std::string msg = "Template is not a ring system: " + smiles;
    throw RDDepict::DepictException(msg);
  }
  auto ri = mol.getRingInfo();
  for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
    if (!ri->numBondRings(i)) {
      std::string msg = "Template is not a ring system: " + smiles;
      throw RDDepict::DepictException(msg);
    }
  }
}

void CoordinateTemplates::loadTemplatesFromPath(
    const std::string &templatePath,
    std::unordered_map<unsigned int, std::vector<std::shared_ptr<RDKit::ROMol>>>
        &templates) {
  std::ifstream cxsmiles(templatePath);
  if (!cxsmiles) {
    std::string msg = "Could not open file " + templatePath;
    throw RDDepict::DepictException(msg);
  }

  // Try loading templates in from this directory, if unsuccessful, keep current
  // templates
  std::string line;
  while (std::getline(cxsmiles, line)) {
    RDKit::ROMol *mol_ptr = RDKit::SmilesToMol(line);
    if (!mol_ptr) {
      std::string msg =
          "Could not load templates from " + templatePath + ": Invalid smiles";
      cxsmiles.close();
      throw RDDepict::DepictException(msg);
    }
    std::shared_ptr<RDKit::ROMol> mol(mol_ptr);
    try {
      assertValidTemplate(*mol, line);
    } catch (RDDepict::DepictException &e) {
      cxsmiles.close();
      throw e;
    }
    templates[mol->getNumAtoms()].push_back(mol);
  }
  cxsmiles.close();
}

void CoordinateTemplates::setRingSystemTemplates(
    const std::string &templatePath) {
  // Try loading templates in from this directory, if unsuccessful, keep current
  // templates
  std::unordered_map<unsigned int, std::vector<std::shared_ptr<RDKit::ROMol>>>
      templates;
  loadTemplatesFromPath(templatePath, templates);
  clearTemplates();
  m_templates = std::move(templates);
}

void CoordinateTemplates::addRingSystemTemplates(
    const std::string &templatePath) {
  // Try loading templates in from this directory, if unsuccessful, keep current
  // templates
  std::unordered_map<unsigned int, std::vector<std::shared_ptr<RDKit::ROMol>>>
      templates;
  loadTemplatesFromPath(templatePath, templates);
  for (auto &kv : templates) {
    m_templates[kv.first].insert(m_templates[kv.first].begin(),
                                 kv.second.begin(), kv.second.end());
  }
}
}  // namespace RDDepict