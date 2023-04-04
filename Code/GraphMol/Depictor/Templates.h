//
//  Copyright (C) 2023 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>

#include "TemplateSmiles.h"

#include <iostream>
#include <fstream>
#include <unordered_map>

class CoordinateTemplates {
 public:
  //! returns a reference to the singleton CoordinateTemplates
  /*
      \return a reference to the singleton CoordinateTemplates

      <b>Notes:</b>
        - if the singleton CoordinateTemplates has already been instantiated
          the singleton will be returned, otherwise the singleton will
          be constructed.

   */
  static CoordinateTemplates& getRingSystemTemplates() {
    static CoordinateTemplates template_mols;
    return template_mols;
  }

  bool hasTemplateOfSize(unsigned int atom_count) {
    if (m_templates.find(atom_count) != m_templates.end()) {
        return true;
    }
    return false;
  }

  const std::vector<std::shared_ptr<RDKit::ROMol>>& getMatchingTemplates(unsigned int atom_count) {
    return m_templates[atom_count];
  }

  bool setRingSystemTemplates(const std::string& template_dir) {
    std::ifstream cxsmiles(template_dir);
    if (!cxsmiles) {
      BOOST_LOG(rdWarningLog)
          << "WARNING: Could not open file " << template_dir << std::endl;
      return false;
    }

    // Try loading templates in from this directory, if unsuccessful, keep current templates
    std::unordered_map<unsigned int, std::vector<std::shared_ptr<RDKit::ROMol>>> templates;
    std::string line;
    bool success = true;
    while (std::getline(cxsmiles, line)) {
      RDKit::ROMol* mol_ptr = RDKit::SmilesToMol(line);
      if (!mol_ptr) {
        BOOST_LOG(rdWarningLog)
            << "WARNING: Could not load templates from " << template_dir
            << ": Invalid smiles" << std::endl;
        success = false;
        break;
      }
      std::shared_ptr<RDKit::ROMol> mol(mol_ptr);
      // check if mol has coordinates -- 3D coordinates are OK
      if (mol->getNumConformers() == 0) {
        BOOST_LOG(rdWarningLog)
            << "WARNING: Could not load templates from " << template_dir
            << ": One or more templates are missing coordinates\n";
        success = false;
        break;
      }
      // template must be a ring system
      if (!isRingSystem(*mol)) {
        BOOST_LOG(rdWarningLog)
            << "WARNING: Could not load templates from " << template_dir
            << ": One or more templates are not a single ring system\n";
        success = false;
        break;
      }
      templates[mol->getNumAtoms()].push_back(mol);
    }
    cxsmiles.close();
    if (success) {
      clearTemplates();
      m_templates = std::move(templates);
    }
    return success;
  }

  void loadDefaultTemplates() {
    clearTemplates();
    // load default templates into m_templates map by atom count
    for (const auto& smiles : TEMPLATE_SMILES) {
        std::shared_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));
        m_templates[mol->getNumAtoms()].push_back(mol);
    }
  }

 private:
  CoordinateTemplates() {
    loadDefaultTemplates();
  }
  CoordinateTemplates(const CoordinateTemplates&) = delete;
  CoordinateTemplates& operator= (const CoordinateTemplates&) = delete;

  void clearTemplates() {
    for (auto& [atom_cout, romols] : m_templates) {
      romols.clear();
    }
    m_templates.clear();
  }

  bool isRingSystem(RDKit::ROMol& mol) {
    // Make sure this template is a ring system (spiro'd ring systems are OK).
    // We can check this by ensuring that every atom is in a ring and that there
    // is only one connected component in the molecular graph.
    if (RDKit::MolOps::getMolFrags(mol).size() != 1) {
      return false;
    }

    // Use symmetrizeSSSR since that is what is used in coordinate generation
    RDKit::MolOps::symmetrizeSSSR(mol);
    auto ri = mol.getRingInfo();
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      if (!ri->numAtomRings(i)) {
        return false;
      }
    }
    return true;
  }

  ~CoordinateTemplates() {
    clearTemplates();
  }

  std::unordered_map<unsigned int, std::vector<std::shared_ptr<RDKit::ROMol>>> m_templates;
};