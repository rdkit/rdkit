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

 private:
  CoordinateTemplates() {
    // load templates into m_templates map by atom count
    for (const auto& smiles : TEMPLATE_SMILES) {
        std::shared_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));
        m_templates[mol->getNumAtoms()].push_back(mol);
    }
  }

  CoordinateTemplates(const CoordinateTemplates&) = delete;
  CoordinateTemplates& operator= (const CoordinateTemplates&) = delete;
  ~CoordinateTemplates() {
    for (auto& [atom_cout, romols] : m_templates) {
      romols.clear();
    }
    m_templates.clear();
  }

  std::unordered_map<unsigned int, std::vector<std::shared_ptr<RDKit::ROMol>>> m_templates;
};