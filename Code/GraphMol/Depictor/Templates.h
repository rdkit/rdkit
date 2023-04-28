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

namespace RDDepict {
class RDKIT_DEPICTOR_EXPORT CoordinateTemplates {
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

  bool hasTemplateOfSize(unsigned int atomCount) {
    if (m_templates.find(atomCount) != m_templates.end()) {
      return true;
    }
    return false;
  }

  const std::vector<std::shared_ptr<RDKit::ROMol>>& getMatchingTemplates(
      unsigned int atomCount) {
    return m_templates[atomCount];
  }

  void setRingSystemTemplates(const std::string& templatePath);
  void addRingSystemTemplates(const std::string& templatePath);

  //! check if a template is considered valid
  /*
      \param template mol
      \param smiles for logging
      \return whether the template is valid

      <b>A template is considered valid if it:</b>
        - has 2D coordinates
        - is a ring system (spiro'd ring systems are OK)
        - consists of only 1 fragment

   */
  static void assertValidTemplate(RDKit::ROMol& mol, const std::string& smiles);

  void loadDefaultTemplates() {
    clearTemplates();
    // load default templates into m_templates map by atom count
    for (const auto& smiles : TEMPLATE_SMILES) {
      std::shared_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));
      m_templates[mol->getNumAtoms()].push_back(mol);
    }
  }

 private:
  CoordinateTemplates() { loadDefaultTemplates(); }
  CoordinateTemplates(const CoordinateTemplates&) = delete;
  CoordinateTemplates& operator=(const CoordinateTemplates&) = delete;

  void clearTemplates() {
    for (auto& [atom_cout, romols] : m_templates) {
      romols.clear();
    }
    m_templates.clear();
  }

  ~CoordinateTemplates() { clearTemplates(); }

  void loadTemplatesFromPath(
      const std::string& templatePath,
      std::unordered_map<
          unsigned int, std::vector<std::shared_ptr<RDKit::ROMol>>>& templates);

  std::unordered_map<unsigned int, std::vector<std::shared_ptr<RDKit::ROMol>>>
      m_templates;
};
}  // namespace RDDepict
