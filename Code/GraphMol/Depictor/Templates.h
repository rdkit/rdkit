
//
//  Copyright ...

#include <RDGeneral/export.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>

#include <iostream>
#include <fstream>
#include <unordered_map>

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

  bool hasTemplateOfSize(unsigned int atom_count) {
    if (m_templates.find(atom_count) != m_templates.end()) {
        return true;
    }
    return false;
  }

  const std::vector<RDKit::ROMol*>& getMatchingTemplates(unsigned int atom_count) {
    return m_templates[atom_count];
  }

 private:
  CoordinateTemplates() {
    // load templates into m_templates map by atom count
    std::string rdbase = getenv("RDBASE");
    std::string fpath = rdbase + "/Code/GraphMol/Depictor/templates.smi";
    std::ifstream cxsmiles_templates (fpath);
    std::string smiles;
    while(std::getline(cxsmiles_templates, smiles)) {
        auto mol = RDKit::SmilesToMol(smiles);
        RDKit::MolOps::symmetrizeSSSR(*mol);
        unsigned int atom_count = mol->getNumAtoms();
        if (m_templates.find(atom_count) == m_templates.end()) {
            m_templates[atom_count] = {mol};
        } else {
            m_templates[atom_count].push_back(mol);
        }
    }
    cxsmiles_templates.close();
  }

  CoordinateTemplates(const CoordinateTemplates&) = delete;
  CoordinateTemplates& operator= (const CoordinateTemplates&) = delete;
  ~CoordinateTemplates() {
    for (auto& [atom_cout, romols] : m_templates) {
        for (auto& romol : romols) {
            delete romol;
        }
    }
  }

  std::unordered_map<unsigned int, std::vector<RDKit::ROMol*>> m_templates;
};