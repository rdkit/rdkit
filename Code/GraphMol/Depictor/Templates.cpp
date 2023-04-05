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

namespace RDDepict {
bool CoordinateTemplates::isValidTemplate(RDKit::ROMol& mol, const std::string& smiles) {
    // template must have 2D coordinates
    if (mol.getNumConformers() == 0) {
        BOOST_LOG(rdWarningLog) << "WARNING: Template is missing coordinates: "
                                << smiles << "\n";
        return false;
    }
    if (mol.getConformer().is3D()) {
        BOOST_LOG(rdWarningLog) << "WARNING: Template has 3D coordinates, 2D "
                                << "coordinates required: " << smiles << "\n";
        return false;
    }

    // Make sure this template is a single ring system (spiro'd ring systems are OK).
    // We can check this by ensuring that every atom is in a ring and that there
    // is only one connected component in the molecular graph.
    if (RDKit::MolOps::getMolFrags(mol).size() != 1) {
        BOOST_LOG(rdWarningLog) << "WARNING: Template consists of multiple fragments, "
                                << "single fragment required: " << smiles << "\n";
        return false;
    }

    // Use symmetrizeSSSR since that is what is used in coordinate generation
    bool includeDativeBonds = true;
    RDKit::MolOps::symmetrizeSSSR(mol, includeDativeBonds);
    auto ri = mol.getRingInfo();
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        if (!ri->numAtomRings(i)) {
            BOOST_LOG(rdWarningLog) << "WARNING: Template is not a ring system: "
                                    << smiles << "\n";
            return false;
        }
    }
    return true;
}


bool CoordinateTemplates::setRingSystemTemplates(const std::string& template_dir) {
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
        if (isValidTemplate(*mol, line)) {
            templates[mol->getNumAtoms()].push_back(mol);
        } else {
            success = false;
            break;
        }
    }
    cxsmiles.close();
    if (success) {
        clearTemplates();
        m_templates = std::move(templates);
    } else {
        BOOST_LOG(rdWarningLog)
            << "WARNING: Could not load templates from " << template_dir << std::endl;
    }
    return success;
}
} // namespace RDDepict