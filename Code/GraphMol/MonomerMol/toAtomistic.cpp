//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "Conversions.h"

#include <queue>
#include <unordered_map>
#include <memory>

#include <boost/filesystem/path.hpp>

#include <GraphMol/RWMol.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "MonomerLibrary.h"
#include "MonomerMol.h"

namespace RDKit
{

namespace
{

using AttachmentMap = std::map<std::pair<unsigned int, unsigned int>,
                               std::pair<unsigned int, unsigned int>>;

static const std::string ATOM_PDB_NAME_PROP{"pdbName"};

std::pair<unsigned int, unsigned int> getAttchpts(const std::string& linkage)
{
    // in form RX-RY, returns {X, Y}
    auto dash = linkage.find('-');
    if (dash == std::string::npos) {
        throw std::runtime_error(
            "Invalid linkage format: " + linkage);
    }
    return {std::stoi(linkage.substr(1, dash - 1)),
            std::stoi(linkage.substr(dash + 2))};
}

void fillAttachmentPointMap(const RDKit::ROMol& new_monomer,
                            AttachmentMap& attachment_points,
                            unsigned int residue_num, unsigned int old_mol_size)
{
    for (const auto& atom : new_monomer.atoms()) {
        unsigned int map_num;
        if (atom->getPropIfPresent(RDKit::common_properties::molAtomMapNumber,
                                   map_num)) {
            for (auto bnd : new_monomer.atomBonds(atom)) {
                // Should be exactly one iteration -- an attachment point can
                // only have one neighbor
                if (attachment_points.find({residue_num, map_num}) !=
                    attachment_points.end()) {
                    throw std::runtime_error(
                        "Invalid attachment point at index " +
                        std::to_string(atom->getIdx()));
                }
                auto atom_to_bond_to =
                    old_mol_size + bnd->getOtherAtomIdx(atom->getIdx());
                auto atom_to_remove = old_mol_size + atom->getIdx();
                attachment_points[{residue_num, map_num}] =
                    std::make_pair(atom_to_bond_to, atom_to_remove);
            }
        }
    }
}

void setResidueInfo(RDKit::RWMol& new_monomer, const std::string& monomer_label,
                    unsigned int residue_number, char chain_id,
                    std::string_view monomer_type,
                    const MonomerLibrary& lib)
{
    std::string residue_name =
        (monomer_type == "PEPTIDE") ? "UNK" : "UNL";

    auto pdb_code = lib.getPdbCode(monomer_label, monomer_type);
    if (pdb_code) {
        residue_name = *pdb_code;
    }

    // Set PDB residue info on the atoms
    for (auto atom : new_monomer.atoms()) {
        auto* res_info = new RDKit::AtomPDBResidueInfo();

        std::string chain_id_str(1, chain_id);
        res_info->setChainId(chain_id_str);
        res_info->setResidueNumber(residue_number);
        res_info->setResidueName(residue_name);
        // to be consistent with RDKit's PDB writer
        res_info->setInsertionCode(" ");

        std::string pdb_name;
        if (atom->getPropIfPresent(ATOM_PDB_NAME_PROP, pdb_name)) {
            res_info->setName(pdb_name);
            // Don't keep property on output molecule
            atom->clearProp(ATOM_PDB_NAME_PROP);
        }

        atom->setMonomerInfo(res_info);
    }
}

std::string getMonomerType(std::string_view polymer_id)
{
    if (polymer_id.find("PEPTIDE") == 0) {
        return "PEPTIDE";
    } else if (polymer_id.find("RNA") == 0) {
        // HELM labels both DNA and RNA as RNA
        return "RNA";
    } else if (polymer_id.find("CHEM") == 0) {
        return "CHEM";
    } else {
        throw std::out_of_range(
            "Invalid polymer id: " + std::string(polymer_id) +
            ". Must be one of PEPTIDE, RNA, CHEM");
    }
}

AttachmentMap addPolymer(RDKit::RWMol& atomistic_mol,
                         const RDKit::RWMol& monomer_mol,
                         const std::string& polymer_id,
                         std::vector<unsigned int>& remove_atoms, char chain_id,
                         const MonomerLibrary& lib)
{
    // Maps residue number and attachment point number to the atom index in
    // atomistic_mol that should be attached to and the atom index of the rgroup
    // that should later be removed
    AttachmentMap attachment_point_map;

    auto chain = getPolymer(monomer_mol, polymer_id);
    auto monomer_type = getMonomerType(polymer_id);
    bool sanitize = false;

    // Add the monomers to the atomistic mol
    for (const auto monomer_idx : chain.atoms) {
        const auto monomer = monomer_mol.getAtomWithIdx(monomer_idx);
        auto monomer_label = monomer->getProp<std::string>(ATOM_LABEL);

        std::string smiles;
        if (monomer->getProp<bool>(SMILES_MONOMER)) {
            smiles = monomer_label;
        } else {
            auto monomer_smiles =
                lib.getMonomerSmiles(monomer_label, monomer_type);
            if (!monomer_smiles) {
                throw std::out_of_range(
                    "Monomer " + monomer_label + " not found in MonomerLibrary");
            }
            smiles = *monomer_smiles;
        }

        std::unique_ptr<RDKit::RWMol> new_monomer(
            RDKit::SmilesToMol(smiles, 0, sanitize));

        if (!new_monomer) {
            // FIXME: I think this is an issue with the HELM parser, see
            // SHARED-11457
            new_monomer.reset(
                RDKit::SmilesToMol("[" + smiles + "]", 0, sanitize));
        }

        if (monomer->getProp<bool>(SMILES_MONOMER)) {
            // SMILES monomers may be in rgroup form like
            // *N[C@H](C(=O)O)S* |$_R1;;;;;;;_R3$| or use atom map numbers like
            // [*:1]N[C@H](C(=O)O)S[*:3]. Translate the RGroup to atom map
            // numbers
            for (auto atom : new_monomer->atoms()) {
                std::string rgroup_label;
                if (atom->getPropIfPresent(RDKit::common_properties::atomLabel,
                                           rgroup_label) &&
                    rgroup_label.find("_R") == 0) {
                    auto rgroup_num = std::stoi(rgroup_label.substr(2));
                    atom->setAtomMapNum(rgroup_num);
                    atom->clearProp(RDKit::common_properties::atomLabel);
                }
            }
        }

        auto residue_number = getResidueNumber(monomer);
        fillAttachmentPointMap(*new_monomer, attachment_point_map,
                               residue_number, atomistic_mol.getNumAtoms());
        setResidueInfo(*new_monomer, monomer_label, residue_number, chain_id,
                       monomer_type, lib);
        atomistic_mol.insertMol(*new_monomer);
    }

    // Add the bonds between monomers and mark the replaced rgroups to be
    // removed
    for (const auto bond_idx : chain.bonds) {
        auto bond = monomer_mol.getBondWithIdx(bond_idx);
        auto [from_rgroup, to_rgroup] =
            getAttchpts(bond->getProp<std::string>(LINKAGE));
        auto from_res = getResidueNumber(bond->getBeginAtom());
        auto to_res = getResidueNumber(bond->getEndAtom());

        if (attachment_point_map.find({from_res, from_rgroup}) ==
                attachment_point_map.end() ||
            attachment_point_map.find({to_res, to_rgroup}) ==
                attachment_point_map.end()) {
            // One of these attachment points is not present
            throw std::runtime_error(
                "Invalid linkage " + bond->getProp<std::string>(LINKAGE) +
                " between monomers " + std::to_string(from_res) + " and " +
                std::to_string(to_res));
        }

        auto [core_aid1, attachment_point1] =
            attachment_point_map.at({from_res, from_rgroup});
        auto [core_aid2, attachment_point2] =
            attachment_point_map.at({to_res, to_rgroup});

        // Dative bonds are only relevant at the monomer mol level for
        // directionality; at the atomistic level they are single bonds
        auto bond_type = bond->getBondType();
        if (bond_type == RDKit::Bond::DATIVE) {
            bond_type = RDKit::Bond::SINGLE;
        }
        [[maybe_unused]] auto atomistic_bond_idx =
            atomistic_mol.addBond(core_aid1, core_aid2, bond_type) - 1;

        remove_atoms.push_back(attachment_point1);
        remove_atoms.push_back(attachment_point2);
    }

    return attachment_point_map;
}

bool isMonomerAtom(const RDKit::Atom* atom)
{
    return atom->hasProp(ATOM_LABEL);
}

std::pair<unsigned int, unsigned int> getMonomerAtomAttchpt(
    const std::string& linkage)
{
    // Parse linkage in form "RX-" (monomer to atom), returns {X, 0}
    // The trailing dash indicates this is a monomer-to-atom connection
    if (linkage.size() >= 3 && linkage[0] == 'R' &&
        linkage.back() == '-') {
        auto rgroup = std::stoi(linkage.substr(1, linkage.size() - 2));
        return {rgroup, 0};
    }
    throw std::runtime_error("Invalid monomer-to-atom linkage format: " +
                             linkage);
}

}  // namespace

std::unique_ptr<RDKit::RWMol> toAtomistic(const RDKit::ROMol& monomer_mol)
{
    auto atomistic_mol = std::make_unique<RDKit::RWMol>();

    // Map from monomer_mol atom index to atomistic_mol atom index for
    // non-monomer atoms
    std::unordered_map<unsigned int, unsigned int> regular_atom_map;

    // First, copy all non-monomer (regular) atoms to the atomistic mol
    for (const auto atom : monomer_mol.atoms()) {
        if (!isMonomerAtom(atom)) {
            auto new_atom = new RDKit::Atom(*atom);
            auto new_idx = atomistic_mol->addAtom(new_atom, true, true);
            regular_atom_map[atom->getIdx()] = new_idx;
        }
    }

    // Add bonds between regular atoms
    for (const auto bond : monomer_mol.bonds()) {
        auto begin_atom = bond->getBeginAtom();
        auto end_atom = bond->getEndAtom();
        if (!isMonomerAtom(begin_atom) && !isMonomerAtom(end_atom)) {
            auto new_begin = regular_atom_map.at(begin_atom->getIdx());
            auto new_end = regular_atom_map.at(end_atom->getIdx());
            atomistic_mol->addBond(new_begin, new_end, bond->getBondType());
        }
    }

    // Get the effective library - try to get from MonomerMol, fall back to global
    std::shared_ptr<MonomerLibrary> lib;
    if (const auto* mm = dynamic_cast<const MonomerMol*>(&monomer_mol)) {
        lib = mm->getEffectiveLibrary();
    } else {
        lib = MonomerLibrary::getGlobalLibrary();
    }

    // Map to track Polymer ID -> attachment point map
    std::unordered_map<std::string, AttachmentMap> polymer_attachment_points;
    std::vector<unsigned int> remove_atoms;
    char chain_id = 'A';
    for (const auto& polymer_id : getPolymerIds(monomer_mol)) {
        polymer_attachment_points[polymer_id] =
            addPolymer(*atomistic_mol, monomer_mol, polymer_id, remove_atoms,
                       chain_id, *lib);
        ++chain_id;
    }

    // Add bonds from interpolymer connections (monomer-to-monomer across chains)
    // and monomer-to-atom connections
    for (const auto bnd : monomer_mol.bonds()) {
        auto begin_atom = bnd->getBeginAtom();
        auto end_atom = bnd->getEndAtom();

        bool begin_is_monomer = isMonomerAtom(begin_atom);
        bool end_is_monomer = isMonomerAtom(end_atom);

        // Skip bonds between regular atoms (already handled above)
        if (!begin_is_monomer && !end_is_monomer) {
            continue;
        }

        // Skip intra-polymer monomer-to-monomer bonds (handled in addPolymer)
        if (begin_is_monomer && end_is_monomer &&
            getPolymerId(begin_atom) == getPolymerId(end_atom)) {
            continue;
        }

        auto linkage = bnd->getProp<std::string>(LINKAGE);

        // Dative bonds are only relevant at the monomer mol level for
        // directionality; at the atomistic level they are single bonds
        auto bond_type = bnd->getBondType();
        if (bond_type == RDKit::Bond::DATIVE) {
            bond_type = RDKit::Bond::SINGLE;
        }

        if (begin_is_monomer && end_is_monomer) {
            // Inter-polymer monomer-to-monomer connection
            auto begin_res = getResidueNumber(begin_atom);
            auto end_res = getResidueNumber(end_atom);
            auto [from_rgroup, to_rgroup] = getAttchpts(linkage);

            const auto& begin_attachment_points =
                polymer_attachment_points.at(getPolymerId(begin_atom));
            const auto& end_attachment_points =
                polymer_attachment_points.at(getPolymerId(end_atom));

            if (begin_attachment_points.find({begin_res, from_rgroup}) ==
                    begin_attachment_points.end() ||
                end_attachment_points.find({end_res, to_rgroup}) ==
                    end_attachment_points.end()) {
                throw std::runtime_error(
                    "Invalid linkage " + linkage + " between monomers " +
                    std::to_string(begin_atom->getIdx()) + " and " +
                    std::to_string(end_atom->getIdx()));
            }

            auto [core_atom1, attachment_point1] =
                begin_attachment_points.at({begin_res, from_rgroup});
            auto [core_atom2, attachment_point2] =
                end_attachment_points.at({end_res, to_rgroup});

            atomistic_mol->addBond(core_atom1, core_atom2, bond_type);
            remove_atoms.push_back(attachment_point1);
            remove_atoms.push_back(attachment_point2);
        } else {
            // Monomer-to-atom connection (hybrid molecule)
            auto [rgroup, _] = getMonomerAtomAttchpt(linkage);

            const RDKit::Atom* monomer =
                begin_is_monomer ? begin_atom : end_atom;
            const RDKit::Atom* regular_atom =
                begin_is_monomer ? end_atom : begin_atom;

            auto monomer_res = getResidueNumber(monomer);
            const auto& attachment_points =
                polymer_attachment_points.at(getPolymerId(monomer));

            if (attachment_points.find({monomer_res, rgroup}) ==
                attachment_points.end()) {
                throw std::runtime_error(
                    "Invalid attachment point R" + std::to_string(rgroup) +
                    " for monomer " + std::to_string(monomer->getIdx()));
            }

            auto [core_atom, attachment_point] =
                attachment_points.at({monomer_res, rgroup});
            auto regular_atom_idx =
                regular_atom_map.at(regular_atom->getIdx());

            atomistic_mol->addBond(core_atom, regular_atom_idx, bond_type);
            remove_atoms.push_back(attachment_point);
        }
    }

    // Remove atoms that represented attachment points and dummy atoms
    atomistic_mol->beginBatchEdit();
    for (auto at_idx : remove_atoms) {
        atomistic_mol->removeAtom(at_idx);
    }
    for (auto at : atomistic_mol->atoms()) {
        if (at->getAtomicNum() == 0) {
            atomistic_mol->removeAtom(at->getIdx());
        }
    }
    atomistic_mol->commitBatchEdit();

    // Let sanitization errors bubble up for now -- it means we did something
    // wrong
    RDKit::MolOps::sanitizeMol(*atomistic_mol);

    // Remove graph hydrogens where some RGroups previously were. This will turn
    // H - NH to NH2, etc
    RDKit::MolOps::removeHs(*atomistic_mol);

    // Remove atom map numbers -- anything remaining is a capping group and the
    // map numbers are no longer meaningful
    for (auto at : atomistic_mol->atoms()) {
        at->setAtomMapNum(0);
    }

    return atomistic_mol;
}

} // namespace RDKit