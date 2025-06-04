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

#include "MonomerDatabase.h"
#include "MonomerMol.h"

namespace RDKit
{

namespace
{
namespace fs = boost::filesystem;

const std::string SGROUP_PROP_CLASS{"CLASS"};
const std::string SGROUP_PROP_LABEL{"LABEL"};
const std::string SGROUP_PROP_RESNUM{"ID"};

using AttachmentMap = std::map<std::pair<unsigned int, unsigned int>,
                               std::pair<unsigned int, unsigned int>>;

static const std::string ATOM_PDB_NAME_PROP{"pdbName"};

const std::unordered_map<std::string, std::string> three_character_codes({
    {"A", "ALA"}, // Alanine
    {"R", "ARG"}, // Arginine
    {"N", "ASN"}, // Asparagine
    {"D", "ASP"}, // Aspartic
    {"C", "CYS"}, // Cysteine
    {"Q", "GLN"}, // Glutamine
    {"E", "GLU"}, // Glutamic
    {"G", "GLY"}, // Glycine
    {"H", "HIS"}, // Histidine
    {"I", "ILE"}, // Isoleucine
    {"L", "LEU"}, // Leucine
    {"K", "LYS"}, // Lysine
    {"M", "MET"}, // Methionine
    {"F", "PHE"}, // Phenylalanine
    {"P", "PRO"}, // Proline
    {"S", "SER"}, // Serine
    {"T", "THR"}, // Threonine
    {"W", "TRP"}, // Tryptophan
    {"Y", "TYR"}, // Tyrosine
    {"V", "VAL"}, // Valine
});

static const std::map<unsigned int, std::string> BIOVIA_ATTCHPT_MAP = {
    {1, "Al"}, // Backbone attachment point
    {2, "Br"}, // Backbone attachment point
    {3, "Cx"}, // Sidechain attachment point
};

std::pair<unsigned int, unsigned int> getAttchpts(const std::string& linkage)
{
    // in form RX-RY, returns {X, Y}
    auto dash = linkage.find('-');
    if (dash == std::string::npos) {
        throw std::runtime_error(
            fmt::format("Invalid linkage format: {}", linkage));
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
                        fmt::format("Invalid attachment point at index {}",
                                    atom->getIdx()));
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
                    ChainType chain_type, unsigned int current_residue,
                    MonomerDatabase& db)
{
    std::string residue_name =
        (chain_type == ChainType::PEPTIDE) ? "UNK" : "UNL";

    auto pdb_code = db.getPdbCode(monomer_label, chain_type);
    if (pdb_code) {
        residue_name = *pdb_code;
    } else if (three_character_codes.find(monomer_label) !=
               three_character_codes.end()) {
        residue_name = three_character_codes.at(monomer_label);
    }

    // Set PDB residue info on the atoms
    for (auto atom : new_monomer.atoms()) {
        auto* res_info = new RDKit::AtomPDBResidueInfo();

        std::string chain_id_str(1, chain_id);
        res_info->setChainId(chain_id_str);
        res_info->setResidueNumber(residue_number);
        res_info->setResidueName(residue_name);
        // to be consistent with the rdkit adapter and RDKit's own PDB writer
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

ChainType getChainType(std::string_view polymer_id)
{
    if (polymer_id.find("PEPTIDE") == 0) {
        return ChainType::PEPTIDE;
    } else if (polymer_id.find("RNA") == 0) {
        // HELM labels both DNA and RNA as RNA
        return ChainType::RNA;
    } else if (polymer_id.find("CHEM") == 0) {
        return ChainType::CHEM;
    } else {
        throw std::out_of_range(fmt::format(
            "Invalid polymer id: {}. Must be one of PEPTIDE, RNA, CHEM",
            polymer_id));
    }
}

AttachmentMap addPolymer(RDKit::RWMol& atomistic_mol,
                         const RDKit::RWMol& monomer_mol,
                         const std::string& polymer_id,
                         std::vector<unsigned int>& remove_atoms, char chain_id,
                         unsigned int& total_residue_count)
{
    // Maps residue number and attachment point number to the atom index in
    // atomistic_mol that should be attached to and the atom index of the rgroup
    // that should later be removed
    AttachmentMap attachment_point_map;

    auto chain = get_polymer(monomer_mol, polymer_id);
    auto chain_type = getChainType(polymer_id);
    bool sanitize = false;

    auto path = getMonomerDbPath();
    if (!path.has_value()) {
        // This shouldn't happen since DEFAULT_MONOMER_DB_PATH_ENV_VAR points to
        // mmshare/data/helm/core_monomerlib.db and that should always exist
        throw std::runtime_error(fmt::format(
            "Could not find monomer database, try setting env variable {}",
            CUSTOM_MONOMER_DB_PATH_ENV_VAR));
    }
    MonomerDatabase db(*path);

    // Add the monomers to the atomistic mol
    for (const auto monomer_idx : chain.atoms) {
        auto monomer = monomer_mol.getAtomWithIdx(monomer_idx);
        auto monomer_label = monomer->getProp<std::string>(ATOM_LABEL);

        std::string smiles;
        if (monomer->getProp<bool>(SMILES_MONOMER)) {
            smiles = monomer_label;
        } else {
            auto monomer_smiles =
                db.getMonomerSmiles(monomer_label, chain_type);
            if (!monomer_smiles) {
                throw std::out_of_range(fmt::format(
                    "Peptide Monomer {} not found in Monomer database",
                    monomer_label));
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
                       chain_type, total_residue_count, db);
        ++total_residue_count;
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
            throw std::runtime_error(fmt::format(
                "Invalid linkage {} between monomers {} and {}",
                bond->getProp<std::string>(LINKAGE), from_res, to_res));
        }

        auto [core_aid1, attachment_point1] =
            attachment_point_map.at({from_res, from_rgroup});
        auto [core_aid2, attachment_point2] =
            attachment_point_map.at({to_res, to_rgroup});

        auto atomistic_bond_idx =
            atomistic_mol.addBond(core_aid1, core_aid2, bond->getBondType()) -
            1;
        remove_atoms.push_back(attachment_point1);
        remove_atoms.push_back(attachment_point2);
    }

    return attachment_point_map;
}
} // namespace

boost::shared_ptr<RDKit::RWMol> mmToAtomistic(const RDKit::ROMol& monomer_mol)
{
    auto atomistic_mol = boost::make_shared<RDKit::RWMol>();

    // Map to track Polymer ID -> attachment point map
    std::unordered_map<std::string, AttachmentMap> polymer_attachment_points;
    std::vector<unsigned int> remove_atoms;
    unsigned int total_residue_count = 1; // 1-based index to label SUP groups
    char chain_id = 'A';
    for (const auto& polymer_id : getPolymerIds(monomer_mol)) {
        polymer_attachment_points[polymer_id] =
            addPolymer(*atomistic_mol, monomer_mol, polymer_id, remove_atoms,
                       chain_id, total_residue_count);
        ++chain_id;
    }

    // Add bonds from interpolymer connections
    for (const auto bnd : monomer_mol.bonds()) {
        auto begin_atom = bnd->getBeginAtom();
        auto end_atom = bnd->getEndAtom();
        if (getPolymerId(begin_atom) == getPolymerId(end_atom)) {
            continue;
        }
        auto begin_res = getResidueNumber(begin_atom);
        auto end_res = getResidueNumber(end_atom);
        auto [from_rgroup, to_rgroup] =
            getAttchpts(bnd->getProp<std::string>(LINKAGE));

        const auto& begin_attachment_points =
            polymer_attachment_points.at(getPolymerId(begin_atom));
        const auto& end_attachment_points =
            polymer_attachment_points.at(getPolymerId(end_atom));

        if (begin_attachment_points.find({begin_res, from_rgroup}) ==
                begin_attachment_points.end() ||
            end_attachment_points.find({end_res, to_rgroup}) ==
                end_attachment_points.end()) {
            // One of these attachment points is not present
            std::string error_msg =
                "Invalid linkage " + bnd->getProp<std::string>(LINKAGE) +
                " between monomers " + std::to_string(begin_atom->getIdx()) + " and " +
                std::to_string(end_atom->getIdx());
            throw std::runtime_error(error_msg);
        }

        auto [core_atom1, attachment_point1] =
            begin_attachment_points.at({begin_res, from_rgroup});
        auto [core_atom2, attachment_point2] =
            end_attachment_points.at({end_res, to_rgroup});
        atomistic_mol->addBond(core_atom1, core_atom2, bnd->getBondType());
        remove_atoms.push_back(attachment_point1);
        remove_atoms.push_back(attachment_point2);
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