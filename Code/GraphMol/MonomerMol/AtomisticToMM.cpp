#include "Conversions.h"

#include <queue>
#include <unordered_map>
#include <memory>

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

const std::string ATTACH_FROM{"attachFrom"};
const std::string MONOMER_IDX1{"monomerIndex1"};
const std::string MONOMER_IDX2{"monomerIndex2"};
const std::string MONOMER_MAP_NUM{"monomerMapNumber"};
const std::string REFERENCE_IDX{"referenceIndex"};

constexpr int SIDECHAIN_IDX = 2;
constexpr int MIN_ATTCHPTS = 2;
constexpr auto NO_ATTACHMENT = std::numeric_limits<unsigned int>::max();

// std::map to allow sequential/ordered iteration
using ChainsAndResidues =
    std::map<std::string,
             std::map<std::pair<int, std::string>, std::vector<unsigned int>>>;

// attachment points 1 and 2 are backbone attachment points
// and 3 is the side chain attachment point
const std::string GENERIC_AMINO_ACID_QUERY =
    "[NX3,NX4+:1][CX4H]([*:3])[CX3](=[OX1])[O,N:2]";
const std::string GLYCINE_AMINO_ACID_QUERY =
    "[N:1][CX4H2][CX3](=[OX1])[O,N:2]"; // no side chain

// SMILES monomer to monomer abbreviation mapping
// temporary for now, for proof of concept
// got most of these from pubchem. include the version with N and O
const std::unordered_map<std::string, std::string> amino_acids = {
    {"CC(N)C(=O)O", "A"},                 // Alanine (Ala)
    {"NC(N)=NCCCC(N)C(=O)O", "R"},        // Arginine (Arg)
    {"NC(=O)CC(N)C(=O)O", "N"},           // Asparagine (Asn)
    {"NC(CC(=O)O)C(=O)O", "D"},           // Aspartic acid (Asp)
    {"NC(CS)C(=O)O", "C"},                // Cysteine (Cys)
    {"NC(=O)CCC(N)C(=O)O", "Q"},          // Glutamine (Gln)
    {"NC(CCC(=O)O)C(=O)O", "E"},          // Glutamic acid (Glu)
    {"NCC(=O)O", "G"},                    // Glycine (Gly)
    {"NC(Cc1cnc[nH]1)C(=O)O", "H"},       // Histidine (His)
    {"CCC(C)C(N)C(=O)O", "I"},            // Isoleucine (Ile)
    {"CC(C)CC(N)C(=O)O", "L"},            // Leucine (Leu)
    {"NCCCCC(N)C(=O)O", "K"},             // Lysine (Lys)
    {"CSCCC(N)C(=O)O", "M"},              // Methionine (Met)
    {"NC(Cc1ccccc1)C(=O)O", "F"},         // Phenylalanine (Phe)
    {"O=C(O)C1CCCN1", "P"},               // Proline (Pro)
    {"NC(CO)C(=O)O", "S"},                // Serine (Ser)
    {"CC(O)C(N)C(=O)O", "T"},             // Threonine (Thr)
    {"NC(Cc1c[nH]c2ccccc12)C(=O)O", "W"}, // Tryptophan (Trp)
    {"NC(Cc1ccc(O)cc1)C(=O)O", "Y"},      // Tyrosine (Tyr)
    {"CC(C)C(N)C(=O)O", "V"},             // Valine (Val)
    {"CC(N)C(N)=O", "A"},
    {"NC(=O)C(N)CCCN=C(N)N", "R"}, // arginine, pubchem version
    {"N=C(N)NCCCC(N)C(N)=O", "R"}, // arginine, from HELM paper examples
                                   // (different double bond placement)
    {"NC(=O)CC(N)C(N)=O", "N"},
    {"NC(=O)C(N)CC(=O)O", "D"},
    {"NC(=O)C(N)CS", "C"},
    {"NC(=O)CCC(N)C(N)=O", "Q"},
    {"NC(=O)C(N)CCC(=O)O", "E"},
    {"NCC(N)=O", "G"},
    {"NC(=O)C(N)Cc1cnc[nH]1", "H"},
    {"CCC(C)C(N)C(N)=O", "I"},
    {"CC(C)CC(N)C(N)=O", "L"},
    {"NCCCCC(N)C(N)=O", "K"},
    {"CSCCC(N)C(N)=O", "M"},
    {"NC(=O)C(N)Cc1ccccc1", "F"},
    {"NC(=O)C1CCCN1", "P"},
    {"NC(=O)C(N)CO", "S"},
    {"CC(O)C(N)C(N)=O", "T"},
    {"NC(=O)C(N)Cc1c[nH]c2ccccc12", "W"},
    {"NC(=O)C(N)Cc1ccc(O)cc1", "Y"},
    {"CC(C)C(N)C(N)=O", "V"}};

static const std::unordered_map<std::string, ChainType> BIOVIA_CHAIN_TYPE_MAP =
    {{"AA", ChainType::PEPTIDE}};

struct Linkage {
    unsigned int monomer_idx1;
    unsigned int monomer_idx2;
    unsigned int attach_from;
    unsigned int attach_to;
    std::string to_string() const
    {
        return "R" + std::to_string(attach_from) + "-R" + std::to_string(attach_to);
    }
};

bool alreadyMatched(const RDKit::ROMol& mol,
                    const std::span<const unsigned int> ids)
{
    // Make sure this match hasn't already been accounted for by a previous
    // match
    for (auto id : ids) {
        auto at = mol.getAtomWithIdx(id);
        unsigned int attch;
        if (at->hasProp(MONOMER_IDX1) &&
            at->getPropIfPresent(ATTACH_FROM, attch) &&
            attch == NO_ATTACHMENT) {
            return false;
        }
    }
    return true;
}

/*
 * Function that takes the SMARTS query and atomistic molecule and adds the atom
 * indices of the matches to the monomers vector
 *
 */
void addMatchesToMonomers(
    const std::string& smarts_query, RDKit::ROMol& atomistic_mol,
    std::vector<std::vector<int>>& monomers,
    std::unordered_map<unsigned int, unsigned int>& attch_pts,
    std::vector<Linkage>& linkages)
{
    std::unique_ptr<RDKit::RWMol> query(RDKit::SmartsToMol(smarts_query));
    // maps SMARTS query index to attachment point #
    std::vector<unsigned int> attch_map(query->getNumAtoms(), NO_ATTACHMENT);
    for (const auto atom : query->atoms()) {
        if (atom->hasProp(RDKit::common_properties::molAtomMapNumber)) {
            attch_map[atom->getIdx()] = atom->getProp<unsigned int>(
                RDKit::common_properties::molAtomMapNumber);
        }
    }

    // Set a final function check that ensures the entire match has not
    // already been accounted for by a previous SMARTS search.
    RDKit::SubstructMatchParameters params;
    params.useChirality = false;
    params.extraFinalCheck = alreadyMatched;
    auto matches = RDKit::SubstructMatch(atomistic_mol, *query, params);

    for (const auto& match : matches) {
        std::vector<int> monomer;
        auto monomer_idx = monomers.size();
        for (const auto& [query_idx, atom_idx] : match) {
            auto atom = atomistic_mol.getAtomWithIdx(atom_idx);
            if (atom->hasProp(MONOMER_IDX1) && atom->hasProp(MONOMER_IDX2)) {
                // This shouldn't happen, sanity check
                throw std::runtime_error(
                    "Atom " + std::to_string(atom->getIdx()) + " belongs to more than 2 monomers");
            }

            if (atom->hasProp(MONOMER_IDX1)) {
                atom->setProp<unsigned int>(MONOMER_IDX2, monomer_idx);
                // ATTACH_FROM property should be set when MONOMER_IDX1 is set
                auto attach_from_idx = atom->getProp<unsigned int>(ATTACH_FROM);

                // Make sure linkages are directionally correct, so R2-R1 or
                // R3-R1 instead of R1-R2 or R1-R3
                if (attach_from_idx >= attch_map[query_idx]) {
                    Linkage link = {atom->getProp<unsigned int>(MONOMER_IDX1),
                                    static_cast<unsigned int>(monomer_idx),
                                    attach_from_idx, attch_map[query_idx]};
                    linkages.push_back(link);
                } else {
                    Linkage link = {static_cast<unsigned int>(monomer_idx),
                                    atom->getProp<unsigned int>(MONOMER_IDX1),
                                    attch_map[query_idx], attach_from_idx};
                    linkages.push_back(link);
                }
            } else {
                atom->setProp<unsigned int>(MONOMER_IDX1, monomer_idx);
                atom->setProp(ATTACH_FROM, attch_map[query_idx]);
            }
            monomer.push_back(atom_idx);

            // if there is a side chain, the attachment point will be at the
            // SIDECHAIN_IDX and will be indicated by the presence of the atom
            // map number. For now, assume there is a single side chain per
            // monomer
            if (query_idx == SIDECHAIN_IDX &&
                query->getAtomWithIdx(query_idx)->hasProp(
                    RDKit::common_properties::molAtomMapNumber)) {
                attch_pts[monomers.size()] = atom_idx;
            }
        }
        monomers.push_back(monomer);
    }
}

void addSidechainToMonomer(const RDKit::ROMol& atomistic_mol,
                           std::vector<int>& monomer, unsigned int monomer_idx,
                           unsigned int attch_at_idx)
{

    // BFS but use MONOMER_IDX as visited marker
    std::queue<unsigned int> q;
    q.push(attch_at_idx);
    while (!q.empty()) {
        auto at_idx = q.front();
        q.pop();
        auto at = atomistic_mol.getAtomWithIdx(at_idx);
        if (!at->hasProp(MONOMER_IDX1)) {
            at->setProp<unsigned int>(MONOMER_IDX1, monomer_idx);
            monomer.push_back(at_idx);
        }
        for (const auto& nbr : atomistic_mol.atomNeighbors(at)) {
            if (!nbr->hasProp(MONOMER_IDX1)) {
                q.push(nbr->getIdx());
            }
        }
    }
}

/*
 * Break an atomistic molecule into monomers
 *
 * Every atom should belong to either 1 or 2 monomers. If an atom belongs to 2
 * monomers, it represents a connection between the two monomers.
 *
 * Input ROMol is labeled as follows:
 * - firstMonomerIdx: index of the first monomer given atom belongs to
 * - secondMonomerIdx: index of the second monomer given atom belongs to (this
 * is optional, means there is a bond between two monomers)
 *
 * Returns a list monomer index sets, which represents the atom indices of each
 * monomer.
 */
void identifyMonomers(RDKit::ROMol& atomistic_mol,
                      std::vector<std::vector<int>>& monomers,
                      std::vector<Linkage>& linkages)
{
    // Approach for identifying monomers:
    // 1. Find all matches with SMARTS queries for amino acids (TODO: Nucleic
    // acids & CHEM)
    // 2. Add side chains to generic matches based on attachment points
    // 3. Identify and group any remaining atoms into 'unclassified' monomers.
    // These are grouped
    //    by connectivity.
    std::unordered_map<unsigned int, unsigned int> attch_pts;
    addMatchesToMonomers(GENERIC_AMINO_ACID_QUERY, atomistic_mol, monomers,
                         attch_pts, linkages);
    addMatchesToMonomers(GLYCINE_AMINO_ACID_QUERY, atomistic_mol, monomers,
                         attch_pts, linkages);
    // TODO: nucleic acids and CHEM monomers

    // now, add sidechains onto each monomer
    for (size_t monomer_idx = 0; monomer_idx < monomers.size(); ++monomer_idx) {
        if (attch_pts.find(monomer_idx) != attch_pts.end()) {
            // there is a sidechain to add!
            addSidechainToMonomer(atomistic_mol, monomers[monomer_idx],
                                  monomer_idx, attch_pts[monomer_idx]);
        }
    }
}

void neutralizeAtoms(RDKit::ROMol& mol)
{
    // Algorithm for neutralizing molecules from
    // https://www.rdkit.org/docs/Cookbook.html#neutralizing-molecules by Noel
    // O’Boyle Will neutralize the molecule by adding or removing hydrogens as
    // needed. This will ensure SMILES can be used to match atomistic structures
    // to the correct monomer.
    static const std::unique_ptr<RDKit::RWMol> neutralize_query(
        RDKit::SmartsToMol(
            "[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]"));
    for (const auto& match : RDKit::SubstructMatch(mol, *neutralize_query)) {
        auto atom = mol.getAtomWithIdx(match[0].second);
        auto chg = atom->getFormalCharge();
        auto hcount = atom->getTotalNumHs();
        atom->setFormalCharge(0);
        atom->setNumExplicitHs(hcount - chg);
        atom->updatePropertyCache();
    }
}

void buildMonomerMol(const RDKit::ROMol& atomistic_mol,
                std::vector<std::vector<int>>& monomers,
                boost::shared_ptr<RDKit::RWMol> monomer_mol,
                std::vector<Linkage>& linkages)
{
    // Start with all atoms in a single peptide chain
    monomer_mol->setProp<bool>(HELM_MODEL, true);

    constexpr bool isomeric_smiles = false;
    int residue_num = 1;
    for (const auto& monomer : monomers) {
        auto monomer_smiles = RDKit::MolFragmentToSmiles(
            atomistic_mol, monomer, nullptr, nullptr, nullptr, isomeric_smiles);
        // We have to roundtrip to canonicalize smiles -- see RDKit issue #7214
        std::unique_ptr<RDKit::RWMol> canon_mol(
            RDKit::SmilesToMol(monomer_smiles));
        neutralizeAtoms(*canon_mol);
        monomer_smiles = RDKit::MolToSmiles(*canon_mol);

        // If the monomer is a known amino acid, use the 1-letter code
        if (amino_acids.find(monomer_smiles) != amino_acids.end()) {
            addMonomer(*monomer_mol, amino_acids.at(monomer_smiles), residue_num,
                       "PEPTIDE1");
        } else {
            addMonomer(*monomer_mol, monomer_smiles, residue_num, "PEPTIDE1",
                       MonomerType::SMILES);
        }
        ++residue_num;

        // TODO: Check for known nucleic acids and CHEM monomers
    }

    for (const auto& link : linkages) {
        // TODO: Non forward linkages
        addConnection(*monomer_mol, link.monomer_idx1, link.monomer_idx2,
                      link.to_string());
    }
}

void removeWaters(RDKit::RWMol& mol)
{
    mol.beginBatchEdit();
    auto is_water = [](const RDKit::Atom* atom) {
        const auto res_info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());

        // TODO: This seems like it shouldn't be the job of this function; there
        // should be some sort of separate preprocessing step that removes
        // waters and other unwanted residues
        if (!res_info) {
            return false;
        }
        // strip whitespace from residue name
        auto res_name = res_info->getResidueName();
        res_name.erase(std::remove(res_name.begin(), res_name.end(), ' '),
                       res_name.end());
        if (res_info && res_name == "HOH") {
            return true;
        }
        return false;
    };
    for (auto atom : mol.atoms()) {
        if (is_water(atom)) {
            mol.removeAtom(atom);
        }
    }
    mol.commitBatchEdit();
}

void findChainsAndResidues(const RDKit::ROMol& mol,
                           ChainsAndResidues& chains_and_residues)
{
    // Find all chains and residues in the molecule
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        const auto atom = mol.getAtomWithIdx(i);
        const auto res_info = static_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());

        if (!res_info) {
            if (atom->getAtomicNum() == 1) {
                // This is a hydrogen that is part of a chain, so add it to the
                // chain
                continue;
            }
            throw std::runtime_error(
                "Atom " + std::to_string(atom->getIdx()) + " does not have residue info");
        }
        const auto chain_id = res_info->getChainId();
        const auto res_num = res_info->getResidueNumber();
        const auto ins_code = res_info->getInsertionCode();
        chains_and_residues[chain_id][std::make_pair(res_num, ins_code)]
            .push_back(i);
    }
}



static void copy_selected_atoms_and_bonds(::RDKit::RWMol& extracted_mol,
                                          const RDKit::ROMol& reference_mol,
                                          const std::vector<bool>& selected_atoms,
                                          const std::vector<bool>& selected_bonds)
{
    std::unordered_map<unsigned int, unsigned int> atom_mapping;
    std::unordered_map<unsigned int, unsigned int> bond_mapping;
    for (const auto& ref_atom : reference_mol.atoms()) {
        if (!selected_atoms[ref_atom->getIdx()]) {
            continue;
        }

        std::unique_ptr<::RDKit::Atom> extracted_atom{ref_atom->copy()};

        static constexpr bool updateLabel = true;
        static constexpr bool takeOwnership = true;
        atom_mapping[ref_atom->getIdx()] = extracted_mol.addAtom(
            extracted_atom.release(), updateLabel, takeOwnership);
    }

    for (const auto& ref_bond : reference_mol.bonds()) {
        if (!selected_bonds[ref_bond->getIdx()]) {
            continue;
        }

        std::unique_ptr<::RDKit::Bond> extracted_bond{ref_bond->copy()};
        extracted_bond->setBeginAtomIdx(
            atom_mapping[ref_bond->getBeginAtomIdx()]);
        extracted_bond->setEndAtomIdx(atom_mapping[ref_bond->getEndAtomIdx()]);

        static constexpr bool takeOwnership = true;
        auto num_bonds =
            extracted_mol.addBond(extracted_bond.release(), takeOwnership);
        bond_mapping[ref_bond->getIdx()] = num_bonds - 1;
    }
}

boost::shared_ptr<RDKit::RWMol>
extractMolFragment(const RDKit::ROMol& mol,
                   const std::vector<unsigned int>& atom_ids)
{

    const auto num_atoms = mol.getNumAtoms();
    std::vector<bool> selected_atoms(num_atoms);
    std::vector<bool> selected_bonds(mol.getNumBonds());
    for (const auto& atom_idx : atom_ids) {
        if (atom_idx < num_atoms) {
            selected_atoms[atom_idx] = true;
        }
    }
    for (const auto& bond : mol.bonds()) {
        if (selected_atoms[bond->getBeginAtomIdx()] &&
            selected_atoms[bond->getEndAtomIdx()]) {
            selected_bonds[bond->getIdx()] = true;
        }
    }

    auto extracted_mol = std::make_unique<::RDKit::RWMol>();
    copy_selected_atoms_and_bonds(*extracted_mol, mol, selected_atoms, selected_bonds);

    // NOTE: Bookmarks are currently not copied
    return extracted_mol;
}


std::string getMonomerSmiles(RDKit::ROMol& mol,
                             const std::vector<unsigned int>& atom_idxs,
                             const std::string& chain_id,
                             const std::pair<int, std::string>& current_key,
                             int res_num, bool end_of_chain)
{
    // Determine the atoms in current_res that connect to adjacent residues
    std::vector<std::pair<int, int>> attch_idxs; // adjacent res num, atom_idx
    for (auto idx : atom_idxs) {
        const auto at = mol.getAtomWithIdx(idx);
        for (const auto neigh : mol.atomNeighbors(at)) {
            if (neigh->getAtomicNum() == 1) {
                // skip Hs
                continue;
            }
            const auto res_info =
                dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                    neigh->getMonomerInfo());
            const auto key = std::make_pair(res_info->getResidueNumber(),
                                            res_info->getInsertionCode());
            if (key != current_key || res_info->getChainId() != chain_id) {
                // neighbor is in different residue, this will be an attachment
                // point
                attch_idxs.push_back(
                    {res_info->getResidueNumber(), at->getIdx()});
            }
        }
    }

    // This part is difficult; we need to figure out the order of the attachment
    // points on this monomer, Ideally we'd find the backbone and ensure to
    // label them directionally correctly using R2-R1
    std::ranges::sort(
        attch_idxs.begin(), attch_idxs.end(),
        [&res_num](const std::pair<int, int>& a, const std::pair<int, int>& b) {
            // we want the most 'adjacent' residue first in terms of residue
            // ordering; this gets more complicated if we have different chains
            if (std::abs(a.first - res_num) == std::abs(b.first - res_num)) {
                return a.first < b.first;
            }
            return std::abs(a.first - res_num) < std::abs(b.first - res_num);
        });

    auto mol_fragment = extractMolFragment(mol, atom_idxs);

    // Add dummy atoms with an attachment point #s
    // For now, all monomers have at least attachment points 1 and 2 so
    // that backbone bonds can be formed (except the beginning monomer)
    int current_attchpt = 1;

    // If this is the beginning of the chain or a branch, start with attachment
    // point 2
    if (!end_of_chain && res_num == 1) {
        current_attchpt = 2;
    }

    static constexpr bool update_label = true;
    static constexpr bool take_ownership = true;
    for (const auto& [_, ref_idx] : attch_idxs) {
        for (auto at : mol_fragment->atoms()) {
            int orig_idx;
            if (at->getPropIfPresent(REFERENCE_IDX, orig_idx) &&
                orig_idx == ref_idx) {
                auto new_at = new RDKit::Atom(0);
                new_at->setProp(RDKit::common_properties::molAtomMapNumber,
                                current_attchpt);
                auto new_at_idx =
                    mol_fragment->addAtom(new_at, update_label, take_ownership);
                mol_fragment->addBond(new_at_idx, at->getIdx(),
                                      RDKit::Bond::SINGLE);
                mol.getAtomWithIdx(orig_idx)->setProp<unsigned int>(
                    MONOMER_MAP_NUM, current_attchpt);
                ++current_attchpt;
            }
        }
    }
    // There should always be enough attachment points so that
    // backbone connections can be made (R1 and R2)
    // TODO: Should this indicate a new chain?
    while (current_attchpt <= MIN_ATTCHPTS) {
        if (end_of_chain && current_attchpt > 1) {
            break;
        }
        auto new_at = new RDKit::Atom(0);
        new_at->setProp(RDKit::common_properties::molAtomMapNumber,
                        current_attchpt);
        mol_fragment->addAtom(new_at, update_label, take_ownership);
        ++current_attchpt;
    }

    // removing hydrogens to keep HELM string readable
    MolOps::removeHs(*mol_fragment);
    return RDKit::MolToSmiles(*mol_fragment);
}

bool sameMonomer(RDKit::RWMol& atomistic_mol,
                 const std::vector<unsigned int>& atom_idxs,
                 const std::string& db_smiles)
{
    // Occasionally SMILES that cannot be kekulized are extracted from atomistic
    // mol, so skip sanitization
    constexpr int debug = 0;
    constexpr bool sanitize = false;

    // Monomer from original atomistic mol and monomer defined by database
    auto monomer_frag = extractMolFragment(atomistic_mol, atom_idxs);
    std::unique_ptr<RDKit::RWMol> db_mol(
        RDKit::SmilesToMol(db_smiles, debug, sanitize));

    // Remove stereochemistry, atom map numbers, and neutralize the atoms
    // leaving groups are removed since they won't always be included in the
    // residue extracted from the atomistic molecule
    auto clean_mol = [](RDKit::RWMol& mol) {
        RDKit::MolOps::removeStereochemistry(mol);
        neutralizeAtoms(mol);
        mol.beginBatchEdit();
        for (auto at : mol.atoms()) {
            unsigned int map_no;
            if (at->getPropIfPresent(RDKit::common_properties::molAtomMapNumber,
                                     map_no)) {
                // Label parent with this map number so we know which map number
                // the linkage uses
                for (const auto& nbr : mol.atomNeighbors(at)) {
                    nbr->setProp(RDKit::common_properties::molAtomMapNumber,
                                 map_no);
                }

                mol.removeAtom(at);
            }
        }
        mol.commitBatchEdit();
        MolOps::removeHs(mol);
        // set aromaticity
        RDKit::MolOps::setAromaticity(mol);
    };
    clean_mol(*monomer_frag);
    clean_mol(*db_mol);

    // The DB monomer has had the leaving groups removed, while the residue
    // extracted from the atomistic mol may still have them present if it is at
    // the beginning or end of a chain. As a result, we need to allow for the DB
    // monomer to have one less atom than the residue extracted from the
    // atomistic mol
    auto match = RDKit::SubstructMatch(*monomer_frag,
                                       *db_mol); // (queryAtomIdx, molAtomIdx)

    if (match.size()) {

        // Any unmapped atom in monomer_frag must map to a leaving group in the
        // db_mol create vector of atom indices in the monomer_frag that are not
        // in the match
        for (auto at : monomer_frag->atoms()) {
            int target_idx = at->getIdx();
            auto it = std::find_if(match[0].begin(), match[0].end(),
                                   [target_idx](const std::pair<int, int>& p) {
                                       return p.second == target_idx;
                                   });
            if (it == match[0].end()) {
                // This atom should be terminal, otherwise the match is
                // definitely wrong
                if (at->getDegree() > 1) {
                    return false;
                }

                // Unmatched atoms should be leaving groups, we can check this
                // by matching
                RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
                boost::tie(nbrIdx, endNbrs) =
                    monomer_frag->getAtomNeighbors(at);

                it = std::find_if(match[0].begin(), match[0].end(),
                                  [nbrIdx](const std::pair<int, int>& p) {
                                      return p.second ==
                                             static_cast<int>(*nbrIdx);
                                  });
                if (it == match[0].end()) {
                    return false;
                }

                auto db_atom = db_mol->getAtomWithIdx(it->first);
                if (!db_atom->hasProp(
                        RDKit::common_properties::molAtomMapNumber)) {
                    return false;
                }
            }
        }

        // we are mapping atomistic atom indices to attachment point numbers
        // (molAtomMapNumber)
        for (auto [db_idx, at_idx] : match[0]) {
            auto db_at = db_mol->getAtomWithIdx(db_idx);
            auto at = monomer_frag->getAtomWithIdx(at_idx);
            unsigned int map_no;
            if (db_at->getPropIfPresent(
                    RDKit::common_properties::molAtomMapNumber, map_no)) {
                auto ref_idx = at->getProp<unsigned int>(REFERENCE_IDX);
                atomistic_mol.getAtomWithIdx(ref_idx)->setProp<unsigned int>(
                    MONOMER_MAP_NUM, map_no);
            }
        }
        return true;
    }

    return false;
}

int getAttchpt(const RDKit::Atom& monomer, const RDKit::Bond& bond,
               const RDKit::RWMol& atomistic_mol)
{
    // This should be the index of the atom in this monomer
    auto target_atom_idx = bond.getBeginAtomIdx();
    if (atomistic_mol.getAtomWithIdx(target_atom_idx)
            ->getProp<unsigned int>(MONOMER_IDX1) != monomer.getIdx()) {
        target_atom_idx = bond.getEndAtomIdx();
    }
    unsigned int attchpt;
    if (atomistic_mol.getAtomWithIdx(target_atom_idx)
            ->getPropIfPresent(MONOMER_MAP_NUM, attchpt)) {
        return attchpt;
    } else {
        return -1;
    }
}

void detectLinkages(RDKit::RWMol& monomer_mol, const RDKit::RWMol& atomistic_mol)
{
    // Find all linkages between monomers
    for (const auto& bond : atomistic_mol.bonds()) {
        const auto begin_atom = bond->getBeginAtom();
        const auto end_atom = bond->getEndAtom();

        unsigned int begin_monomer_idx;
        unsigned int end_monomer_idx;
        if (!begin_atom->getPropIfPresent(MONOMER_IDX1, begin_monomer_idx) ||
            !end_atom->getPropIfPresent(MONOMER_IDX1, end_monomer_idx)) {
            // One of these is a hydrogen or unmapped for some reason
            continue;
        }

        // Check if this is already present as a backbone linkage
        if (begin_monomer_idx == end_monomer_idx ||
            monomer_mol.getBondBetweenAtoms(begin_monomer_idx, end_monomer_idx) !=
                nullptr) {
            continue;
        }

        const auto begin_monomer = monomer_mol.getAtomWithIdx(begin_monomer_idx);
        const auto end_monomer = monomer_mol.getAtomWithIdx(end_monomer_idx);

        auto begin_attchpt = getAttchpt(*begin_monomer, *bond, atomistic_mol);
        auto end_attchpt = getAttchpt(*end_monomer, *bond, atomistic_mol);

        if (begin_attchpt == -1 || end_attchpt == -1) {
            // This happens when the input atomistic mol has bonds between
            // residues on atoms that are not marked as attachment points in the
            // monomer DB.
            auto begin_res_info =
                dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                    begin_atom->getMonomerInfo());
            auto end_res_info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                end_atom->getMonomerInfo());
            if (begin_res_info && end_res_info) {
                std::cerr << "Bond between residue " << begin_res_info->getResidueName()
                          << begin_res_info->getResidueNumber() << " in chain "
                          << begin_res_info->getChainId() << " and residue "
                          << end_res_info->getResidueName() << end_res_info->getResidueNumber()
                          << " in chain " << end_res_info->getChainId()
                          << " do not correspond to attachment points in their residues.\n";
            } else {
                std::cerr << "Bond between atoms " << begin_atom->getIdx()
                          << " and " << end_atom->getIdx() << " in monomers "
                          << begin_monomer->getProp<std::string>(ATOM_LABEL)
                          << " and " << end_monomer->getProp<std::string>(ATOM_LABEL)
                          << " do not correspond to attachment points in their residues.\n";
            }
            continue;
        }

        // Backbone connections (R2-R1) should be added in that order when
        // possible.
        if (begin_attchpt == 2 && end_attchpt == 1) {
            addConnection(monomer_mol, begin_monomer_idx, end_monomer_idx, "R2-R1");
        } else if (begin_attchpt == 1 && end_attchpt == 2) {
            addConnection(monomer_mol, end_monomer_idx, begin_monomer_idx, "R2-R1");
        } else if (begin_monomer_idx < end_monomer_idx) {
            addConnection(monomer_mol, begin_monomer_idx, end_monomer_idx,
                          "R" + std::to_string(begin_attchpt) + "-R" + std::to_string(end_attchpt));
        } else {
            addConnection(monomer_mol, end_monomer_idx, begin_monomer_idx,
                          "R" + std::to_string(end_attchpt) + "-R" + std::to_string(begin_attchpt));
        }
    }
}

std::optional<std::tuple<std::string, std::string, ChainType>>
getHelmInfo(MonomerDatabase& db, const RDKit::Atom* atom)
{
    // This comes from something like a PDB or MAE file
    auto res_info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
        atom->getMonomerInfo());
    auto res_name = res_info->getResidueName();
    res_name.erase(
        std::remove_if(res_name.begin(), res_name.end(), ::isspace),
        res_name.end());
    return db.getHelmInfo(res_name);
}

boost::shared_ptr<RDKit::RWMol>
pdbInfoAtomisticToMM(const RDKit::ROMol& input_mol)
{
    // Make RWMol and remove waters
    RDKit::RWMol mol(input_mol);
    removeWaters(mol);

    // Set reference index for SMILES fragments
    for (auto at : mol.atoms()) {
        at->setProp(REFERENCE_IDX, at->getIdx());
    }

    // Map chain_id -> {residue mols}
    ChainsAndResidues chains_and_residues;
    findChainsAndResidues(mol, chains_and_residues);

    // Eventually, this will be connecting to a database of monomers or accessing an
    // in-memory datastructure
    MonomerDatabase db;
    std::map<ChainType, unsigned int> chain_counts = {{ChainType::PEPTIDE, 0},
                                                      {ChainType::RNA, 0},
                                                      {ChainType::DNA, 0},
                                                      {ChainType::CHEM, 0}};
    boost::shared_ptr<RDKit::RWMol> monomer_mol = boost::make_shared<RDKit::RWMol>();
    for (const auto& [chain_id, residues] : chains_and_residues) {
        // Use first residue to determine chain type. We assume that PDB data
        // is correct and there aren't multiple chain types in a single chain.
        // Default chain type is PEPTIDE if not specified.
        auto helm_info = getHelmInfo(
            db, mol.getAtomWithIdx(residues.begin()->second[0]));
        auto chain_type =
            helm_info ? std::get<2>(*helm_info) : ChainType::PEPTIDE;
        std::string helm_chain_id = toString(chain_type) +
                                     std::to_string(++chain_counts[chain_type]);
        // Assuming residues are ordered correctly
        size_t res_num = 1;
        for (const auto& [key, atom_idxs] : residues) {
            helm_info = getHelmInfo(db, mol.getAtomWithIdx(atom_idxs[0]));
            bool end_of_chain = res_num == residues.size();
            size_t this_monomer;
            if (helm_info &&
                sameMonomer(mol, atom_idxs, std::get<1>(*helm_info))) {
                // Standard residue in monomer DB, Verify that the fragment
                // labeled as the residue matches what is in the monomer
                // database
                this_monomer = addMonomer(*monomer_mol, std::get<0>(*helm_info),
                                          res_num, helm_chain_id);
            } else {
                auto smiles = getMonomerSmiles(mol, atom_idxs, chain_id, key,
                                               res_num, end_of_chain);
                this_monomer = addMonomer(*monomer_mol, smiles, res_num,
                                          helm_chain_id, MonomerType::SMILES);
            }

            // Track which atoms are in which monomer
            for (auto idx : atom_idxs) {
                mol.getAtomWithIdx(idx)->setProp<unsigned int>(MONOMER_IDX1,
                                                               this_monomer);
            }
            ++res_num;
        }
    }
    detectLinkages(*monomer_mol, mol);
    return monomer_mol;
}

bool hasPdbResidueInfo(const RDKit::ROMol& mol)
{
    return std::any_of(mol.atoms().begin(), mol.atoms().end(),
                       [](const auto& atom) {
                           return atom->getMonomerInfo() &&
                                  atom->getMonomerInfo()->getMonomerType() ==
                                      RDKit::AtomMonomerInfo::PDBRESIDUE;
                       });
}
} // unnamed namespace

boost::shared_ptr<RDKit::RWMol> atomisticToMonomerMol(const RDKit::ROMol& mol,
                                                      bool use_residue_info)
{
    // Use residue information to build the monomeristic molecule. Assumes that the
    // residue information is correct, and throws if any residue information
    // is missing.
    RDKit::ROMol atomistic_mol(mol);
    if (use_residue_info) {
        if (hasPdbResidueInfo(atomistic_mol)) {
            auto monomer_mol = pdbInfoAtomisticToMM(atomistic_mol);
            assignChains(*monomer_mol);
            return monomer_mol;
        } else {
            throw std::runtime_error(
                "No residue information found in molecule");
        }
    }

    // Work on copy, for now
    std::vector<std::vector<int>> monomers;
    std::vector<Linkage> linkages;
    identifyMonomers(atomistic_mol, monomers, linkages);

    boost::shared_ptr<RDKit::RWMol> monomer_mol = boost::make_shared<RDKit::RWMol>();
    buildMonomerMol(atomistic_mol, monomers, monomer_mol, linkages);
    assignChains(*monomer_mol);

    // TODO
    // Now that we have the monomeristic mol, we need to set the properties needed by the
    // HELM writer and other functions that work with monomeristic mols created by the
    // HELM parser. I think this will include a few steps
    // 1. Break the monomeristic Mol into polymers -- this would be done by connectivity
    // and monomer type (peptide, rna, dna, chem)
    // 2. Insure that the linkage information is correct -- backbone vs not, etc
    // 3. Set the polymers as substance groups on the molecule, and set
    // monomer-specific properties
    // 4. Maybe: make sure monomeristic monomer indices are in connectivity order

    return monomer_mol;
}

} // namespace RDKit