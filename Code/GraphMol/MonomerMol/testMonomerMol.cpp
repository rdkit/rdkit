//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <sstream>
#include "RDGeneral/test.h"
#include <catch2/catch_all.hpp>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/FileParseException.h>
#include <boost/algorithm/string.hpp>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/Chirality.h>

#include <GraphMol/FileParsers/SequenceParsers.h>

#include <GraphMol/MonomerMol/Conversions.h>
#include <GraphMol/MonomerMol/MonomerMol.h>
#include <GraphMol/MonomerInfo.h>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <string>

using namespace RDKit;

/*
 * Temporary, simple FASTA parser to show how to use the MonomerMol
*/
static std::string to_fasta(const ROMol& mol)
{
    // make a set of all the chains in the molecule
    std::set<std::string> chains;
    for (auto atom : mol.atoms()) {
        chains.insert(getPolymerId(atom));
    }
    std::stringstream fasta;
    bool first = true;
    for (auto chain_id : chains) {
        if (!first) {
            fasta << "\n";
        }
        first = false;
        fasta << ">Chain " << chain_id << "\n"; // add title
        auto chain = getPolymer(mol, chain_id);
        for (auto atom : chain.atoms) {
            // Really there should be a lookup to ensure that these are all 1 letter
            // And a sort by residue number & connectivity
            const auto* res_info = static_cast<const RDKit::AtomPDBResidueInfo*>(mol.getAtomWithIdx(atom)->getMonomerInfo());
            fasta << res_info->getResidueName();
        }
    }
    return fasta.str();
}

static void neutralize_atoms(RDKit::ROMol& mol)
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

static std::unique_ptr<RDKit::ROMol> resolve_his(const RDKit::ROMol& mol)
{
    // Some structures may contain different protonation states for histidine,
    // but we currently map all of them to the same single letter code 'H' in
    // the monomeric representation. Since we want to test the roundtrip
    // conversion against the original, we need to resolve the histidine
    // protonation state to what is in the monomer database.
    std::string smiles = RDKit::MolToSmiles(mol);

    std::vector<std::string> targets = {"Cc1c[nH]cn1"};
    std::string replace_with = "Cc1cnc[nH]1";
    for (const auto& target : targets) {
        size_t pos = 0;
        while ((pos = smiles.find(target, pos)) != std::string::npos) {
            smiles.replace(pos, target.length(), replace_with);
            pos += replace_with.length();
        }
    }
    return std::unique_ptr<RDKit::ROMol>(RDKit::SmilesToMol(smiles));
}

static void remove_solvents(RDKit::RWMol& rwmol)
{
    std::vector<unsigned int> atoms_to_remove;
    for (const auto& atom : rwmol.atoms()) {
        const auto res_info = static_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());

        if (res_info != nullptr) {
            std::string residue_name = res_info->getResidueName();
            if (residue_name == "HOH" || residue_name == "S04") {
                atoms_to_remove.push_back(atom->getIdx());
            }
        }
    }
    rwmol.beginBatchEdit();
    for (unsigned int atom_idx : atoms_to_remove) {
        rwmol.removeAtom(atom_idx);
    }
    rwmol.commitBatchEdit();
}

static bool same_roundtrip_mol(const RDKit::ROMol& original,
                               const RDKit::ROMol& roundtrip,
                               unsigned int num_chains) {
    // Verify the roundtrip produces a molecule with expected number of atoms,
    // roundtripped structure has additional oxygen atoms at the end of the
    // chains so there may be some extra atoms in the roundtrip (up to 2 *
    // number of chains)
    RDKit::SubstructMatchParameters params;
    params.maxMatches = 1;
    auto match = RDKit::SubstructMatch(roundtrip,
                                       original, params);
    if (match.empty()) {
      return false;
    }
    if (match[0].size() != original.getNumAtoms()) {
      return false;
    }
    return (roundtrip.getNumAtoms() - original.getNumAtoms()) <= (2 * num_chains);
}

TEST_CASE("FASTAConversions") {
  SECTION("SIMPLE") {
    // Build MonomerMol with a single chain
    MonomerMol monomer_mol;
    monomer_mol.addMonomer("R", 1, "PEPTIDE1");
    monomer_mol.addMonomer("D");
    monomer_mol.addMonomer("K");
    monomer_mol.addMonomer("I");
    monomer_mol.addMonomer("T");

    monomer_mol.addConnection(0, 1, ConnectionType::FORWARD);
    monomer_mol.addConnection(1, 2, ConnectionType::FORWARD);
    monomer_mol.addConnection(2, 3, ConnectionType::FORWARD);
    monomer_mol.addConnection(3, 4, ConnectionType::FORWARD);

    CHECK(monomer_mol.getNumAtoms() == 5);
    CHECK(monomer_mol.getNumBonds() == 4);
    CHECK(std::string(">Chain PEPTIDE1\nRDKIT") == to_fasta(monomer_mol));
  }

  SECTION("MultipleChains") {
    // Build MonomerMol with two chains
    MonomerMol monomer_mol;
    auto midx1 = monomer_mol.addMonomer("R", 1, "A");
    auto midx2 = monomer_mol.addMonomer("D");

    auto midx3 = monomer_mol.addMonomer("K", 1, "B");
    auto midx4 = monomer_mol.addMonomer("I");
    auto midx5 = monomer_mol.addMonomer("T");

    monomer_mol.addConnection(midx1, midx2, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx3, midx4, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx4, midx5, ConnectionType::FORWARD);

    CHECK(monomer_mol.getNumAtoms() == 5);
    CHECK(monomer_mol.getNumBonds() == 3);

    CHECK(std::string(">Chain A\nRD\n>Chain B\nKIT") == to_fasta(monomer_mol));
  }

  SECTION("UsingSequenceReader") {
    std::string seq = "CGCGAATTACCGCG";
    // the sequence parser creates an atomistic molecule
    auto atomistic_mol = SequenceToMol(seq);
    CHECK(atomistic_mol);

    // PDB info is used to convert the atomistic sturcture into a MonomerMol
    auto monomer_mol = toMonomeric(*atomistic_mol);
    CHECK(std::string(">Chain PEPTIDE1\nCGCGAATTACCGCG") == to_fasta(*monomer_mol));
    delete atomistic_mol;
  }
}

TEST_CASE("Conversions") {
  SECTION("toAtomistic") {
    std::string seq = "CGCGA";
    auto atomistic_mol = SequenceToMol(seq);

    MonomerMol monomer_mol;
    auto midx1 = monomer_mol.addMonomer("C", 1, "PEPTIDE1");
    auto midx2 = monomer_mol.addMonomer("G");
    auto midx3 = monomer_mol.addMonomer("C");
    auto midx4 = monomer_mol.addMonomer("G");
    auto midx5 = monomer_mol.addMonomer("A");

    monomer_mol.addConnection(midx1, midx2, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx2, midx3, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx3, midx4, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx4, midx5, ConnectionType::FORWARD);
    auto atomistic_mol2 = toAtomistic(monomer_mol);

    // atomistic structure is same as using sequence parser
    std::string smi1 = MolToSmiles(*atomistic_mol);
    std::string smi2 = MolToSmiles(*atomistic_mol2);
    CHECK(smi1 == smi2);

    delete atomistic_mol;
  }

  SECTION("toAtomisticWithBranch") {
    // This is equivalent to HELM string "PEPTIDE1{A.D(C)P}$$$$"
    MonomerMol monomer_mol;
    auto midx1 = monomer_mol.addMonomer("A", 1, "PEPTIDE1");
    auto midx2 = monomer_mol.addMonomer("D");
    auto midx3 = monomer_mol.addMonomer("C");
    auto midx4 = monomer_mol.addMonomer("P");

    monomer_mol.addConnection(midx1, midx2, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx2, midx3, ConnectionType::SIDECHAIN);
    monomer_mol.addConnection(midx2, midx4, ConnectionType::FORWARD);

    std::string smi = MolToSmiles(*toAtomistic(monomer_mol));
    CHECK(smi == "C[C@H](N)C(=O)N[C@@H](CC(=O)N[C@@H](CS)C(=O)O)C(=O)N1CCC[C@H]1C(=O)O");
  }


  SECTION("toAtomisticWithDisulfide") {
    std::string helm = "PEPTIDE1{C.A.A.A.C}$PEPTIDE1,PEPTIDE1,1:R3-5:R3$$$V2.0";
    auto atomistic_mol = HELMToMol(helm);

    MonomerMol monomer_mol;
    auto midx1 = monomer_mol.addMonomer("C", 1, "PEPTIDE1");
    auto midx2 = monomer_mol.addMonomer("A");
    auto midx3 = monomer_mol.addMonomer("A");
    auto midx4 = monomer_mol.addMonomer("A");
    auto midx5 = monomer_mol.addMonomer("C");

    monomer_mol.addConnection(midx1, midx2, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx2, midx3, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx3, midx4, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx4, midx5, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx1, midx5, ConnectionType::CROSSLINK);

    std::string smi1 = MolToSmiles(*atomistic_mol);
    std::string smi2 = MolToSmiles(*toAtomistic(monomer_mol));
    CHECK(smi1 == smi2);
    delete atomistic_mol;
  }

  SECTION("toAtomisticCyclicPeptide") {
    std::string helm = "PEPTIDE1{F.Y.K.A.R.L}$PEPTIDE1,PEPTIDE1,6:R2-1:R1$$$V2.0";
    auto atomistic_mol = HELMToMol(helm);

    MonomerMol monomer_mol;
    auto midx1 = monomer_mol.addMonomer("F", 1, "PEPTIDE1");
    auto midx2 = monomer_mol.addMonomer("Y");
    auto midx3 = monomer_mol.addMonomer("K");
    auto midx4 = monomer_mol.addMonomer("A");
    auto midx5 = monomer_mol.addMonomer("R");
    auto midx6 = monomer_mol.addMonomer("L");
    monomer_mol.addConnection(midx1, midx2, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx2, midx3, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx3, midx4, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx4, midx5, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx5, midx6, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx6, midx1, ConnectionType::FORWARD);

    CHECK(monomer_mol.getNumAtoms() == 6);
    CHECK(monomer_mol.getNumBonds() == 6);

    std::string smi1 = MolToSmiles(*atomistic_mol);
    std::string smi2 = MolToSmiles(*toAtomistic(monomer_mol));
    CHECK(smi1 == smi2);
    delete atomistic_mol;
  }

  SECTION("toAtomisticCyclicPeptideOutOfOrder") {
    std::string helm = "PEPTIDE1{F.Y.K.A.R.L}$PEPTIDE1,PEPTIDE1,6:R2-1:R1$$$V2.0";
    auto atomistic_mol = HELMToMol(helm);


    // toAtomistic should still work even when monomers are added out of order
    MonomerMol monomer_mol;
    auto midx1 = monomer_mol.addMonomer("A", 1, "PEPTIDE1");
    auto midx2 = monomer_mol.addMonomer("K");
    auto midx3 = monomer_mol.addMonomer("Y");
    auto midx4 = monomer_mol.addMonomer("F");
    auto midx5 = monomer_mol.addMonomer("L");
    auto midx6 = monomer_mol.addMonomer("R");
    monomer_mol.addConnection(midx4, midx3, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx3, midx2, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx2, midx1, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx1, midx6, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx6, midx5, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx5, midx4, ConnectionType::FORWARD);
    assignChains(monomer_mol);

    CHECK(monomer_mol.getNumAtoms() == 6);
    CHECK(monomer_mol.getNumBonds() == 6);

    std::string smi1 = MolToSmiles(*atomistic_mol);
    std::string smi2 = MolToSmiles(*toAtomistic(monomer_mol));
    CHECK(smi1 == smi2);
    delete atomistic_mol;
  }

  SECTION("toAtomisticSmilesMonomer") {
    // This is equivelent to HELM string "PEPTIDE1{A.[O=C([C@H]1CCCN1[*:1])[*:2]].P}$$$$"
    MonomerMol monomer_mol;
    auto midx1 = monomer_mol.addMonomer("A", 1, "PEPTIDE1");
    auto midx2 = monomer_mol.addMonomer("O=C([C@H]1CCCN1[*:1])[*:2]", MonomerType::SMILES);
    auto midx3 = monomer_mol.addMonomer("P");
    monomer_mol.addConnection(midx1, midx2, ConnectionType::FORWARD);
    monomer_mol.addConnection(midx2, midx3, ConnectionType::FORWARD);
    assignChains(monomer_mol);
    CHECK(monomer_mol.getNumAtoms() == 3);
    CHECK(monomer_mol.getNumBonds() == 2);
    auto atomistic_mol = toAtomistic(monomer_mol);
    std::string smi1 = MolToSmiles(*atomistic_mol);
    CHECK(smi1 == "C[C@H](N)C(=O)N1CCC[C@@H]1C(=O)N1CCC[C@H]1C(=O)O");
  }

  SECTION("toMonomeric1DNG") {
    std::string pdbfile = getenv("RDBASE");
    pdbfile += "/Code/GraphMol/MonomerMol/test_data/1dng.pdb";
    auto mol = PDBFileToMol(pdbfile);
    auto monomer_mol = toMonomeric(*mol);
    CHECK(std::string(">Chain PEPTIDE1\nQAPAYEEAAEELAKS") == to_fasta(*monomer_mol));
    CHECK(monomer_mol->getNumAtoms() == 15);
    CHECK(monomer_mol->getNumBonds() == 14);

    auto roundtrip = toAtomistic(*monomer_mol);
    neutralize_atoms(*mol);
    neutralize_atoms(*roundtrip);
    CHECK(same_roundtrip_mol(*mol, *roundtrip, getPolymerIds(*monomer_mol).size()));
  }

  SECTION("toMonomeric2N65") {
    // Example with multiple chains and a disulfide bond between chains
    std::string pdbfile = getenv("RDBASE");
    pdbfile += "/Code/GraphMol/MonomerMol/test_data/2n65.pdb";
    auto mol = PDBFileToMol(pdbfile);
    auto monomer_mol = toMonomeric(*mol);
    CHECK(std::string(">Chain PEPTIDE1\nVARGWKRKCPLFGKGG\n>Chain PEPTIDE2\nVARGWKRKCPLFGKGG") == to_fasta(*monomer_mol));
    CHECK(monomer_mol->getNumAtoms() == 32);
    CHECK(monomer_mol->getNumBonds() == 31);

    auto roundtrip = toAtomistic(*monomer_mol);
    neutralize_atoms(*mol);
    auto resolved_his = resolve_his(*mol);
    neutralize_atoms(*roundtrip);
    CHECK(same_roundtrip_mol(*resolved_his, *roundtrip, getPolymerIds(*monomer_mol).size()));
  }

  SECTION("toMonomeric4QAF") {
    // Example with disulfide bond between chains and within chains
    std::string pdbfile = getenv("RDBASE");
    pdbfile += "/Code/GraphMol/MonomerMol/test_data/4qaf.pdb";
    auto mol = PDBFileToMol(pdbfile);
    remove_solvents(*mol); // we don't care about water or S04

    auto monomer_mol = toMonomeric(*mol);

    // There should be 4 chains
    auto polymer_ids = getPolymerIds(*monomer_mol);
    CHECK(polymer_ids.size() == 4);

    CHECK(monomer_mol->getNumAtoms() == 395);
    CHECK(monomer_mol->getNumBonds() == 392);

    auto roundtrip = toAtomistic(*monomer_mol);
    neutralize_atoms(*mol);
    neutralize_atoms(*roundtrip);
    CHECK(same_roundtrip_mol(*mol, *roundtrip, polymer_ids.size()));
  }

  SECTION("toMonomeric5VAV") {
    // Example with disulfide bond and R2-R1 cycle closure
    std::string pdbfile = getenv("RDBASE");
    pdbfile += "/Code/GraphMol/MonomerMol/test_data/5vav.pdb";
    auto mol = PDBFileToMol(pdbfile);
    auto monomer_mol = toMonomeric(*mol);
    remove_solvents(*mol); // we don't care about water or S04

    // Checking sequence doesn't make sense due to SMILES monomer
    // Two extra bonds than atoms because this structure has two cycle closures, one
    // by a disulfide bond and another from an R2-R1 closure
    CHECK(monomer_mol->getNumAtoms() == 14);
    CHECK(monomer_mol->getNumBonds() == 15);

    auto roundtrip = toAtomistic(*monomer_mol);
    neutralize_atoms(*mol);
    neutralize_atoms(*roundtrip);
    CHECK(same_roundtrip_mol(*mol, *roundtrip, getPolymerIds(*monomer_mol).size()));
  }
}

TEST_CASE("DativeBonds") {
  SECTION("BackboneConnectionsUseDativeBonds") {
    // Test that backbone connections (R2-R1) use dative bonds
    MonomerMol mol;
    auto m1 = mol.addMonomer("A", 1, "PEPTIDE1");
    auto m2 = mol.addMonomer("G");

    mol.addConnection(m1, m2, ConnectionType::FORWARD);

    CHECK(mol.getNumBonds() == 1);
    auto bond = mol.getBondWithIdx(0);

    // Should be a dative bond since R2 > R1
    CHECK(bond->getBondType() == Bond::DATIVE);
    // Linkage should be R2-R1 (higher to lower)
    CHECK(bond->getProp<std::string>(LINKAGE) == "R2-R1");
    // Direction: begin atom should have the higher attachment point (R2)
    CHECK(bond->getBeginAtomIdx() == m1);
    CHECK(bond->getEndAtomIdx() == m2);
  }

  SECTION("SidechainConnectionsUseDativeBonds") {
    // Test that sidechain connections (R3-R1) use dative bonds
    MonomerMol mol;
    auto m1 = mol.addMonomer("D", 1, "PEPTIDE1");  // Aspartic acid has R3
    auto m2 = mol.addMonomer("C");                  // Cysteine

    mol.addConnection(m1, m2, ConnectionType::SIDECHAIN);

    CHECK(mol.getNumBonds() == 1);
    auto bond = mol.getBondWithIdx(0);

    // Should be a dative bond since R3 > R1
    CHECK(bond->getBondType() == Bond::DATIVE);
    CHECK(bond->getProp<std::string>(LINKAGE) == "R3-R1");
  }

  SECTION("CrosslinkConnectionsUseSingleBonds") {
    // Test that crosslink connections (R3-R3) use single bonds (same attachment)
    MonomerMol mol;
    auto m1 = mol.addMonomer("C", 1, "PEPTIDE1");
    auto m2 = mol.addMonomer("C", 2, "PEPTIDE1");

    mol.addConnection(m1, m2, ConnectionType::CROSSLINK);

    CHECK(mol.getNumBonds() == 1);
    auto bond = mol.getBondWithIdx(0);

    // Should be a single bond since R3 == R3
    CHECK(bond->getBondType() == Bond::SINGLE);
    CHECK(bond->getProp<std::string>(LINKAGE) == "R3-R3");
  }

  SECTION("DirectionalityIsNormalized") {
    // Test that linkage is always normalized to higher-to-lower
    MonomerMol mol;
    auto m1 = mol.addMonomer("A", 1, "PEPTIDE1");
    auto m2 = mol.addMonomer("G");

    // Add connection with R1-R2 (lower to higher) - should be normalized
    mol.addConnection(m1, m2, "R1-R2");

    auto bond = mol.getBondWithIdx(0);

    // Should be normalized to R2-R1 with swapped atom order
    CHECK(bond->getBondType() == Bond::DATIVE);
    CHECK(bond->getProp<std::string>(LINKAGE) == "R2-R1");
    // m2 should be the begin atom (has R2), m1 should be end (has R1)
    CHECK(bond->getBeginAtomIdx() == m2);
    CHECK(bond->getEndAtomIdx() == m1);
  }

  SECTION("ToAtomisticConvertsDativeToSingle") {
    // Test that toAtomistic converts dative bonds back to single bonds
    MonomerMol monomer_mol;
    auto m1 = monomer_mol.addMonomer("A", 1, "PEPTIDE1");
    auto m2 = monomer_mol.addMonomer("G");
    monomer_mol.addConnection(m1, m2, ConnectionType::FORWARD);

    // Verify it's dative at the monomer level
    CHECK(monomer_mol.getBondWithIdx(0)->getBondType() == Bond::DATIVE);

    // Convert to atomistic
    auto atomistic = toAtomistic(monomer_mol);

    // All bonds should be single (or aromatic, etc.) - no dative bonds
    for (auto bond : atomistic->bonds()) {
      CHECK(bond->getBondType() != Bond::DATIVE);
    }
  }
}

TEST_CASE("HybridMonomerMol") {
  SECTION("isMonomer") {
    // Test that isMonomer correctly distinguishes monomers from atoms
    MonomerMol mol;

    // Add a monomer
    auto monomerIdx = mol.addMonomer("A", 1, "PEPTIDE1");

    // Add a regular atom (carbon)
    Atom carbon(6);
    auto atomIdx = mol.addAtom(&carbon);

    CHECK(mol.isMonomer(monomerIdx) == true);
    CHECK(mol.isMonomer(atomIdx) == false);
  }

  SECTION("addBondAtomToAtom") {
    // Test that addBond works for atom-to-atom bonds
    MonomerMol mol;

    // Add two regular atoms
    Atom carbon1(6);
    Atom carbon2(6);
    auto idx1 = mol.addAtom(&carbon1);
    auto idx2 = mol.addAtom(&carbon2);

    // Should work - both are regular atoms
    auto numBonds = mol.addBond(idx1, idx2, Bond::SINGLE);
    CHECK(numBonds == 1);
    CHECK(mol.getNumBonds() == 1);
  }

  SECTION("addBondMonomerToAtom") {
    // Test addBond between a monomer and an atom with attachment point
    MonomerMol mol;

    // Add a monomer (Alanine)
    auto monomerIdx = mol.addMonomer("A", 1, "PEPTIDE1");

    // Add a regular atom (nitrogen, simulating an N-terminal modification)
    Atom nitrogen(7);
    auto atomIdx = mol.addAtom(&nitrogen);

    // Add bond from monomer's R1 attachment point to the atom
    auto numBonds = mol.addBond(monomerIdx, atomIdx, 1, Bond::SINGLE);
    CHECK(numBonds == 1);
    CHECK(mol.getNumBonds() == 1);

    // Check that the linkage property is set correctly
    auto bond = mol.getBondWithIdx(0);
    CHECK(bond->getProp<std::string>(LINKAGE) == "R1-");
  }

  SECTION("addBondMonomerToAtomReversedOrder") {
    // Test that addBond works even if atom comes before monomer in args
    MonomerMol mol;

    auto monomerIdx = mol.addMonomer("A", 1, "PEPTIDE1");
    Atom oxygen(8);
    auto atomIdx = mol.addAtom(&oxygen);

    // Pass atom first, monomer second - should still work
    auto numBonds = mol.addBond(atomIdx, monomerIdx, 2, Bond::SINGLE);
    CHECK(numBonds == 1);

    // Linkage should still be set correctly
    auto bond = mol.getBondWithIdx(0);
    CHECK(bond->getProp<std::string>(LINKAGE) == "R2-");
  }

  SECTION("addBondMonomerToMonomerThrows") {
    // Test that addBond throws for monomer-to-monomer bonds
    MonomerMol mol;

    auto m1 = mol.addMonomer("A", 1, "PEPTIDE1");
    auto m2 = mol.addMonomer("G");

    // Should throw - must use addConnection for monomer-to-monomer
    CHECK_THROWS_AS(mol.addBond(m1, m2, Bond::SINGLE), std::invalid_argument);
  }

  SECTION("addBondMixedWithoutAttachmentThrows") {
    // Test that addBond without attachment point throws for mixed bonds
    MonomerMol mol;

    auto monomerIdx = mol.addMonomer("A", 1, "PEPTIDE1");
    Atom carbon(6);
    auto atomIdx = mol.addAtom(&carbon);

    // Should throw - need to specify attachment point
    CHECK_THROWS_AS(mol.addBond(monomerIdx, atomIdx, Bond::SINGLE),
                    std::invalid_argument);
  }

  SECTION("hybridMoleculeConstruction") {
    // Build a hybrid molecule: monomer - atom - atom - monomer
    // This could represent a peptide with a non-standard linker
    MonomerMol mol;

    // Add monomers
    auto m1 = mol.addMonomer("A", 1, "PEPTIDE1");
    auto m2 = mol.addMonomer("G", 2, "PEPTIDE1");

    // Add linker atoms (e.g., -CH2-CH2-)
    Atom c1(6);
    Atom c2(6);
    auto a1 = mol.addAtom(&c1);
    auto a2 = mol.addAtom(&c2);

    // Connect: monomer1 -R2-> atom1 - atom2 <-R1- monomer2
    mol.addBond(m1, a1, 2, Bond::SINGLE);  // Alanine R2 to first carbon
    mol.addBond(a1, a2, Bond::SINGLE);      // Carbon to carbon
    mol.addBond(m2, a2, 1, Bond::SINGLE);  // Glycine R1 to second carbon

    CHECK(mol.getNumAtoms() == 4);
    CHECK(mol.getNumBonds() == 3);

    // Verify monomer detection
    CHECK(mol.isMonomer(m1) == true);
    CHECK(mol.isMonomer(m2) == true);
    CHECK(mol.isMonomer(a1) == false);
    CHECK(mol.isMonomer(a2) == false);

    // Verify linkage properties on monomer-atom bonds
    auto bond1 = mol.getBondBetweenAtoms(m1, a1);
    auto bond2 = mol.getBondBetweenAtoms(a1, a2);
    auto bond3 = mol.getBondBetweenAtoms(m2, a2);

    CHECK(bond1->getProp<std::string>(LINKAGE) == "R2-");
    CHECK_FALSE(bond2->hasProp(LINKAGE));  // atom-atom bond has no linkage
    CHECK(bond3->getProp<std::string>(LINKAGE) == "R1-");
  }
}