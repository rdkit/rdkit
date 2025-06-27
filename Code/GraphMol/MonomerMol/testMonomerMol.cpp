//
//  Copyright (C) 2002-2024 Rachel Walker and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
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
std::string to_fasta(const ROMol& mol)
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

void neutralize_atoms(RDKit::ROMol& mol)
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

TEST_CASE("FASTAConversions") {
  SECTION("SIMPLE") {
    // Build MonomerMol with a single chain
    RWMol monomer_mol;
    addMonomer(monomer_mol, "R", 1, "PEPTIDE1");
    addMonomer(monomer_mol, "D");
    addMonomer(monomer_mol, "K");
    addMonomer(monomer_mol, "I");
    addMonomer(monomer_mol, "T");

    addConnection(monomer_mol, 0, 1, ConnectionType::FORWARD);
    addConnection(monomer_mol, 1, 2, ConnectionType::FORWARD);
    addConnection(monomer_mol, 2, 3, ConnectionType::FORWARD);
    addConnection(monomer_mol, 3, 4, ConnectionType::FORWARD);

    CHECK(std::string(">Chain PEPTIDE1\nRDKIT") == to_fasta(monomer_mol));
  }

  SECTION("MultipleChains") {
    // Build MonomerMol with two chains
    RWMol monomer_mol;
    auto midx1 = addMonomer(monomer_mol, "R", 1, "A");
    auto midx2 = addMonomer(monomer_mol, "D");

    auto midx3 = addMonomer(monomer_mol, "K", 1, "B");
    auto midx4 = addMonomer(monomer_mol, "I");
    auto midx5 = addMonomer(monomer_mol, "T");

    addConnection(monomer_mol, midx1, midx2, ConnectionType::FORWARD);
    addConnection(monomer_mol, midx3, midx4, ConnectionType::FORWARD);
    addConnection(monomer_mol, midx4, midx5, ConnectionType::FORWARD);

    CHECK(std::string(">Chain A\nRD\n>Chain B\nKIT") == to_fasta(monomer_mol));
  }

  SECTION("UsingSequenceReader") {
    std::string seq = "CGCGAATTACCGCG";
    // the sequence parser creates an atomistic molecule
    auto atomistic_mol = SequenceToMol(seq);
    CHECK(atomistic_mol);

    // PDB info is used to convert the atomistic sturcture into a MonomerMol
    auto monomer_mol = atomisticToMonomerMol(*atomistic_mol);
    CHECK(std::string(">Chain PEPTIDE1\nCGCGAATTACCGCG") == to_fasta(*monomer_mol));
  }
}

TEST_CASE("Conversions") {
  SECTION("MonomerMolToAtomistic") {
    std::string seq = "CGCGA";
    auto atomistic_mol = SequenceToMol(seq);

    RWMol monomer_mol;
    auto midx1 = addMonomer(monomer_mol, "C", 1, "PEPTIDE1");
    auto midx2 = addMonomer(monomer_mol, "G");
    auto midx3 = addMonomer(monomer_mol, "C");
    auto midx4 = addMonomer(monomer_mol, "G");
    auto midx5 = addMonomer(monomer_mol, "A");

    addConnection(monomer_mol, midx1, midx2, ConnectionType::FORWARD);
    addConnection(monomer_mol, midx2, midx3, ConnectionType::FORWARD);
    addConnection(monomer_mol, midx3, midx4, ConnectionType::FORWARD);
    addConnection(monomer_mol, midx4, midx5, ConnectionType::FORWARD);
    auto atomistic_mol2 = monomerMolToAtomsitic(monomer_mol);

    // atomistic structure is same as using sequence parser
    std::string smi1 = MolToSmiles(*atomistic_mol);
    std::string smi2 = MolToSmiles(*atomistic_mol2);
    CHECK(smi1 == smi2);
  }

  SECTION("AtomisticToMonomerMol") {
    std::string pdbfile = getenv("RDBASE");
    pdbfile += "/Code/GraphMol/MonomerMol/test_data/1dng.pdb";
    auto mol = PDBFileToMol(pdbfile);
    auto monomer_mol = atomisticToMonomerMol(*mol);
    CHECK(std::string(">Chain PEPTIDE1\nQAPAYEEAAEELAKS") == to_fasta(*monomer_mol));

    auto atomistic_mol2 = monomerMolToAtomsitic(*monomer_mol);
    neutralize_atoms(*mol);
    neutralize_atoms(*atomistic_mol2);
    std::string smi1 = MolToSmiles(*mol);
    std::string smi2 = MolToSmiles(*atomistic_mol2);
    CHECK(smi1 == smi2);
  }
}