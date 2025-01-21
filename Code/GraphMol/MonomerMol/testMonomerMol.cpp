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
        chains.insert(get_polymer_id(atom));
    }
    std::stringstream fasta;
    bool first = true;
    for (auto chain_id : chains) {
        if (!first) {
            fasta << "\n";
        }
        first = false;
        fasta << ">Chain " << chain_id << "\n"; // add title
        auto chain = get_polymer(mol, chain_id);
        for (auto atom : chain.atoms) {
            // Really there should be a lookup to ensure that these are all 1 letter
            // And a sort by residue number & connectivity
            const auto* res_info = static_cast<const RDKit::AtomPDBResidueInfo*>(mol.getAtomWithIdx(atom)->getMonomerInfo());
            fasta << res_info->getResidueName();
        }
    }
    return fasta.str();
}

TEST_CASE("FASTAConversions") {
  SECTION("SIMPLE") {
    // Build MonomerMol
    RWMol monomer_mol;
    add_monomer(monomer_mol, "R", 1, "A");
    add_monomer(monomer_mol, "D");
    add_monomer(monomer_mol, "K");
    add_monomer(monomer_mol, "I");
    add_monomer(monomer_mol, "T");

    add_connection(monomer_mol, 0, 1, ConnectionType::FORWARD);
    add_connection(monomer_mol, 1, 2, ConnectionType::FORWARD);
    add_connection(monomer_mol, 2, 3, ConnectionType::FORWARD);
    add_connection(monomer_mol, 3, 4, ConnectionType::FORWARD);

    std::cerr << to_fasta(monomer_mol);
    CHECK(std::string(">Chain A\nRDKIT") == to_fasta(monomer_mol));
  }

  SECTION("MultipleChains") {
    RWMol monomer_mol;
    auto midx1 = add_monomer(monomer_mol, "R", 1, "A");
    auto midx2 = add_monomer(monomer_mol, "D");

    auto midx3 = add_monomer(monomer_mol, "K", 1, "B");
    auto midx4 = add_monomer(monomer_mol, "I");
    auto midx5 = add_monomer(monomer_mol, "T");

    add_connection(monomer_mol, midx1, midx2, ConnectionType::FORWARD);
    add_connection(monomer_mol, midx3, midx4, ConnectionType::FORWARD);
    add_connection(monomer_mol, midx4, midx5, ConnectionType::FORWARD);

    CHECK(std::string(">>Chain A\nRD\n>Chain B\nKIT") == to_fasta(monomer_mol));
  }

  SECTION("UsingFASTAReader") {
    std::string seq = "CGCGAATTACCGCG";
    // the sequence parser creates an atomistic molecules
    auto atomistic_mol = SequenceToMol(seq);
    CHECK(atomistic_mol);

    auto monomer_mol = atomisticToMonomerMol(*atomistic_mol);
    CHECK(std::string(">>Chain A\nRD\n>Chain B\nKIT") == to_fasta(*monomer_mol));
  }
}
