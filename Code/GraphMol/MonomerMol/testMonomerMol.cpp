//
//  Copyright (C) 2002-2024 Rachel Walker and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
//  Copyright (C) 2002-2021 Greg Landrum and other RDKit contributors
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


#include <GraphMol/MonomerMol/MonomerMol.h>
#include <GraphMol/MonomerInfo.h>


#include <string>

using namespace RDKit;

std::string to_fasta(const ROMol& mol)
{
    // make a set of all the chains in the molecule
    std::set<std::string> chains;
    for (auto atom : mol.atoms()) {
        chains.insert(get_polymer_id(atom));
    }
    std::stringstream fasta;
    for (auto chain_id : chains) {
        fasta << ">Chain " << chain_id << "\n"; // add title
        auto chain = get_polymer(mol, chain_id);
        for (auto atom : chain.atoms) {
            // Really there should be a lookup to ensure that these are all 1 letter
            // And a sort by residue number & connectivity
            const auto* res_info = static_cast<const RDKit::AtomPDBResidueInfo*>(mol.getAtomWithIdx(atom)->getMonomerInfo());
            fasta << res_info->getResidueName();
        }
    }
    // fasta << ">Chain " << polymer_id << "\n"; // add title
    // // Really there should be a lookup to ensure that these are all 1 letter
    // // And a sort by residue number & connectivity
    // for (auto a: mol.atoms()) {
    //     auto* res_info = static_cast<RDKit::AtomPDBResidueInfo*>(a->getMonomerInfo());
    //     fasta << res_info->getResidueName();
    // }
    return fasta.str();
}

RWMol from_fasta(std::string fasta)
{
    std::stringstream fasta_stream(fasta);
    std::string name;
    std::getline(fasta_stream, name);
    RWMol monomer_mol;
    if (name.size() > 2) {
        monomer_mol.setProp("_Name", name.substr(2, std::string::npos));
    }
    char c;
    while (fasta_stream >> c) {
        if (monomer_mol.getNumAtoms() == 0) {
            add_monomer(monomer_mol, {&c, 1}, 1, "A");
        } else {
            add_monomer(monomer_mol, {&c, 1});
        }
    }
    return monomer_mol;
}


TEST_CASE("MonomerMolToFASTA") {
  SECTION("SIMPLE") {

    // Build the sequence RDKIT
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

  SECTION("AtomisticToMonomer") {
    
  }
}
