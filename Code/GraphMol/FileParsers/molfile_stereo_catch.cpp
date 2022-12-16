//
//  Copyright (C) 2022 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RDGeneral/test.h"
#include "catch.hpp"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolWriters.h>

using namespace RDKit;

TEST_CASE("Github #5863: failure in WedgeMolBonds") {
  SECTION("as reported") {
    auto mol =
        "C[C@]12C=CC(=O)C=C1CC[C@@H]1[C@@H]2[C@@H](O)C[C@@]2(C)[C@H]1CC[C@]2(O)C(=O)CO |(1.94354,1.43772,;2.70098,0.14301,;3.44351,1.44633,;4.94349,1.45494,;5.70093,0.160228,;7.20091,0.168838,;4.9584,-1.14309,;3.45842,-1.1517,;2.71589,-2.45502,;1.21592,-2.46363,;0.458474,-1.16892,;1.20101,0.1344,;0.443562,1.42911,;1.18609,2.73243,;-1.05641,1.4205,;-1.79895,0.117181,;-2.53364,1.42493,;-1.0415,-1.17753,;-2.03878,-2.29799,;-3.41258,-1.69576,;-3.26435,-0.203103,;-3.6102,1.25648,;-4.76433,-0.211712,;-5.50686,-1.51503,;-5.52177,1.083,;-7.02175,1.07439,)|"_smiles;
    REQUIRE(mol);
    unsigned int radius = 2;
    unsigned int atomId = 2;
    auto env = findAtomEnvironmentOfRadiusN(*mol, radius + 1, atomId);
    std::unique_ptr<ROMol> frag(Subgraphs::pathToSubmol(*mol, env));
    REQUIRE(frag);
    WedgeMolBonds(*frag, &frag->getConformer());
    CHECK(frag->getBondBetweenAtoms(9, 10)->getBondDir() !=
          Bond::BondDir::NONE);
    std::cerr << MolToV3KMolBlock(*frag) << std::endl;
  }
}