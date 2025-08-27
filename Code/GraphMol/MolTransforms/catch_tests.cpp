//
//  Copyright (c) 2024 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
///
#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <algorithm>
#include <regex>

using namespace RDKit;
using std::unique_ptr;

TEST_CASE("chirality flip in canonicalizeMolConformer") {
  SECTION("first discovery") {
    auto m1 = R"CTAB(
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 11 10 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.533948 2.465972 -1.220733 0
M  V30 2 C -1.361222 1.165424 -1.939720 0 CFG=2
M  V30 3 C -2.694145 0.606596 -2.153754 0
M  V30 4 N -3.766061 0.178135 -2.319155 0
M  V30 5 C -0.566772 0.289201 -1.072783 0
M  V30 6 C 0.084234 -0.420637 -0.345400 0
M  V30 7 C 0.867082 -1.237691 0.508367 0
M  V30 8 C 1.549171 -1.919699 1.232613 0
M  V30 9 C 2.339643 -2.731169 2.106347 0
M  V30 10 C 2.993349 -3.411935 2.854815 0
M  V30 11 O 3.747921 -4.178375 3.716269 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=1
M  V30 2 1 2 3
M  V30 3 3 3 4
M  V30 4 1 2 5
M  V30 5 3 5 6
M  V30 6 1 6 7
M  V30 7 3 7 8
M  V30 8 1 8 9
M  V30 9 3 9 10
M  V30 10 1 10 11
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m1);
    auto m2 = R"CTAB(
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 11 10 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -3.278200 0.712725 -0.066592 0
M  V30 2 C -2.682025 -0.648154 -0.271665 0 CFG=1
M  V30 3 C -2.915612 -1.059871 -1.656844 0
M  V30 4 N -3.080403 -1.390524 -2.739848 0
M  V30 5 C -1.232937 -0.539032 -0.023304 0
M  V30 6 C -0.029030 -0.463778 0.167869 0
M  V30 7 C 1.366183 -0.398546 0.377840 0
M  V30 8 C 2.562372 -0.361417 0.548269 0
M  V30 9 C 3.973502 -0.322846 0.772691 0
M  V30 10 C 5.164230 -0.269523 0.960648 0
M  V30 11 O 6.516239 -0.218493 1.152837 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=3
M  V30 2 1 2 3
M  V30 3 3 3 4
M  V30 4 1 2 5
M  V30 5 3 5 6
M  V30 6 1 6 7
M  V30 7 3 7 8
M  V30 8 1 8 9
M  V30 9 3 9 10
M  V30 10 1 10 11
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m2);

    auto m1_cp = ROMol(*m1);
    auto m2_cp = ROMol(*m2);
    MolTransforms::canonicalizeConformer(m1_cp.getConformer());
    MolTransforms::canonicalizeConformer(m2_cp.getConformer());

    MolOps::assignChiralTypesFrom3D(m1_cp);
    MolOps::assignChiralTypesFrom3D(m2_cp);
    CHECK(m1_cp.getAtomWithIdx(1)->getChiralTag() ==
          m1->getAtomWithIdx(1)->getChiralTag());
    CHECK(m2_cp.getAtomWithIdx(1)->getChiralTag() ==
          m2->getAtomWithIdx(1)->getChiralTag());
  }
}

TEST_CASE("Other flips") {
  // These came from PubChem's 3D conformer collection:
  // https://ftp.ncbi.nlm.nih.gov/pubchem/Compound_3D/10_conf_per_cmpd
  // They are: esomeprazole, 64478
  std::vector<ROMOL_SPTR> mols{
      R"(FC1(F)[C@]2(F)C(F)(F)[C@]3(F)C(F)(F)[C@@]1(F)C(F)(F)[C@](F)(C2(F)F)C3(F)F |(0.9902,1.5037,-2.1461;1.0672,0.4066,-1.3547;2.1203,-0.3186,-1.8028;1.2929,0.8217,0.1055;2.4333,1.5466,0.1986;0.1063,1.6692,0.5847;0.0063,2.7967,-0.1599;0.3037,2.0687,1.8642;-1.1864,0.8486,0.4796;-2.233,1.5972,0.9026;-1.4104,0.4334,-0.9809;-1.5471,1.5312,-1.7633;-2.5642,-0.2677,-1.096;-0.225,-0.4149,-1.4612;-0.4234,-0.7809,-2.7501;-0.1063,-1.6693,-0.5847;0.9185,-2.4444,-1.0141;-1.2285,-2.4211,-0.6902;0.1185,-1.2554,0.8761;0.2231,-2.3629,1.6489;1.4104,-0.4335,0.9809;2.472,-1.1788,0.5893;1.6392,-0.0847,2.27;-1.0672,-0.4066,1.3548;-2.2125,-1.1279,1.2961;-0.8981,-0.0571,2.6528),wD:3.3,13.14,wU:8.8,18.19|)"_smiles,
      R"(COc1ccc2nc([S@@](=O)Cc3ncc(C)c(OC)c3C)[nH]c2c1 |(0.641362,7.44566,-1.90542;1.74659,6.85341,-1.22739;1.60806,5.55453,-0.835526;2.65944,4.92299,-0.166852;2.55231,3.59036,0.250442;1.35968,2.89747,-0.0183527;1.00287,1.60246,0.279491;-0.215269,1.46666,-0.194606;-1.16791,-0.0427214,-0.0819799;-1.95615,-0.218445,-1.36497;0.169454,-1.27136,-0.106371;-0.362868,-2.6725,-0.0211404;-0.63428,-3.26081,-1.20402;-1.11383,-4.52338,-1.14388;-1.3311,-5.21502,0.0368067;-1.86296,-6.60553,0.0061196;-1.0376,-4.56935,1.23282;-1.23419,-5.20992,2.42029;-0.141111,-5.9658,2.93629;-0.542453,-3.2702,1.21184;-0.218994,-2.55639,2.4807;-0.66478,2.61454,-0.785806;0.329678,3.55345,-0.686829;0.414062,4.87992,-1.11164),wD:8.9|)"_smiles,
  };
  for (auto &mol : mols) {
    REQUIRE(mol);
    auto mol_cp = ROMol(*mol);
    auto rms1 = MolAlign::alignMol(mol_cp, *mol);
    MolTransforms::canonicalizeMol(mol_cp);
    auto rms2 = MolAlign::alignMol(mol_cp, *mol);
    CHECK_THAT(rms1 - rms2, Catch::Matchers::WithinAbs(0.0, 0.001));
  }
}