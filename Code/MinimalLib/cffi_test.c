/*
  Copyright (C) 2021 Greg Landrum

  @@ All Rights Reserved @@
  This file is part of the RDKit.
  The contents are covered by the terms of the BSD license
  which is included in the file license.txt, found at the root
  of the RDKit source tree.
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "cffiwrapper.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>

void test_io() {
  char *pkl;
  size_t pkl_size;
  size_t pkl2_size;
  char *pkl2;

  printf("--------------------------\n");
  printf("  test_io\n");

  pkl = get_mol("c1cc(O)ccc1", &pkl_size, "");
  assert(pkl);
  assert(pkl_size > 0);

  char *smiles = get_smiles(pkl, pkl_size, NULL);
  assert(!strcmp(smiles, "Oc1ccccc1"));
  free(smiles);
  smiles = NULL;

  smiles = get_cxsmiles(pkl, pkl_size, NULL);
  assert(!strcmp(smiles, "Oc1ccccc1"));
  free(smiles);
  smiles = NULL;

  char *json = get_json(pkl, pkl_size, NULL);
  assert(strstr(json, "commonchem"));

  pkl2 = get_mol(json, &pkl2_size, "");
  assert(pkl2);
  assert(pkl2_size > 0);
  smiles = get_smiles(pkl2, pkl2_size, NULL);
  assert(!strcmp(smiles, "Oc1ccccc1"));
  free(smiles);
  smiles = NULL;
  free(pkl2);
  pkl2 = NULL;
  free(json);
  json = NULL;

  //---------
  // failures

  // kekulization
  pkl2 = get_mol("c1cccc1", &pkl2_size, "");
  assert(pkl2 == NULL);
  assert(!pkl2_size);

  // bad input
  pkl2 = get_mol("foo", &pkl2_size, "");
  assert(pkl2 == NULL);
  assert(!pkl2_size);

  // valence
  pkl2 = get_mol("CO(C)C", &pkl2_size, "");
  assert(pkl2 == NULL);
  assert(!pkl2_size);

  // bad molblock
  pkl2 = get_mol(
      "  Mrv1921 05042106432D\n\
\n\
  2  1  0  0  0  0            999 V2000\n\
   -7.3214    3.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -6.6070    4.1625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
  1  2  1  0  0  0  0\n\
M  END",
      &pkl2_size, "");
  assert(pkl2 == NULL);
  assert(!pkl2_size);

  //---------
  // options
  pkl2 = get_mol("[H]C", &pkl2_size, "{\"removeHs\":false}");
  assert(pkl2 != NULL);
  assert(pkl2_size);
  smiles = get_smiles(pkl2, pkl2_size, NULL);
  assert(!strcmp(smiles, "[H]C"));
  free(smiles);
  smiles = NULL;
  free(pkl2);
  pkl2 = NULL;

  pkl2 = get_mol("[H]C", &pkl2_size, NULL);
  assert(pkl2 != NULL);
  assert(pkl2_size);
  free(pkl2);
  pkl2 = NULL;

  smiles = get_smiles(pkl, pkl_size, "{\"canonical\":false}");
  assert(!strcmp(smiles, "c1cc(O)ccc1"));
  free(smiles);
  smiles = NULL;

  //---------
  // mol block
  char *molblock = get_molblock(pkl, pkl_size, NULL);
  pkl2 = get_mol(molblock, &pkl2_size, "");
  assert(pkl2);
  assert(pkl2_size > 0);
  smiles = get_smiles(pkl2, pkl2_size, NULL);
  assert(!strcmp(smiles, "Oc1ccccc1"));
  free(smiles);
  smiles = NULL;
  free(pkl2);
  pkl2 = NULL;
  free(molblock);
  molblock = NULL;

  molblock = get_v3kmolblock(pkl, pkl_size, NULL);
  pkl2 = get_mol(molblock, &pkl2_size, "");
  assert(pkl2);
  assert(pkl2_size > 0);
  smiles = get_smiles(pkl2, pkl2_size, NULL);
  assert(!strcmp(smiles, "Oc1ccccc1"));
  free(smiles);
  smiles = NULL;
  free(pkl2);
  pkl2 = NULL;
  free(molblock);
  molblock = NULL;

  molblock = get_v3kmolblock(pkl, pkl_size, "{\"kekulize\":false}");
  assert(strstr(molblock, "M  V30 1 4 1 2"));
  free(molblock);
  molblock = NULL;

  //---------
  // InChI
  char *inchi = get_inchi(pkl, pkl_size, NULL);
  assert(!strcmp(inchi, "InChI=1S/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H"));
  free(inchi);
  inchi = get_inchikey_for_inchi("InChI=1S/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H");
  assert(!strcmp(inchi, "ISWSIDIOOBJBQZ-UHFFFAOYSA-N"));
  free(inchi);

  inchi = get_inchi(pkl, pkl_size, "{\"options\":\"/FixedH\"}");
  assert(!strcmp(inchi, "InChI=1/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H"));
  free(inchi);

  molblock = get_molblock(pkl, pkl_size, NULL);
  inchi = get_inchi_for_molblock(molblock, NULL);
  assert(!strcmp(inchi, "InChI=1S/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H"));
  free(inchi);
  inchi = get_inchi_for_molblock(molblock, "{\"options\":\"/FixedH\"}");
  assert(!strcmp(inchi, "InChI=1/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H"));
  free(inchi);
  free(molblock);

  //---------
  // queries
  char *smarts = get_smarts(pkl, pkl_size, NULL);
  assert(!strcmp(smarts, "[#6]1:[#6]:[#6](-[#8]):[#6]:[#6]:[#6]:1"));

  pkl2 = get_qmol(smarts, &pkl2_size, "");
  assert(pkl2);
  free(smarts);
  smarts = get_smarts(pkl2, pkl2_size, NULL);
  assert(!strcmp(smarts, "[#6]1:[#6]:[#6](-[#8]):[#6]:[#6]:[#6]:1"));
  free(smarts);
  free(pkl2);

  free(pkl);
  pkl = NULL;
  printf("  done\n");
  printf("--------------------------\n");
}

void test_svg() {
  char *pkl;
  size_t pkl_size;

  printf("--------------------------\n");
  printf("  test_svg\n");

  pkl = get_mol("c1cc(O)ccc1", &pkl_size, "");
  assert(pkl);
  assert(pkl_size > 0);

  char *svg = get_svg(pkl, pkl_size, "{\"width\":350,\"height\":300}");
  assert(strstr(svg, "width='350px'"));
  assert(strstr(svg, "height='300px'"));
  assert(strstr(svg, "</svg>"));
  free(svg);

  svg = get_svg(pkl, pkl_size, "{\"atoms\": [0, 1, 2], \"width\":127}");
  assert(strstr(svg, "fill:#FF7F7F"));
  assert(strstr(svg, "width='127px'"));
  assert(strstr(svg, "</svg>"));
  free(svg);

  free(pkl);
  pkl = NULL;
  printf("  done\n");
  printf("--------------------------\n");
}

void test_flexicanvas() {
  char *pkl;
  size_t pkl_size;

  printf("--------------------------\n");
  printf("  test_flexicanvas\n");

  pkl = get_mol("CCCC", &pkl_size, "");
  assert(pkl);
  assert(pkl_size > 0);

  char *svg = get_svg(pkl, pkl_size, "{\"width\":-1,\"height\":-1}");
  assert(strstr(svg, "width='95px'"));
  assert(strstr(svg, "height='21px'"));
  assert(strstr(svg, "</svg>"));
  free(svg);

  free(pkl);
  pkl = NULL;
  printf("  done\n");
  printf("--------------------------\n");
}

void test_substruct() {
  printf("--------------------------\n");
  printf("  test_substruct\n");

  char *mpkl;
  size_t mpkl_size;
  mpkl = get_mol("Cl[C@H](F)C[C@H](F)Cl", &mpkl_size, "");
  char *qpkl;
  size_t qpkl_size;
  qpkl = get_qmol("Cl[C@@H](F)C", &qpkl_size, "");

  char *json = get_substruct_match(mpkl, mpkl_size, qpkl, qpkl_size, "");
  assert(!strcmp(json, "{\"atoms\":[0,1,2,3],\"bonds\":[0,1,2]}"));
  free(json);
  json = get_substruct_matches(mpkl, mpkl_size, qpkl, qpkl_size, "");
  assert(!strcmp(json,
                 "[{\"atoms\":[0,1,2,3],\"bonds\":[0,1,2]},{\"atoms\":[6,4,5,3]"
                 ",\"bonds\":[5,4,3]}]"));
  free(json);

  // make sure using parameters works
  json = get_substruct_match(mpkl, mpkl_size, qpkl, qpkl_size,
                             "{\"useChirality\":true}");
  assert(!strcmp(json, "{\"atoms\":[6,4,5,3],\"bonds\":[5,4,3]}"));
  free(json);
  json = get_substruct_matches(mpkl, mpkl_size, qpkl, qpkl_size,
                               "{\"useChirality\":true}");
  assert(!strcmp(json, "[{\"atoms\":[6,4,5,3],\"bonds\":[5,4,3]}]"));
  free(json);

  free(mpkl);
  mpkl = NULL;
  free(qpkl);
  qpkl = NULL;
  printf("  done\n");
  printf("--------------------------\n");
}

void test_descriptors() {
  printf("--------------------------\n");
  printf("  test_descriptors\n");

  char *mpkl;
  size_t mpkl_size;
  mpkl = get_mol("c1nccc(O)c1", &mpkl_size, "");

  char *descrs = get_descriptors(mpkl, mpkl_size);
  assert(strstr(descrs, "lipinskiHBA"));
  assert(strstr(descrs, "NumAliphaticRings"));
  assert(strstr(descrs, "chi3v"));
  assert(strstr(descrs, "amw"));
  free(descrs);

  free(mpkl);
  mpkl = NULL;
  printf("  done\n");
  printf("--------------------------\n");
}

void test_fingerprints() {
  printf("--------------------------\n");
  printf("  test_fingerprints\n");

  char *mpkl;
  size_t mpkl_size;
  mpkl = get_mol("c1nccc(O)c1", &mpkl_size, "");
  const char *mfp_json = "{\"radius\":2,\"nBits\":1024}";
  char *fp = get_morgan_fp(mpkl, mpkl_size, mfp_json);
  assert(strlen(fp) == 1024);
  free(fp);
  fp = get_morgan_fp(mpkl, mpkl_size, "{\"radius\":1,\"nBits\":64}");
  assert(!strcmp(
      fp, "0011000000100000010000100000000000001001010000000000000000100000"));
  free(fp);
  fp = get_morgan_fp(mpkl, mpkl_size, "{\"radius\":2,\"nBits\":64}");
  assert(!strcmp(
      fp, "0011000000100000010000100000000000001001010000000010000010100001"));
  free(fp);
  size_t nbytes;
  fp = get_morgan_fp_as_bytes(mpkl, mpkl_size, &nbytes,
                              "{\"radius\":2,\"nBits\":64}");
  assert(nbytes == 8);
  free(fp);

  fp = get_rdkit_fp(mpkl, mpkl_size, "{\"nBits\":64}");
  assert(!strcmp(
      fp, "1111011000111100011011011111011001111111110010000111000011111111"));
  free(fp);
  fp = get_rdkit_fp(mpkl, mpkl_size, "{\"nBits\":64,\"maxPath\":4}");
  assert(!strcmp(
      fp, "1111011000110000010001010111000001101011100010000001000011100011"));
  free(fp);
  fp = get_rdkit_fp_as_bytes(mpkl, mpkl_size, &nbytes, "{\"nBits\":64}");
  assert(nbytes == 8);
  free(fp);

  fp = get_pattern_fp(mpkl, mpkl_size, "{\"nBits\":64}");
  assert(!strcmp(
      fp, "1011111111111111110011011111011001011110111101111110101111010011"));
  free(fp);
  fp = get_pattern_fp_as_bytes(mpkl, mpkl_size, &nbytes, "{\"nBits\":64}");
  assert(nbytes == 8);
  free(fp);

  fp = get_pattern_fp(mpkl, mpkl_size,
                      "{\"nBits\":64,\"tautomericFingerprint\":true}");
  assert(!strcmp(
      fp, "1011111111111111110011111111011101011111111101111111111111011011"));
  free(fp);

  free(mpkl);
  mpkl = NULL;
  printf("  done\n");
  printf("--------------------------\n");
}

void test_modifications() {
  printf("--------------------------\n");
  printf("  test_modifications\n");
  char *mpkl;
  size_t mpkl_size;
  mpkl = get_mol("CCC", &mpkl_size, "");

  assert(add_hs(&mpkl, &mpkl_size) > 0);
  char *ctab = get_molblock(mpkl, mpkl_size, NULL);
  assert(strstr(ctab, " H "));
  free(ctab);

  assert(remove_all_hs(&mpkl, &mpkl_size) > 0);
  ctab = get_molblock(mpkl, mpkl_size, NULL);
  assert(!strstr(ctab, " H "));
  free(ctab);

  free(mpkl);
  mpkl = NULL;
  printf("  done\n");
  printf("--------------------------\n");
}

void test_coords() {
  printf("--------------------------\n");
  printf("  test_coords\n");
  char *mpkl;
  size_t mpkl_size;
  mpkl = get_mol("C1CNC1CC", &mpkl_size, "");

  char *cxsmi = get_cxsmiles(mpkl, mpkl_size, NULL);
  // no cxsmiles yet
  assert(!strstr(cxsmi, "|"));

  prefer_coordgen(0);
  set_2d_coords(&mpkl, &mpkl_size);
  free(cxsmi);
  cxsmi = get_cxsmiles(mpkl, mpkl_size, NULL);
  // since we have coords there's something there:
  assert(strstr(cxsmi, "|"));
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  prefer_coordgen(1);
  set_2d_coords(&mpkl, &mpkl_size);
  char *cxsmi2 = get_cxsmiles(mpkl, mpkl_size, NULL);
  assert(strstr(cxsmi2, "|"));
  assert(strcmp(cxsmi, cxsmi2));
  free(cxsmi2);
#endif
  free(cxsmi);

  // aligned
  char *tpkl;
  size_t tpkl_size;
  tpkl = get_mol(
      "\n\
  Mrv2102 04062106432D\n\
\n\
  0  0  0     0  0            999 V3000\n\
M  V30 BEGIN CTAB\n\
M  V30 COUNTS 6 6 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 R# 0 2.1032 0 0 RGROUPS=(1 1)\n\
M  V30 2 C 0 0.5632 0 0\n\
M  V30 3 C 1.0889 -0.5258 0 0\n\
M  V30 4 C 0 -1.6147 0 0\n\
M  V30 5 N -1.0889 -0.5258 0 0\n\
M  V30 6 R# 2.6289 -0.5258 0 0 RGROUPS=(1 2)\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 1 2\n\
M  V30 2 1 2 3\n\
M  V30 3 1 3 4\n\
M  V30 4 1 4 5\n\
M  V30 5 1 2 5\n\
M  V30 6 1 3 6\n\
M  V30 END BOND\n\
M  V30 END CTAB\n\
M  END\n",
      &tpkl_size, "");
  assert(!set_2d_coords_aligned(&mpkl, &mpkl_size, tpkl, tpkl_size,
                                "{\"acceptFailure\":false}"));
  assert(
      set_2d_coords_aligned(&mpkl, &mpkl_size, tpkl, tpkl_size,
                            "{\"allowRGroups\":true,\"acceptFailure\":false}"));
  free(tpkl);

  // Github #4121: set_2d_coords_aligned crashes if template mol has no
  // conformer
  tpkl = get_mol("C1CNC1", &tpkl_size, NULL);
  assert(!set_2d_coords_aligned(&mpkl, &mpkl_size, tpkl, tpkl_size, ""));
  free(tpkl);

  // 3D
  assert(add_hs(&mpkl, &mpkl_size));
  assert(set_3d_coords(&mpkl, &mpkl_size, "") > 0);
  char *cxsmi3 = get_cxsmiles(mpkl, mpkl_size, NULL);
  assert(set_3d_coords(&mpkl, &mpkl_size, "{\"randomSeed\":123}") > 0);
  assert(set_3d_coords(&mpkl, &mpkl_size,
                       "{\"randomSeed\":123,\"coordMap\":{\"3\":[0,0,0],\"4\":["
                       "0,0,1.5],\"5\":[0,1.5,1.5]}}") > 0);
  cxsmi = get_cxsmiles(mpkl, mpkl_size, NULL);
  // since we have coords there's something there:
  assert(strstr(cxsmi, "|"));
  // coords generated with two different seeds differ:
  assert(strcmp(cxsmi, cxsmi3));
  free(cxsmi3);
  free(cxsmi);

  free(mpkl);
  mpkl = NULL;
  printf("  done\n");
  printf("--------------------------\n");
}

void test_standardize() {
  printf("--------------------------\n");
  printf("  test_standardize\n");
  disable_logging();
  char *mpkl;
  size_t mpkl_size;
  mpkl = get_mol("[Pt]CCN(=O)=O", &mpkl_size, "{\"sanitize\":false}");
  char *smi = get_smiles(mpkl, mpkl_size, "");
  assert(!strcmp(smi, "O=N(=O)CC[Pt]"));
  free(smi);
  assert(cleanup(&mpkl, &mpkl_size, "") > 0);
  smi = get_smiles(mpkl, mpkl_size, "");
  assert(!strcmp(smi, "[CH2-]C[N+](=O)[O-].[Pt+]"));
  free(smi);

  assert(fragment_parent(&mpkl, &mpkl_size, "") > 0);
  smi = get_smiles(mpkl, mpkl_size, "");
  assert(!strcmp(smi, "[CH2-]C[N+](=O)[O-]"));
  free(smi);

  assert(charge_parent(&mpkl, &mpkl_size, "") > 0);
  smi = get_smiles(mpkl, mpkl_size, "");
  assert(!strcmp(smi, "CC[N+](=O)[O-]"));
  free(smi);
  free(mpkl);

  mpkl = get_mol("[Pt]CCN(=O)=O", &mpkl_size, "{\"sanitize\":false}");
  assert(charge_parent(&mpkl, &mpkl_size, "") > 0);
  smi = get_smiles(mpkl, mpkl_size, "");
  assert(!strcmp(smi, "CC[N+](=O)[O-]"));
  free(smi);
  free(mpkl);

  mpkl = get_mol("[Pt]CCN(=O)=O", &mpkl_size, "{\"sanitize\":false}");
  assert(charge_parent(&mpkl, &mpkl_size, "{\"skipStandardize\":true}") > 0);
  smi = get_smiles(mpkl, mpkl_size, "");
  assert(!strcmp(smi, "[CH2-]C[N+](=O)[O-].[Pt+]"));
  free(smi);
  free(mpkl);

  mpkl = get_mol("[CH2-]CN(=O)=O", &mpkl_size, "{\"sanitize\":false}");
  assert(neutralize(&mpkl, &mpkl_size, "") > 0);
  smi = get_smiles(mpkl, mpkl_size, "");
  assert(!strcmp(smi, "CCN(=O)=O"));
  free(smi);
  free(mpkl);

  mpkl = get_mol("[O-]c1cc(C(=O)O)ccc1", &mpkl_size, "");
  assert(reionize(&mpkl, &mpkl_size, "") > 0);
  smi = get_smiles(mpkl, mpkl_size, "");
  assert(!strcmp(smi, "O=C([O-])c1cccc(O)c1"));
  free(smi);
  free(mpkl);

  mpkl = get_mol("OC(O)C(=N)CO", &mpkl_size, "");
  assert(canonical_tautomer(&mpkl, &mpkl_size, "") > 0);
  smi = get_smiles(mpkl, mpkl_size, "");
  assert(!strcmp(smi, "NC(CO)C(=O)O"));
  free(smi);
  free(mpkl);

  enable_logging();
  printf("  done\n");
  printf("--------------------------\n");
}

int main() {
  enable_logging();
  char *vers = version();
  printf("hello %s\n", vers);
  free(vers);

  test_io();
  test_svg();
  test_flexicanvas();
  test_substruct();
  test_descriptors();
  test_fingerprints();
  test_modifications();
  test_coords();
  test_standardize();
  return 0;
}