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
#include <math.h>
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

  svg = get_svg(
      pkl, pkl_size,
      "{\"highlightAtomColors\":{"
      "\"0\": [1.0, 0.0, 0.0],"
      "\"1\": [1.0, 0.0, 0.0],"
      "\"2\": [1.0, 0.0, 0.0],"
      "\"3\": [0.0, 1.0, 0.0],"
      "\"4\": [1.0, 0.0, 0.0],"
      "\"5\": [1.0, 0.0, 0.0],"
      "\"6\": [1.0, 0.0, 0.0]"
      "}, \"highlightBondColors\":{"
      "\"2\": [0.0, 0.7, 0.9]"
      "}, \"highlightAtomRadii\":{"
      "\"0\": 0.1,"
      "\"1\": 0.1,"
      "\"2\": 0.1,"
      "\"3\": 0.8,"
      "\"4\": 0.1,"
      "\"5\": 0.1,"
      "\"6\": 0.1"
      "}, \"atoms\": [0, 1, 2, 3, 4, 5, 6],  \"bonds\": [2], \"width\":127}");
  assert(!strstr(svg, "fill:#FF7F7F"));
  assert(strstr(svg, "fill:#FF0000"));
  assert(strstr(svg, "fill:#00FF00"));
  assert(strstr(svg, "fill:#00B2E5"));
  assert(strstr(svg, "width='127px'"));
  assert(strstr(svg, "</svg>"));
  free(svg);

  free(pkl);
  pkl = NULL;
  printf("  done\n");
  printf("--------------------------\n");
}

unsigned int count_occurrences(const char *str, const char *substr) {
  unsigned int count = 0;
  const char *from = str;
  while (from = strstr(from, substr)) {
    ++count;
    ++from;
  }
  return count;
}

void test_rxn_svg() {
  char *pkl;
  size_t pkl_size;

  printf("--------------------------\n");
  printf("  test_rxn_svg\n");

  pkl = get_rxn("[CH3:1][OH:2]>>[CH2:1]=[OH0:2]", &pkl_size, "");
  assert(pkl);
  assert(pkl_size > 0);

  char *svg = get_rxn_svg(pkl, pkl_size, "{\"width\":350,\"height\":300}");
  assert(strstr(svg, "width='350px'"));
  assert(strstr(svg, "height='300px'"));
  assert(strstr(svg, "</svg>"));
  free(pkl);

  pkl = get_rxn(
      "$RXN\n\
\n\
      RDKit\n\
\n\
  1  1\n\
$MOL\n\
\n\
     RDKit          2D\n\
\n\
  2  1  0  0  0  0  0  0  0  0999 V2000\n\
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0\n\
    1.2990    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  2  0  0\n\
  1  2  6  0\n\
V    1 [C&H3:1]\n\
V    2 [O&H1:2]\n\
M  END\n\
$MOL\n\
\n\
     RDKit          2D\n\
\n\
  2  1  0  0  0  0  0  0  0  0999 V2000\n\
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0\n\
    1.2990    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  2  0  0\n\
  1  2  2  0\n\
V    1 [C&H2:1]\n\
V    2 [O&H0:2]\n\
M  END",
      &pkl_size, "");
  char *svg_from_block = get_rxn_svg(pkl, pkl_size, "");
  assert(strstr(svg, "</svg>"));
  assert(count_occurrences(svg, "<path") > 0);
  assert(count_occurrences(svg, "<path") ==
         count_occurrences(svg_from_block, "<path"));
  free(svg);
  free(svg_from_block);
  free(pkl);

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

  assert(!get_morgan_fp(NULL, 0, NULL));
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
  assert(!get_morgan_fp_as_bytes(NULL, 0, &nbytes, NULL));
  fp = get_morgan_fp_as_bytes(mpkl, mpkl_size, &nbytes,
                              "{\"radius\":2,\"nBits\":64}");
  assert(nbytes == 8);
  free(fp);

  assert(!get_rdkit_fp(NULL, 0, NULL));
  fp = get_rdkit_fp(mpkl, mpkl_size, "{\"nBits\":64}");
  assert(!strcmp(
      fp, "1111011000111100011011011111011001111111110010000111000011111111"));
  free(fp);
  fp = get_rdkit_fp(mpkl, mpkl_size, "{\"nBits\":64,\"maxPath\":4}");
  assert(!strcmp(
      fp, "1111011000110000010001010111000001101011100010000001000011100011"));
  free(fp);
  assert(!get_rdkit_fp_as_bytes(NULL, 0, &nbytes, NULL));
  fp = get_rdkit_fp_as_bytes(mpkl, mpkl_size, &nbytes, "{\"nBits\":64}");
  assert(nbytes == 8);
  free(fp);

  assert(!get_pattern_fp(NULL, 0, NULL));
  fp = get_pattern_fp(mpkl, mpkl_size, "{\"nBits\":64}");
  assert(!strcmp(
      fp, "1011111111111111110011011111011001011110111101111110101111010011"));
  free(fp);
  assert(!get_pattern_fp_as_bytes(NULL, 0, &nbytes, NULL));
  fp = get_pattern_fp_as_bytes(mpkl, mpkl_size, &nbytes, "{\"nBits\":64}");
  assert(nbytes == 8);
  free(fp);

  fp = get_pattern_fp(mpkl, mpkl_size,
                      "{\"nBits\":64,\"tautomericFingerprint\":true}");
  assert(!strcmp(
      fp, "1011111111111111110011111111011101011111111101111111111111011011"));
  free(fp);

  assert(!get_topological_torsion_fp(NULL, 0, NULL));
  fp = get_topological_torsion_fp(mpkl, mpkl_size, "{\"nBits\":64}");
  assert(!strcmp(
      fp, "1100000000000000000000000000000000000000000011000000000011001100"));
  free(fp);
  assert(!get_topological_torsion_fp_as_bytes(NULL, 0, &nbytes, NULL));
  fp = get_topological_torsion_fp_as_bytes(mpkl, mpkl_size, &nbytes,
                                           "{\"nBits\":64}");
  assert(nbytes == 8);
  free(fp);

  assert(!get_atom_pair_fp(NULL, 0, NULL));
  fp = get_atom_pair_fp(mpkl, mpkl_size, "{\"nBits\":64}");
  assert(!strcmp(
      fp, "0000000000000000000000001110111011101100000000000000000011001100"));
  free(fp);
  fp = get_atom_pair_fp(mpkl, mpkl_size, "{\"nBits\":64,\"minLength\":2}");
  assert(!strcmp(
      fp, "0000000000000000000000001000111011100000000000000000000011001100"));
  free(fp);
  fp = get_atom_pair_fp(mpkl, mpkl_size, "{\"nBits\":64,\"maxLength\":2}");
  assert(!strcmp(
      fp, "0000000000000000000000001110111011001100000000000000000011000000"));
  free(fp);
  assert(!get_atom_pair_fp_as_bytes(NULL, 0, &nbytes, NULL));
  fp = get_atom_pair_fp_as_bytes(mpkl, mpkl_size, &nbytes, "{\"nBits\":64}");
  assert(nbytes == 8);
  free(fp);

#ifdef RDK_BUILD_AVALON_SUPPORT
  assert(!get_avalon_fp(NULL, 0, NULL));
  fp = get_avalon_fp(mpkl, mpkl_size, "{\"nBits\":64}");
  assert(!strcmp(
      fp, "0001000001100011110001011100011100000110100010001110110001110011"));
  free(fp);
  assert(!get_avalon_fp_as_bytes(NULL, 0, nbytes, NULL));
  fp = get_avalon_fp_as_bytes(mpkl, mpkl_size, &nbytes, "{\"nBits\":64}");
  assert(nbytes == 8);
  free(fp);
#endif

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

float **_get_coord_array(char *mpkl, size_t mpkl_size) {
  char *molblock = get_molblock(mpkl, mpkl_size, NULL);
  assert(molblock);
  size_t molblock_len = strlen(molblock);
  char *line_start = molblock;
  char *line_end = strpbrk(line_start, "\n");
  int line_num = 0;
  float **res = NULL;
  unsigned int i;
  unsigned int j;
  unsigned int s;
  unsigned int e;
  unsigned int atom_num = 0;
  int n_atoms = 0;
  while (line_end) {
    *line_end = '\0';
    if (line_num == 3) {
      assert(strlen(line_start) > 3);
      line_start[3] = '\0';
      sscanf(line_start, "%d", &n_atoms);
      if (n_atoms) {
        res = (float **)malloc((n_atoms + 1) * sizeof(float *));
        assert(res);
        res[n_atoms] = NULL;
        for (i = 0; i < n_atoms; ++i) {
          res[i] = (float *)malloc(3 * sizeof(float));
          assert(res[i]);
        }
      }
    } else if (line_num > 3 && line_num < 4 + n_atoms && res) {
      for (i = 0; i < 3; ++i) {
        j = 2 - i;
        s = j * 10;
        e = (j + 1) * 10;
        line_start[e] = '\0';
        sscanf(&line_start[s], "%f", &res[atom_num][j]);
      }
      ++atom_num;
    }
    line_start = line_end + 1;
    if (line_start >= molblock + molblock_len) {
      break;
    }
    line_end = strpbrk(line_start, "\n");
    if (!line_end) {
      line_end = molblock + molblock_len;
    }
    ++line_num;
  }
  free(molblock);
  return res;
}

void _free_coord_array(float **coord_array) {
  size_t i = 0;
  if (coord_array) {
    while (coord_array[i]) {
      free(coord_array[i++]);
    }
    free(coord_array);
  }
}

float _sq_dist(float *xyz1, float *xyz2) {
  float sqd = 0.;
  unsigned int i;
  for (i = 0; i < 3; ++i) {
    sqd += (xyz1[i] - xyz2[i]) * (xyz1[i] - xyz2[i]);
  }
  return sqd;
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
  prefer_coordgen(0);
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
                                "{\"acceptFailure\":false}", NULL));
  char *match_json;
  assert(set_2d_coords_aligned(
      &mpkl, &mpkl_size, tpkl, tpkl_size,
      "{\"allowRGroups\":true,\"acceptFailure\":false}", &match_json));
  assert(!strcmp(match_json, "{\"atoms\":[4,3,0,1,2],\"bonds\":[3,5,0,1,2]}"));
  free(match_json);
  assert(set_2d_coords_aligned(
      &mpkl, &mpkl_size, tpkl, tpkl_size,
      "{\"allowRGroups\":true,\"acceptFailure\":false,\"alignOnly\":true}", &match_json));
  assert(!strcmp(match_json, "{\"atoms\":[4,3,0,1,2],\"bonds\":[3,5,0,1,2]}"));
  free(match_json);
  free(tpkl);

  // Github #4121: set_2d_coords_aligned crashes if template mol has no
  // conformer
  tpkl = get_mol("C1CNC1", &tpkl_size, NULL);
  assert(!set_2d_coords_aligned(&mpkl, &mpkl_size, tpkl, tpkl_size, "", NULL));
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

  // test set_2d_coords_aligned
  tpkl = get_mol(
      "\n\
     RDKit          2D\n\
\n\
  6  6  0  0  0  0  0  0  0  0999 V2000\n\
  -13.7477    6.3036    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n\
  -13.7477    4.7567    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
  -12.6540    3.6628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
  -13.7477    2.5691    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
  -14.8414    3.6628    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n\
  -11.1071    3.6628    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n\
  1  2  1  0\n\
  2  3  1  0\n\
  3  4  1  0\n\
  4  5  1  0\n\
  2  5  1  0\n\
  3  6  1  0\n\
M  RGP  2   1   1   6   2\n\
M  END\n",
      &tpkl_size, "");
  mpkl = get_mol(
      "\n\
     RDKit          2D\n\
\n\
 18 22  0  0  0  0  0  0  0  0999 V2000\n\
    4.3922   -1.5699    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    2.9211   -2.0479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    2.5995   -0.5349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    3.3731    0.8046    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    4.8441    1.2825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    4.0704   -0.0568    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    1.8666    0.7748    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.7736   -0.3197    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.7749   -1.8666    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.7718   -1.8679    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.7731   -0.3208    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.8679    0.7718    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -4.0718   -0.0598    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -4.3933   -1.5729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -2.9222   -2.0509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -2.6008   -0.5379    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -3.3744    0.8016    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -4.8454    1.2795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
  9 10  1  0\n\
 11 10  1  0\n\
 11  8  1  0\n\
  8  9  1  0\n\
  4  5  1  0\n\
  6  5  1  0\n\
  7  6  1  0\n\
  3  4  1  0\n\
  3  7  1  0\n\
  1  6  1  0\n\
  2  3  1  0\n\
  2  1  1  0\n\
 17 18  1  0\n\
 13 18  1  0\n\
 12 13  1  0\n\
 16 17  1  0\n\
 16 12  1  0\n\
 14 13  1  0\n\
 15 16  1  0\n\
 15 14  1  0\n\
 12 11  1  0\n\
  8  7  1  0\n\
M  END\n",
      &mpkl_size, "");
  float **mol_coords = _get_coord_array(mpkl, mpkl_size);
  assert(mol_coords);
  float bond_length11_12 = sqrt(_sq_dist(mol_coords[11], mol_coords[12]));
  float bond_length5_6 = sqrt(_sq_dist(mol_coords[5], mol_coords[6]));
  assert(fabs(bond_length11_12 - bond_length5_6) < 1.e-4);
  assert(bond_length11_12 > 2.3);
  float bond_length_ali11_12;
  float bond_length_ali5_6;
  float **tpl_coords = _get_coord_array(tpkl, tpkl_size);
  assert(tpl_coords);
  float **mol_ali_coords = NULL;
  unsigned int mol_indices[] = {11, 10, 7, 8, 9, 6};
  unsigned int mol_idx;
  unsigned int tpl_idx;
  unsigned int i;
  char details[200];
  char *align_only_choices[] = { "false", "true" };
  char *mpkl2_molblock_before = NULL;
  char *mpkl2_molblock_after = NULL;
  char *mpkl_smi = NULL;
  size_t mpkl_smi_size;
  char *mpkl2 = NULL;
  size_t mpkl2_size;

  for (i = 0; i < 2; ++i) {
    // this has no initial coordinates and matches the template
    mpkl_smi = get_mol("C1CC2CCC1N2C1CNC1N1C2CCC1CC2", &mpkl_smi_size, "");
    assert(!has_coords(mpkl_smi, mpkl_smi_size));
    memset(details, 0, 200);
    sprintf(details, "{\"acceptFailure\":false,\"allowRGroups\":true,\"alignOnly\":%s}", align_only_choices[i]);
    assert(set_2d_coords_aligned(
        &mpkl_smi, &mpkl_smi_size, tpkl, tpkl_size, details, &match_json));
    assert(!strcmp(match_json, "{\"atoms\":[11,10,7,8,9,6],\"bonds\":[10,18,7,8,9,6]}"));
    free(match_json);
    // coordinates should be present as alignment has taken place anyway
    assert(has_coords(mpkl_smi, mpkl_smi_size));
    free(mpkl_smi);

    mpkl2_size = mpkl_size;
    mpkl2 = malloc(mpkl2_size);
    assert(mpkl2);
    memcpy(mpkl2, mpkl, mpkl2_size);
    memset(details, 0, 200);
    sprintf(details, "{\"allowRGroups\":true,\"alignOnly\":%s}", align_only_choices[i]);
    assert(set_2d_coords_aligned(
        &mpkl2, &mpkl2_size, tpkl, tpkl_size, details, &match_json));
    assert(!strcmp(match_json,
                  "{\"atoms\":[11,10,7,8,9,6],\"bonds\":[20,2,3,0,1,21]}"));
    free(match_json);
    mol_ali_coords = _get_coord_array(mpkl2, mpkl2_size);
    for (tpl_idx = 0; tpl_idx < sizeof(mol_indices) / sizeof(mol_indices[0]);
        ++tpl_idx) {
      mol_idx = mol_indices[tpl_idx];
      assert(_sq_dist(tpl_coords[tpl_idx], mol_ali_coords[mol_idx]) < 1.e-4);
    }
    bond_length_ali11_12 = sqrt(_sq_dist(mol_ali_coords[11], mol_ali_coords[12]));
    bond_length_ali5_6 = sqrt(_sq_dist(mol_ali_coords[5], mol_ali_coords[6]));
    assert(fabs(bond_length_ali11_12 - bond_length_ali5_6) < 1.e-4);
    if (i) {
      assert(bond_length_ali11_12 > 2.3);
    } else {
      assert(bond_length_ali11_12 < 1.6);
    }
    _free_coord_array(mol_ali_coords);
    free(mpkl2);

    memset(details, 0, 200);
    sprintf(details, "{\"useCoordGen\":true,\"acceptFailure\":false,\"allowRGroups\":true,\"alignOnly\":%s}", align_only_choices[i]);
    // this has no initial coordinates and does not match the template
    mpkl_smi = get_mol("C1CC2CCC1N2C1CCNC1N1C2CCC1CC2", &mpkl_smi_size, "");
    assert(!has_coords(mpkl_smi, mpkl_smi_size));
    // This should fail
    assert(!set_2d_coords_aligned(
        &mpkl_smi, &mpkl_smi_size, tpkl, tpkl_size, details, &match_json));
    assert(!match_json);
    // coordinates should be absent since alignment has not taken place
    assert(!has_coords(mpkl_smi, mpkl_smi_size));

    // this has initial coordinates and does not match the template
    mpkl2_size = mpkl_smi_size;
    mpkl2 = malloc(mpkl2_size);
    assert(mpkl2);
    memcpy(mpkl2, mpkl_smi, mpkl2_size);
    assert(!has_coords(mpkl2, mpkl2_size));
    set_2d_coords(&mpkl2, &mpkl2_size);
    assert(has_coords(mpkl2, mpkl2_size));
    mpkl2_molblock_before = get_molblock(mpkl2, mpkl2_size, "");
    // This should fail
    assert(!set_2d_coords_aligned(
        &mpkl2, &mpkl2_size, tpkl, tpkl_size, details, &match_json));
    assert(!match_json);
    // coordinates should be unchanged since alignment has not taken place
    assert(has_coords(mpkl2, mpkl2_size));
    mpkl2_molblock_after = get_molblock(mpkl2, mpkl2_size, "");
    assert(!strcmp(mpkl2_molblock_before, mpkl2_molblock_after));
    free(mpkl2_molblock_after);

    memset(details, 0, 200);
    sprintf(details, "{\"useCoordGen\":true,\"acceptFailure\":true,\"allowRGroups\":true,\"alignOnly\":%s}", align_only_choices[i]);
    // This should do a simple coordinate generation, no alignment
    assert(set_2d_coords_aligned(
        &mpkl2, &mpkl2_size, tpkl, tpkl_size, details, &match_json));
    assert(!strcmp(match_json, "{}"));
    free(match_json);
    // coordinates should have changed since coordinate generation has taken place anyway using CoordGen
    assert(has_coords(mpkl2, mpkl2_size));
    mpkl2_molblock_after = get_molblock(mpkl2, mpkl2_size, "");
    assert(strcmp(mpkl2_molblock_before, mpkl2_molblock_after));
    free(mpkl2_molblock_before);
    free(mpkl2_molblock_after);
    free(mpkl2);

    // this has no initial coordinates and does not match the template
    assert(set_2d_coords_aligned(
        &mpkl_smi, &mpkl_smi_size, tpkl, tpkl_size, details, &match_json));
    assert(!strcmp(match_json, "{}"));
    free(match_json);
    // coordinates should be present since coordinate generation has taken place anyway using CoordGen
    assert(has_coords(mpkl_smi, mpkl_smi_size));
    free(mpkl_smi);
  }

  _free_coord_array(mol_coords);
  _free_coord_array(tpl_coords);
  free(mpkl);
  free(tpkl);

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
  test_rxn_svg();
  test_flexicanvas();
  test_substruct();
  test_descriptors();
  test_fingerprints();
  test_modifications();
  test_coords();
  test_standardize();
  return 0;
}
