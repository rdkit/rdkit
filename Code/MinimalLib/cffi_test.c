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
#ifdef WIN32
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#define _DEFINED_USE_MATH_DEFINES
#endif
#endif
#include <math.h>
#ifdef _DEFINED_USE_MATH_DEFINES
#undef _DEFINED_USE_MATH_DEFINES
#undef _USE_MATH_DEFINES
#endif
#include "cffiwrapper.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>

static const char molblock_native_wedging[] =
    "\n\
  MJ201100                      \n\
\n\
 18 21  0  0  1  0  0  0  0  0999 V2000\n\
   -0.8540   -1.4441    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.3019   -0.8310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.5185   -0.9172    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.8540   -0.1635    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.6825    0.6434    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.1379    0.7296    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.5504    1.4441    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.4734   -0.0239    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.2409    0.3885    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.6609   -1.2726    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -2.2130   -1.8857    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.9580   -2.6703    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.1511   -2.8419    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.5990   -2.2287    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -3.0201   -1.7143    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -3.5720   -2.3275    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -3.3171   -3.1121    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -2.5100   -3.2835    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n\
  2  3  1  0  0  0  0\n\
  4  3  1  0  0  0  0\n\
  4  5  1  0  0  0  0\n\
  6  5  1  0  0  0  0\n\
  6  7  1  1  0  0  0\n\
  6  8  1  0  0  0  0\n\
  8  9  1  1  0  0  0\n\
  8  2  1  0  0  0  0\n\
  4  9  1  1  0  0  0\n\
  2  1  1  1  0  0  0\n\
 10 11  1  0  0  0  0\n\
 11 12  2  0  0  0  0\n\
 12 13  1  0  0  0  0\n\
 13 14  2  0  0  0  0\n\
  1 10  2  0  0  0  0\n\
  1 14  1  0  0  0  0\n\
 15 16  2  0  0  0  0\n\
 16 17  1  0  0  0  0\n\
 11 15  1  0  0  0  0\n\
 17 18  2  0  0  0  0\n\
 12 18  1  0  0  0  0\n\
M  END\n";

static const char quinoline_scaffold[] =
    "\n\
  MJ201100                      \n\
\n\
 10 11  0  0  1  0  0  0  0  0999 V2000\n\
   -8.1001    2.8219    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -8.8145    2.4094    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -8.8145    1.5843    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -8.1001    1.1718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -7.3856    1.5843    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -7.3856    2.4094    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -6.6711    1.1718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -5.9566    1.5842    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -5.9566    2.4092    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -6.6711    2.8218    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
  2  3  1  0  0  0  0\n\
  3  4  2  0  0  0  0\n\
  4  5  1  0  0  0  0\n\
  5  6  2  0  0  0  0\n\
  7  8  2  0  0  0  0\n\
  8  9  1  0  0  0  0\n\
  9 10  2  0  0  0  0\n\
  5  7  1  0  0  0  0\n\
 10  6  1  0  0  0  0\n\
  1  2  2  0  0  0  0\n\
  6  1  1  0  0  0  0\n\
M  END\n";

void find_wedged_bonds(char *molblock, int *have1, int *have6) {
  molblock = strdup(molblock);
  assert(molblock);
  size_t molblock_len = strlen(molblock);
  char *line_start = molblock;
  char *line_end = strpbrk(line_start, "\n");
  int line_num = 0;
  unsigned int i;
  unsigned int j;
  unsigned int s;
  unsigned int e;
  int n_atoms = -1;
  int n_bonds = -1;
  int b[4];
  *have1 = 0;
  *have6 = 0;
  while (line_end) {
    *line_end = '\0';
    if (line_num == 3) {
      assert(strlen(line_start) > 6);
      line_start[6] = '\0';
      sscanf(&line_start[3], "%d", &n_bonds);
      line_start[3] = '\0';
      sscanf(line_start, "%d", &n_atoms);
      assert(n_atoms >= 0 && n_bonds >= 0);
    } else if (line_num > 3 + n_atoms && line_num < 4 + n_atoms + n_bonds) {
      for (i = 0; i < 4; ++i) {
        j = 3 - i;
        s = j * 3;
        e = (j + 1) * 3;
        line_start[e] = '\0';
        sscanf(&line_start[s], "%d", &b[j]);
      }
      assert(b[0] >= 1 && b[0] <= n_atoms);
      assert(b[1] >= 1 && b[1] <= n_atoms);
      assert(b[2] == 1 || b[2] == 2);
      assert(b[3] == 0 || b[3] == 1 || b[3] == 6);
      if (b[3] == 1) {
        *have1 = 1;
      }
      if (b[3] == 6) {
        *have6 = 1;
      }
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
}

int extract_bond_coords(char *svg, char *bond, double *coord1, double *coord2) {
  svg = strdup(svg);
  assert(svg);
  char *line = strtok(svg, "\n");
  char *str = NULL;
  double dummy[2];
  coord1 = coord1 ? coord1 : dummy;
  coord2 = coord2 ? coord2 : dummy;
  while (line) {
    str = strstr(line, bond);
    if (str) {
      str = strstr(str, "M ");
    }
    if (str) {
      assert(sscanf(str, "M %lf,%lf L %lf,%lf", &coord1[0], &coord1[1],
                    &coord2[0], &coord2[1]) == 4);
      break;
    }
    line = strtok(NULL, "\n");
  }
  free(svg);
  return (str ? 1 : 0);
}

double angle_deg_between_vectors(double *v1, double *v2) {
  return 180 / M_PI *
         acos((v1[0] * v2[0] + v1[1] * v2[1]) /
              sqrt((v1[0] * v1[0] + v1[1] * v1[1]) *
                   (v2[0] * v2[0] + v2[1] * v2[1])));
}

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
  assert(strstr(json, "rdkitjson"));

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
  char *molblock = get_molblock(NULL, pkl_size, NULL);
  assert(!molblock);
  molblock = get_molblock(pkl, 0, NULL);
  assert(!molblock);
  molblock = get_molblock(pkl, pkl_size, NULL);
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

  pkl2 = get_mol(molblock_native_wedging, &pkl2_size, "");
  assert(pkl2);
  assert(pkl2_size > 0);
  molblock = get_molblock(pkl2, pkl2_size, NULL);
  assert(strstr(molblock, "4  3  1  6"));
  assert(!strstr(molblock, "H  "));
  free(molblock);
  molblock = get_molblock(pkl2, pkl2_size, "{\"useMolBlockWedging\":true}");
  assert(!strstr(molblock, "4  3  1  6"));
  assert(strstr(molblock, "6  7  1  1"));
  assert(!strstr(molblock, "H  "));
  free(molblock);
  molblock = get_molblock(pkl2, pkl2_size, "{\"addChiralHs\":true}");
  assert(strstr(molblock, "H  "));
  free(molblock);
  // Here we want to test that the original molblock wedging is preserved and
  // inverted as the coordinates are rigid-body rotated
  size_t scaffold_pkl_size;
  int have1;
  int have6;
  char *scaffold = get_mol(quinoline_scaffold, &scaffold_pkl_size, NULL);
  assert(set_2d_coords_aligned(&pkl2, &pkl2_size, scaffold, scaffold_pkl_size,
                               "{\"acceptFailure\":false,\"alignOnly\":true}",
                               NULL));
  molblock = get_molblock(pkl2, pkl2_size, "{\"useMolBlockWedging\":true}");
  find_wedged_bonds(molblock, &have1, &have6);
  assert(!have1 && have6);
  assert(!strstr(molblock, "4  3  1  6"));
  assert(strstr(molblock, "6  7  1  6"));
  assert(!strstr(molblock, "H  "));
  free(molblock);
  free(pkl2);
  // Here we want to test that the original molblock wedging gets cleared
  // and hence wedging is recomputed as the coordinates are re-generated
  pkl2 = get_mol(molblock_native_wedging, &pkl2_size, "");
  assert(pkl2);
  assert(pkl2_size > 0);
  assert(set_2d_coords_aligned(&pkl2, &pkl2_size, scaffold, scaffold_pkl_size,
                               "{\"acceptFailure\":false}", NULL));
  molblock = get_molblock(pkl2, pkl2_size, "{\"useMolBlockWedging\":true}");
  find_wedged_bonds(molblock, &have1, &have6);
  assert(have1 && have6);
  free(molblock);
  free(pkl2);
  free(scaffold);

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

  pkl = get_mol(molblock_native_wedging, &pkl_size, "");
  assert(pkl);
  assert(pkl_size > 0);
  char *svg1 = get_svg(pkl, pkl_size, "{\"width\":350,\"height\":300}");
  assert(strstr(svg1, "width='350px'"));
  assert(strstr(svg1, "height='300px'"));
  assert(strstr(svg1, "</svg>"));
  assert(strstr(svg1, "atom-17"));
  assert(strstr(svg1, "atom-18"));
  assert(strstr(svg1, "atom-19"));
  char *svg2 = get_svg(
      pkl, pkl_size,
      "{\"width\":350,\"height\":300,\"useMolBlockWedging\":true,\"wedgeBonds\":false,\"addChiralHs\":false}");
  assert(strstr(svg2, "width='350px'"));
  assert(strstr(svg2, "height='300px'"));
  assert(strstr(svg2, "</svg>"));
  assert(strstr(svg2, "atom-17"));
  assert(!strstr(svg2, "atom-18"));
  assert(!strstr(svg2, "atom-19"));
  free(svg1);
  free(svg2);

  free(pkl);
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
  assert(strstr(svg, "width='87px'"));
  assert(strstr(svg, "height='19px'"));
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

  assert(!get_maccs_fp(NULL, 0));
  fp = get_maccs_fp(mpkl, mpkl_size);
  assert(!strcmp(
      fp,
      "00000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000001100000000000000100000001000001000000000101000100000000100001000111110"));
  free(fp);
  assert(!get_maccs_fp_as_bytes(NULL, 0, &nbytes));
  fp = get_maccs_fp_as_bytes(mpkl, mpkl_size, &nbytes);
  assert(nbytes == 21);
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
  char *mpkl2;
  size_t mpkl2_size;
  char *mpkl3;
  size_t mpkl3_size;
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
  char *tpkl2;
  size_t tpkl2_size;
  char *tpkl3;
  size_t tpkl3_size;
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
      "{\"allowRGroups\":true,\"acceptFailure\":false,\"alignOnly\":true}",
      &match_json));
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
  assert(has_coords(mpkl, mpkl_size) == 3);
  assert(set_3d_coords(&mpkl, &mpkl_size,
                       "{\"randomSeed\":123,\"coordMap\":{\"3\":[0,0,0],\"4\":["
                       "0,0,1.5],\"5\":[0,1.5,1.5]}}") > 0);
  assert(has_coords(mpkl, mpkl_size) == 3);
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
  tpkl3 = get_mol(
      "\n\
  MJ201100                      \n\
\n\
 12 13  0  0  0  0  0  0  0  0999 V2000\n\
   -0.5398    0.0400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.3648    0.0400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.7773   -0.6745    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.3649   -1.3889    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.5399   -1.3889    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.1273   -0.6744    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.6976   -0.6744    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.9167    0.6531    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -2.6704    0.3176    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -2.5842   -0.5028    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -3.3849    0.7302    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.7451    1.4600    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n\
  1  2  2  0  0  0  0\n\
  2  3  1  0  0  0  0\n\
  3  4  2  0  0  0  0\n\
  4  5  1  0  0  0  0\n\
  5  6  2  0  0  0  0\n\
  6  1  1  0  0  0  0\n\
  6  7  1  0  0  0  0\n\
  8  9  2  0  0  0  0\n\
  2  8  1  0  0  0  0\n\
  9 10  1  0  0  0  0\n\
  3 10  1  0  0  0  0\n\
  9 11  1  0  0  0  0\n\
  8 12  1  0  0  0  0\n\
M  ALS   7 10 F H   C   N   O   F   P   S   Cl  Br  I   \n\
M  ALS  11 10 F H   C   N   O   F   P   S   Cl  Br  I   \n\
M  ALS  12 10 F H   C   N   O   F   P   S   Cl  Br  I   \n\
M  END\n",
      &tpkl3_size, "");
  mpkl3 = get_mol(
      "\n\
  MJ201100                      \n\
\n\
 13 14  0  0  0  0  0  0  0  0999 V2000\n\
   -0.6112    0.3665    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.3648    0.0310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.4510   -0.7895    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.7836   -1.2744    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.0299   -0.9389    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.0562   -0.1183    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.8099    0.2172    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -2.1184    0.3666    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -2.6705   -0.2464    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -2.2580   -0.9608    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.6374   -1.4238    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.8961    1.0377    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.5512   -2.2443    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
  1  2  2  0  0  0  0\n\
  2  3  1  0  0  0  0\n\
  3  4  2  0  0  0  0\n\
  4  5  1  0  0  0  0\n\
  5  6  2  0  0  0  0\n\
  6  1  1  0  0  0  0\n\
  8  9  2  0  0  0  0\n\
  2  8  1  0  0  0  0\n\
  9 10  1  0  0  0  0\n\
  3 10  1  0  0  0  0\n\
  6  7  1  0  0  0  0\n\
  5 11  1  0  0  0  0\n\
  7 12  1  0  0  0  0\n\
 11 13  1  0  0  0  0\n\
M  END\n",
      &mpkl3_size, "");
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
  char *align_only_choices[] = {"false", "true"};
  char *mpkl2_molblock_before = NULL;
  char *mpkl2_molblock_after = NULL;
  char *mpkl_smi = NULL;
  size_t mpkl_smi_size;

  for (i = 0; i < 2; ++i) {
    // this has no initial coordinates and matches the template
    mpkl_smi = get_mol("C1CC2CCC1N2C1CNC1N1C2CCC1CC2", &mpkl_smi_size, "");
    assert(!has_coords(mpkl_smi, mpkl_smi_size));
    memset(details, 0, 200);
    sprintf(details,
            "{\"acceptFailure\":false,\"allowRGroups\":true,\"alignOnly\":%s}",
            align_only_choices[i]);
    assert(set_2d_coords_aligned(&mpkl_smi, &mpkl_smi_size, tpkl, tpkl_size,
                                 details, &match_json));
    assert(!strcmp(match_json,
                   "{\"atoms\":[11,10,7,8,9,6],\"bonds\":[10,18,7,8,9,6]}"));
    free(match_json);
    // coordinates should be present as alignment has taken place anyway
    assert(has_coords(mpkl_smi, mpkl_smi_size) == 2);
    free(mpkl_smi);

    mpkl2_size = mpkl_size;
    mpkl2 = malloc(mpkl2_size);
    assert(mpkl2);
    memcpy(mpkl2, mpkl, mpkl2_size);
    memset(details, 0, 200);
    sprintf(details, "{\"allowRGroups\":true,\"alignOnly\":%s}",
            align_only_choices[i]);
    assert(set_2d_coords_aligned(&mpkl2, &mpkl2_size, tpkl, tpkl_size, details,
                                 &match_json));
    assert(!strcmp(match_json,
                   "{\"atoms\":[11,10,7,8,9,6],\"bonds\":[20,2,3,0,1,21]}"));
    free(match_json);
    mol_ali_coords = _get_coord_array(mpkl2, mpkl2_size);
    for (tpl_idx = 0; tpl_idx < sizeof(mol_indices) / sizeof(mol_indices[0]);
         ++tpl_idx) {
      mol_idx = mol_indices[tpl_idx];
      assert(_sq_dist(tpl_coords[tpl_idx], mol_ali_coords[mol_idx]) < 1.e-4);
    }
    bond_length_ali11_12 =
        sqrt(_sq_dist(mol_ali_coords[11], mol_ali_coords[12]));
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
    sprintf(
        details,
        "{\"useCoordGen\":true,\"acceptFailure\":false,\"allowRGroups\":true,\"alignOnly\":%s}",
        align_only_choices[i]);
    // this has no initial coordinates and does not match the template
    mpkl_smi = get_mol("C1CC2CCC1N2C1CCNC1N1C2CCC1CC2", &mpkl_smi_size, "");
    assert(!has_coords(mpkl_smi, mpkl_smi_size));
    // This should fail
    assert(!set_2d_coords_aligned(&mpkl_smi, &mpkl_smi_size, tpkl, tpkl_size,
                                  details, &match_json));
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
    assert(has_coords(mpkl2, mpkl2_size) == 2);
    mpkl2_molblock_before = get_molblock(mpkl2, mpkl2_size, "");
    // This should fail
    assert(!set_2d_coords_aligned(&mpkl2, &mpkl2_size, tpkl, tpkl_size, details,
                                  &match_json));
    assert(!match_json);
    // coordinates should be unchanged since alignment has not taken place
    assert(has_coords(mpkl2, mpkl2_size) == 2);
    mpkl2_molblock_after = get_molblock(mpkl2, mpkl2_size, "");
    assert(!strcmp(mpkl2_molblock_before, mpkl2_molblock_after));
    free(mpkl2_molblock_after);

    memset(details, 0, 200);
    sprintf(
        details,
        "{\"useCoordGen\":true,\"acceptFailure\":true,\"allowRGroups\":true,\"alignOnly\":%s}",
        align_only_choices[i]);
    // This should do a simple coordinate generation, no alignment
    assert(set_2d_coords_aligned(&mpkl2, &mpkl2_size, tpkl, tpkl_size, details,
                                 &match_json));
    assert(!strcmp(match_json, "{}"));
    free(match_json);
    // coordinates should have changed since coordinate generation has taken
    // place anyway using CoordGen
    assert(has_coords(mpkl2, mpkl2_size) == 2);
    mpkl2_molblock_after = get_molblock(mpkl2, mpkl2_size, "");
    assert(strcmp(mpkl2_molblock_before, mpkl2_molblock_after));
    free(mpkl2_molblock_before);
    free(mpkl2_molblock_after);
    free(mpkl2);

    // this has no initial coordinates and does not match the template
    assert(set_2d_coords_aligned(&mpkl_smi, &mpkl_smi_size, tpkl, tpkl_size,
                                 details, &match_json));
    assert(!strcmp(match_json, "{}"));
    free(match_json);
    // coordinates should be present since coordinate generation has taken place
    // anyway using CoordGen
    assert(has_coords(mpkl_smi, mpkl_smi_size) == 2);
    free(mpkl_smi);

    memset(details, 0, 200);
    sprintf(details,
            "{\"acceptFailure\":false,\"allowRGroups\":true,\"alignOnly\":%s}",
            align_only_choices[i]);
    assert(set_2d_coords_aligned(&mpkl3, &mpkl3_size, tpkl3, tpkl3_size,
                                 details, &match_json));
    assert(!strcmp(
        match_json,
        "{\"atoms\":[0,1,2,3,4,5,6,7,8,9],\"bonds\":[0,1,2,3,4,5,10,6,7,8,9]}"));
    free(match_json);
  }

  _free_coord_array(mol_coords);
  _free_coord_array(tpl_coords);
  free(mpkl);
  free(tpkl);
  free(mpkl3);
  free(tpkl3);

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

void test_get_mol_frags() {
  printf("--------------------------\n");
  printf("  test_get_mol_frags\n");
  char *mpkl;
  char *smi;
  size_t mpkl_size;
  size_t *frags_pkl_sz_array = NULL;
  size_t num_frags = 0;
  char **frags_mpkl_array = NULL;
  char *mappings_json = NULL;
  size_t i;

  mpkl = get_mol("n1ccccc1.CC(C)C.OCCCN", &mpkl_size, "");
  const char *expected_frag_smiles[] = {"c1ccncc1", "CC(C)C", "NCCCO"};
  const char *expected_frag_smiles_non_sanitized[] = {"CN(C)(C)C", "c1ccc1"};
  const char *expected_mappings =
      "{\"frags\":[0,0,0,0,0,0,1,1,1,1,2,2,2,2,2],\"fragsMolAtomMapping\":[[0,1,2,3,4,5],[6,7,8,9],[10,11,12,13,14]]}";

  frags_mpkl_array =
      get_mol_frags(mpkl, mpkl_size, &frags_pkl_sz_array, &num_frags, "", NULL);
  assert(frags_mpkl_array);
  assert(frags_pkl_sz_array);
  assert(num_frags == 3);
  for (i = 0; i < num_frags; ++i) {
    assert(frags_pkl_sz_array[i]);
    smi = get_smiles(frags_mpkl_array[i], frags_pkl_sz_array[i], NULL);
    assert(smi);
    assert(!strcmp(smi, expected_frag_smiles[i]));
    free(smi);
    free(frags_mpkl_array[i]);
    frags_mpkl_array[i] = NULL;
  }
  free(frags_mpkl_array);
  frags_mpkl_array = NULL;
  free(frags_pkl_sz_array);
  frags_pkl_sz_array = NULL;

  frags_mpkl_array = get_mol_frags(mpkl, mpkl_size, &frags_pkl_sz_array,
                                   &num_frags, "", &mappings_json);
  assert(frags_mpkl_array);
  assert(frags_pkl_sz_array);
  assert(mappings_json);
  assert(num_frags == 3);
  for (i = 0; i < num_frags; ++i) {
    assert(frags_pkl_sz_array[i]);
    smi = get_smiles(frags_mpkl_array[i], frags_pkl_sz_array[i], NULL);
    assert(smi);
    assert(!strcmp(smi, expected_frag_smiles[i]));
    free(smi);
    free(frags_mpkl_array[i]);
    frags_mpkl_array[i] = NULL;
  }
  free(frags_mpkl_array);
  frags_mpkl_array = NULL;
  free(frags_pkl_sz_array);
  frags_pkl_sz_array = NULL;
  assert(!strcmp(mappings_json, expected_mappings));
  free(mappings_json);
  mappings_json = NULL;
  free(mpkl);
  mpkl = NULL;

  mpkl = get_mol("N(C)(C)(C)C.c1ccc1", &mpkl_size, "{\"sanitize\":false}");
  frags_mpkl_array =
      get_mol_frags(mpkl, mpkl_size, &frags_pkl_sz_array, &num_frags, "", NULL);
  assert(!frags_mpkl_array);
  assert(!frags_pkl_sz_array);
  assert(num_frags == 0);
  frags_mpkl_array =
      get_mol_frags(mpkl, mpkl_size, &frags_pkl_sz_array, &num_frags,
                    "{\"sanitizeFrags\":false}", NULL);
  assert(frags_mpkl_array);
  assert(frags_pkl_sz_array);
  assert(num_frags == 2);
  for (i = 0; i < num_frags; ++i) {
    assert(frags_pkl_sz_array[i]);
    smi = get_smiles(frags_mpkl_array[i], frags_pkl_sz_array[i], NULL);
    assert(smi);
    assert(!strcmp(smi, expected_frag_smiles_non_sanitized[i]));
    free(smi);
    free(frags_mpkl_array[i]);
    frags_mpkl_array[i] = NULL;
  }
  free(frags_mpkl_array);
  frags_mpkl_array = NULL;
  free(frags_pkl_sz_array);
  frags_pkl_sz_array = NULL;
  free(mpkl);
  mpkl = NULL;

  printf("  done\n");
  printf("--------------------------\n");
}

void get_wedged_mol_and_inverted_wedges(char **wedged_pkl,
                                        size_t *wedged_pkl_size,
                                        char **inverted_wedges) {
  *wedged_pkl = get_mol(
      "\n\
     RDKit          2D\n\
\n\
 29 34  0  0  1  0  0  0  0  0999 V2000\n\
    1.3719    5.1304    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.5985    3.7907    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.9482    3.7907    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.7216    5.1304    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -3.2685    5.1304    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -3.8994    3.5835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -2.5597    4.3569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -2.5597    5.9038    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -3.8994    6.6771    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -5.2389    5.9038    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -6.5784    6.6771    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -5.2389    4.3569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    1.3719    2.4510    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.5985    1.1115    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    1.3719   -0.2276    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    2.9188   -0.2276    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    3.6921    1.1115    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    2.9188    2.4510    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    5.2389    1.1115    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    6.0124   -0.2276    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    5.2389   -1.5673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    3.6921   -1.5673    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n\
    3.8996   -5.0201    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    5.2391   -4.2467    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    3.5777   -6.5331    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    4.9909   -5.9040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    6.0124   -2.9070    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    6.3306   -6.6772    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    6.5784   -5.0201    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
  2  1  1  1\n\
  2  3  1  0\n\
  3  4  1  0\n\
  5  4  1  6\n\
  5  6  1  0\n\
  6  7  1  0\n\
  7  8  1  0\n\
  9  8  1  1\n\
  5  9  1  0\n\
  9 10  1  0\n\
 10 11  1  1\n\
 10 12  1  0\n\
  6 12  1  1\n\
  2 13  1  0\n\
 13 14  2  0\n\
 14 15  1  0\n\
 15 16  2  0\n\
 16 17  1  0\n\
 17 18  2  0\n\
 13 18  1  0\n\
 17 19  1  0\n\
 19 20  1  0\n\
 20 21  1  0\n\
 21 22  1  0\n\
 16 22  1  0\n\
 23 24  1  0\n\
 23 25  1  0\n\
 25 26  1  0\n\
 24 27  1  0\n\
 27 26  1  0\n\
 26 28  1  0\n\
 24 29  1  0\n\
 28 29  1  0\n\
 21 27  1  0\n\
M  END\n",
      wedged_pkl_size, "");
  assert(*wedged_pkl);
  *inverted_wedges = strdup(
      "  2  1  1  6\n\
  2  3  1  0\n\
  3  4  1  0\n\
  5  4  1  1\n\
  5  6  1  0\n\
  6  7  1  0\n\
  7  8  1  0\n\
  9  8  1  6\n\
  5  9  1  0\n\
  9 10  1  0\n\
 10 11  1  6\n\
 10 12  1  0\n\
  6 12  1  6\n\
  2 13  1  0\n\
 13 14  2  0\n\
 14 15  1  0\n\
 15 16  2  0\n\
 16 17  1  0\n\
 17 18  2  0\n\
 13 18  1  0\n\
 17 19  1  0\n\
 19 20  1  0\n\
 20 21  1  0\n\
 21 22  1  0\n\
 16 22  1  0\n\
 23 24  1  0\n\
 23 25  1  0\n\
 25 26  1  0\n\
 24 27  1  0\n\
 27 26  1  0\n\
 26 28  1  0\n\
 24 29  1  0\n\
 28 29  1  0\n\
 21 27  1  0\n");
  assert(*inverted_wedges);
}

void test_wedging_all_within_scaffold() {
  printf("--------------------------\n");
  printf("  test_wedging_all_within_scaffold\n");
  char *mpkl;
  size_t mpkl_size;
  char *inverted_wedges;
  get_wedged_mol_and_inverted_wedges(&mpkl, &mpkl_size, &inverted_wedges);
  size_t tpkl_size;
  char *tpkl = get_mol(
      "\n\
     RDKit          2D\n\
\n\
 13 14  0  0  1  0  0  0  0  0999 V2000\n\
   -1.6549    2.5755    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.8814    1.2358    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.6653    1.2358    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    1.4385    2.5755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    2.9854    2.5755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    3.6161    1.0286    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    2.2766    1.8019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    2.2766    3.3487    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    3.6161    4.1222    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    4.9558    3.3487    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    6.2953    4.1222    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n\
    4.9558    1.8019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.6549   -0.1037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
  2  1  1  0\n\
  2  3  1  0\n\
  3  4  1  0\n\
  5  4  1  1\n\
  5  6  1  0\n\
  6  7  1  6\n\
  7  8  1  0\n\
  9  8  1  6\n\
  5  9  1  0\n\
  9 10  1  0\n\
 10 11  1  6\n\
 10 12  1  0\n\
  6 12  1  0\n\
  2 13  1  6\n\
M  END\n",
      &tpkl_size, "");
  // the "alignOnly" alignment should succeed and preserve molblock wedging
  // (inverted with respect to the original molecule)
  // it should feature a narrow angle between the bridge bonds
  // as the original geometry of the bridge is preserved
  size_t mpkl_copy_size;
  char *mpkl_copy;
  char *molblock;
  char *svg;
  mpkl_copy_size = mpkl_size;
  mpkl_copy = malloc(mpkl_size);
  assert(mpkl_copy);
  memcpy(mpkl_copy, mpkl, mpkl_size);
  assert(set_2d_coords_aligned(&mpkl_copy, &mpkl_copy_size, tpkl, tpkl_size,
                               "{\"acceptFailure\":false,\"alignOnly\":true}",
                               NULL));
  molblock =
      get_molblock(mpkl_copy, mpkl_copy_size, "{\"useMolBlockWedging\":true}");
  svg = get_svg(
      mpkl_copy, mpkl_copy_size,
      "{\"width\":350,\"height\":300,\"useMolBlockWedging\":true,\"wedgeBonds\":true,\"addChiralHs\":false}");
  double xy23[2];
  double xy26[2];
  double xy25[2];
  double v1[2];
  double v2[2];
  assert(extract_bond_coords(svg, "atom-23 atom-26", xy23, xy26));
  assert(extract_bond_coords(svg, "atom-26 atom-25", NULL, xy25));
  v1[0] = xy23[0] - xy26[0];
  v1[1] = xy23[1] - xy26[1];
  v2[0] = xy25[0] - xy26[0];
  v2[1] = xy25[1] - xy26[1];
  double v1v2Theta = angle_deg_between_vectors(v1, v2);
  assert(v1v2Theta > 10.0 && v1v2Theta < 15.0);
  assert(strstr(molblock, inverted_wedges));
  free(mpkl_copy);
  free(molblock);
  free(svg);
  // the "rebuild" alignment should succeed and preserve molblock wedging
  // (inverted with respect to the original molecule)
  // it should feature a much wider angle between the bridge bonds as the
  // bridged system is entirely rebuilt since it is not part of the scaffold
  mpkl_copy_size = mpkl_size;
  mpkl_copy = malloc(mpkl_size);
  assert(mpkl_copy);
  memcpy(mpkl_copy, mpkl, mpkl_size);
  assert(set_2d_coords_aligned(&mpkl_copy, &mpkl_copy_size, tpkl, tpkl_size,
                               "{\"acceptFailure\":false}", NULL));
  molblock =
      get_molblock(mpkl_copy, mpkl_copy_size, "{\"useMolBlockWedging\":true}");
  svg = get_svg(
      mpkl_copy, mpkl_copy_size,
      "{\"width\":350,\"height\":300,\"useMolBlockWedging\":true,\"wedgeBonds\":true,\"addChiralHs\":false}");
  assert(extract_bond_coords(svg, "atom-23 atom-26", xy23, xy26));
  assert(extract_bond_coords(svg, "atom-26 atom-25", NULL, xy25));
  v1[0] = xy23[0] - xy26[0];
  v1[1] = xy23[1] - xy26[1];
  v2[0] = xy25[0] - xy26[0];
  v2[1] = xy25[1] - xy26[1];
  v1v2Theta = angle_deg_between_vectors(v1, v2);
  assert(v1v2Theta > 105.0 && v1v2Theta < 110.0);
  assert(strstr(molblock, inverted_wedges));
  free(mpkl_copy);
  free(molblock);
  free(svg);
  // the "rebuildCoordGen" alignment should succeed and clear original wedging
  // it should feature an even wider angle between the bridge bonds as CoordGen
  // has a template for the bridged system.
  // Additionally, CoordGen also rebuilds the scaffold, therefore original
  // wedging should be cleared
  mpkl_copy_size = mpkl_size;
  mpkl_copy = malloc(mpkl_size);
  assert(mpkl_copy);
  memcpy(mpkl_copy, mpkl, mpkl_size);
  assert(set_2d_coords_aligned(&mpkl_copy, &mpkl_copy_size, tpkl, tpkl_size,
                               "{\"acceptFailure\":false,\"useCoordGen\":true}",
                               NULL));
  molblock =
      get_molblock(mpkl_copy, mpkl_copy_size, "{\"useMolBlockWedging\":true}");
  svg = get_svg(
      mpkl_copy, mpkl_copy_size,
      "{\"width\":350,\"height\":300,\"useMolBlockWedging\":true,\"wedgeBonds\":true,\"addChiralHs\":false}");
  assert(extract_bond_coords(svg, "atom-23 atom-26", xy23, xy26));
  assert(extract_bond_coords(svg, "atom-26 atom-25", NULL, xy25));
  v1[0] = xy23[0] - xy26[0];
  v1[1] = xy23[1] - xy26[1];
  v2[0] = xy25[0] - xy26[0];
  v2[1] = xy25[1] - xy26[1];
  v1v2Theta = angle_deg_between_vectors(v1, v2);
  assert(v1v2Theta > 145.0 && v1v2Theta < 150.0);
  assert(!strstr(molblock, inverted_wedges));
  free(mpkl_copy);
  free(molblock);
  free(svg);
  free(mpkl);
  free(inverted_wedges);
  free(tpkl);
}

void test_wedging_outside_scaffold() {
  printf("--------------------------\n");
  printf("  test_wedging_outside_scaffold\n");
  char *mpkl;
  size_t mpkl_size;
  char *inverted_wedges;
  get_wedged_mol_and_inverted_wedges(&mpkl, &mpkl_size, &inverted_wedges);
  size_t tpkl_size;
  char *tpkl = get_mol(
      "\n\
     RDKit          2D\n\
\n\
  9 10  0  0  1  0  0  0  0  0999 V2000\n\
   -0.8816    0.5663    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.6651    0.5663    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    1.2958   -0.9804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.0435   -0.2072    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.0435    1.3395    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    1.2958    2.1129    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    2.6355    1.3395    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    3.9750    2.1129    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n\
    2.6355   -0.2072    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
  2  1  1  1\n\
  2  3  1  0\n\
  3  4  1  6\n\
  4  5  1  0\n\
  6  5  1  6\n\
  2  6  1  0\n\
  6  7  1  0\n\
  7  8  1  6\n\
  7  9  1  0\n\
  3  9  1  0\n\
M  END\n",
      &tpkl_size, "");
  // the "alignOnly" alignment should succeed and preserve molblock wedging
  // (inverted with respect to the original molecule)
  // it should feature a narrow angle between the bridge bonds
  // as the original geometry of the bridge is preserved
  size_t mpkl_copy_size;
  char *mpkl_copy;
  char *molblock;
  char *svg;
  mpkl_copy_size = mpkl_size;
  mpkl_copy = malloc(mpkl_size);
  assert(mpkl_copy);
  memcpy(mpkl_copy, mpkl, mpkl_size);
  assert(set_2d_coords_aligned(&mpkl_copy, &mpkl_copy_size, tpkl, tpkl_size,
                               "{\"acceptFailure\":false,\"alignOnly\":true}",
                               NULL));
  molblock =
      get_molblock(mpkl_copy, mpkl_copy_size, "{\"useMolBlockWedging\":true}");
  svg = get_svg(
      mpkl_copy, mpkl_copy_size,
      "{\"width\":350,\"height\":300,\"useMolBlockWedging\":true,\"wedgeBonds\":true,\"addChiralHs\":false}");
  double xy23[2];
  double xy26[2];
  double xy25[2];
  double v1[2];
  double v2[2];
  assert(extract_bond_coords(svg, "atom-23 atom-26", xy23, xy26));
  assert(extract_bond_coords(svg, "atom-26 atom-25", NULL, xy25));
  v1[0] = xy23[0] - xy26[0];
  v1[1] = xy23[1] - xy26[1];
  v2[0] = xy25[0] - xy26[0];
  v2[1] = xy25[1] - xy26[1];
  double v1v2Theta = angle_deg_between_vectors(v1, v2);
  assert(v1v2Theta > 10.0 && v1v2Theta < 15.0);
  assert(strstr(molblock, inverted_wedges));
  free(mpkl_copy);
  free(molblock);
  free(svg);
  // the "rebuild" alignment should succeed and clear molblock wedging
  // it should feature a much wider angle between the bridge bonds as the
  // bridged system is entirely rebuilt since it is not part of the scaffold
  mpkl_copy_size = mpkl_size;
  mpkl_copy = malloc(mpkl_size);
  assert(mpkl_copy);
  memcpy(mpkl_copy, mpkl, mpkl_size);
  assert(set_2d_coords_aligned(&mpkl_copy, &mpkl_copy_size, tpkl, tpkl_size,
                               "{\"acceptFailure\":false}", NULL));
  molblock =
      get_molblock(mpkl_copy, mpkl_copy_size, "{\"useMolBlockWedging\":true}");
  svg = get_svg(
      mpkl_copy, mpkl_copy_size,
      "{\"width\":350,\"height\":300,\"useMolBlockWedging\":true,\"wedgeBonds\":true,\"addChiralHs\":false}");
  assert(extract_bond_coords(svg, "atom-23 atom-26", xy23, xy26));
  assert(extract_bond_coords(svg, "atom-26 atom-25", NULL, xy25));
  v1[0] = xy23[0] - xy26[0];
  v1[1] = xy23[1] - xy26[1];
  v2[0] = xy25[0] - xy26[0];
  v2[1] = xy25[1] - xy26[1];
  v1v2Theta = angle_deg_between_vectors(v1, v2);
  assert(v1v2Theta > 105.0 && v1v2Theta < 110.0);
  assert(!strstr(molblock, inverted_wedges));
  free(mpkl_copy);
  free(molblock);
  free(svg);
  // the "rebuildCoordGen" alignment should succeed and clear original wedging
  // it should feature an even wider angle between the bridge bonds as CoordGen
  // has a template for the bridged system.
  // Additionally, CoordGen also rebuilds the scaffold, therefore original
  // wedging should be cleared
  mpkl_copy_size = mpkl_size;
  mpkl_copy = malloc(mpkl_size);
  assert(mpkl_copy);
  memcpy(mpkl_copy, mpkl, mpkl_size);
  assert(set_2d_coords_aligned(&mpkl_copy, &mpkl_copy_size, tpkl, tpkl_size,
                               "{\"acceptFailure\":false,\"useCoordGen\":true}",
                               NULL));
  molblock =
      get_molblock(mpkl_copy, mpkl_copy_size, "{\"useMolBlockWedging\":true}");
  svg = get_svg(
      mpkl_copy, mpkl_copy_size,
      "{\"width\":350,\"height\":300,\"useMolBlockWedging\":true,\"wedgeBonds\":true,\"addChiralHs\":false}");
  assert(extract_bond_coords(svg, "atom-23 atom-26", xy23, xy26));
  assert(extract_bond_coords(svg, "atom-26 atom-25", NULL, xy25));
  v1[0] = xy23[0] - xy26[0];
  v1[1] = xy23[1] - xy26[1];
  v2[0] = xy25[0] - xy26[0];
  v2[1] = xy25[1] - xy26[1];
  v1v2Theta = angle_deg_between_vectors(v1, v2);
  assert(v1v2Theta > 145.0 && v1v2Theta < 150.0);
  assert(!strstr(molblock, inverted_wedges));
  free(mpkl_copy);
  free(molblock);
  free(svg);
  free(mpkl);
  free(inverted_wedges);
  free(tpkl);
}

void test_wedging_if_no_match() {
  printf("--------------------------\n");
  printf("  test_wedging_if_no_match\n");
  char *mpkl;
  size_t mpkl_size;
  char *inverted_wedges;
  get_wedged_mol_and_inverted_wedges(&mpkl, &mpkl_size, &inverted_wedges);
  size_t tpkl_size;
  char *tpkl = get_mol(
      "\n\
     RDKit          2D\n\
\n\
 13 14  0  0  1  0  0  0  0  0999 V2000\n\
   -1.6549    2.5755    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.8814    1.2358    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.6653    1.2358    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    1.4385    2.5755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    2.9854    2.5755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    3.6161    1.0286    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    2.2766    1.8019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    2.2766    3.3487    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    3.6161    4.1222    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    4.9558    3.3487    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    6.2953    4.1222    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n\
    4.9558    1.8019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.6549   -0.1037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
  2  1  1  0\n\
  2  3  1  0\n\
  3  4  1  0\n\
  5  4  1  1\n\
  5  6  1  0\n\
  6  7  1  6\n\
  7  8  1  0\n\
  9  8  1  6\n\
  5  9  1  0\n\
  9 10  1  0\n\
 10 11  1  6\n\
 10 12  1  0\n\
  6 12  1  0\n\
  2 13  1  6\n\
M  END\n",
      &tpkl_size, "");
  char *orig_molblock =
      get_molblock(mpkl, mpkl_size, "{\"useMolBlockWedging\":true}");
  // the "alignOnly" alignment should return "" if acceptFailure is false
  // and preserve the original coordinates
  char *mpkl_copy;
  size_t mpkl_copy_size;
  char *molblock;
  char *match;
  mpkl_copy_size = mpkl_size;
  mpkl_copy = malloc(mpkl_size);
  assert(mpkl_copy);
  memcpy(mpkl_copy, mpkl, mpkl_size);
  assert(!set_2d_coords_aligned(&mpkl_copy, &mpkl_copy_size, tpkl, tpkl_size,
                                "{\"acceptFailure\":false,\"alignOnly\":true}",
                                &match));
  assert(!match);
  molblock =
      get_molblock(mpkl_copy, mpkl_copy_size, "{\"useMolBlockWedging\":true}");
  assert(!strcmp(molblock, orig_molblock));
  assert(!strstr(molblock, inverted_wedges));
  free(mpkl_copy);
  free(molblock);
  // the "alignOnly" alignment should return "{}" if acceptFailure is true
  // and generate new coordinates, hence wedging should be cleared
  mpkl_copy_size = mpkl_size;
  mpkl_copy = malloc(mpkl_size);
  assert(mpkl_copy);
  memcpy(mpkl_copy, mpkl, mpkl_size);
  assert(set_2d_coords_aligned(&mpkl_copy, &mpkl_copy_size, tpkl, tpkl_size,
                               "{\"acceptFailure\":true,\"alignOnly\":true}",
                               &match));
  assert(!strcmp(match, "{}"));
  free(match);
  molblock =
      get_molblock(mpkl_copy, mpkl_copy_size, "{\"useMolBlockWedging\":true}");
  assert(strcmp(molblock, orig_molblock));
  assert(!strstr(molblock, inverted_wedges));
  free(mpkl_copy);
  free(molblock);
  // the "rebuild" alignment should return "" if acceptFailure is false
  // and preserve the original coordinates
  mpkl_copy_size = mpkl_size;
  mpkl_copy = malloc(mpkl_size);
  assert(mpkl_copy);
  memcpy(mpkl_copy, mpkl, mpkl_size);
  assert(!set_2d_coords_aligned(&mpkl_copy, &mpkl_copy_size, tpkl, tpkl_size,
                                "{\"acceptFailure\":false}", &match));
  assert(!match);
  molblock =
      get_molblock(mpkl_copy, mpkl_copy_size, "{\"useMolBlockWedging\":true}");
  assert(!strcmp(molblock, orig_molblock));
  assert(!strstr(molblock, inverted_wedges));
  free(mpkl_copy);
  free(molblock);
  // the "rebuild" alignment should return "{}" if acceptFailure is true
  // and generate new coordinates, hence wedging should be cleared
  mpkl_copy_size = mpkl_size;
  mpkl_copy = malloc(mpkl_size);
  assert(mpkl_copy);
  memcpy(mpkl_copy, mpkl, mpkl_size);
  assert(set_2d_coords_aligned(&mpkl_copy, &mpkl_copy_size, tpkl, tpkl_size,
                               "{\"acceptFailure\":true}", &match));
  assert(!strcmp(match, "{}"));
  free(match);
  molblock =
      get_molblock(mpkl_copy, mpkl_copy_size, "{\"useMolBlockWedging\":true}");
  assert(strcmp(molblock, orig_molblock));
  assert(!strstr(molblock, inverted_wedges));
  free(mpkl_copy);
  free(molblock);
  // the "rebuildCoordGen" alignment should return "" if acceptFailure is false
  // and preserve the original coordinates
  mpkl_copy_size = mpkl_size;
  mpkl_copy = malloc(mpkl_size);
  assert(mpkl_copy);
  memcpy(mpkl_copy, mpkl, mpkl_size);
  assert(!set_2d_coords_aligned(
      &mpkl_copy, &mpkl_copy_size, tpkl, tpkl_size,
      "{\"acceptFailure\":false,\"useCoordGen\":true}", &match));
  assert(!match);
  molblock =
      get_molblock(mpkl_copy, mpkl_copy_size, "{\"useMolBlockWedging\":true}");
  assert(!strcmp(molblock, orig_molblock));
  assert(!strstr(molblock, inverted_wedges));
  free(mpkl_copy);
  free(molblock);
  // the "rebuildCoordGen" alignment should return "{}" if acceptFailure is true
  // and generate new coordinates, hence wedging should be cleared
  mpkl_copy_size = mpkl_size;
  mpkl_copy = malloc(mpkl_size);
  assert(mpkl_copy);
  memcpy(mpkl_copy, mpkl, mpkl_size);
  assert(set_2d_coords_aligned(&mpkl_copy, &mpkl_copy_size, tpkl, tpkl_size,
                               "{\"acceptFailure\":true,\"useCoordGen\":true}",
                               &match));
  assert(!strcmp(match, "{}"));
  free(match);
  molblock =
      get_molblock(mpkl_copy, mpkl_copy_size, "{\"useMolBlockWedging\":true}");
  assert(strcmp(molblock, orig_molblock));
  assert(!strstr(molblock, inverted_wedges));
  free(mpkl_copy);
  free(molblock);
  free(mpkl);
  free(inverted_wedges);
  free(orig_molblock);
  free(tpkl);
}

void test_removehs() {
  printf("--------------------------\n");
  printf("  test_removehs\n");
  size_t mpkl_size;
  char *mpkl = get_mol("N1C=CC(=O)c2ccc(N(C)(C)(C)(C)C)cc12", &mpkl_size,
                       "{\"sanitize\":false,\"removeHs\":false}");
  char *smi = get_smiles(mpkl, mpkl_size, "");
  assert(mpkl);
  assert(!strcmp(smi, "CN(C)(C)(C)(C)c1ccc2c(c1)NC=CC2=O"));
  free(smi);
  free(mpkl);
}

void test_use_legacy_stereo() {
  printf("--------------------------\n");
  printf("  test_use_legacy_stereo\n");
  short orig_setting = use_legacy_stereo_perception(1);
  char *mpkl;
  size_t mpkl_size;
  mpkl = get_mol("O[C@@]1(C)C/C(/C1)=C(/C)\\CC", &mpkl_size, "");
  assert(mpkl);
  assert(mpkl_size > 0);
  char *smiles = get_smiles(mpkl, mpkl_size, NULL);
  assert(!strcmp(smiles, "CCC(C)=C1CC(C)(O)C1"));
  free(smiles);
  free(mpkl);
  use_legacy_stereo_perception(0);
  mpkl = get_mol("O[C@@]1(C)C/C(/C1)=C(/C)\\CC", &mpkl_size, "");
  assert(mpkl);
  assert(mpkl_size > 0);
  smiles = get_smiles(mpkl, mpkl_size, NULL);
  assert(!strcmp(smiles, "CC/C(C)=C1\\C[C@](C)(O)C1"));
  free(smiles);
  free(mpkl);
  use_legacy_stereo_perception(orig_setting);
}

void test_allow_non_tetrahedral_chirality() {
  printf("--------------------------\n");
  printf("  test_allow_non_tetrahedral_chirality\n");
  char ctab[] =
      "\n\
  Mrv2108 09132105183D          \n\
\n\
  5  4  0  0  0  0            999 V2000\n\
   -1.2500    1.4518    0.0000 Pt  0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.2500    2.2768    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.4250    1.4518    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n\
   -2.0750    1.4518    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.2500    0.6268    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n\
  1  2  1  0  0  0  0\n\
  1  3  1  0  0  0  0\n\
  1  4  1  0  0  0  0\n\
  1  5  1  0  0  0  0\n\
M  END\n";
  short orig_setting = allow_non_tetrahedral_chirality(1);
  char *mpkl;
  size_t mpkl_size;
  mpkl = get_mol(ctab, &mpkl_size, "");
  assert(mpkl);
  assert(mpkl_size > 0);
  char *smiles = get_smiles(mpkl, mpkl_size, NULL);
  assert(!strcmp(smiles, "F[Pt@SP3](F)(Cl)Cl"));
  free(smiles);
  free(mpkl);
  allow_non_tetrahedral_chirality(0);
  mpkl = get_mol(ctab, &mpkl_size, "");
  assert(mpkl);
  assert(mpkl_size > 0);
  smiles = get_smiles(mpkl, mpkl_size, NULL);
  assert(!strcmp(smiles, "F[Pt](F)(Cl)Cl"));
  free(smiles);
  free(mpkl);
  allow_non_tetrahedral_chirality(orig_setting);
}

void test_query_colour() {
  printf("--------------------------\n");
  printf("  test_queryColour\n");
  char smarts[] = "c1ccc2nc([*:1])nc([*:2])c2c1";
  char *pkl;
  size_t pkl_size;
  pkl = get_qmol(smarts, &pkl_size, "");
  assert(pkl);
  assert(pkl_size > 0);
  char *svg = get_svg(pkl, pkl_size, "{\"width\":350,\"height\":300}");
  assert(strstr(svg, "#7F7F7F"));
  assert(strstr(svg, "</svg>"));
  free(svg);
  svg = get_svg(pkl, pkl_size,
                "{\"width\":350,\"height\":300,\"queryColour\":[0.0,0.0,0.0]}");
  assert(!strstr(svg, "#7F7F7F"));
  assert(strstr(svg, "</svg>"));
  free(svg);
  free(pkl);
}

void test_alignment_r_groups_aromatic_ring() {
  printf("--------------------------\n");
  printf("  test_alignment_r_groups_aromatic_ring\n");
  char *mpkl;
  size_t mpkl_size;
  char *tpkl;
  size_t tpkl_size;
  char *match_json = NULL;
  mpkl = get_mol("c1ccc2nccnc2c1", &mpkl_size, "");
  assert(mpkl);
  assert(mpkl_size > 0);
  tpkl = get_mol(
      "\n\
  MJ201100                      \n\
\n\
  8  8  0  0  0  0  0  0  0  0999 V2000\n\
   -1.0263   -0.3133    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n\
   -2.4553    0.5116    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.7408   -0.7258    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.7408   -1.5509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -2.4553   -1.9633    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -3.1698   -1.5509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -3.1698   -0.7258    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -2.4553   -0.3133    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
  3  1  1  0  0  0  0\n\
  8  2  1  0  0  0  0\n\
  4  3  2  0  0  0  0\n\
  5  4  1  0  0  0  0\n\
  6  5  2  0  0  0  0\n\
  7  6  1  0  0  0  0\n\
  8  3  1  0  0  0  0\n\
  8  7  2  0  0  0  0\n\
M  RGP  2   1   2   2   1\n\
M  END\n",
      &tpkl_size, "");
  assert(tpkl);
  assert(tpkl_size > 0);
  assert(set_2d_coords_aligned(
      &mpkl, &mpkl_size, tpkl, tpkl_size,
      "{\"acceptFailure\":false,\"allowRGroups\":true}", &match_json));
  assert(match_json);
  assert(!strcmp(match_json,
                 "{\"atoms\":[4,7,3,2,1,0,9,8],\"bonds\":[3,7,2,1,0,9,10,8]}"));
  free(match_json);
  free(mpkl);
  mpkl = get_mol(
      "\n\
  MJ201100                      \n\
\n\
 10 11  0  0  0  0  0  0  0  0999 V2000\n\
    3.6937    2.5671    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    2.8687    2.5671    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    2.4561    1.8526    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    2.8687    1.1382    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    3.6937    1.1381    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    4.1062    1.8526    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    4.9313    1.8527    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n\
    5.3438    2.5671    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    4.9313    3.2816    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    4.1062    3.2816    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n\
  5  6  1  0  0  0  0\n\
  4  5  2  0  0  0  0\n\
  3  4  1  0  0  0  0\n\
  2  3  2  0  0  0  0\n\
  1  6  2  0  0  0  0\n\
  1  2  1  0  0  0  0\n\
  8  9  1  0  0  0  0\n\
  9 10  2  0  0  0  0\n\
 10  1  1  0  0  0  0\n\
  7  8  2  0  0  0  0\n\
  6  7  1  0  0  0  0\n\
M  END\n",
      &mpkl_size, "");
  assert(mpkl);
  assert(mpkl_size > 0);
  assert(set_2d_coords_aligned(
      &mpkl, &mpkl_size, tpkl, tpkl_size,
      "{\"acceptFailure\":false,\"allowRGroups\":true,\"alignOnly\":true}",
      &match_json));
  assert(match_json);
  assert(!strcmp(match_json,
                 "{\"atoms\":[6,9,5,4,3,2,1,0],\"bonds\":[10,8,0,1,2,3,4,5]}"));
  free(match_json);
  free(mpkl);
  free(tpkl);
}

void test_partial_sanitization() {
  printf("--------------------------\n");
  printf("  test_partial_sanitization\n");
  char *mpkl;
  char *fp;
  size_t mpkl_size;
  const char *mfp_json = "{\"radius\":2,\"nBits\":32}";
  const char *otherfp_json = "{\"nBits\":32}";
  mpkl =
      get_mol("C1CCC2CCCC2C1", &mpkl_size,
              "{\"sanitize\":false,\"removeHs\":false,\"assignStereo\":false}");
  assert(mpkl);
  assert(mpkl_size > 0);
  fp = get_morgan_fp(mpkl, mpkl_size, mfp_json);
  assert(fp);
  assert(strlen(fp) == 32);
  free(fp);
  free(mpkl);
  mpkl = get_mol(
      "C1CCC2CCCC2C1", &mpkl_size,
      "{\"sanitize\":false,\"removeHs\":false,\"assignStereo\":false,\"fastFindRings\":false}");
  assert(mpkl);
  assert(mpkl_size > 0);
  fp = get_morgan_fp(mpkl, mpkl_size, mfp_json);
  assert(!fp);
  fp = get_rdkit_fp(mpkl, mpkl_size, otherfp_json);
  assert(fp);
  free(fp);
  fp = get_pattern_fp(mpkl, mpkl_size, otherfp_json);
  assert(fp);
  free(fp);
  fp = get_atom_pair_fp(mpkl, mpkl_size, otherfp_json);
  assert(fp);
  free(fp);
  fp = get_maccs_fp(mpkl, mpkl_size);
  assert(!fp);
#ifdef RDK_BUILD_AVALON_SUPPORT
  fp = get_avalon_fp(mpkl, mpkl_size, otherfp_json);
  assert(fp);
  free(fp);
#endif
  free(mpkl);
}

void test_capture_logs() {
  printf("--------------------------\n");
  printf("  test_capture_logs\n");
  char *mpkl;
  char *log_buffer;
  void *null_handle = NULL;
  size_t mpkl_size;
  void *log_handle;
  typedef struct {
    const char *type;
    void *(*func)(const char *);
  } capture_test;
  capture_test tests[] = {{"tee", set_log_tee}, {"capture", set_log_capture}};
  for (size_t i = 0; i < sizeof(tests) / sizeof(capture_test); ++i) {
    printf("%zu. %s\n", i + 1, tests[i].type);
    log_handle = tests[i].func("dummy");
    assert(!log_handle);
    log_handle = tests[i].func("rdApp.*");
    assert(log_handle);
    assert(!get_log_buffer(null_handle));
    log_buffer = get_log_buffer(log_handle);
    assert(log_buffer);
    assert(!strlen(log_buffer));
    free(log_buffer);
    mpkl = get_mol("CN(C)(C)C", &mpkl_size, "");
    assert(!mpkl);
    log_buffer = get_log_buffer(log_handle);
    assert(log_buffer);
    assert(strstr(
        log_buffer,
        "Explicit valence for atom # 1 N, 4, is greater than permitted"));
    free(log_buffer);
    assert(!clear_log_buffer(null_handle));
    assert(clear_log_buffer(log_handle));
    log_buffer = get_log_buffer(log_handle);
    assert(log_buffer);
    assert(!strlen(log_buffer));
    free(log_buffer);
    assert(!destroy_log_handle(null_handle));
    assert(!destroy_log_handle(&null_handle));
    assert(destroy_log_handle(&log_handle));
    assert(!log_handle);
  }
}

void test_relabel_mapped_dummies() {
  printf("--------------------------\n");
  printf("  test_relabel_mapped_dummies\n");
  char *mpkl;
  size_t mpkl_size;
  char *smiles;
  mpkl = get_mol("c1cc([4*:2])c([3*:1])cn1", &mpkl_size, "");
  smiles = get_cxsmiles(mpkl, mpkl_size, NULL);
  assert(!strcmp(
      smiles,
      "c1cc([4*:2])c([3*:1])cn1 |atomProp:3.molAtomMapNumber.2:3.dummyLabel.*:5.molAtomMapNumber.1:5.dummyLabel.*|"));
  free(smiles);
  free(mpkl);
  mpkl = get_mol("c1cc([4*:2])c([3*:1])cn1", &mpkl_size,
                 "{\"mappedDummiesAreRGroups\":true}");
  smiles = get_cxsmiles(mpkl, mpkl_size, NULL);
  assert(
      !strcmp(smiles, "*c1ccncc1* |atomProp:0.dummyLabel.R2:7.dummyLabel.R1|"));
  free(smiles);
  free(mpkl);
}

unsigned int count_matches(const char *svg, const char **stereo_array,
                           size_t stereo_array_len) {
  char *svg_copy = strdup(svg);
  unsigned int i = 0;
  char *line = strtok(svg_copy, "\n");
  while (line && i < stereo_array_len) {
    if (strstr(line, stereo_array[i])) {
      ++i;
    } else if (i) {
      break;
    }
    line = strtok(NULL, "\n");
  }
  free(svg_copy);
  return i;
}

void test_assign_cip_labels() {
  printf("--------------------------\n");
  printf("  test_assign_cip_labels\n");
  char *mpkl;
  size_t mpkl_size;
  char *svg;
  static const char *STEREO_SMI = "C/C=C/c1ccccc1[S@@](C)=O";
  static const char *S_STEREO[3] = {">(<", ">S<", ">)<"};
  static const char *R_STEREO[3] = {">(<", ">R<", ">)<"};
  short orig_setting = use_legacy_stereo_perception(1);
  mpkl = get_mol(STEREO_SMI, &mpkl_size, "");
  svg = get_svg(mpkl, mpkl_size,
                "{\"noFreetype\":true,\"addStereoAnnotation\":true}");
  assert(count_matches(svg, S_STEREO, 3) == 3);
  assert(count_matches(svg, R_STEREO, 3) < 3);
  free(svg);
  free(mpkl);
  use_legacy_stereo_perception(0);
  mpkl = get_mol(STEREO_SMI, &mpkl_size, "");
  svg = get_svg(mpkl, mpkl_size,
                "{\"noFreetype\":true,\"addStereoAnnotation\":true}");
  assert(count_matches(svg, S_STEREO, 3) < 3);
  assert(count_matches(svg, R_STEREO, 3) < 3);
  free(svg);
  free(mpkl);
  mpkl = get_mol(STEREO_SMI, &mpkl_size, "{\"assignCIPLabels\":true}");
  svg = get_svg(mpkl, mpkl_size,
                "{\"noFreetype\":true,\"addStereoAnnotation\":true}");
  assert(count_matches(svg, S_STEREO, 3) < 3);
  assert(count_matches(svg, R_STEREO, 3) == 3);
  free(svg);
  free(mpkl);
  use_legacy_stereo_perception(orig_setting);
}

void test_smiles_smarts_params() {
  printf("--------------------------\n");
  printf("  test_smiles_smarts_params\n");
  const char *amoxicillin_pub_chem =
      "CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)[C@@H](C3=CC=C(C=C3)O)N)C(=O)O)C";
  const char *bicyclo221heptane =
      "\n\
     RDKit          2D\n\
\n\
  9 10  0  0  1  0  0  0  0  0999 V2000\n\
   -2.8237   -1.3088    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.5723   -0.3996    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.5723    1.1473    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.1011    1.6253    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    1.3701    1.1474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    1.3701   -0.3995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    2.6217   -1.3087    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.1009   -0.8775    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.8083    0.3739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
  2  1  1  1\n\
  2  3  1  0\n\
  4  3  1  0\n\
  4  5  1  0\n\
  6  5  1  0\n\
  6  7  1  1\n\
  6  8  1  0\n\
  8  9  1  1\n\
  8  2  1  0\n\
  4  9  1  1\n\
M  END\n\
";
  const char *atom_prop = " |atomProp:1.atomProp.1&#46;234|";
  const char *chiral_smarts = "N-[C@H](-C(-O)=O)-C(-C)-C";
  const char *chiral_smarts_with_atom_prop =
      "N-[C@H](-C(-O)=O)-C(-C)-C |atomProp:1.atomProp.1&#46;234|";
  char *mpkl;
  size_t mpkl_size;
  char *mpkl_atom_prop;
  size_t mpkl_atom_prop_size;
  char *canonical_smiles;
  char *non_canonical_smiles;
  char *canonical_smiles_no_stereo;
  char *non_canonical_smiles_no_stereo;
  char *canonical_cxsmiles;
  char *non_canonical_cxsmiles_no_stereo;
  char *non_canonical_cxsmiles_no_stereo_atom_prop;
  char *cxsmiles_with_atom_prop;
  char *smarts;
  char *ptr;
  char *ptr_end;
  unsigned int i;
  mpkl = get_mol(amoxicillin_pub_chem, &mpkl_size, "");
  assert(mpkl);
  const char *empty_json[3] = {NULL, "", "{}"};
  for (i = 0; i < 3; ++i) {
    canonical_smiles = get_smiles(mpkl, mpkl_size, empty_json[i]);
    assert(canonical_smiles);
    assert(!strcmp(
        canonical_smiles,
        "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O"));
    free(canonical_smiles);
  }
  non_canonical_smiles =
      get_smiles(mpkl, mpkl_size, "{\"canonical\":\"false\"}");
  assert(non_canonical_smiles);
  assert(!strcmp(
      non_canonical_smiles,
      "CC1(C)[C@H](C(=O)O)N2[C@H](S1)[C@H](NC(=O)[C@@H](c1ccc(O)cc1)N)C2=O"));
  free(non_canonical_smiles);
  canonical_smiles_no_stereo =
      get_smiles(mpkl, mpkl_size, "{\"doIsomericSmiles\":false}");
  assert(canonical_smiles_no_stereo);
  assert(!strcmp(canonical_smiles_no_stereo,
                 "CC1(C)SC2C(NC(=O)C(N)c3ccc(O)cc3)C(=O)N2C1C(=O)O"));
  free(canonical_smiles_no_stereo);
  non_canonical_smiles_no_stereo = get_smiles(
      mpkl, mpkl_size, "{\"doIsomericSmiles\":false,\"canonical\":false}");
  assert(non_canonical_smiles_no_stereo);
  assert(!strcmp(non_canonical_smiles_no_stereo,
                 "CC1(C)C(C(=O)O)N2C(S1)C(NC(=O)C(c1ccc(O)cc1)N)C2=O"));
  free(non_canonical_smiles_no_stereo);
  free(mpkl);
  mpkl =
      get_mol(bicyclo221heptane, &mpkl_size, "{\"useMolBlockWedging\":true}");
  assert(mpkl);
  for (i = 0; i < 3; ++i) {
    canonical_cxsmiles = get_cxsmiles(mpkl, mpkl_size, empty_json[i]);
    assert(canonical_cxsmiles);
    ptr = strstr(canonical_cxsmiles, " |");
    assert(ptr);
    *ptr = '\0';
    assert(!strcmp(canonical_cxsmiles, "N[C@@H]1C[C@@H]2C[C@H]1[C@@H](O)C2"));
    ++ptr;
    assert(*ptr == '|');
    ptr = strstr(ptr, "),");
    assert(ptr);
    ptr += 2;
    ptr_end = strstr(ptr, "|");
    assert(ptr_end);
    *ptr_end = '\0';
    assert(!strcmp(ptr, "wD:3.9,wU:1.0,5.4,6.7"));
    free(canonical_cxsmiles);
  }
  canonical_cxsmiles =
      get_cxsmiles(mpkl, mpkl_size,
                   "{\"restoreBondDirOption\":\"RestoreBondDirOptionTrue\"}");
  assert(canonical_cxsmiles);
  ptr = strstr(canonical_cxsmiles, " |");
  assert(ptr);
  *ptr = '\0';
  assert(!strcmp(canonical_cxsmiles, "N[C@@H]1C[C@@H]2C[C@H]1[C@@H](O)C2"));
  ++ptr;
  assert(*ptr == '|');
  ptr = strstr(ptr, "),");
  assert(ptr);
  ptr += 2;
  ptr_end = strstr(ptr, "|");
  assert(ptr_end);
  *ptr_end = '\0';
  assert(!strcmp(ptr, "wU:1.0,3.3,5.4,6.7"));
  free(canonical_cxsmiles);
  non_canonical_cxsmiles_no_stereo = get_cxsmiles(
      mpkl, mpkl_size,
      "{\"doIsomericSmiles\":false,\"canonical\":false,\"CX_ALL_BUT_COORDS\":true}");
  assert(non_canonical_cxsmiles_no_stereo);
  assert(!strcmp(non_canonical_cxsmiles_no_stereo, "OC1CC2CC(N)C1C2"));
  non_canonical_cxsmiles_no_stereo_atom_prop =
      realloc(non_canonical_cxsmiles_no_stereo,
              strlen(non_canonical_cxsmiles_no_stereo) + strlen(atom_prop) + 1);
  assert(non_canonical_cxsmiles_no_stereo_atom_prop);
  assert(strcat(non_canonical_cxsmiles_no_stereo_atom_prop, atom_prop) ==
         non_canonical_cxsmiles_no_stereo_atom_prop);
  mpkl_atom_prop = get_mol(non_canonical_cxsmiles_no_stereo_atom_prop,
                           &mpkl_atom_prop_size, "");
  assert(mpkl_atom_prop);
  free(non_canonical_cxsmiles_no_stereo_atom_prop);
  cxsmiles_with_atom_prop = get_cxsmiles(mpkl_atom_prop, mpkl_atom_prop_size,
                                         "{\"CX_ALL_BUT_COORDS\":true}");
  assert(cxsmiles_with_atom_prop);
  assert(!strcmp(cxsmiles_with_atom_prop,
                 "NC1CC2CC(O)C1C2 |atomProp:5.atomProp.1&#46;234|"));
  free(cxsmiles_with_atom_prop);
  free(mpkl_atom_prop);
  free(mpkl);
  mpkl = get_qmol(chiral_smarts, &mpkl_size, "");
  assert(mpkl);
  for (i = 0; i < 3; ++i) {
    smarts = get_smarts(mpkl, mpkl_size, empty_json[i]);
    assert(smarts);
    assert(!strcmp(smarts, "N-[C@&H1](-C(-O)=O)-C(-C)-C"));
    free(smarts);
  }
  smarts = get_smarts(mpkl, mpkl_size, "{\"doIsomericSmiles\":false}");
  assert(smarts);
  assert(!strcmp(smarts, "N-[C&H1](-C(-O)=O)-C(-C)-C"));
  free(smarts);
  free(mpkl);
  mpkl = get_qmol(chiral_smarts_with_atom_prop, &mpkl_size, "");
  assert(mpkl);
  for (i = 0; i < 3; ++i) {
    smarts = get_cxsmarts(mpkl, mpkl_size, empty_json[i]);
    assert(smarts);
    assert(!strcmp(
        smarts, "N-[C@&H1](-C(-O)=O)-C(-C)-C |atomProp:1.atomProp.1&#46;234|"));
    free(smarts);
  }
  smarts = get_cxsmarts(mpkl, mpkl_size, "{\"doIsomericSmiles\":false}");
  assert(smarts);
  assert(!strcmp(smarts,
                 "N-[C&H1](-C(-O)=O)-C(-C)-C |atomProp:1.atomProp.1&#46;234|"));
  free(smarts);
  free(mpkl);
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
  test_get_mol_frags();
  test_wedging_all_within_scaffold();
  test_wedging_outside_scaffold();
  test_wedging_if_no_match();
  test_removehs();
  test_use_legacy_stereo();
  test_allow_non_tetrahedral_chirality();
  test_query_colour();
  test_alignment_r_groups_aromatic_ring();
  test_partial_sanitization();
  test_capture_logs();
  test_relabel_mapped_dummies();
  test_assign_cip_labels();
  test_smiles_smarts_params();
  return 0;
}
