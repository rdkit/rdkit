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
#include "cffiwrapper.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>

void test_io(){
  char *pkl;
  size_t pkl_size;
  size_t pkl2_size;
  char *pkl2;
  
  printf("--------------------------\n");
  printf("  test_io\n");
  
  pkl = get_mol("c1cc(O)ccc1",&pkl_size);
  assert(pkl);
  assert(pkl_size>0);

  char *smiles=get_smiles(pkl,pkl_size);
  assert(!strcmp(smiles,"Oc1ccccc1"));
  free(smiles);
  smiles=NULL;

  smiles=get_cxsmiles(pkl,pkl_size);
  assert(!strcmp(smiles,"Oc1ccccc1"));
  free(smiles);
  smiles=NULL;

  char *json=get_json(pkl,pkl_size);
  assert(strstr(json,"commonchem"));

  pkl2=get_mol(json,&pkl2_size);
  assert(pkl2);
  assert(pkl2_size>0);
  smiles=get_smiles(pkl2,pkl2_size);
  assert(!strcmp(smiles,"Oc1ccccc1"));
  free(smiles);
  smiles=NULL;
  free(pkl2);
  pkl2 = NULL;
  free(json);
  json=NULL;

  //---------
  // failures

  // kekulization
  pkl2 = get_mol("c1cccc1",&pkl2_size);
  assert(pkl2==NULL);
  assert(!pkl2_size);

  // bad input
  pkl2 = get_mol("foo",&pkl2_size);
  assert(pkl2==NULL);
  assert(!pkl2_size);
  
  // valence
  pkl2 = get_mol("CO(C)C",&pkl2_size);
  assert(pkl2==NULL);
  assert(!pkl2_size);





  char *molblock = get_molblock(pkl,pkl_size);
  pkl2=get_mol(molblock,&pkl2_size);
  assert(pkl2);
  assert(pkl2_size>0);
  smiles=get_smiles(pkl2,pkl2_size);
  assert(!strcmp(smiles,"Oc1ccccc1"));
  free(smiles);
  smiles=NULL;
  free(pkl2);
  pkl2 = NULL;
  free(molblock);
  molblock=NULL;

  molblock = get_v3kmolblock(pkl,pkl_size);
  pkl2=get_mol(molblock,&pkl2_size);
  assert(pkl2);
  assert(pkl2_size>0);
  smiles=get_smiles(pkl2,pkl2_size);
  assert(!strcmp(smiles,"Oc1ccccc1"));
  free(smiles);
  smiles=NULL;
  free(pkl2);
  pkl2 = NULL;
  free(molblock);
  molblock=NULL;

  char *inchi = get_inchi(pkl,pkl_size);
  assert(!strcmp(inchi,"InChI=1S/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H"));
  free(inchi);
  inchi = get_inchikey_for_inchi("InChI=1S/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H");
  assert(!strcmp(inchi,"ISWSIDIOOBJBQZ-UHFFFAOYSA-N"));
  free(inchi);

  molblock = get_molblock(pkl,pkl_size);
  inchi = get_inchi_for_molblock(molblock);
  assert(!strcmp(inchi,"InChI=1S/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H"));
  free(inchi);
  free(molblock);

  char *smarts = get_smarts(pkl,pkl_size);
  assert(!strcmp(smarts,"[#6]1:[#6]:[#6](-[#8]):[#6]:[#6]:[#6]:1"));
  
  pkl2 = get_qmol(smarts,&pkl2_size);
  assert(pkl2);
  free(smarts);
  smarts = get_smarts(pkl2,pkl2_size);
  assert(!strcmp(smarts,"[#6]1:[#6]:[#6](-[#8]):[#6]:[#6]:[#6]:1"));
  free(smarts);
  free(pkl2);


  free(pkl);
  pkl=NULL;
  printf("  done\n");
  printf("--------------------------\n");
  
}



void test_svg(){
  char *pkl;
  size_t pkl_size;
  
  printf("--------------------------\n");
  printf("  test_svg\n");
  
  pkl = get_mol("c1cc(O)ccc1",&pkl_size);
  assert(pkl);
  assert(pkl_size>0);

  char *svg = get_svg(pkl,pkl_size,"{\"width\":350,\"height\":300}");
  assert(strstr(svg,"width='350px'"));
  assert(strstr(svg,"height='300px'"));
  assert(strstr(svg,"</svg>"));
  free(svg);

  svg = get_svg(pkl,pkl_size,"{\"atoms\": [0, 1, 2], \"width\":127}");
  assert(strstr(svg,"fill:#FF7F7F"));
  assert(strstr(svg,"width='127px'"));
  assert(strstr(svg,"</svg>"));
  free(svg);

  free(pkl);
  pkl=NULL;
  printf("  done\n");
  printf("--------------------------\n"); 
}

void test_substruct(){
  printf("--------------------------\n");
  printf("  test_substruct\n");
  
  char *mpkl;
  size_t mpkl_size;
  mpkl = get_mol("Cl[C@H](F)C[C@H](F)Cl",&mpkl_size);
  char *qpkl;
  size_t qpkl_size;
  qpkl = get_qmol("Cl[C@@H](F)C",&qpkl_size);

  char *json = get_substruct_match(mpkl,mpkl_size,qpkl,qpkl_size,"");
  assert(!strcmp(json,"{\"atoms\":[0,1,2,3],\"bonds\":[0,1,2]}"));
  free(json);
  json = get_substruct_matches(mpkl,mpkl_size,qpkl,qpkl_size,"");
  assert(!strcmp(json,"[{\"atoms\":[0,1,2,3],\"bonds\":[0,1,2]},{\"atoms\":[6,4,5,3],\"bonds\":[5,4,3]}]"));
  free(json);

  // make sure using parameters works
  json = get_substruct_match(mpkl,mpkl_size,qpkl,qpkl_size,"{\"useChirality\":true}");
  assert(!strcmp(json,"{\"atoms\":[6,4,5,3],\"bonds\":[5,4,3]}"));
  free(json);
  json = get_substruct_matches(mpkl,mpkl_size,qpkl,qpkl_size,"{\"useChirality\":true}");
  assert(!strcmp(json,"[{\"atoms\":[6,4,5,3],\"bonds\":[5,4,3]}]"));
  free(json);
  
  free(mpkl);
  mpkl=NULL;
  free(qpkl);
  qpkl=NULL;
  printf("  done\n");
  printf("--------------------------\n"); 
}

void test_descriptors(){
  printf("--------------------------\n");
  printf("  test_descriptors\n");
  
  char *mpkl;
  size_t mpkl_size;
  mpkl = get_mol("c1nccc(O)c1",&mpkl_size);

  char *descrs = get_descriptors(mpkl,mpkl_size);
  assert(strstr(descrs,"lipinskiHBA"));
  assert(strstr(descrs,"NumAliphaticRings"));
  assert(strstr(descrs,"chi3v"));
  free(descrs);

  free(mpkl);
  mpkl=NULL;
  printf("  done\n");
  printf("--------------------------\n"); 
}

void test_fingerprints(){
  printf("--------------------------\n");
  printf("  test_fingerprints\n");
  
  char *mpkl;
  size_t mpkl_size;
  mpkl = get_mol("c1nccc(O)c1",&mpkl_size);
  const char *mfp_json="{\"radius\":2,\"nBits\":1024}";
  char *fp = get_morgan_fp(mpkl,mpkl_size,mfp_json);
  assert(strlen(fp)==1024);
  free(fp);
  fp = get_morgan_fp(mpkl,mpkl_size,"{\"radius\":1,\"nBits\":64}");
  assert(!strcmp(fp,"0011000000100000010000100000000000001001010000000000000000100000"));
  free(fp);
  fp = get_morgan_fp(mpkl,mpkl_size,"{\"radius\":2,\"nBits\":64}");
  assert(!strcmp(fp,"0011000000100000010000100000000000001001010000000010000010100001"));
  free(fp);
  size_t nbytes;
  fp = get_morgan_fp_as_bytes(mpkl,mpkl_size,&nbytes,"{\"radius\":2,\"nBits\":64}");
  assert(nbytes==8);
  free(fp);

  fp = get_rdkit_fp(mpkl,mpkl_size,64);
  assert(!strcmp(fp,"1111011000111100011011011111011001111111110010000111000011111111"));
  free(fp);
  fp = get_rdkit_fp_as_bytes(mpkl,mpkl_size,&nbytes,64);
  assert(nbytes==8);
  free(fp);

  fp = get_pattern_fp(mpkl,mpkl_size,0,64);
  assert(!strcmp(fp,"1011111111111111110011011111011001011110111101111110101111010011"));
  free(fp);
  fp = get_pattern_fp_as_bytes(mpkl,mpkl_size,&nbytes,0,64);
  assert(nbytes==8);
  free(fp);

  fp = get_pattern_fp(mpkl,mpkl_size,1,64);
  assert(!strcmp(fp,"1011111111111111110011111111011101011111111101111111111111011011"));
  free(fp);


  free(mpkl);
  mpkl=NULL;
  printf("  done\n");
  printf("--------------------------\n"); 
}

void test_modifications(){
  printf("--------------------------\n");
  printf("  test_modifications\n");
  char *mpkl;
  size_t mpkl_size;
  mpkl = get_mol("CCC",&mpkl_size);
  
  assert(add_hs(&mpkl,&mpkl_size)>0);
  char *ctab = get_molblock(mpkl,mpkl_size);
  assert(strstr(ctab," H "));
  free(ctab);

  assert(remove_hs(&mpkl,&mpkl_size)>0);
  ctab = get_molblock(mpkl,mpkl_size);
  assert(!strstr(ctab," H "));
  free(ctab);

  free(mpkl);
  mpkl=NULL;
  printf("  done\n");
  printf("--------------------------\n"); 
}


void test_coords(){
  printf("--------------------------\n");
  printf("  test_coords\n");
  char *mpkl;
  size_t mpkl_size;
  mpkl = get_mol("C1CCC1CC",&mpkl_size);
  
  char *cxsmi = get_cxsmiles(mpkl,mpkl_size);
  // no cxsmiles yet
  assert(!strstr(cxsmi,"|"));

  prefer_coordgen(0);
  set_2d_coords(&mpkl,&mpkl_size);
  free(cxsmi);
  cxsmi = get_cxsmiles(mpkl,mpkl_size);
  // since we have coords there's something there:
  assert(strstr(cxsmi,"|"));
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  prefer_coordgen(1);
  set_2d_coords(&mpkl,&mpkl_size);
  char *cxsmi2 = get_cxsmiles(mpkl,mpkl_size);
  assert(strstr(cxsmi2,"|"));
  assert(strcmp(cxsmi,cxsmi2));
  free(cxsmi2);
#endif
  free(cxsmi);

  // 3D
  assert(add_hs(&mpkl,&mpkl_size));
  assert(set_3d_coords(&mpkl,&mpkl_size,"")>0);
  char *cxsmi3 = get_cxsmiles(mpkl,mpkl_size); 
  assert(set_3d_coords(&mpkl,&mpkl_size,"{\"randomSeed\":123}")>0);
  cxsmi = get_cxsmiles(mpkl,mpkl_size);
  // since we have coords there's something there:
  assert(strstr(cxsmi,"|"));
  // coords generated with two different seeds differ:
  assert(strcmp(cxsmi,cxsmi3));
  free(cxsmi3);
  free(cxsmi);


  free(mpkl);
  mpkl=NULL;
  printf("  done\n");
  printf("--------------------------\n"); 
}

int main(){
  enable_logging();
  printf("hello %s\n",version()); 
  test_io();
  test_svg();
  test_substruct();
  test_descriptors();
  test_fingerprints();
  test_modifications();
  test_coords();
  return 0;
}