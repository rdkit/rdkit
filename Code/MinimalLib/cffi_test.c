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

  char *svg = get_svg(pkl,pkl_size,350,300);
  assert(strstr(svg,"width='350px'"));
  assert(strstr(svg,"height='300px'"));
  assert(strstr(svg,"</svg>"));
  free(svg);

  svg = get_svg_with_highlights(pkl,pkl_size,"{\"atoms\": [0, 1, 2], \"width\":127}");
  assert(strstr(svg,"fill:#FF7F7F"));
  assert(strstr(svg,"width='127px'"));
  assert(strstr(svg,"</svg>"));
  free(svg);

  free(pkl);
  pkl=NULL;
  printf("  done\n");
  printf("--------------------------\n");
  
}

int main(){
  printf("hello %s\n",version()); 
  test_io();
  test_svg();
  return 0;
}