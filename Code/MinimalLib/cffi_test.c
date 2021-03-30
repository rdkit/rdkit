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


  free(pkl);
  pkl=NULL;
  printf("  done\n");
  printf("--------------------------\n");
  
}

int main(){
  printf("hello %s\n",version()); 
  test_io();
  return 0;
}