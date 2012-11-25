// $Id$
//
//  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met: 
//
//     * Redistributions of source code must retain the above copyright 
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following 
//       disclaimer in the documentation and/or other materials provided 
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
//       nor the names of its contributors may be used to endorse or promote 
//       products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include "rdkit.h"
#include "fmgr.h"

/***************** Mol operations ***********************/

#define MOLCMPFUNC( type, action, ret )                                 \
  PG_FUNCTION_INFO_V1(mol_##type);                                      \
  Datum           mol_##type(PG_FUNCTION_ARGS);                         \
  Datum                                                                 \
  mol_##type(PG_FUNCTION_ARGS)                                          \
  {                                                                     \
    CROMol    a, b;                                                     \
    int             res;                                                \
                                                                        \
    fcinfo->flinfo->fn_extra = SearchMolCache(                          \
                                              fcinfo->flinfo->fn_extra, \
                                              fcinfo->flinfo->fn_mcxt,  \
                                              PG_GETARG_DATUM(0),       \
                                              NULL, &a, NULL);          \
    fcinfo->flinfo->fn_extra = SearchMolCache(                          \
                                              fcinfo->flinfo->fn_extra, \
                                              fcinfo->flinfo->fn_mcxt,  \
                                              PG_GETARG_DATUM(1),       \
                                              NULL, &b, NULL);          \
    res = molcmp(a, b);                                                 \
    PG_RETURN_##ret( res action 0 );                                    \
  }                                                                     \
  /* keep compiler quiet - no extra ; */                                \
  extern int no_such_variable

MOLCMPFUNC(lt, <, BOOL);
MOLCMPFUNC(le, <=, BOOL);
MOLCMPFUNC(eq, ==, BOOL);
MOLCMPFUNC(ge, >=, BOOL);
MOLCMPFUNC(gt, >, BOOL);
MOLCMPFUNC(ne, !=, BOOL);
MOLCMPFUNC(cmp, +, INT32);

PG_FUNCTION_INFO_V1(is_valid_smiles);
Datum           is_valid_smiles(PG_FUNCTION_ARGS);
Datum
is_valid_smiles(PG_FUNCTION_ARGS) {
  char  *data = PG_GETARG_CSTRING(0);
  PG_RETURN_BOOL(isValidSmiles(data));
}

PG_FUNCTION_INFO_V1(is_valid_smarts);
Datum           is_valid_smarts(PG_FUNCTION_ARGS);
Datum
is_valid_smarts(PG_FUNCTION_ARGS) {
  char  *data = PG_GETARG_CSTRING(0);
  PG_RETURN_BOOL(isValidSmarts(data));
}

PG_FUNCTION_INFO_V1(is_valid_ctab);
Datum           is_valid_ctab(PG_FUNCTION_ARGS);
Datum
is_valid_ctab(PG_FUNCTION_ARGS) {
  char  *data = PG_GETARG_CSTRING(0);
  PG_RETURN_BOOL(isValidCTAB(data));
}

PG_FUNCTION_INFO_V1(is_valid_mol_pkl);
Datum           is_valid_mol_pkl(PG_FUNCTION_ARGS);
Datum
is_valid_mol_pkl(PG_FUNCTION_ARGS) {
  bytea    *data = PG_GETARG_BYTEA_P(0);
  int len=VARSIZE(data)-VARHDRSZ;
  bool res=isValidMolBlob(VARDATA(data),len);
  PG_FREE_IF_COPY(data, 0);
  PG_RETURN_BOOL(res);           
}


PG_FUNCTION_INFO_V1(mol_substruct);
Datum           mol_substruct(PG_FUNCTION_ARGS);
Datum
mol_substruct(PG_FUNCTION_ARGS) {
  CROMol  i,
    a;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0), 
                                            NULL, &i, NULL);
  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(1), 
                                            NULL, &a, NULL);

  PG_RETURN_BOOL(MolSubstruct(i, a));             
}

PG_FUNCTION_INFO_V1(mol_rsubstruct);
Datum           mol_rsubstruct(PG_FUNCTION_ARGS);
Datum
mol_rsubstruct(PG_FUNCTION_ARGS) {
  CROMol  i,
    a;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0), 
                                            NULL, &a, NULL);
  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(1), 
                                            NULL, &i, NULL);

  PG_RETURN_BOOL(MolSubstruct(i, a));             
}


PG_FUNCTION_INFO_V1(mol_amw);
Datum           mol_amw(PG_FUNCTION_ARGS);
Datum
mol_amw(PG_FUNCTION_ARGS) {
  CROMol        i;
  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0), 
                                            NULL, &i, NULL);
  PG_RETURN_FLOAT4(MolAMW(i));
}
PG_FUNCTION_INFO_V1(mol_logp);
Datum           mol_logp(PG_FUNCTION_ARGS);
Datum
mol_logp(PG_FUNCTION_ARGS) {
  CROMol        i;
  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0), 
                                            NULL, &i, NULL);
  PG_RETURN_FLOAT4(MolLogP(i));
}
PG_FUNCTION_INFO_V1(mol_hba);
Datum           mol_hba(PG_FUNCTION_ARGS);
Datum
mol_hba(PG_FUNCTION_ARGS) {
  CROMol        i;
  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0), 
                                            NULL, &i, NULL);
  PG_RETURN_INT32(MolHBA(i));
}
PG_FUNCTION_INFO_V1(mol_hbd);
Datum           mol_hbd(PG_FUNCTION_ARGS);
Datum
mol_hbd(PG_FUNCTION_ARGS) {
  CROMol        i;
  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0), 
                                            NULL, &i, NULL);
  PG_RETURN_INT32(MolHBD(i));
}
PG_FUNCTION_INFO_V1(mol_numatoms);
Datum           mol_numatoms(PG_FUNCTION_ARGS);
Datum
mol_numatoms(PG_FUNCTION_ARGS) {
  CROMol        i;
  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0), 
                                            NULL, &i, NULL);
  PG_RETURN_INT32(MolNumAtoms(i));
}
PG_FUNCTION_INFO_V1(mol_numheavyatoms);
Datum           mol_numheavyatoms(PG_FUNCTION_ARGS);
Datum
mol_numheavyatoms(PG_FUNCTION_ARGS) {
  CROMol        i;
  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0), 
                                            NULL, &i, NULL);
  PG_RETURN_INT32(MolNumHeavyAtoms(i));
}
PG_FUNCTION_INFO_V1(mol_numrotatablebonds);
Datum           mol_numrotatablebonds(PG_FUNCTION_ARGS);
Datum
mol_numrotatablebonds(PG_FUNCTION_ARGS) {
  CROMol        i;
  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0), 
                                            NULL, &i, NULL);
  PG_RETURN_INT32(MolNumRotatableBonds(i));
}
PG_FUNCTION_INFO_V1(mol_numheteroatoms);
Datum           mol_numheteroatoms(PG_FUNCTION_ARGS);
Datum
mol_numheteroatoms(PG_FUNCTION_ARGS) {
  CROMol        i;
  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0), 
                                            NULL, &i, NULL);
  PG_RETURN_INT32(MolNumHeteroatoms(i));
}
PG_FUNCTION_INFO_V1(mol_numrings);
Datum           mol_numrings(PG_FUNCTION_ARGS);
Datum
mol_numrings(PG_FUNCTION_ARGS) {
  CROMol        i;
  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0), 
                                            NULL, &i, NULL);
  PG_RETURN_INT32(MolNumRings(i));
}
PG_FUNCTION_INFO_V1(mol_tpsa);
Datum           mol_tpsa(PG_FUNCTION_ARGS);
Datum
mol_tpsa(PG_FUNCTION_ARGS) {
  CROMol        i;
  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0), 
                                            NULL, &i, NULL);
  PG_RETURN_FLOAT4(MolTPSA(i));
}

PG_FUNCTION_INFO_V1(mol_inchi);
Datum           mol_inchi(PG_FUNCTION_ARGS);
Datum
mol_inchi(PG_FUNCTION_ARGS) {
  CROMol  mol;
  const char    *str;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);
  str = MolInchi(mol);
  char *res=pnstrdup(str, strlen(str));
  free((void *)str);
  PG_RETURN_CSTRING( res );
}

PG_FUNCTION_INFO_V1(mol_inchikey);
Datum           mol_inchikey(PG_FUNCTION_ARGS);
Datum
mol_inchikey(PG_FUNCTION_ARGS) {
  CROMol  mol;
  const char    *str;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);
  str = MolInchiKey(mol);
  char *res=pnstrdup(str, strlen(str));
  free((void *)str);
  PG_RETURN_CSTRING( res );
}
