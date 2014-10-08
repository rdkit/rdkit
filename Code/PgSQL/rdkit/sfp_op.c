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

static int 
sfpcmp(SparseFingerPrint *a, SparseFingerPrint *b) {
  int res;

  res = memcmp(VARDATA(a), VARDATA(b), Min(VARSIZE(a), VARSIZE(b)) - VARHDRSZ);
  if ( res )
    return res;

  if (VARSIZE(a) == VARSIZE(b))
    return 0;
  return (VARSIZE(a) > VARSIZE(b)) ? 1 : -1; 
}


#define sfpCMPFUNC( type, action, ret )                                 \
  PG_FUNCTION_INFO_V1(sfp_##type);                                      \
  Datum           sfp_##type(PG_FUNCTION_ARGS);                         \
  Datum                                                                 \
  sfp_##type(PG_FUNCTION_ARGS)                                          \
  {                                                                     \
    SparseFingerPrint    *a, *b;                                        \
    int             res;                                                \
                                                                        \
    fcinfo->flinfo->fn_extra = SearchSparseFPCache(                     \
                                                   fcinfo->flinfo->fn_extra, \
                                                   fcinfo->flinfo->fn_mcxt, \
                                                   PG_GETARG_DATUM(0),  \
                                                   &a, NULL, NULL);     \
    fcinfo->flinfo->fn_extra = SearchSparseFPCache(                     \
                                                   fcinfo->flinfo->fn_extra, \
                                                   fcinfo->flinfo->fn_mcxt, \
                                                   PG_GETARG_DATUM(1),  \
                                                   &b, NULL, NULL);     \
    res = sfpcmp(a, b);                                                 \
    PG_RETURN_##ret( res action 0 );                                    \
  }                                                                     \
  /* keep compiler quiet - no extra ; */                                \
  extern int no_such_variable

sfpCMPFUNC(lt, <, BOOL);
sfpCMPFUNC(le, <=, BOOL);
sfpCMPFUNC(eq, ==, BOOL);
sfpCMPFUNC(ge, >=, BOOL);
sfpCMPFUNC(gt, >, BOOL);
sfpCMPFUNC(ne, !=, BOOL);
sfpCMPFUNC(cmp, +, INT32);

PG_FUNCTION_INFO_V1(sfp_tanimoto_sml);
Datum           sfp_tanimoto_sml(PG_FUNCTION_ARGS);
Datum
sfp_tanimoto_sml(PG_FUNCTION_ARGS) {
  MolSparseFingerPrint    asfp,
    bsfp;
  double                  res;

  fcinfo->flinfo->fn_extra = SearchSparseFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(0), 
                                                 NULL, &asfp, NULL);
  fcinfo->flinfo->fn_extra = SearchSparseFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(1), 
                                                 NULL, &bsfp, NULL);

  res = calcSparseTanimotoSml(asfp, bsfp); 

  PG_RETURN_FLOAT8(res);          
}

PG_FUNCTION_INFO_V1(sfp_tanimoto_sml_op);
Datum           sfp_tanimoto_sml_op(PG_FUNCTION_ARGS);
Datum
sfp_tanimoto_sml_op(PG_FUNCTION_ARGS) {
  MolSparseFingerPrint    asfp,
    bsfp;
  double                  res;

  fcinfo->flinfo->fn_extra = SearchSparseFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(0), 
                                                 NULL, &asfp, NULL);
  fcinfo->flinfo->fn_extra = SearchSparseFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(1), 
                                                 NULL, &bsfp, NULL);

  res = calcSparseTanimotoSml(asfp, bsfp); 
  PG_RETURN_BOOL( res >= getTanimotoLimit() );
}


#ifdef USE_SFP_OBJECTS
PG_FUNCTION_INFO_V1(sfp_dice_sml);
Datum           sfp_dice_sml(PG_FUNCTION_ARGS);
Datum
sfp_dice_sml(PG_FUNCTION_ARGS) {
  MolSparseFingerPrint    asfp,
    bsfp;
  double                  res;

  fcinfo->flinfo->fn_extra = SearchSparseFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(0), 
                                                 NULL, &asfp, NULL);
  fcinfo->flinfo->fn_extra = SearchSparseFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(1), 
                                                 NULL, &bsfp, NULL);

  res = calcSparseDiceSml(asfp, bsfp); 

  PG_RETURN_FLOAT8(res);          
}

PG_FUNCTION_INFO_V1(sfp_dice_sml_op);
Datum           sfp_dice_sml_op(PG_FUNCTION_ARGS);
Datum
sfp_dice_sml_op(PG_FUNCTION_ARGS) {
  MolSparseFingerPrint    asfp,
    bsfp;
  double                  res;

  fcinfo->flinfo->fn_extra = SearchSparseFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(0), 
                                                 NULL, &asfp, NULL);
  fcinfo->flinfo->fn_extra = SearchSparseFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(1), 
                                                 NULL, &bsfp, NULL);

  res = calcSparseDiceSml(asfp, bsfp); 
  PG_RETURN_BOOL( res >= getDiceLimit() );
}
#else
PG_FUNCTION_INFO_V1(sfp_dice_sml);
Datum           sfp_dice_sml(PG_FUNCTION_ARGS);
Datum
sfp_dice_sml(PG_FUNCTION_ARGS) {
  const char *a,*b;
  unsigned int sza,szb;
  double                        res;

  bytea *ba=PG_GETARG_BYTEA_P(0);
  bytea *bb=PG_GETARG_BYTEA_P(1);
  a = VARDATA(ba);
  sza=VARSIZE(ba)-VARHDRSZ;
  b = VARDATA(bb);
  szb=VARSIZE(bb)-VARHDRSZ;

  res = calcSparseStringDiceSml(a,sza,b,szb);

  PG_RETURN_FLOAT8(res);                
}

PG_FUNCTION_INFO_V1(sfp_dice_sml_op);
Datum           sfp_dice_sml_op(PG_FUNCTION_ARGS);
Datum
sfp_dice_sml_op(PG_FUNCTION_ARGS) {
  const char *a,*b;
  unsigned int sza,szb;
  double                        res;

  bytea *ba=PG_GETARG_BYTEA_P(0);
  bytea *bb=PG_GETARG_BYTEA_P(1);
  a = VARDATA(ba);
  sza=VARSIZE(ba)-VARHDRSZ;
  b = VARDATA(bb);
  szb=VARSIZE(bb)-VARHDRSZ;

  res = calcSparseStringDiceSml(a,sza,b,szb);

  PG_RETURN_BOOL(res >= getDiceLimit() );               
}

PG_FUNCTION_INFO_V1(sfp_allvals_gt);
Datum           sfp_allvals_gt(PG_FUNCTION_ARGS);
Datum
sfp_allvals_gt(PG_FUNCTION_ARGS) {
  const char *a;
  unsigned int sza;
  bool res;

  bytea *ba=PG_GETARG_BYTEA_P(0);
  int tgt=PG_GETARG_INT32(1);
  a = VARDATA(ba);
  sza=VARSIZE(ba)-VARHDRSZ;

  res = calcSparseStringAllValsGT(a,sza,tgt);

  PG_RETURN_BOOL(res);
}
PG_FUNCTION_INFO_V1(sfp_allvals_lt);
Datum           sfp_allvals_lt(PG_FUNCTION_ARGS);
Datum
sfp_allvals_lt(PG_FUNCTION_ARGS) {
  const char *a;
  unsigned int sza;
  bool res;

  bytea *ba=PG_GETARG_BYTEA_P(0);
  int tgt=PG_GETARG_INT32(1);
  a = VARDATA(ba);
  sza=VARSIZE(ba)-VARHDRSZ;

  res = calcSparseStringAllValsLT(a,sza,tgt);

  PG_RETURN_BOOL(res);
}
#endif // USE_SFP_OBJECTS



PG_FUNCTION_INFO_V1(sfp_add);
Datum           sfp_add(PG_FUNCTION_ARGS);
Datum
sfp_add(PG_FUNCTION_ARGS) {
  MolSparseFingerPrint    asfp,
    bsfp;
  MolSparseFingerPrint    fp;
  SparseFingerPrint       *sfp;

  fcinfo->flinfo->fn_extra = SearchSparseFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(0), 
                                                 NULL, &asfp, NULL);
  fcinfo->flinfo->fn_extra = SearchSparseFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(1), 
                                                 NULL, &bsfp, NULL);

  fp = addSFP(asfp, bsfp); 
  sfp = deconstructMolSparseFingerPrint(fp);
  freeMolSparseFingerPrint(fp);

  PG_RETURN_SPARSEFINGERPRINT_P(sfp);
}

PG_FUNCTION_INFO_V1(sfp_subtract);
Datum           sfp_subtract(PG_FUNCTION_ARGS);
Datum
sfp_subtract(PG_FUNCTION_ARGS) {
  MolSparseFingerPrint    asfp,
    bsfp;
  MolSparseFingerPrint    fp;
  SparseFingerPrint       *sfp;

  fcinfo->flinfo->fn_extra = SearchSparseFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(0), 
                                                 NULL, &asfp, NULL);
  fcinfo->flinfo->fn_extra = SearchSparseFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(1), 
                                                 NULL, &bsfp, NULL);

  fp = subtractSFP(asfp, bsfp); 
  sfp = deconstructMolSparseFingerPrint(fp);
  freeMolSparseFingerPrint(fp);

  PG_RETURN_SPARSEFINGERPRINT_P(sfp);
}


PG_FUNCTION_INFO_V1(morgan_fp);
Datum       morgan_fp(PG_FUNCTION_ARGS);
Datum
morgan_fp(PG_FUNCTION_ARGS) {
  CROMol  mol;
  MolSparseFingerPrint    fp;
  SparseFingerPrint       *sfp;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);

  fp = makeMorganSFP(mol, PG_GETARG_INT32(1) /* radius */ );
  sfp = deconstructMolSparseFingerPrint(fp);
  freeMolSparseFingerPrint(fp);

  PG_RETURN_SPARSEFINGERPRINT_P(sfp);
}
PG_FUNCTION_INFO_V1(featmorgan_fp);
Datum       featmorgan_fp(PG_FUNCTION_ARGS);
Datum
featmorgan_fp(PG_FUNCTION_ARGS) {
  CROMol  mol;
  MolSparseFingerPrint    fp;
  SparseFingerPrint       *sfp;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);

  fp = makeFeatMorganSFP(mol, PG_GETARG_INT32(1) /* radius */ );
  sfp = deconstructMolSparseFingerPrint(fp);
  freeMolSparseFingerPrint(fp);

  PG_RETURN_SPARSEFINGERPRINT_P(sfp);
}
PG_FUNCTION_INFO_V1(atompair_fp);
Datum       atompair_fp(PG_FUNCTION_ARGS);
Datum
atompair_fp(PG_FUNCTION_ARGS) {
  CROMol  mol;
  MolSparseFingerPrint    fp;
  SparseFingerPrint       *sfp;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);

  fp = makeAtomPairSFP(mol);
  sfp = deconstructMolSparseFingerPrint(fp);
  freeMolSparseFingerPrint(fp);

  PG_RETURN_SPARSEFINGERPRINT_P(sfp);
}
PG_FUNCTION_INFO_V1(torsion_fp);
Datum       torsion_fp(PG_FUNCTION_ARGS);
Datum
torsion_fp(PG_FUNCTION_ARGS) {
  CROMol  mol;
  MolSparseFingerPrint    fp;
  SparseFingerPrint       *sfp;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);

  fp = makeTopologicalTorsionSFP(mol);
  sfp = deconstructMolSparseFingerPrint(fp);
  freeMolSparseFingerPrint(fp);

  PG_RETURN_SPARSEFINGERPRINT_P(sfp);
}

PG_FUNCTION_INFO_V1(reaction_difference_fp);
Datum       reaction_difference_fp(PG_FUNCTION_ARGS);
Datum
reaction_difference_fp(PG_FUNCTION_ARGS) {
  CChemicalReaction  rxn;
  MolSparseFingerPrint    fp;
  SparseFingerPrint       *sfp;

  fcinfo->flinfo->fn_extra = SearchChemReactionCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &rxn, NULL);

  fp = makeReactionDifferenceSFP(rxn, getReactionDifferenceFpSize(),
		  PG_GETARG_INT32(1) /* Fingerprinttype */ );
  sfp = deconstructMolSparseFingerPrint(fp);
  freeMolSparseFingerPrint(fp);

  PG_RETURN_SPARSEFINGERPRINT_P(sfp);
}

