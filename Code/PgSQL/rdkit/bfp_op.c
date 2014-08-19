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
bfpcmp(BitmapFingerPrint *a, BitmapFingerPrint *b) {
  int res;

  res = memcmp(VARDATA(a), VARDATA(b), Min(VARSIZE(a), VARSIZE(b)) - VARHDRSZ);
  if ( res )
    return res;

  if (VARSIZE(a) == VARSIZE(b))
    return 0;
  return (VARSIZE(a) > VARSIZE(b)) ? 1 : -1; 
}


#define bfpCMPFUNC( type, action, ret )                                 \
  PG_FUNCTION_INFO_V1(bfp_##type);                                      \
  Datum           bfp_##type(PG_FUNCTION_ARGS);                         \
  Datum                                                                 \
  bfp_##type(PG_FUNCTION_ARGS)                                          \
  {                                                                     \
    BitmapFingerPrint    *a, *b;                                        \
    int             res;                                                \
                                                                        \
    fcinfo->flinfo->fn_extra = SearchBitmapFPCache(                     \
                                                   fcinfo->flinfo->fn_extra, \
                                                   fcinfo->flinfo->fn_mcxt, \
                                                   PG_GETARG_DATUM(0),  \
                                                   &a, NULL, NULL);     \
    fcinfo->flinfo->fn_extra = SearchBitmapFPCache(                     \
                                                   fcinfo->flinfo->fn_extra, \
                                                   fcinfo->flinfo->fn_mcxt, \
                                                   PG_GETARG_DATUM(1),  \
                                                   &b, NULL, NULL);     \
    res = bfpcmp(a, b);                                                 \
    PG_RETURN_##ret( res action 0 );                                    \
  }                                                                     \
  /* keep compiler quiet - no extra ; */                                \
  extern int no_such_variable

bfpCMPFUNC(lt, <, BOOL);
bfpCMPFUNC(le, <=, BOOL);
bfpCMPFUNC(eq, ==, BOOL);
bfpCMPFUNC(ge, >=, BOOL);
bfpCMPFUNC(gt, >, BOOL);
bfpCMPFUNC(ne, !=, BOOL);
bfpCMPFUNC(cmp, +, INT32);

PG_FUNCTION_INFO_V1(bfp_tanimoto_sml);
Datum           bfp_tanimoto_sml(PG_FUNCTION_ARGS);
Datum
bfp_tanimoto_sml(PG_FUNCTION_ARGS) {
  MolBitmapFingerPrint    abfp,
    bbfp;
  double                  res;

  fcinfo->flinfo->fn_extra = SearchBitmapFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(0), 
                                                 NULL, &abfp, NULL);
  fcinfo->flinfo->fn_extra = SearchBitmapFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(1), 
                                                 NULL, &bbfp, NULL);

  res = calcBitmapTanimotoSml(abfp, bbfp); 

  PG_RETURN_FLOAT8(res);          
}
PG_FUNCTION_INFO_V1(bfp_tversky_sml);
Datum           bfp_tversky_sml(PG_FUNCTION_ARGS);
Datum
bfp_tversky_sml(PG_FUNCTION_ARGS) {
  MolBitmapFingerPrint    abfp,
    bbfp;
  double                  res;

  fcinfo->flinfo->fn_extra = SearchBitmapFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(0), 
                                                 NULL, &abfp, NULL);
  fcinfo->flinfo->fn_extra = SearchBitmapFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(1), 
                                                 NULL, &bbfp, NULL);

  res = calcBitmapTverskySml(abfp, bbfp, PG_GETARG_FLOAT4(2),PG_GETARG_FLOAT4(3) ); 

  PG_RETURN_FLOAT8(res);          
}
#if PG_VERSION_NUM >= 90100
PG_FUNCTION_INFO_V1(bfp_tanimoto_dist);
Datum           bfp_tanimoto_dist(PG_FUNCTION_ARGS);
Datum
bfp_tanimoto_dist(PG_FUNCTION_ARGS) {
  MolBitmapFingerPrint    abfp,
    bbfp;
  double                  res;

  fcinfo->flinfo->fn_extra = SearchBitmapFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(0),
                                                 NULL, &abfp, NULL);
  fcinfo->flinfo->fn_extra = SearchBitmapFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(1),
                                                 NULL, &bbfp, NULL);

  res = 1.0 - calcBitmapTanimotoSml(abfp, bbfp);

  PG_RETURN_FLOAT8(res);
}
#endif

PG_FUNCTION_INFO_V1(bfp_tanimoto_sml_op);
Datum           bfp_tanimoto_sml_op(PG_FUNCTION_ARGS);
Datum
bfp_tanimoto_sml_op(PG_FUNCTION_ARGS) {
  MolBitmapFingerPrint    abfp,
    bbfp;
  double                  res;

  fcinfo->flinfo->fn_extra = SearchBitmapFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(0), 
                                                 NULL, &abfp, NULL);
  fcinfo->flinfo->fn_extra = SearchBitmapFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(1), 
                                                 NULL, &bbfp, NULL);

  res = calcBitmapTanimotoSml(abfp, bbfp); 
  PG_RETURN_BOOL( res >= getTanimotoLimit() );
}

PG_FUNCTION_INFO_V1(bfp_dice_sml);
Datum           bfp_dice_sml(PG_FUNCTION_ARGS);
Datum
bfp_dice_sml(PG_FUNCTION_ARGS) {
  MolBitmapFingerPrint    abfp,
    bbfp;
  double                  res;

  fcinfo->flinfo->fn_extra = SearchBitmapFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(0), 
                                                 NULL, &abfp, NULL);
  fcinfo->flinfo->fn_extra = SearchBitmapFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(1), 
                                                 NULL, &bbfp, NULL);

  res = calcBitmapDiceSml(abfp, bbfp); 

  PG_RETURN_FLOAT8(res);          
}

#if PG_VERSION_NUM >= 90100
PG_FUNCTION_INFO_V1(bfp_dice_dist);
Datum           bfp_dice_dist(PG_FUNCTION_ARGS);
Datum
bfp_dice_dist(PG_FUNCTION_ARGS) {
  MolBitmapFingerPrint    abfp,
    bbfp;
  double                  res;

  fcinfo->flinfo->fn_extra = SearchBitmapFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(0),
                                                 NULL, &abfp, NULL);
  fcinfo->flinfo->fn_extra = SearchBitmapFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(1),
                                                 NULL, &bbfp, NULL);

  res = 1.0 - calcBitmapDiceSml(abfp, bbfp);

  PG_RETURN_FLOAT8(res);
}
#endif

PG_FUNCTION_INFO_V1(bfp_dice_sml_op);
Datum           bfp_dice_sml_op(PG_FUNCTION_ARGS);
Datum
bfp_dice_sml_op(PG_FUNCTION_ARGS) {
  MolBitmapFingerPrint    abfp,
    bbfp;
  double                  res;

  fcinfo->flinfo->fn_extra = SearchBitmapFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(0), 
                                                 NULL, &abfp, NULL);
  fcinfo->flinfo->fn_extra = SearchBitmapFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(1), 
                                                 NULL, &bbfp, NULL);

  res = calcBitmapDiceSml(abfp, bbfp); 
  PG_RETURN_BOOL( res >= getDiceLimit() );
}

PG_FUNCTION_INFO_V1(bfp_size);
Datum           bfp_size(PG_FUNCTION_ARGS);
Datum
bfp_size(PG_FUNCTION_ARGS) {
  MolBitmapFingerPrint    bfp;

  fcinfo->flinfo->fn_extra = SearchBitmapFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(0), 
                                                 NULL, &bfp, NULL);

  PG_RETURN_INT32(MolBitmapFingerPrintSize(bfp));
}

PG_FUNCTION_INFO_V1(layered_fp);
Datum       layered_fp(PG_FUNCTION_ARGS);
Datum
layered_fp(PG_FUNCTION_ARGS) {
  CROMol  mol;
  MolBitmapFingerPrint    fp;
  BitmapFingerPrint               *sfp;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);

  fp = makeLayeredBFP(mol);
  sfp = deconstructMolBitmapFingerPrint(fp);
  freeMolBitmapFingerPrint(fp);

  PG_RETURN_BITMAPFINGERPRINT_P(sfp);
}

PG_FUNCTION_INFO_V1(rdkit_fp);
Datum       rdkit_fp(PG_FUNCTION_ARGS);
Datum
rdkit_fp(PG_FUNCTION_ARGS) {
  CROMol  mol;
  MolBitmapFingerPrint    fp;
  BitmapFingerPrint               *sfp;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);

  fp = makeRDKitBFP(mol);
  sfp = deconstructMolBitmapFingerPrint(fp);
  freeMolBitmapFingerPrint(fp);

  PG_RETURN_BITMAPFINGERPRINT_P(sfp);
}


PG_FUNCTION_INFO_V1(morganbv_fp);
Datum       morganbv_fp(PG_FUNCTION_ARGS);
Datum
morganbv_fp(PG_FUNCTION_ARGS) {
  CROMol  mol;
  MolBitmapFingerPrint    fp;
  BitmapFingerPrint               *sfp;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);

  fp = makeMorganBFP(mol, PG_GETARG_INT32(1) /* radius */ );
  sfp = deconstructMolBitmapFingerPrint(fp);
  freeMolBitmapFingerPrint(fp);

  PG_RETURN_BITMAPFINGERPRINT_P(sfp);
}
PG_FUNCTION_INFO_V1(featmorganbv_fp);
Datum       featmorganbv_fp(PG_FUNCTION_ARGS);
Datum
featmorganbv_fp(PG_FUNCTION_ARGS) {
  CROMol  mol;
  MolBitmapFingerPrint    fp;
  BitmapFingerPrint               *sfp;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);

  fp = makeFeatMorganBFP(mol, PG_GETARG_INT32(1) /* radius */ );
  sfp = deconstructMolBitmapFingerPrint(fp);
  freeMolBitmapFingerPrint(fp);

  PG_RETURN_BITMAPFINGERPRINT_P(sfp);
}
PG_FUNCTION_INFO_V1(atompairbv_fp);
Datum       atompairbv_fp(PG_FUNCTION_ARGS);
Datum
atompairbv_fp(PG_FUNCTION_ARGS) {
  CROMol  mol;
  MolBitmapFingerPrint    fp;
  BitmapFingerPrint               *sfp;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);

  fp = makeAtomPairBFP(mol);
  sfp = deconstructMolBitmapFingerPrint(fp);
  freeMolBitmapFingerPrint(fp);

  PG_RETURN_BITMAPFINGERPRINT_P(sfp);
}
PG_FUNCTION_INFO_V1(torsionbv_fp);
Datum       torsionbv_fp(PG_FUNCTION_ARGS);
Datum
torsionbv_fp(PG_FUNCTION_ARGS) {
  CROMol  mol;
  MolBitmapFingerPrint    fp;
  BitmapFingerPrint               *sfp;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);

  fp = makeTopologicalTorsionBFP(mol);
  sfp = deconstructMolBitmapFingerPrint(fp);
  freeMolBitmapFingerPrint(fp);

  PG_RETURN_BITMAPFINGERPRINT_P(sfp);
}
PG_FUNCTION_INFO_V1(maccs_fp);
Datum       maccs_fp(PG_FUNCTION_ARGS);
Datum
maccs_fp(PG_FUNCTION_ARGS) {
  CROMol  mol;
  MolBitmapFingerPrint    fp;
  BitmapFingerPrint               *sfp;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);

  fp = makeMACCSBFP(mol);
  sfp = deconstructMolBitmapFingerPrint(fp);
  freeMolBitmapFingerPrint(fp);

  PG_RETURN_BITMAPFINGERPRINT_P(sfp);
}

PG_FUNCTION_INFO_V1(avalon_fp);
Datum       avalon_fp(PG_FUNCTION_ARGS);
Datum
avalon_fp(PG_FUNCTION_ARGS) {
  CROMol  mol;
  MolBitmapFingerPrint    fp;
  BitmapFingerPrint               *sfp;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);

  fp = makeAvalonBFP(mol,
                     PG_GETARG_BOOL(1), /* isQuery */
                     PG_GETARG_UINT32(2) /* flags */
                     );
  sfp = deconstructMolBitmapFingerPrint(fp);
  freeMolBitmapFingerPrint(fp);

  PG_RETURN_BITMAPFINGERPRINT_P(sfp);
}

PG_FUNCTION_INFO_V1(reaction_structural_bfp);
Datum       reaction_structural_bfp(PG_FUNCTION_ARGS);
Datum
reaction_structural_bfp(PG_FUNCTION_ARGS) {
  CChemicalReaction rxn;
  MolBitmapFingerPrint fp;
  BitmapFingerPrint *sfp;

  fcinfo->flinfo->fn_extra = SearchChemReactionCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &rxn, NULL);

  fp = makeReactionBFP(rxn, getReactionSubstructFpSize(), PG_GETARG_INT32(1) /* fpType */);
  sfp = deconstructMolBitmapFingerPrint(fp);
  freeMolBitmapFingerPrint(fp);

  PG_RETURN_BITMAPFINGERPRINT_P(sfp);
}
