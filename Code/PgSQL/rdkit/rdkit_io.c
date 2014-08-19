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
#if PG_VERSION_NUM>=90000
#include "utils/bytea.h"
#endif
#include "utils/builtins.h"

PG_MODULE_MAGIC;

PG_FUNCTION_INFO_V1(mol_in);
Datum           mol_in(PG_FUNCTION_ARGS);
Datum
mol_in(PG_FUNCTION_ARGS) {
  char    *data = PG_GETARG_CSTRING(0);
  CROMol  mol;
  Mol     *res;

  mol = parseMolText(data,false,false);
  if(!mol){
    ereport(ERROR,
            (errcode(ERRCODE_DATA_EXCEPTION),
             errmsg("could not construct molecule")));
  }
  res = deconstructROMol(mol);
  freeCROMol(mol);

  PG_RETURN_MOL_P(res);           
}

PG_FUNCTION_INFO_V1(mol_recv);
Datum           mol_recv(PG_FUNCTION_ARGS);
Datum
mol_recv(PG_FUNCTION_ARGS) {
  bytea    *data = PG_GETARG_BYTEA_P(0);
  int len=VARSIZE(data)-VARHDRSZ;
  CROMol  mol;
  Mol     *res;
  mol = parseMolBlob(VARDATA(data),len);
  res = deconstructROMol(mol);
  freeCROMol(mol);

  PG_FREE_IF_COPY(data, 0);

  PG_RETURN_MOL_P(res);           
}


PG_FUNCTION_INFO_V1(mol_out);
Datum           mol_out(PG_FUNCTION_ARGS);
Datum
mol_out(PG_FUNCTION_ARGS) {
  CROMol  mol;
  char    *str;
  int     len;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);
  str = makeMolText(mol, &len,false);

  PG_RETURN_CSTRING( pnstrdup(str, len) );
}

PG_FUNCTION_INFO_V1(mol_send);
Datum           mol_send(PG_FUNCTION_ARGS);
Datum
mol_send(PG_FUNCTION_ARGS) {
  CROMol  mol;
  bytea    *res;
  char *str;
  int     len;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);
  str = makeMolBlob(mol, &len);
  res=(bytea *)palloc(len+VARHDRSZ);
  SET_VARSIZE(res,len+VARHDRSZ);
  memcpy(VARDATA(res),str,len);
  PG_RETURN_BYTEA_P( res );
}

PG_FUNCTION_INFO_V1(mol_from_ctab);
Datum           mol_from_ctab(PG_FUNCTION_ARGS);
Datum
mol_from_ctab(PG_FUNCTION_ARGS) {
  char    *data = PG_GETARG_CSTRING(0);
  bool keepConformer = PG_GETARG_BOOL(1);
  CROMol  mol;
  Mol     *res;

  mol = parseMolCTAB(data,keepConformer,true);
  if(!mol) PG_RETURN_NULL();
  res = deconstructROMol(mol);
  freeCROMol(mol);

  PG_RETURN_MOL_P(res);           
}

PG_FUNCTION_INFO_V1(mol_from_smarts);
Datum           mol_from_smarts(PG_FUNCTION_ARGS);
Datum
mol_from_smarts(PG_FUNCTION_ARGS) {
  char    *data = PG_GETARG_CSTRING(0);
  CROMol  mol;
  Mol     *res;

  mol = parseMolText(data,true,true);
  if(!mol) PG_RETURN_NULL();
  res = deconstructROMol(mol);
  freeCROMol(mol);

  PG_RETURN_MOL_P(res);           
}

PG_FUNCTION_INFO_V1(mol_from_smiles);
Datum           mol_from_smiles(PG_FUNCTION_ARGS);
Datum
mol_from_smiles(PG_FUNCTION_ARGS) {
  char    *data = PG_GETARG_CSTRING(0);
  CROMol  mol;
  Mol     *res;

  mol = parseMolText(data,false,true);
  if(!mol) PG_RETURN_NULL();
  res = deconstructROMol(mol);
  freeCROMol(mol);

  PG_RETURN_MOL_P(res);           
}

PG_FUNCTION_INFO_V1(mol_to_ctab);
Datum           mol_to_ctab(PG_FUNCTION_ARGS);
Datum
mol_to_ctab(PG_FUNCTION_ARGS) {
  CROMol  mol;
  char    *str;
  int     len;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);

  bool createDepictionIfMissing = PG_GETARG_BOOL(1);
  str = makeCtabText(mol, &len, createDepictionIfMissing);

  PG_RETURN_CSTRING( pnstrdup(str, len) );
}

PG_FUNCTION_INFO_V1(mol_to_smiles);
Datum           mol_to_smiles(PG_FUNCTION_ARGS);
Datum
mol_to_smiles(PG_FUNCTION_ARGS) {
  CROMol  mol;
  char    *str;
  int     len;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);
  str = makeMolText(mol, &len,false);

  PG_RETURN_CSTRING( pnstrdup(str, len) );
}

PG_FUNCTION_INFO_V1(mol_to_smarts);
Datum           mol_to_smarts(PG_FUNCTION_ARGS);
Datum
mol_to_smarts(PG_FUNCTION_ARGS) {
  CROMol  mol;
  char    *str;
  int     len;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);
  str = makeMolText(mol, &len,true);

  PG_RETURN_CSTRING( pnstrdup(str, len) );
}

PG_FUNCTION_INFO_V1(mol_from_pkl);
Datum           mol_from_pkl(PG_FUNCTION_ARGS);
Datum
mol_from_pkl(PG_FUNCTION_ARGS) {
  bytea    *data = PG_GETARG_BYTEA_P(0);
  int len=VARSIZE(data)-VARHDRSZ;
  CROMol  mol;
  Mol     *res;
  mol = parseMolBlob(VARDATA(data),len);
  res = deconstructROMol(mol);
  freeCROMol(mol);

  PG_FREE_IF_COPY(data, 0);

  PG_RETURN_MOL_P(res);           
}

PG_FUNCTION_INFO_V1(mol_to_pkl);
Datum           mol_to_pkl(PG_FUNCTION_ARGS);
Datum
mol_to_pkl(PG_FUNCTION_ARGS) {
  CROMol  mol;
  bytea    *res;
  char *str;
  int     len;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);
  str = makeMolBlob(mol, &len);
  res=(bytea *)palloc(len+VARHDRSZ);
  SET_VARSIZE(res,len+VARHDRSZ);
  memcpy(VARDATA(res),str,len);
  PG_RETURN_BYTEA_P( res );
}

PG_FUNCTION_INFO_V1(qmol_in);
Datum           qmol_in(PG_FUNCTION_ARGS);
Datum
qmol_in(PG_FUNCTION_ARGS) {
  char    *data = PG_GETARG_CSTRING(0);
  CROMol  mol;
  Mol     *res;

  mol = parseMolText(data,true,false);
  if(!mol){
    ereport(ERROR,
            (errcode(ERRCODE_DATA_EXCEPTION),
             errmsg("could not construct molecule")));
  }
  res = deconstructROMol(mol);
  freeCROMol(mol);

  PG_RETURN_MOL_P(res);           
}


PG_FUNCTION_INFO_V1(qmol_out);
Datum           qmol_out(PG_FUNCTION_ARGS);
Datum
qmol_out(PG_FUNCTION_ARGS) {
  CROMol  mol;
  char    *str;
  int     len;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &mol, NULL);
  str = makeMolText(mol, &len,true);

  PG_RETURN_CSTRING( pnstrdup(str, len) );
}


PG_FUNCTION_INFO_V1(bfp_in);
Datum           bfp_in(PG_FUNCTION_ARGS);
Datum
bfp_in(PG_FUNCTION_ARGS) {
  MolBitmapFingerPrint    fp;
  BitmapFingerPrint       *b = DatumGetBitmapFingerPrintP(DirectFunctionCall1(
                                                                              byteain,
                                                                              PG_GETARG_DATUM(0)
                                                                              ));

  /* check correctness */
  fp = constructMolBitmapFingerPrint(b);
  freeMolBitmapFingerPrint(fp);

  PG_RETURN_BITMAPFINGERPRINT_P(b);
}

PG_FUNCTION_INFO_V1(bfp_out);
Datum           bfp_out(PG_FUNCTION_ARGS);
Datum
bfp_out(PG_FUNCTION_ARGS) {
  PG_RETURN_DATUM( DirectFunctionCall1( byteaout, PG_GETARG_DATUM(0) ) );
}


PG_FUNCTION_INFO_V1(bfp_from_binary_text);
Datum           bfp_from_binary_text(PG_FUNCTION_ARGS);
Datum
bfp_from_binary_text(PG_FUNCTION_ARGS) {
  MolBitmapFingerPrint    fp;
  BitmapFingerPrint       *b =PG_GETARG_BYTEA_P(0);

  fp = constructMolBitmapFingerPrint(b);
  freeMolBitmapFingerPrint(fp);

  PG_RETURN_BITMAPFINGERPRINT_P(b);
}

PG_FUNCTION_INFO_V1(bfp_to_binary_text);
Datum           bfp_to_binary_text(PG_FUNCTION_ARGS);
Datum
bfp_to_binary_text(PG_FUNCTION_ARGS) {
  MolBitmapFingerPrint    abfp;
  fcinfo->flinfo->fn_extra = SearchBitmapFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(0), 
                                                 NULL, &abfp, NULL);
  
  BitmapFingerPrint *b=deconstructMolBitmapFingerPrint(abfp);
  PG_RETURN_BYTEA_P( b );
}


PG_FUNCTION_INFO_V1(sfp_in);
Datum           sfp_in(PG_FUNCTION_ARGS);
Datum
sfp_in(PG_FUNCTION_ARGS) {
  MolSparseFingerPrint    fp;
  SparseFingerPrint       *b = DatumGetSparseFingerPrintP(DirectFunctionCall1(
                                                                              byteain,
                                                                              PG_GETARG_DATUM(0)
                                                                              ));

  /* check correctness */
  fp = constructMolSparseFingerPrint(b);
  freeMolSparseFingerPrint(fp);

  PG_RETURN_SPARSEFINGERPRINT_P(b);
}

PG_FUNCTION_INFO_V1(sfp_out);
Datum           sfp_out(PG_FUNCTION_ARGS);
Datum
sfp_out(PG_FUNCTION_ARGS) {
  PG_RETURN_DATUM( DirectFunctionCall1( byteaout, PG_GETARG_DATUM(0) ) );
}

PG_FUNCTION_INFO_V1(rdkit_version);
Datum           rdkit_version(PG_FUNCTION_ARGS);
Datum
rdkit_version(PG_FUNCTION_ARGS) {
  char    *ver = "" RDKITVER;
  char    buf[1024];
  Assert(strlen(ver) == 6);
  snprintf(buf, sizeof(buf), "%d.%d.%d", 
           atoi( pnstrdup(ver, 2) ),  
           atoi( pnstrdup(ver + 2 , 2) ),  
           atoi( pnstrdup(ver + 4, 2) ));

  PG_RETURN_TEXT_P(cstring_to_text(buf));
}


/* chemical reactions */

PG_FUNCTION_INFO_V1(reaction_in);
Datum           reaction_in(PG_FUNCTION_ARGS);
Datum
reaction_in(PG_FUNCTION_ARGS) {
  char    *data = PG_GETARG_CSTRING(0);
  CChemicalReaction  rxn;
  ChemReactionBA     *rxnBA;

  rxn = parseChemReactText(data,false,false);

  if(!rxn){
    ereport(ERROR,
            (errcode(ERRCODE_DATA_EXCEPTION),
             errmsg("could not construct chemical reaction")));
  }
  rxnBA = deconstructChemReact(rxn);
  freeChemReaction(rxn);

  PG_RETURN_CHEMREACTION_P(rxnBA);           
}

PG_FUNCTION_INFO_V1(reaction_recv);
Datum           reaction_recv(PG_FUNCTION_ARGS);
Datum
reaction_recv(PG_FUNCTION_ARGS) {
  bytea    *data = PG_GETARG_BYTEA_P(0);
  int len=VARSIZE(data)-VARHDRSZ;
  CChemicalReaction  rxn;
  ChemReactionBA     *rxnBA;

  rxn = parseChemReactBlob(VARDATA(data),len);

  rxnBA = deconstructChemReact(rxn);
  freeChemReaction(rxn);

  PG_FREE_IF_COPY(data, 0);

  PG_RETURN_CHEMREACTION_P(rxnBA);           
}


PG_FUNCTION_INFO_V1(reaction_out);
Datum           reaction_out(PG_FUNCTION_ARGS);
Datum
reaction_out(PG_FUNCTION_ARGS) {
	CChemicalReaction  rxn;
  char    *str;
  int     len;

  fcinfo->flinfo->fn_extra = SearchChemReactionCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &rxn, NULL);
  str = makeChemReactText(rxn, &len,false);

  PG_RETURN_CSTRING( pnstrdup(str, len) );
}

PG_FUNCTION_INFO_V1(reaction_send);
Datum           reaction_send(PG_FUNCTION_ARGS);
Datum
reaction_send(PG_FUNCTION_ARGS) {
	CChemicalReaction  rxn;
  bytea    *res;
  char *str;
  int     len;

  fcinfo->flinfo->fn_extra = SearchChemReactionCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &rxn, NULL);
  str = makeChemReactBlob(rxn, &len);

  res=(bytea *)palloc(len+VARHDRSZ);
  SET_VARSIZE(res,len+VARHDRSZ);
  memcpy(VARDATA(res),str,len);
  PG_RETURN_BYTEA_P( res );
}

PG_FUNCTION_INFO_V1(reaction_from_ctab);
Datum           reaction_from_ctab(PG_FUNCTION_ARGS);
Datum
reaction_from_ctab(PG_FUNCTION_ARGS) {
  char    *data = PG_GETARG_CSTRING(0);
  CChemicalReaction  rxn;
  ChemReactionBA     *rxnBA;

  rxn = parseChemReactCTAB(data,true);
  if(!rxn) PG_RETURN_NULL();
  rxnBA = deconstructChemReact(rxn);
  freeChemReaction(rxn);

  PG_RETURN_CHEMREACTION_P(rxnBA);
}

PG_FUNCTION_INFO_V1(reaction_from_smarts);
Datum           reaction_from_smarts(PG_FUNCTION_ARGS);
Datum
reaction_from_smarts(PG_FUNCTION_ARGS) {
  char    *data = PG_GETARG_CSTRING(0);
  CChemicalReaction  rxn;
  ChemReactionBA     *rxnBA;

  rxn = parseChemReactText(data,true,true);
  if(!rxn) PG_RETURN_NULL();
  rxnBA = deconstructChemReact(rxn);
  freeChemReaction(rxn);

  PG_RETURN_CHEMREACTION_P(rxnBA);
}

PG_FUNCTION_INFO_V1(reaction_from_smiles);
Datum           reaction_from_smiles(PG_FUNCTION_ARGS);
Datum
reaction_from_smiles(PG_FUNCTION_ARGS) {
  char    *data = PG_GETARG_CSTRING(0);
  CChemicalReaction  rxn;
  ChemReactionBA     *rxnBA;

  rxn = parseChemReactText(data,false,true);
  if(!rxn) PG_RETURN_NULL();
  rxnBA = deconstructChemReact(rxn);
  freeChemReaction(rxn);

  PG_RETURN_CHEMREACTION_P(rxnBA);
}

PG_FUNCTION_INFO_V1(reaction_to_ctab);
Datum           reaction_to_ctab(PG_FUNCTION_ARGS);
Datum
reaction_to_ctab(PG_FUNCTION_ARGS) {
  CChemicalReaction  rxn;
  char    *str;
  int     len;

  fcinfo->flinfo->fn_extra = SearchChemReactionCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &rxn, NULL);

  str = makeCTABChemReact(rxn, &len);

  PG_RETURN_CSTRING( pnstrdup(str, len) );
}

PG_FUNCTION_INFO_V1(reaction_to_smiles);
Datum           reaction_to_smiles(PG_FUNCTION_ARGS);
Datum
reaction_to_smiles(PG_FUNCTION_ARGS) {
  CChemicalReaction  rxn;
  char    *str;
  int     len;

  fcinfo->flinfo->fn_extra = SearchChemReactionCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &rxn, NULL);
  str = makeChemReactText(rxn, &len,false);

  PG_RETURN_CSTRING( pnstrdup(str, len) );
}


PG_FUNCTION_INFO_V1(reaction_to_smarts);
Datum           reaction_to_smarts(PG_FUNCTION_ARGS);
Datum
reaction_to_smarts(PG_FUNCTION_ARGS) {
  CChemicalReaction  rxn;
  char    *str;
  int     len;

  fcinfo->flinfo->fn_extra = SearchChemReactionCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &rxn, NULL);
  str = makeChemReactText(rxn, &len,true);

  PG_RETURN_CSTRING( pnstrdup(str, len) );
}


