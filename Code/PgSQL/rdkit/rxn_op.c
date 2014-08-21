// $Id$
//
//  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
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
/***************** chem reaction operations ***********************/


#define CHEMREACTDESCR( name, func, ret )                               \
  PG_FUNCTION_INFO_V1(reaction_##name);                                      \
  Datum           reaction_##name(PG_FUNCTION_ARGS);                         \
  Datum                                                                 \
  reaction_##name(PG_FUNCTION_ARGS){                                         \
  CChemicalReaction        rxn;                                         \
  fcinfo->flinfo->fn_extra = SearchChemReactionCache(                   \
                                            fcinfo->flinfo->fn_extra,   \
                                            fcinfo->flinfo->fn_mcxt,    \
                                            PG_GETARG_DATUM(0),         \
                                            NULL, &rxn, NULL);          \
  PG_RETURN_##ret( func(rxn) );                                           \
}
  

CHEMREACTDESCR(numreactants,ChemReactNumReactants,INT32)
CHEMREACTDESCR(numproducts,ChemReactNumProducts,INT32)
CHEMREACTDESCR(numagents,ChemReactNumAgents,INT32)

PG_FUNCTION_INFO_V1(reaction_substruct);
Datum           reaction_substruct(PG_FUNCTION_ARGS);
Datum
reaction_substruct(PG_FUNCTION_ARGS) {
  CChemicalReaction  rxn, rxn2;

  fcinfo->flinfo->fn_extra = SearchChemReactionCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &rxn, NULL);
  fcinfo->flinfo->fn_extra = SearchChemReactionCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(1),
                                            NULL, &rxn2, NULL);

  PG_RETURN_BOOL(ReactionSubstruct(rxn, rxn2));
}

PG_FUNCTION_INFO_V1(reaction_rsubstruct);
Datum           reaction_rsubstruct(PG_FUNCTION_ARGS);
Datum
reaction_rsubstruct(PG_FUNCTION_ARGS) {
  CChemicalReaction rxn, rxn2;

  fcinfo->flinfo->fn_extra = SearchChemReactionCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(0),
                                            NULL, &rxn, NULL);
  fcinfo->flinfo->fn_extra = SearchChemReactionCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(1),
                                            NULL, &rxn2, NULL);

  PG_RETURN_BOOL(ReactionSubstruct(rxn2, rxn));
}

PG_FUNCTION_INFO_V1(reaction_substructFP);
Datum           reaction_substructFP(PG_FUNCTION_ARGS);
Datum
reaction_substructFP(PG_FUNCTION_ARGS) {
  CChemicalReaction rxn, rxn2;

  fcinfo->flinfo->fn_extra = SearchChemReactionCache(
	                                            fcinfo->flinfo->fn_extra,
	                                            fcinfo->flinfo->fn_mcxt,
	                                            PG_GETARG_DATUM(0),
	                                            NULL, &rxn, NULL);
  fcinfo->flinfo->fn_extra = SearchChemReactionCache(
	                                            fcinfo->flinfo->fn_extra,
	                                            fcinfo->flinfo->fn_mcxt,
	                                            PG_GETARG_DATUM(1),
	                                            NULL, &rxn2, NULL);

  PG_RETURN_BOOL(ReactionSubstructFP(rxn, rxn2));
}

PG_FUNCTION_INFO_V1(reaction_rsubstructFP);
Datum           reaction_rsubstructFP(PG_FUNCTION_ARGS);
Datum
reaction_rsubstructFP(PG_FUNCTION_ARGS) {
  CChemicalReaction rxn, rxn2;

  fcinfo->flinfo->fn_extra = SearchChemReactionCache(
	                                            fcinfo->flinfo->fn_extra,
	                                            fcinfo->flinfo->fn_mcxt,
	                                            PG_GETARG_DATUM(0),
	                                            NULL, &rxn, NULL);
  fcinfo->flinfo->fn_extra = SearchChemReactionCache(
	                                            fcinfo->flinfo->fn_extra,
	                                            fcinfo->flinfo->fn_mcxt,
	                                            PG_GETARG_DATUM(1),
	                                            NULL, &rxn2, NULL);

  PG_RETURN_BOOL(ReactionSubstructFP(rxn2, rxn));
}

#define REACTIONCMPFUNC( type, action, ret )                            \
  PG_FUNCTION_INFO_V1(reaction_##type);                                 \
  Datum           reaction_##type(PG_FUNCTION_ARGS);                    \
  Datum                                                                 \
  reaction_##type(PG_FUNCTION_ARGS)                                     \
  {                                                                     \
	CChemicalReaction rxn, rxn2;                                        \
    int res;                                                            \
                                                                        \
    fcinfo->flinfo->fn_extra = SearchChemReactionCache(                 \
                                              fcinfo->flinfo->fn_extra, \
                                              fcinfo->flinfo->fn_mcxt,  \
                                              PG_GETARG_DATUM(0),       \
                                              NULL, &rxn, NULL);        \
    fcinfo->flinfo->fn_extra = SearchChemReactionCache(                 \
                                              fcinfo->flinfo->fn_extra, \
                                              fcinfo->flinfo->fn_mcxt,  \
                                              PG_GETARG_DATUM(1),       \
                                              NULL, &rxn2, NULL);       \
    res = reactioncmp(rxn, rxn2);                                       \
    PG_RETURN_##ret( res action 0 );                                    \
  }                                                                     \
//  /* keep compiler quiet - no extra ; */                                \
//  extern int no_such_variable

REACTIONCMPFUNC(eq, ==, BOOL);
REACTIONCMPFUNC(ne, !=, BOOL);
