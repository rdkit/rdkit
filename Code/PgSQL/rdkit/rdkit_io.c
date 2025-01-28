//
//  Copyright (c) 2010-2021, Novartis Institutes for BioMedical Research Inc.
//    and other RDKit contributors
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
//       products derived from this software without specific prior written
//       permission.
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
#include <postgres.h>
#include <fmgr.h>
#include <utils/bytea.h>
#include <utils/builtins.h>

#include "rdkit.h"
#include "cache.h"

PG_MODULE_MAGIC;

PGDLLEXPORT Datum mol_in(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_in);
Datum mol_in(PG_FUNCTION_ARGS) {
  char *data = PG_GETARG_CSTRING(0);
  CROMol mol;
  Mol *res;

  mol = parseMolText(data, false, false, false, false);
  if (!mol) {
    ereport(ERROR, (errcode(ERRCODE_DATA_EXCEPTION),
                    errmsg("could not construct molecule")));
  }
  res = deconstructROMolWithQueryProperties(mol);
  freeCROMol(mol);

  PG_RETURN_MOL_P(res);
}

PGDLLEXPORT Datum mol_recv(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_recv);
Datum mol_recv(PG_FUNCTION_ARGS) {
  bytea *data = PG_GETARG_BYTEA_P(0);
  int len = VARSIZE(data) - VARHDRSZ;
  CROMol mol;
  Mol *res;
  mol = parseMolBlob(VARDATA(data), len);
  res = deconstructROMolWithQueryProperties(mol);
  freeCROMol(mol);

  PG_FREE_IF_COPY(data, 0);

  PG_RETURN_MOL_P(res);
}

PGDLLEXPORT Datum mol_out(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_out);
Datum mol_out(PG_FUNCTION_ARGS) {
  CROMol mol;
  char *str;
  int len;

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);
  str = makeMolText(mol, &len, false, true, true);

  PG_RETURN_CSTRING(pnstrdup(str, len));
}

PGDLLEXPORT Datum mol_send(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_send);
Datum mol_send(PG_FUNCTION_ARGS) {
  CROMol mol;
  bytea *res;
  char *str;
  int len;

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);
  str = makeMolBlob(mol, &len);
  res = (bytea *)palloc(len + VARHDRSZ);
  SET_VARSIZE(res, len + VARHDRSZ);
  memcpy(VARDATA(res), str, len);
  PG_RETURN_BYTEA_P(res);
}

PGDLLEXPORT Datum mol_from_ctab(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_from_ctab);
Datum mol_from_ctab(PG_FUNCTION_ARGS) {
  char *data = PG_GETARG_CSTRING(0);
  bool keepConformer = PG_GETARG_BOOL(1);
  CROMol mol;
  Mol *res;

  mol = parseMolCTAB(data, keepConformer, true, false);
  if (!mol) {
    PG_RETURN_NULL();
  }
  res = deconstructROMol(mol);
  freeCROMol(mol);

  PG_RETURN_MOL_P(res);
}

PGDLLEXPORT Datum qmol_from_ctab(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(qmol_from_ctab);
Datum qmol_from_ctab(PG_FUNCTION_ARGS) {
  char *data = PG_GETARG_CSTRING(0);
  bool keepConformer = PG_GETARG_BOOL(1);
  CROMol mol;
  Mol *res;

  mol = parseMolCTAB(data, keepConformer, true, true);
  if (!mol) {
    PG_RETURN_NULL();
  }
  res = deconstructROMol(mol);
  freeCROMol(mol);

  PG_RETURN_MOL_P(res);
}

PGDLLEXPORT Datum mol_from_smarts(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_from_smarts);
Datum mol_from_smarts(PG_FUNCTION_ARGS) {
  char *data = PG_GETARG_CSTRING(0);
  CROMol mol;
  Mol *res;

  mol = parseMolText(data, true, true, false, false);
  if (!mol) {
    PG_RETURN_NULL();
  }
  res = deconstructROMol(mol);
  freeCROMol(mol);

  PG_RETURN_MOL_P(res);
}

PGDLLEXPORT Datum mol_from_smiles(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_from_smiles);
Datum mol_from_smiles(PG_FUNCTION_ARGS) {
  char *data = PG_GETARG_CSTRING(0);
  CROMol mol;
  Mol *res;

  mol = parseMolText(data, false, true, false, true);
  if (!mol) {
    PG_RETURN_NULL();
  }
  res = deconstructROMolWithQueryProperties(mol);
  freeCROMol(mol);

  PG_RETURN_MOL_P(res);
}

PGDLLEXPORT Datum qmol_from_smiles(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(qmol_from_smiles);
Datum qmol_from_smiles(PG_FUNCTION_ARGS) {
  char *data = PG_GETARG_CSTRING(0);
  CROMol mol;
  Mol *res;

  mol = parseMolText(data, false, true, true, false);
  if (!mol) {
    PG_RETURN_NULL();
  }
  res = deconstructROMol(mol);
  freeCROMol(mol);

  PG_RETURN_MOL_P(res);
}

PGDLLEXPORT Datum mol_to_ctab(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_to_ctab);
Datum mol_to_ctab(PG_FUNCTION_ARGS) {
  CROMol mol;
  char *str;
  int len;

  bool createDepictionIfMissing = PG_GETARG_BOOL(1);
  bool usev3000 = PG_GETARG_BOOL(2);

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);

  str = makeCtabText(mol, &len, createDepictionIfMissing, usev3000);

  PG_RETURN_CSTRING(pnstrdup(str, len));
}

PGDLLEXPORT Datum mol_to_v3kctab(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_to_v3kctab);
Datum mol_to_v3kctab(PG_FUNCTION_ARGS) {
  CROMol mol;
  char *str;
  int len;

  bool createDepictionIfMissing = PG_GETARG_BOOL(1);
  bool usev3000 = 1;

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);

  str = makeCtabText(mol, &len, createDepictionIfMissing, usev3000);

  PG_RETURN_CSTRING(pnstrdup(str, len));
}

PGDLLEXPORT Datum mol_to_smiles(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_to_smiles);
Datum mol_to_smiles(PG_FUNCTION_ARGS) {
  CROMol mol;
  char *str;
  int len;

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);
  bool doIsomeric = PG_GETARG_BOOL(1);
  str = makeMolText(mol, &len, false, false, doIsomeric);

  PG_RETURN_CSTRING(pnstrdup(str, len));
}

PGDLLEXPORT Datum mol_to_cxsmiles(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_to_cxsmiles);
Datum mol_to_cxsmiles(PG_FUNCTION_ARGS) {
  CROMol mol;
  char *str;
  int len;

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);
  bool doIsomeric = PG_GETARG_BOOL(1);
  str = makeMolText(mol, &len, false, true, doIsomeric);

  PG_RETURN_CSTRING(pnstrdup(str, len));
}

PGDLLEXPORT Datum mol_to_smarts(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_to_smarts);
Datum mol_to_smarts(PG_FUNCTION_ARGS) {
  CROMol mol;
  char *str;
  int len;

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);
  bool dummy = false; // arg is ignored by makeMolText for smarts output
  str = makeMolText(mol, &len, true, false, dummy);

  PG_RETURN_CSTRING(pnstrdup(str, len));
}

PGDLLEXPORT Datum mol_to_cxsmarts(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_to_cxsmarts);
Datum mol_to_cxsmarts(PG_FUNCTION_ARGS) {
  CROMol mol;
  char *str;
  int len;

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);
  bool dummy = true; // arg is ignored by makeMolText for cxsmarts output
  str = makeMolText(mol, &len, true, true, dummy);

  PG_RETURN_CSTRING(pnstrdup(str, len));
}

PGDLLEXPORT Datum mol_from_pkl(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_from_pkl);
Datum mol_from_pkl(PG_FUNCTION_ARGS) {
  bytea *data = PG_GETARG_BYTEA_P(0);
  int len = VARSIZE(data) - VARHDRSZ;
  CROMol mol;
  Mol *res;
  mol = parseMolBlob(VARDATA(data), len);
  res = deconstructROMol(mol);
  freeCROMol(mol);

  PG_FREE_IF_COPY(data, 0);

  PG_RETURN_MOL_P(res);
}

PGDLLEXPORT Datum mol_to_pkl(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_to_pkl);
Datum mol_to_pkl(PG_FUNCTION_ARGS) {
  CROMol mol;
  bytea *res;
  char *str;
  int len;

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);
  str = makeMolBlob(mol, &len);
  res = (bytea *)palloc(len + VARHDRSZ);
  SET_VARSIZE(res, len + VARHDRSZ);
  memcpy(VARDATA(res), str, len);
  PG_RETURN_BYTEA_P(res);
}

PGDLLEXPORT Datum mol_to_json(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_to_json);
Datum mol_to_json(PG_FUNCTION_ARGS) {
  CROMol mol;
  const char *str;
  char *res;

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);
  str = makeMolJSON(mol);
  if (*str == 0) {
    free((void *)str);
    PG_RETURN_NULL();
  } else {
    res = pnstrdup(str, strlen(str));
    free((void *)str);
    PG_RETURN_CSTRING(res);
  }
}

PGDLLEXPORT Datum mol_from_json(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_from_json);
Datum mol_from_json(PG_FUNCTION_ARGS) {
  char *data = PG_GETARG_CSTRING(0);
  CROMol mol;
  Mol *res;

  mol = parseMolJSON(data, true);
  if (!mol) {
    PG_RETURN_NULL();
  }
  res = deconstructROMol(mol);
  freeCROMol(mol);

  PG_RETURN_MOL_P(res);
}

PGDLLEXPORT Datum qmol_in(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(qmol_in);
Datum qmol_in(PG_FUNCTION_ARGS) {
  char *data = PG_GETARG_CSTRING(0);
  CROMol mol;
  Mol *res;

  mol = parseMolText(data, true, false, false, false);
  if (!mol) {
    ereport(ERROR, (errcode(ERRCODE_DATA_EXCEPTION),
                    errmsg("could not construct molecule")));
  }
  res = deconstructROMol(mol);
  freeCROMol(mol);

  PG_RETURN_MOL_P(res);
}

PGDLLEXPORT Datum qmol_out(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(qmol_out);
Datum qmol_out(PG_FUNCTION_ARGS) {
  CROMol mol;
  char *str;
  int len;

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);
  bool dummy = false; // arg is ignored by makeMolText for smarts output
  str = makeMolText(mol, &len, true, false, dummy);

  PG_RETURN_CSTRING(pnstrdup(str, len));
}


/* xqmols */
PGDLLEXPORT Datum xqmol_in(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(xqmol_in);
Datum xqmol_in(PG_FUNCTION_ARGS) {
  char *data = PG_GETARG_CSTRING(0);
  CXQMol mol;
  XQMol *res;

  mol = parseXQMolText(data);
  if (!mol) {
    ereport(ERROR, (errcode(ERRCODE_DATA_EXCEPTION),
                    errmsg("could not construct extended query molecule")));
  }
  res = deconstructXQMol(mol);
  freeCXQMol(mol);

  PG_RETURN_XQMOL_P(res);
}

PGDLLEXPORT Datum xqmol_recv(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(xqmol_recv);
Datum xqmol_recv(PG_FUNCTION_ARGS) {
  bytea *data = PG_GETARG_BYTEA_P(0);
  int len = VARSIZE(data) - VARHDRSZ;
  CXQMol mol;
  XQMol *res;
  mol = parseXQMolBlob(VARDATA(data), len);
  res = deconstructXQMol(mol);
  freeCXQMol(mol);

  PG_FREE_IF_COPY(data, 0);

  PG_RETURN_XQMOL_P(res);
}

PGDLLEXPORT Datum xqmol_out(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(xqmol_out);
Datum xqmol_out(PG_FUNCTION_ARGS) {
  CXQMol mol;
  char *str;
  int len;

  fcinfo->flinfo->fn_extra =
      searchXQMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);
  str = makeXQMolText(mol, &len);

  PG_RETURN_CSTRING(pnstrdup(str, len));
}

PGDLLEXPORT Datum xqmol_send(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(xqmol_send);
Datum xqmol_send(PG_FUNCTION_ARGS) {
  CXQMol mol;
  bytea *res;
  char *str;
  int len;

  fcinfo->flinfo->fn_extra =
      searchXQMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);
  str = makeXQMolBlob(mol, &len);
  res = (bytea *)palloc(len + VARHDRSZ);
  SET_VARSIZE(res, len + VARHDRSZ);
  memcpy(VARDATA(res), str, len);
  PG_RETURN_BYTEA_P(res);
}



PGDLLEXPORT Datum bfp_in(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(bfp_in);
Datum bfp_in(PG_FUNCTION_ARGS) {
  CBfp fp;
  Bfp *b = DatumGetBfpP(DirectFunctionCall1(byteain, PG_GETARG_DATUM(0)));

  /* check correctness */
  fp = constructCBfp(b);
  freeCBfp(fp);

  PG_RETURN_BFP_P(b);
}

PGDLLEXPORT Datum bfp_out(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(bfp_out);
Datum bfp_out(PG_FUNCTION_ARGS) {
  PG_RETURN_DATUM(DirectFunctionCall1(byteaout, PG_GETARG_DATUM(0)));
}

PGDLLEXPORT Datum bfp_from_binary_text(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(bfp_from_binary_text);
Datum bfp_from_binary_text(PG_FUNCTION_ARGS) {
  CBfp fp;
  Bfp *b = PG_GETARG_BYTEA_P(0);

  fp = constructCBfp(b);
  freeCBfp(fp);

  PG_RETURN_BFP_P(b);
}

PGDLLEXPORT Datum bfp_to_binary_text(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(bfp_to_binary_text);
Datum bfp_to_binary_text(PG_FUNCTION_ARGS) {
  CBfp abfp;
  fcinfo->flinfo->fn_extra =
      searchBfpCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &abfp, NULL);

  PG_RETURN_BYTEA_P(deconstructCBfp(abfp));
}

PGDLLEXPORT Datum sfp_in(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(sfp_in);
Datum sfp_in(PG_FUNCTION_ARGS) {
  CSfp fp;
  Sfp *b = DatumGetSfpP(DirectFunctionCall1(byteain, PG_GETARG_DATUM(0)));

  /* check correctness */
  fp = constructCSfp(b);
  freeCSfp(fp);

  PG_RETURN_SFP_P(b);
}

PGDLLEXPORT Datum sfp_out(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(sfp_out);
Datum sfp_out(PG_FUNCTION_ARGS) {
  PG_RETURN_DATUM(DirectFunctionCall1(byteaout, PG_GETARG_DATUM(0)));
}

PGDLLEXPORT Datum rdkit_version(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(rdkit_version);
Datum rdkit_version(PG_FUNCTION_ARGS) {
  char *ver = "" RDKITVER;
  char buf[1024];
  Assert(strlen(ver) == 6);
  snprintf(buf, sizeof(buf), "%d.%d.%d", atoi(pnstrdup(ver, 2)),
           atoi(pnstrdup(ver + 2, 2)), atoi(pnstrdup(ver + 4, 2)));

  PG_RETURN_TEXT_P(cstring_to_text(buf));
}

PGDLLEXPORT Datum rdkit_toolkit_version(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(rdkit_toolkit_version);
Datum rdkit_toolkit_version(PG_FUNCTION_ARGS) {
  const char *ver = "" RDK_TOOLKIT_VERSION;
  PG_RETURN_TEXT_P(cstring_to_text(ver));
}

/* chemical reactions */

PGDLLEXPORT Datum reaction_in(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(reaction_in);
Datum reaction_in(PG_FUNCTION_ARGS) {
  char *data = PG_GETARG_CSTRING(0);
  CChemicalReaction crxn;
  Reaction *rxn;

  crxn = parseChemReactText(data, false, false);

  if (!crxn) {
    ereport(ERROR, (errcode(ERRCODE_DATA_EXCEPTION),
                    errmsg("could not construct chemical reaction")));
  }
  rxn = deconstructChemReact(crxn);
  freeChemReaction(crxn);

  PG_RETURN_REACTION_P(rxn);
}

PGDLLEXPORT Datum reaction_recv(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(reaction_recv);
Datum reaction_recv(PG_FUNCTION_ARGS) {
  bytea *data = PG_GETARG_BYTEA_P(0);
  int len = VARSIZE(data) - VARHDRSZ;
  CChemicalReaction crxn;
  Reaction *rxn;

  crxn = parseChemReactBlob(VARDATA(data), len);

  rxn = deconstructChemReact(crxn);
  freeChemReaction(crxn);

  PG_FREE_IF_COPY(data, 0);

  PG_RETURN_REACTION_P(rxn);
}

PGDLLEXPORT Datum reaction_out(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(reaction_out);
Datum reaction_out(PG_FUNCTION_ARGS) {
  CChemicalReaction rxn;
  char *str;
  int len;

  fcinfo->flinfo->fn_extra =
      searchReactionCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                          PG_GETARG_DATUM(0), NULL, &rxn, NULL);
  str = makeChemReactText(rxn, &len, false);

  PG_RETURN_CSTRING(pnstrdup(str, len));
}

PGDLLEXPORT Datum reaction_send(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(reaction_send);
Datum reaction_send(PG_FUNCTION_ARGS) {
  CChemicalReaction rxn;
  bytea *res;
  char *str;
  int len;

  fcinfo->flinfo->fn_extra =
      searchReactionCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                          PG_GETARG_DATUM(0), NULL, &rxn, NULL);
  str = makeChemReactBlob(rxn, &len);

  res = (bytea *)palloc(len + VARHDRSZ);
  SET_VARSIZE(res, len + VARHDRSZ);
  memcpy(VARDATA(res), str, len);
  PG_RETURN_BYTEA_P(res);
}

PGDLLEXPORT Datum reaction_from_ctab(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(reaction_from_ctab);
Datum reaction_from_ctab(PG_FUNCTION_ARGS) {
  char *data = PG_GETARG_CSTRING(0);
  CChemicalReaction crxn;
  Reaction *rxn;

  crxn = parseChemReactCTAB(data, true);
  if (!crxn) {
    PG_RETURN_NULL();
  }
  rxn = deconstructChemReact(crxn);
  freeChemReaction(crxn);

  PG_RETURN_REACTION_P(rxn);
}

PGDLLEXPORT Datum reaction_from_smarts(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(reaction_from_smarts);
Datum reaction_from_smarts(PG_FUNCTION_ARGS) {
  char *data = PG_GETARG_CSTRING(0);
  CChemicalReaction crxn;
  Reaction *rxn;

  crxn = parseChemReactText(data, true, true);
  if (!crxn) {
    PG_RETURN_NULL();
  }
  rxn = deconstructChemReact(crxn);
  freeChemReaction(crxn);

  PG_RETURN_REACTION_P(rxn);
}

PGDLLEXPORT Datum reaction_from_smiles(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(reaction_from_smiles);
Datum reaction_from_smiles(PG_FUNCTION_ARGS) {
  char *data = PG_GETARG_CSTRING(0);
  CChemicalReaction crxn;
  Reaction *rxn;

  crxn = parseChemReactText(data, false, true);
  if (!crxn) {
    PG_RETURN_NULL();
  }
  rxn = deconstructChemReact(crxn);
  freeChemReaction(crxn);

  PG_RETURN_REACTION_P(rxn);
}

PGDLLEXPORT Datum reaction_to_ctab(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(reaction_to_ctab);
Datum reaction_to_ctab(PG_FUNCTION_ARGS) {
  CChemicalReaction rxn;
  char *str;
  int len;

  fcinfo->flinfo->fn_extra =
      searchReactionCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                          PG_GETARG_DATUM(0), NULL, &rxn, NULL);

  str = makeCTABChemReact(rxn, &len);

  PG_RETURN_CSTRING(pnstrdup(str, len));
}

PGDLLEXPORT Datum reaction_to_smiles(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(reaction_to_smiles);
Datum reaction_to_smiles(PG_FUNCTION_ARGS) {
  CChemicalReaction rxn;
  char *str;
  int len;

  fcinfo->flinfo->fn_extra =
      searchReactionCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                          PG_GETARG_DATUM(0), NULL, &rxn, NULL);
  str = makeChemReactText(rxn, &len, false);

  PG_RETURN_CSTRING(pnstrdup(str, len));
}

PGDLLEXPORT Datum reaction_to_smarts(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(reaction_to_smarts);
Datum reaction_to_smarts(PG_FUNCTION_ARGS) {
  CChemicalReaction rxn;
  char *str;
  int len;

  fcinfo->flinfo->fn_extra =
      searchReactionCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                          PG_GETARG_DATUM(0), NULL, &rxn, NULL);
  str = makeChemReactText(rxn, &len, true);

  PG_RETURN_CSTRING(pnstrdup(str, len));
}
