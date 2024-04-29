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

#include "rdkit.h"
#include "cache.h"

/***************** Mol operations ***********************/

#define MOLCMPFUNC(type, action, ret)                                     \
  PGDLLEXPORT Datum mol_##type(PG_FUNCTION_ARGS);                         \
  PG_FUNCTION_INFO_V1(mol_##type);                                        \
  Datum mol_##type(PG_FUNCTION_ARGS) {                                    \
    CROMol a, b;                                                          \
    int res;                                                              \
                                                                          \
    fcinfo->flinfo->fn_extra =                                            \
        searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt, \
                       PG_GETARG_DATUM(0), NULL, &a, NULL);               \
    fcinfo->flinfo->fn_extra =                                            \
        searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt, \
                       PG_GETARG_DATUM(1), NULL, &b, NULL);               \
    res = molcmp(a, b);                                                   \
    PG_RETURN_##ret(res action 0);                                        \
  }                                                                       \
  /* keep compiler quiet - no extra ; */                                  \
  extern int no_such_variable

MOLCMPFUNC(lt, <, BOOL);
MOLCMPFUNC(le, <=, BOOL);
MOLCMPFUNC(eq, ==, BOOL);
MOLCMPFUNC(ge, >=, BOOL);
MOLCMPFUNC(gt, >, BOOL);
MOLCMPFUNC(ne, !=, BOOL);
MOLCMPFUNC(cmp, +, INT32);

PGDLLEXPORT Datum is_valid_smiles(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(is_valid_smiles);
Datum is_valid_smiles(PG_FUNCTION_ARGS) {
  char *data = PG_GETARG_CSTRING(0);
  PG_RETURN_BOOL(isValidSmiles(data));
}

PGDLLEXPORT Datum is_valid_smarts(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(is_valid_smarts);
Datum is_valid_smarts(PG_FUNCTION_ARGS) {
  char *data = PG_GETARG_CSTRING(0);
  PG_RETURN_BOOL(isValidSmarts(data));
}

PGDLLEXPORT Datum is_valid_ctab(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(is_valid_ctab);
Datum is_valid_ctab(PG_FUNCTION_ARGS) {
  char *data = PG_GETARG_CSTRING(0);
  PG_RETURN_BOOL(isValidCTAB(data));
}

PGDLLEXPORT Datum is_valid_mol_pkl(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(is_valid_mol_pkl);
Datum is_valid_mol_pkl(PG_FUNCTION_ARGS) {
  bytea *data = PG_GETARG_BYTEA_P(0);
  int len = VARSIZE(data) - VARHDRSZ;
  bool res = isValidMolBlob(VARDATA(data), len);
  PG_FREE_IF_COPY(data, 0);
  PG_RETURN_BOOL(res);
}

PGDLLEXPORT Datum mol_substruct(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_substruct);
Datum mol_substruct(PG_FUNCTION_ARGS) {
  CROMol i, a;

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &i, NULL);
  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(1), NULL, &a, NULL);

  PG_RETURN_BOOL(MolSubstruct(i, a, false, false));
}

PGDLLEXPORT Datum mol_substruct_query(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_substruct_query);
Datum mol_substruct_query(PG_FUNCTION_ARGS) {
  CROMol i, a;

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &i, NULL);
  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(1), NULL, &a, NULL);

  PG_RETURN_BOOL(MolSubstruct(i, a, false, true));
}

PGDLLEXPORT Datum mol_substruct_chiral(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_substruct_chiral);
Datum mol_substruct_chiral(PG_FUNCTION_ARGS) {
  CROMol i, a;

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &i, NULL);
  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(1), NULL, &a, NULL);

  PG_RETURN_BOOL(MolSubstruct(i, a, true, false));
}

PGDLLEXPORT Datum mol_rsubstruct(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_rsubstruct);
Datum mol_rsubstruct(PG_FUNCTION_ARGS) {
  CROMol i, a;

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &i, NULL);
  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(1), NULL, &a, NULL);

  PG_RETURN_BOOL(MolSubstruct(a, i, false, false));
}

PGDLLEXPORT Datum mol_rsubstruct_query(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_rsubstruct_query);
Datum mol_rsubstruct_query(PG_FUNCTION_ARGS) {
  CROMol i, a;

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &i, NULL);
  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(1), NULL, &a, NULL);

  PG_RETURN_BOOL(MolSubstruct(a, i, false, true));
}

PGDLLEXPORT Datum mol_rsubstruct_chiral(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_rsubstruct_chiral);
Datum mol_rsubstruct_chiral(PG_FUNCTION_ARGS) {
  CROMol i, a;

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &i, NULL);
  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(1), NULL, &a, NULL);

  PG_RETURN_BOOL(MolSubstruct(a, i, true, false));
}

PGDLLEXPORT Datum mol_substruct_count(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_substruct_count);
Datum mol_substruct_count(PG_FUNCTION_ARGS) {
  CROMol i, a;
  bool uniquify = PG_GETARG_BOOL(2);

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &i, NULL);
  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(1), NULL, &a, NULL);

  PG_RETURN_INT32(MolSubstructCount(i, a, uniquify, false));
}
PGDLLEXPORT Datum mol_substruct_count_chiral(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_substruct_count_chiral);
Datum mol_substruct_count_chiral(PG_FUNCTION_ARGS) {
  CROMol i, a;
  bool uniquify = PG_GETARG_BOOL(2);

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &i, NULL);
  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(1), NULL, &a, NULL);

  PG_RETURN_INT32(MolSubstructCount(i, a, uniquify, true));
}

PGDLLEXPORT Datum mol_xq_substruct(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_xq_substruct);
Datum mol_xq_substruct(PG_FUNCTION_ARGS) {
  CROMol i;
  CXQMol a;

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &i, NULL);
  fcinfo->flinfo->fn_extra =
      searchXQMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                       PG_GETARG_DATUM(1), NULL, &a, NULL);

  PG_RETURN_BOOL(XQMolSubstruct(i, a, false, false));
}

PGDLLEXPORT Datum mol_xq_rsubstruct(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_xq_rsubstruct);
Datum mol_xq_rsubstruct(PG_FUNCTION_ARGS) {
  CROMol i;
  CXQMol a;

  fcinfo->flinfo->fn_extra =
      searchXQMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                       PG_GETARG_DATUM(0), NULL, &a, NULL);
  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(1), NULL, &i, NULL);

  PG_RETURN_BOOL(XQMolSubstruct(i, a, false, false));
}

PGDLLEXPORT Datum mol_xq_substruct_query(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_xq_substruct_query);
Datum mol_xq_substruct_query(PG_FUNCTION_ARGS) {
  CROMol i;
  CXQMol a;

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &i, NULL);
  fcinfo->flinfo->fn_extra =
      searchXQMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                       PG_GETARG_DATUM(1), NULL, &a, NULL);

  PG_RETURN_BOOL(XQMolSubstruct(i, a, false, true));
}
PGDLLEXPORT Datum mol_xq_rsubstruct_query(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_xq_rsubstruct_query);
Datum mol_xq_rsubstruct_query(PG_FUNCTION_ARGS) {
  CROMol i;
  CXQMol a;

  fcinfo->flinfo->fn_extra =
      searchXQMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                       PG_GETARG_DATUM(0), NULL, &a, NULL);
  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(1), NULL, &i, NULL);

  PG_RETURN_BOOL(XQMolSubstruct(i, a, false, true));
}

#define MOLDESCR(name, func, ret)                                         \
  PGDLLEXPORT Datum mol_##name(PG_FUNCTION_ARGS);                         \
  PG_FUNCTION_INFO_V1(mol_##name);                                        \
  Datum mol_##name(PG_FUNCTION_ARGS) {                                    \
    CROMol i;                                                             \
    fcinfo->flinfo->fn_extra =                                            \
        searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt, \
                       PG_GETARG_DATUM(0), NULL, &i, NULL);               \
    PG_RETURN_##ret(func(i));                                             \
  }

MOLDESCR(amw, MolAMW, FLOAT4)
MOLDESCR(exactmw, MolExactMW, FLOAT4)
MOLDESCR(logp, MolLogP, FLOAT4)
MOLDESCR(tpsa, MolTPSA, FLOAT4)
MOLDESCR(labuteasa, MolLabuteASA, FLOAT4)
MOLDESCR(hba, MolHBA, INT32)
MOLDESCR(hbd, MolHBD, INT32)
MOLDESCR(numatoms, MolNumAtoms, INT32)
MOLDESCR(numheavyatoms, MolNumHeavyAtoms, INT32)
MOLDESCR(numrotatablebonds, MolNumRotatableBonds, INT32)
MOLDESCR(numheteroatoms, MolNumHeteroatoms, INT32)
MOLDESCR(numrings, MolNumRings, INT32)
MOLDESCR(numaromaticrings, MolNumAromaticRings, INT32)
MOLDESCR(numaliphaticrings, MolNumAliphaticRings, INT32)
MOLDESCR(numsaturatedrings, MolNumSaturatedRings, INT32)
MOLDESCR(numaromaticheterocycles, MolNumAromaticHeterocycles, INT32)
MOLDESCR(numaliphaticheterocycles, MolNumAliphaticHeterocycles, INT32)
MOLDESCR(numsaturatedheterocycles, MolNumSaturatedHeterocycles, INT32)
MOLDESCR(numaromaticcarbocycles, MolNumAromaticCarbocycles, INT32)
MOLDESCR(numaliphaticcarbocycles, MolNumAliphaticCarbocycles, INT32)
MOLDESCR(numsaturatedcarbocycles, MolNumSaturatedCarbocycles, INT32)
MOLDESCR(numheterocycles, MolNumHeterocycles, INT32)
MOLDESCR(numamidebonds, MolNumAmideBonds, INT32)

MOLDESCR(fractioncsp3, MolFractionCSP3, FLOAT4)
MOLDESCR(chi0n, MolChi0n, FLOAT4)
MOLDESCR(chi1n, MolChi1n, FLOAT4)
MOLDESCR(chi2n, MolChi2n, FLOAT4)
MOLDESCR(chi3n, MolChi3n, FLOAT4)
MOLDESCR(chi4n, MolChi4n, FLOAT4)
MOLDESCR(chi0v, MolChi0v, FLOAT4)
MOLDESCR(chi1v, MolChi1v, FLOAT4)
MOLDESCR(chi2v, MolChi2v, FLOAT4)
MOLDESCR(chi3v, MolChi3v, FLOAT4)
MOLDESCR(chi4v, MolChi4v, FLOAT4)
MOLDESCR(kappa1, MolKappa1, FLOAT4)
MOLDESCR(kappa2, MolKappa2, FLOAT4)
MOLDESCR(kappa3, MolKappa3, FLOAT4)
MOLDESCR(hallkieralpha, MolHallKierAlpha, FLOAT4)
MOLDESCR(phi, MolPhi, FLOAT4)

MOLDESCR(numspiroatoms, MolNumSpiroAtoms, INT32)
MOLDESCR(numbridgeheadatoms, MolNumBridgeheadAtoms, INT32)

PGDLLEXPORT Datum mol_formula(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_formula);
Datum mol_formula(PG_FUNCTION_ARGS) {
  CROMol mol;
  char *str;
  int len;

  bool separateIsotopes = PG_GETARG_BOOL(1);
  bool abbreviateHIsotopes = PG_GETARG_BOOL(2);

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);

  str = makeMolFormulaText(mol, &len, separateIsotopes, abbreviateHIsotopes);

  PG_RETURN_CSTRING(pnstrdup(str, len));
}

PGDLLEXPORT Datum mol_inchi(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_inchi);
Datum mol_inchi(PG_FUNCTION_ARGS) {
  CROMol mol;
  const char *str;
  char *res, *opts = PG_GETARG_CSTRING(1);

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);
  str = MolInchi(mol, opts);
  if (*str == 0) {
    free((void *)str);
    PG_RETURN_NULL();
  } else {
    res = pnstrdup(str, strlen(str));
    free((void *)str);
    PG_RETURN_CSTRING(res);
  }
}

PGDLLEXPORT Datum mol_inchikey(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_inchikey);
Datum mol_inchikey(PG_FUNCTION_ARGS) {
  CROMol mol;
  const char *str;
  char *res, *opts = PG_GETARG_CSTRING(1);

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);
  str = MolInchiKey(mol, opts);
  if (*str == 0) {
    free((void *)str);
    PG_RETURN_NULL();
  } else {
    res = pnstrdup(str, strlen(str));
    free((void *)str);
    PG_RETURN_CSTRING(res);
  }
}
PGDLLEXPORT Datum mol_murckoscaffold(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_murckoscaffold);
Datum mol_murckoscaffold(PG_FUNCTION_ARGS) {
  CROMol mol;
  Mol *res;
  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);
  CROMol scaffold = MolMurckoScaffold(mol);
  if (!scaffold) {
    PG_RETURN_NULL();
  }
  res = deconstructROMol(scaffold);
  freeCROMol(scaffold);

  PG_RETURN_MOL_P(res);
}

PGDLLEXPORT Datum mol_nm_hash(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_nm_hash);
Datum mol_nm_hash(PG_FUNCTION_ARGS) {
  CROMol mol;
  char *opts = PG_GETARG_CSTRING(1);
  int len;
  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);
  Assert(mol != 0);
  char *str = computeNMMolHash(mol, opts);
  Assert(str != 0 && strlen(str) != 0);
  char *res = pnstrdup(str, strlen(str));
  free((void *)str);
  PG_RETURN_CSTRING(res);
}

PGDLLEXPORT Datum mol_adjust_query_properties(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_adjust_query_properties);
Datum mol_adjust_query_properties(PG_FUNCTION_ARGS) {
  CROMol mol;
  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);
  Assert(mol != 0);
  char *data = PG_GETARG_CSTRING(1);

  CROMol adj = MolAdjustQueryProperties(mol, data);
  if (!adj) {
    PG_RETURN_NULL();
  }
  Mol *res = deconstructROMolWithQueryProperties(adj);
  freeCROMol(adj);

  PG_RETURN_MOL_P(res);
}

PGDLLEXPORT Datum mol_to_svg(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_to_svg);
Datum mol_to_svg(PG_FUNCTION_ARGS) {
  CROMol mol;
  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);
  Assert(mol != 0);

  char *legend = PG_GETARG_CSTRING(1);
  unsigned int w = PG_GETARG_UINT32(2);
  unsigned int h = PG_GETARG_UINT32(3);
  char *params = PG_GETARG_CSTRING(4);

  char *str = MolGetSVG(mol, w, h, legend, params);
  char *res = pnstrdup(str, strlen(str));
  free((void *)str);
  PG_RETURN_CSTRING(res);
}

PGDLLEXPORT Datum mol_to_xqmol(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(mol_to_xqmol);
Datum mol_to_xqmol(PG_FUNCTION_ARGS) {
  CROMol mol;
  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(0), NULL, &mol, NULL);
  Assert(mol != 0);

  bool doEnumeration = PG_GETARG_BOOL(1);
  bool doTautomers = PG_GETARG_BOOL(2);
  bool adjustQueryProperties = PG_GETARG_BOOL(3);
  char *params = PG_GETARG_CSTRING(4);

  CXQMol xqm = MolToXQMol(mol, doEnumeration, doTautomers,
                          adjustQueryProperties, params);

  if (!xqm) {
    PG_RETURN_NULL();
  }
  XQMol *res = deconstructXQMol(xqm);
  freeCXQMol(xqm);

  PG_RETURN_MOL_P(res);
}

/*** fmcs ***/

PGDLLEXPORT Datum fmcs_smiles(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(fmcs_smiles);
Datum fmcs_smiles(PG_FUNCTION_ARGS) {
  char *str = PG_GETARG_CSTRING(0);
  char *params = PG_GETARG_CSTRING(1);
  // elog(WARNING, str);

  str = findMCSsmiles(str, params);
  // elog(WARNING, str);
  Assert(str != 0);

  char *res = pnstrdup(str, strlen(str));
  free((void *)str);
  PG_RETURN_CSTRING(res);
}

PGDLLEXPORT Datum fmcs_smiles_transition(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(fmcs_smiles_transition);
Datum fmcs_smiles_transition(PG_FUNCTION_ARGS) {
  if (!AggCheckCallContext(fcinfo, NULL) || PG_ARGISNULL(0)) {
    ereport(
        ERROR,
        (errmsg(
            "fmcs_smiles_transition() called in out of aggregate context")));
  } else if (!PG_ARGISNULL(0)) {  // Called in aggregate context...
    text *t0 = PG_GETARG_TEXT_P(0);
    text *tsmiles = PG_GETARG_TEXT_P(1);
    // char    *s0 = VARDATA(t0);
    // char    *smiles = VARDATA(tsmiles);
    // elog(WARNING, "fmcs_trans: next iteration in the same run");
    // elog(WARNING, s0);
    // elog(WARNING, smiles);

    int32 ts_size = VARSIZE(t0) + 1 + VARSIZE(tsmiles) - VARHDRSZ;
    text *ts = (text *)palloc(ts_size);  // new return value
    SET_VARSIZE(ts, ts_size);
    memcpy(VARDATA(ts), VARDATA(t0), VARSIZE(t0) - VARHDRSZ);
    *(char *)(VARDATA(ts) + VARSIZE(t0) - VARHDRSZ) = ' ';
    memcpy(VARDATA(ts) + VARSIZE(t0) - VARHDRSZ + 1, VARDATA(tsmiles),
           VARSIZE(tsmiles) - VARHDRSZ);
    // elog(WARNING, VARDATA(ts));
    PG_RETURN_TEXT_P(ts);
  }
}
//------------------------

char *Mol2Smiles(CROMol data);

PGDLLEXPORT Datum fmcs_mol2s_transition(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(fmcs_mol2s_transition);
Datum fmcs_mol2s_transition(PG_FUNCTION_ARGS) {
  // elog(WARNING, (PG_ARGISNULL(0)) ? "arg 0 is NULL" : "arg 0 is NOT NULL");
  // elog(WARNING, (PG_ARGISNULL(1)) ? "arg 1 is NULL" : "arg 1 is NOT NULL");

  if (!AggCheckCallContext(fcinfo, NULL)) {
    // elog(WARNING, "fmcs_mol2s_transition() called in out of aggregate
    // context");
    ereport(
        ERROR,
        (errmsg("fmcs_mol2s_transition() called in out of aggregate context")));
  }
  if (PG_ARGISNULL(0) && !PG_ARGISNULL(1)) {  // first call
    /// elog(WARNING, "fmcs_mol2s_transition() called first time");
    CROMol mol;
    int len, ts_size;
    char *smiles;
    // elog(WARNING, "mol=%p, fcinfo: %p, %p", mol, fcinfo->flinfo->fn_extra,
    //      fcinfo->flinfo->fn_mcxt);
    fcinfo->flinfo->fn_extra =
        searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                       PG_GETARG_DATUM(1), NULL, &mol, NULL);

    smiles = makeMolText(mol, &len, false, false, true);

    //        char *smiles= Mol2Smiles(mol);
    //        int   len   = strlen(smiles);
    /// elog(WARNING, smiles);

    ts_size = len + VARHDRSZ;
    text *ts = (text *)palloc(ts_size);  // new return value
    SET_VARSIZE(ts, ts_size);
    memcpy(VARDATA(ts), smiles, len);
    PG_RETURN_TEXT_P(ts);
  } else if (!PG_ARGISNULL(0) &&
             !PG_ARGISNULL(1)) {  // Called in aggregate context...
    text *t0 = PG_GETARG_TEXT_P(0);
    /// elog(WARNING, "fmcs_mol2s_transition(): next iteration in the same
    /// run");
    // elog(WARNING, VARDATA(t0));

    // mol_to_smiles():
    CROMol mol;
    int len;
    // elog(WARNING, "mol=%p, fcinfo: %p, %p", mol, fcinfo->flinfo->fn_extra,
    //      fcinfo->flinfo->fn_mcxt);
    fcinfo->flinfo->fn_extra =
        searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                       PG_GETARG_DATUM(1), NULL, &mol, NULL);

    char *smiles = makeMolText(mol, &len, false, false, true);
    //        char *smiles= Mol2Smiles(mol);
    //        int   len   = strlen(smiles);
    /// elog(WARNING, smiles);

    int32 ts_size = VARSIZE(t0) + 1 + len;
    text *ts = (text *)palloc(ts_size);  // new return value
    SET_VARSIZE(ts, ts_size);
    memcpy(VARDATA(ts), VARDATA(t0), VARSIZE(t0) - VARHDRSZ);
    *(char *)(VARDATA(ts) + VARSIZE(t0) - VARHDRSZ) = ' ';
    memcpy(VARDATA(ts) + VARSIZE(t0) - VARHDRSZ + 1, smiles, len);
    // elog(WARNING, VARDATA(ts));
    PG_RETURN_TEXT_P(ts);
  }
  //------
  /// elog(WARNING, "fmcs_mol2s_transition(): return empty text block");
  {
    int32 ts_size = VARHDRSZ;
    text *ts = (text *)palloc(ts_size);  // return empty text block
    SET_VARSIZE(ts, ts_size);
    PG_RETURN_TEXT_P(ts);
  }
}
//------------------------

// fmcs_mol:
PGDLLEXPORT Datum fmcs_mols(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(fmcs_mols);
Datum fmcs_mols(PG_FUNCTION_ARGS) {
  // elog(WARNING, "fmcs_mols() called. FINALFUNC");
  void *lst = PG_GETARG_POINTER(0);
  // char t[256];
  // sprintf(t,"lst=%p, fcinfo: %p, %p", lst, fcinfo->flinfo->fn_extra,
  // fcinfo->flinfo->fn_mcxt);
  // elog(WARNING, t);
  //  char *params = PG_GETARG_CSTRING(1);
  char *str = findMCS(lst, NULL);  // params);
  Assert(str != 0);
  int32 ts_size = VARHDRSZ + strlen(str);
  text *ts = (text *)palloc(ts_size);
  SET_VARSIZE(ts, ts_size);
  memcpy(VARDATA(ts), str, strlen(str));
  free((void *)str);
  PG_RETURN_TEXT_P(ts);
}

PGDLLEXPORT Datum fmcs_mol_transition(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(fmcs_mol_transition);
Datum fmcs_mol_transition(PG_FUNCTION_ARGS) {
  // elog(WARNING, (PG_ARGISNULL(0)) ? "arg 0 is NULL" : "arg 0 is NOT NULL");
  // elog(WARNING, (PG_ARGISNULL(1)) ? "arg 1 is NULL" : "arg 1 is NOT NULL");

  if (!AggCheckCallContext(fcinfo, NULL)) {
    ereport(
        ERROR,
        (errmsg("fmcs_mol_transition() called in out of aggregate context")));
  }
  if (PG_ARGISNULL(0) && !PG_ARGISNULL(1)) {  // first call
    // elog(WARNING, "fmcs_mol_transition() called first time");
    void *lst = NULL;
    Mol *mol = PG_GETARG_MOL_P(1);
    lst = addMol2list(NULL, mol);
    // char t[256];
    // sprintf(t,"mol=%p, lst=%p, fcinfo: %p, %p", mol, lst,
    // fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt);
    // elog(WARNING, t);
    // fcinfo->flinfo->fn_extra = lst;
    PG_RETURN_POINTER(lst);
  } else if (!PG_ARGISNULL(0) &&
             !PG_ARGISNULL(1)) {  // Called in aggregate context...
    // elog(WARNING, "fmcs_mol_transition(): next iteration in the same run");
    Mol *mol = PG_GETARG_MOL_P(1);
    void *lst = addMol2list(PG_GETARG_POINTER(0), mol);
    // char t[256];
    // sprintf(t,"mol=%p, lst=%p, fcinfo: %p, %p", mol, lst,
    // fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt);
    // elog(WARNING, t);
    PG_RETURN_POINTER(lst);
  }
}
//------------------------

//==================================
/*
#include <utils/array.h>
PGDLLEXPORT Datum           fmcs_mol(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(fmcs_mol);
Datum
fmcs_mol(PG_FUNCTION_ARGS) {
elog(WARNING, "fmcs_mol(): FINAL function in the same run.");
  char *str = "";
  ArrayType *ma = (ArrayType*)PG_GETARG_ARRAYTYPE_P(0);
  Assert(ma != 0);
  int   len = VARSIZE(ma) - VARHDRSZ;
  if(len>1){
    CROMol  *mols = (CROMol*) palloc(len*sizeof(CROMol));
    int i;
    for(i=0; i < len; i++){
        mols[i] = ((CROMol*)ARR_DATA_PTR(ma))[i];
elog(WARNING, "fmcs_mol(): call searchMolCache(...).");
//--
      fcinfo->flinfo->fn_extra = searchMolCache(
                                                fcinfo->flinfo->fn_extra,
                                                fcinfo->flinfo->fn_mcxt,
                                                ((CROMol*)ARR_DATA_PTR(ma))[i],
                                                NULL, mols[i], NULL);
//--
      if(mols[i]==0){
        pfree(mols);
        mols = NULL;
      }
      Assert(mols[i] != 0);
    }
elog(WARNING, "fmcs_mol(): call findMCS(...).");
    char *params = "";
    str = findMCS(mols, len, params);
    pfree(mols); mols = NULL;
  }
  else
    elog(WARNING, "fmcs_mol(): mols is empty.");
  Assert(str != 0);
  len = strlen(str);
  text* tmcs = (text*) palloc(len+1+VARHDRSZ);
  SET_VARSIZE(tmcs, len+1+VARHDRSZ);
  memcpy(VARDATA(tmcs), str, len+1);
  PG_RETURN_TEXT_P(tmcs);
}
PGDLLEXPORT Datum           fmcs_mol_transition(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(fmcs_mol_transition);
Datum
fmcs_mol_transition(PG_FUNCTION_ARGS) {
    if ( ! AggCheckCallContext(fcinfo, NULL)){
      ereport(ERROR, (errmsg("fmcs_mol_transition() called in out of aggregate
context")));
      PG_RETURN_NULL();
    }
    else{ // Called in aggregate context...
      ArrayType *mols = NULL;
      if(PG_ARGISNULL(0)){
elog(WARNING, "fmcs_mol_transition(): first call. Allocate state type array of
CROMol pointers");
if(PG_ARGISNULL(1))
  elog(WARNING, "fmcs_mol_transition(): first call. Argument [1] is null.");
         Oid elmtype = get_fn_expr_argtype(fcinfo->flinfo, 1);	// CROMol
         mols = construct_empty_array(elmtype);
         Assert(mols != 0);
       }
       else{
elog(WARNING, "fmcs_mol_transition(): next iteration in the same run. Append
CROMol pointer");
         mols = PG_GETARG_ARRAYTYPE_P(0);
         CROMol mol = PG_GETARG_DATUM(1);
         mols;
       }
       PG_RETURN_ARRAYTYPE_P(mols);
    }
}
*/

/*
PGDLLEXPORT Datum           fmcs_transition(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(fmcs_transition);
Datum
fmcs_transition(PG_FUNCTION_ARGS) {
  text    *tsmiles = PG_GETARG_TEXT(1);
//  int      tsmiles_len = VARSIZE(tsmiles);
  char    *smiles = pgstrdup(VARDATA(tsmiles), VARSIZE(tsmiles));
  char    *params = PG_GETARG_CSTRING(2);
elog(WARNING, smiles);
//  Assert(smiles != 0);
    if (AggCheckCallContext(fcinfo, NULL))
    {
        // Called in aggregate context...
        if (PG_ARGISNULL(0)) {
            // ... for the first time in a run, so the state in the 1st
            // argument is null. Create a state-holder array by copying the
            // second input array and return it.
elog(WARNING, "fmcs_trans: the first time in a run");
  PG_RETURN_TEXT(tsmiles);
//            PG_RETURN_POINTER(copy_intArrayType(b));
        }
        else {
            // ... for a later invocation in the same run, so we'll modify the
state array directly.
elog(WARNING, "fmcs_trans: next iteration in the same run");
  PG_RETURN_TEXT(tsmiles);
        }
    }
    else
    {
elog(WARNING, "fmcs_trans: Not in aggregate context");
        // Not in aggregate context
        if (PG_ARGISNULL(0))
            ereport(ERROR, (errmsg("First operand must be non-null")));
        else
;;;;;;
  PG_RETURN_TEXT(tsmiles);
    }
}
*/
