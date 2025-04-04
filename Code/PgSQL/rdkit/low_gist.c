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
#include <access/gist.h>
#include <access/skey.h>
#include <utils/memutils.h>

#include "rdkit.h"
#include "cache.h"

#define GETENTRY(vec, pos) ((bytea *)DatumGetPointer((vec)->vector[(pos)].key))

PGDLLEXPORT Datum gslfp_compress(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gslfp_compress);
Datum gslfp_compress(PG_FUNCTION_ARGS) {
  GISTENTRY *entry = (GISTENTRY *)PG_GETARG_POINTER(0);
  GISTENTRY *retval = entry;

  if (entry->leafkey) {
    CSfp fp = constructCSfp(DatumGetSfpP(entry->key));

    retval = (GISTENTRY *)palloc(sizeof(GISTENTRY));

    gistentryinit(*retval,
                  PointerGetDatum(makeLowSparseFingerPrint(fp, NUMRANGE)),
                  entry->rel, entry->page, entry->offset, false);
    freeCSfp(fp);
  }

  PG_RETURN_POINTER(retval);
}

PGDLLEXPORT Datum gslfp_decompress(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gslfp_decompress);
Datum gslfp_decompress(PG_FUNCTION_ARGS) {
  GISTENTRY *entry = (GISTENTRY *)PG_GETARG_POINTER(0);
  bytea *key = (bytea *)PG_DETOAST_DATUM(entry->key);

  if (key != (bytea *)DatumGetPointer(entry->key)) {
    GISTENTRY *retval = (GISTENTRY *)palloc(sizeof(GISTENTRY));

    gistentryinit(*retval, PointerGetDatum(key), entry->rel, entry->page,
                  entry->offset, false);

    PG_RETURN_POINTER(retval);
  }

  PG_RETURN_POINTER(entry);
}

static void adjustKey(IntRange *s, IntRange *k) {
  int j;

  for (j = 0; j < NUMRANGE; j++) {
    /* set minimal non-zero value */
    if (k[j].low > 0 && (s[j].low == 0 || k[j].low < s[j].low))
      s[j].low = k[j].low;
    /* set maximum value */
    if (k[j].high > s[j].high) s[j].high = k[j].high;
  }
}

PGDLLEXPORT Datum gslfp_union(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gslfp_union);
Datum gslfp_union(PG_FUNCTION_ARGS) {
  GistEntryVector *entryvec = (GistEntryVector *)PG_GETARG_POINTER(0);
  int *size = (int *)PG_GETARG_POINTER(1);
  int32 i;
  bytea *result, *key;
  IntRange *s, *k;

  *size = VARHDRSZ + NUMRANGE * sizeof(IntRange);
  result = palloc0(*size);
  SET_VARSIZE(result, *size);
  key = GETENTRY(entryvec, 0);
  memcpy(VARDATA(result), VARDATA(key), NUMRANGE * sizeof(IntRange));

  s = (IntRange *)VARDATA(result);
  for (i = 1; i < entryvec->n; i++) {
    key = GETENTRY(entryvec, i);
    k = (IntRange *)VARDATA(key);

    adjustKey(s, k);
  }

  PG_RETURN_POINTER(result);
}

/*
 * Same method
 */

PGDLLEXPORT Datum gslfp_same(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gslfp_same);
Datum gslfp_same(PG_FUNCTION_ARGS) {
  bytea *a = (bytea *)PG_GETARG_POINTER(0);
  bytea *b = (bytea *)PG_GETARG_POINTER(1);
  bool *result = (bool *)PG_GETARG_POINTER(2);

  *result = (memcmp(VARDATA(a), VARDATA(b), VARSIZE(a) - VARHDRSZ) == 0)
                ? true
                : false;

  PG_RETURN_POINTER(result);
}

static uint32 distance(bytea *a, bytea *b) {
  int i;
  uint32 dist = 0;
  IntRange *as = (IntRange *)VARDATA(a), *bs = (IntRange *)VARDATA(b);

  if (VARSIZE(a) != VARSIZE(b))
    elog(ERROR, "All fingerprints should be the same length");
  for (i = 0; i < NUMRANGE; i++) {
    if (as[i].low > bs[i].low)
      dist += as[i].low - bs[i].low;
    else if (as[i].low < bs[i].low)
      dist += bs[i].low - as[i].low;

    if (as[i].high > bs[i].high)
      dist += as[i].high - bs[i].high;
    else if (as[i].high < bs[i].high)
      dist += bs[i].high - as[i].high;
  }
  return dist;
}

static uint32 penalty(bytea *origval, bytea *newval) {
  int i;
  uint32 dist = 0;
  IntRange *as = (IntRange *)VARDATA(origval),
           *bs = (IntRange *)VARDATA(newval);

  if (VARSIZE(origval) != VARSIZE(newval))
    elog(ERROR, "All fingerprints should be the same length");

  for (i = 0; i < NUMRANGE; i++) {
    if (bs[i].low > 0) {
      if (as[i].low == 0)
        dist += bs[i].low;
      else if (bs[i].low < as[i].low)
        dist += as[i].low - bs[i].low;
    }

    if (bs[i].high > as[i].high) dist += bs[i].high - as[i].high;
  }

  return dist;
}

PGDLLEXPORT Datum gslfp_penalty(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gslfp_penalty);
Datum gslfp_penalty(PG_FUNCTION_ARGS) {
  GISTENTRY *origentry =
      (GISTENTRY *)PG_GETARG_POINTER(0); /* always ISSIGNKEY */
  GISTENTRY *newentry = (GISTENTRY *)PG_GETARG_POINTER(1);
  float *p = (float *)PG_GETARG_POINTER(2);
  bytea *origval = (bytea *)DatumGetPointer(origentry->key);
  bytea *newval = (bytea *)DatumGetPointer(newentry->key);

  *p = (float)penalty(origval, newval);

  PG_RETURN_POINTER(p);
}

/*
 * Picksplit method
 */

#define WISH_F(a, b, c) \
  (double)(-(double)(((a) - (b)) * ((a) - (b)) * ((a) - (b))) * (c))

typedef struct {
  OffsetNumber pos;
  uint32 cost;
} SPLITCOST;

static int comparecost(const void *va, const void *vb) {
  SPLITCOST *a = (SPLITCOST *)va;
  SPLITCOST *b = (SPLITCOST *)vb;

  if (a->cost == b->cost)
    return 0;
  else
    return (a->cost > b->cost) ? 1 : -1;
}

PGDLLEXPORT Datum gslfp_picksplit(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gslfp_picksplit);
Datum gslfp_picksplit(PG_FUNCTION_ARGS) {
  GistEntryVector *entryvec = (GistEntryVector *)PG_GETARG_POINTER(0);
  GIST_SPLITVEC *v = (GIST_SPLITVEC *)PG_GETARG_POINTER(1);
  OffsetNumber k, j;
  bytea *datum_l, *datum_r;
  int32 size_alpha, size_beta;
  int32 size_waste, waste = -1;
  int32 nbytes;
  OffsetNumber seed_1 = 0, seed_2 = 0;
  OffsetNumber *left, *right;
  OffsetNumber maxoff;
  SPLITCOST *costvector;

  maxoff = entryvec->n - 1;
  nbytes = (maxoff + 2) * sizeof(OffsetNumber);
  v->spl_left = (OffsetNumber *)palloc(nbytes);
  v->spl_right = (OffsetNumber *)palloc(nbytes);

  for (k = FirstOffsetNumber; k < maxoff; k = OffsetNumberNext(k)) {
    for (j = OffsetNumberNext(k); j <= maxoff; j = OffsetNumberNext(j)) {
      size_waste = distance(GETENTRY(entryvec, j), GETENTRY(entryvec, k));
      if (size_waste > waste) {
        waste = size_waste;
        seed_1 = k;
        seed_2 = j;
      }
    }
  }

  left = v->spl_left;
  v->spl_nleft = 0;
  right = v->spl_right;
  v->spl_nright = 0;

  if (waste == 0) {
    int len = VARHDRSZ + sizeof(IntRange) * NUMRANGE;
    /* all entries are all the same */

    for (k = FirstOffsetNumber; k <= maxoff; k = OffsetNumberNext(k)) {
      if (k <= (maxoff - FirstOffsetNumber + 1) / 2) {
        v->spl_left[v->spl_nleft] = k;
        v->spl_nleft++;
      } else {
        v->spl_right[v->spl_nright] = k;
        v->spl_nright++;
      }
    }

    datum_l = palloc(len);
    memcpy(datum_l, GETENTRY(entryvec, FirstOffsetNumber), len);
    v->spl_ldatum = PointerGetDatum(datum_l);
    datum_r = palloc(len);
    memcpy(datum_r, GETENTRY(entryvec, FirstOffsetNumber), len);
    v->spl_rdatum = PointerGetDatum(datum_r);

    Assert(v->spl_nleft + v->spl_nright == maxoff);
    PG_RETURN_POINTER(v);
  }

  if (seed_1 == 0 || seed_2 == 0) {
    seed_1 = 1;
    seed_2 = 2;
  }

  /* form initial .. */
  datum_l = palloc(sizeof(IntRange) * NUMRANGE + VARHDRSZ);
  memcpy(datum_l, GETENTRY(entryvec, seed_1),
         sizeof(IntRange) * NUMRANGE + VARHDRSZ);
  datum_r = palloc(sizeof(IntRange) * NUMRANGE + VARHDRSZ);
  memcpy(datum_r, GETENTRY(entryvec, seed_2),
         sizeof(IntRange) * NUMRANGE + VARHDRSZ);

  /* sort before ... */
  costvector = (SPLITCOST *)palloc(sizeof(SPLITCOST) * maxoff);
  for (j = FirstOffsetNumber; j <= maxoff; j = OffsetNumberNext(j)) {
    costvector[j - 1].pos = j;
    size_alpha = distance(datum_l, GETENTRY(entryvec, j));
    size_beta = distance(datum_r, GETENTRY(entryvec, j));
    costvector[j - 1].cost = (size_alpha > size_beta)
                                 ? (size_alpha - size_beta)
                                 : (size_beta - size_alpha);
  }
  qsort((void *)costvector, maxoff, sizeof(SPLITCOST), comparecost);

  for (k = 0; k < maxoff; k++) {
    j = costvector[k].pos;
    if (j == seed_1) {
      *left++ = j;
      v->spl_nleft++;
      continue;
    } else if (j == seed_2) {
      *right++ = j;
      v->spl_nright++;
      continue;
    }

    size_alpha = distance(GETENTRY(entryvec, j), datum_l);
    size_beta = distance(GETENTRY(entryvec, j), datum_r);

    if (size_alpha < size_beta + WISH_F(v->spl_nleft, v->spl_nright, 0.01)) {
      IntRange *as = (IntRange *)VARDATA(datum_l),
               *bs = (IntRange *)VARDATA(GETENTRY(entryvec, j));

      adjustKey(as, bs);
      *left++ = j;
      v->spl_nleft++;
    } else {
      IntRange *as = (IntRange *)VARDATA(datum_r),
               *bs = (IntRange *)VARDATA(GETENTRY(entryvec, j));

      adjustKey(as, bs);
      *right++ = j;
      v->spl_nright++;
    }
  }
  *right = *left = FirstOffsetNumber;
  v->spl_ldatum = PointerGetDatum(datum_l);
  v->spl_rdatum = PointerGetDatum(datum_r);

  Assert(v->spl_nleft + v->spl_nright == maxoff);

  PG_RETURN_POINTER(v);
}

/*
 * Consistent function
 */

PGDLLEXPORT Datum gslfp_consistent(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gslfp_consistent);
Datum gslfp_consistent(PG_FUNCTION_ARGS) {
  GISTENTRY *entry = (GISTENTRY *)PG_GETARG_POINTER(0);
  StrategyNumber strategy = (StrategyNumber)PG_GETARG_UINT16(2);
  bool *recheck = (bool *)PG_GETARG_POINTER(4);
  bytea *key = (bytea *)DatumGetPointer(entry->key);
  CSfp data;
  int querySum, keySum, overlapUp, overlapDown;

  fcinfo->flinfo->fn_extra =
      searchSfpCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(1), NULL, &data, NULL);

  *recheck = true; /* we use signature, so it's needed to recheck */

  countLowOverlapValues(key, data, NUMRANGE, &querySum, &keySum, &overlapUp,
                        &overlapDown);

  PG_RETURN_BOOL(calcConsistency(GIST_LEAF(entry), strategy,
                                 overlapUp,   /* nCommonUp */
                                 overlapDown, /* nCommonDown */
                                 keySum,      /* nKey */
                                 querySum     /* nQuery */
                                 ));
}
