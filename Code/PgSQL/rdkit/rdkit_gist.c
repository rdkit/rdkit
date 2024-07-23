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
#include <utils/sortsupport.h>
#include <utils/memutils.h>
#include <math.h>

#include "rdkit.h"
#include "guc.h"
#include "cache.h"
#include "bitstring.h"

#define SIGLEN(x) (VARSIZE(x) - VARHDRSZ)
#define SIGLENBIT(x) (SIGLEN(x) * 8) /* see makesign */

#define ISALLTRUE(x) (VARSIZE(x) <= VARHDRSZ)

#define GETENTRY(vec, pos) ((bytea *)DatumGetPointer((vec)->vector[(pos)].key))

/*
 * Compress/decompress
 */
static GISTENTRY *compressAllTrue(GISTENTRY *entry) {
  GISTENTRY *retval = entry;

  bytea *b = (bytea *)DatumGetPointer(entry->key);

  bool allTrue = bitstringAllTrue(SIGLEN(b), (uint8 *)VARDATA(b));

  if (!allTrue) {
    return retval;
  }

  retval = (GISTENTRY *)palloc(sizeof(GISTENTRY));
  b = palloc(VARHDRSZ);
  SET_VARSIZE(b, VARHDRSZ);

  gistentryinit(*retval, PointerGetDatum(b), entry->rel, entry->page,
                entry->offset, false);

  return retval;
}

PGDLLEXPORT Datum gmol_compress(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gmol_compress);
Datum gmol_compress(PG_FUNCTION_ARGS) {
  GISTENTRY *entry = (GISTENTRY *)PG_GETARG_POINTER(0);
  GISTENTRY *retval = entry;

  if (entry->leafkey) {
    CROMol m = constructROMol(DatumGetMolP(entry->key));

    retval = (GISTENTRY *)palloc(sizeof(GISTENTRY));

    gistentryinit(*retval, PointerGetDatum(makeMolSignature(m)), entry->rel,
                  entry->page, entry->offset, false);
    freeCROMol(m);
  } else if (!ISALLTRUE(DatumGetPointer(entry->key))) {
    retval = compressAllTrue(entry);
  }

  PG_RETURN_POINTER(retval);
}

PGDLLEXPORT Datum gsfp_compress(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gsfp_compress);
Datum gsfp_compress(PG_FUNCTION_ARGS) {
  GISTENTRY *entry = (GISTENTRY *)PG_GETARG_POINTER(0);
  GISTENTRY *retval = entry;

  if (entry->leafkey) {
    CSfp fp = constructCSfp(DatumGetSfpP(entry->key));

    retval = (GISTENTRY *)palloc(sizeof(GISTENTRY));

    gistentryinit(*retval, PointerGetDatum(makeSfpSignature(fp, NUMBITS)),
                  entry->rel, entry->page, entry->offset, false);
    freeCSfp(fp);
  } else if (!ISALLTRUE(DatumGetPointer(entry->key))) {
    retval = compressAllTrue(entry);
  }

  PG_RETURN_POINTER(retval);
}

PGDLLEXPORT Datum gmol_decompress(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gmol_decompress);
Datum gmol_decompress(PG_FUNCTION_ARGS) {
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

PGDLLEXPORT Datum gmol_union(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gmol_union);
Datum gmol_union(PG_FUNCTION_ARGS) {
  GistEntryVector *entryvec = (GistEntryVector *)PG_GETARG_POINTER(0);
  int *size = (int *)PG_GETARG_POINTER(1);

  int i, signlen;
  bytea *result, *key;
  unsigned char *s, *k;

  int numentries = entryvec->n;

  for (i = 0; i < numentries; ++i) {
    key = GETENTRY(entryvec, i);
    if (ISALLTRUE(key)) {
      *size = VARHDRSZ;
      result = palloc(VARHDRSZ);
      SET_VARSIZE(result, VARHDRSZ);
      PG_RETURN_POINTER(result);
    }
  }

  key = GETENTRY(entryvec, 0);
  signlen = SIGLEN(key);
  *size = VARHDRSZ + signlen;
  result = palloc(*size);
  SET_VARSIZE(result, *size);
  memcpy(VARDATA(result), VARDATA(key), signlen);

  s = (uint8 *)VARDATA(result);
  for (i = 1; i < entryvec->n; ++i) {
    key = GETENTRY(entryvec, i);
    k = (uint8 *)VARDATA(key);

    if (SIGLEN(key) != signlen) {
      elog(ERROR, "All fingerprints should be the same length");
    }

    bitstringUnion(signlen, s, k);
  }

  PG_RETURN_POINTER(result);
}

/*
 * Same method
 */

PGDLLEXPORT Datum gmol_same(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gmol_same);
Datum gmol_same(PG_FUNCTION_ARGS) {
  bytea *a = (bytea *)PG_GETARG_POINTER(0);
  bytea *b = (bytea *)PG_GETARG_POINTER(1);
  bool *result = (bool *)PG_GETARG_POINTER(2);

  *result = (VARSIZE(a) == VARSIZE(b)) /* alltrue ot not */
            && (memcmp(VARDATA(a), VARDATA(b), VARSIZE(a) - VARHDRSZ) == 0);

  PG_RETURN_POINTER(result);
}

/*
 * Penalty method
 */
static int hemdistsign(bytea *a, bytea *b) {
  int siglen = SIGLEN(a);
  uint8 *as = (uint8 *)VARDATA(a);
  uint8 *bs = (uint8 *)VARDATA(b);

  if (siglen != SIGLEN(b)) {
    elog(ERROR, "All fingerprints should be the same length");
  }

  return bitstringHemDistance(siglen, as, bs);
}

static int soergeldistsign(bytea *a, bytea *b) {
  unsigned int siglen = SIGLEN(a);

  if (siglen != SIGLEN(b)) {
    elog(ERROR, "All fingerprints should be the same length");
  }

  uint8 *as = (uint8 *)VARDATA(a);
  uint8 *bs = (uint8 *)VARDATA(b);

  return bitstringTanimotoDistance(siglen, as, bs);
}

static int hemdist(bytea *a, bytea *b) {
  if (ISALLTRUE(a)) {
    if (ISALLTRUE(b))
      return 0;
    else
      return SIGLENBIT(b) - bitstringWeight(SIGLEN(b), (uint8 *)VARDATA(b));
  } else if (ISALLTRUE(b)) {
    return SIGLENBIT(a) - bitstringWeight(SIGLEN(a), (uint8 *)VARDATA(a));
  }
  return hemdistsign(a, b);
}

static int soergeldist(bytea *a, bytea *b) {
  double d;

  if (ISALLTRUE(a)) {
    if (ISALLTRUE(b))
      return 0;
    else
      // FIXME shouldn't it be double(sizebitvec(b))/SIGLENBIT(b); ?
      return SIGLENBIT(b) - bitstringWeight(SIGLEN(b), (uint8 *)VARDATA(b));
  } else if (ISALLTRUE(b)) {
    // FIXME shouldn't it be double(sizebitvec(a))/SIGLENBIT(a); ?
    return SIGLENBIT(a) - bitstringWeight(SIGLEN(a), (uint8 *)VARDATA(a));
  }
  return (int)floor(10000 * soergeldistsign(a, b));
}

PGDLLEXPORT Datum gmol_penalty(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gmol_penalty);
Datum gmol_penalty(PG_FUNCTION_ARGS) {
  GISTENTRY *origentry =
      (GISTENTRY *)PG_GETARG_POINTER(0); /* always ISSIGNKEY */
  GISTENTRY *newentry = (GISTENTRY *)PG_GETARG_POINTER(1);
  float *penalty = (float *)PG_GETARG_POINTER(2);
  bytea *origval = (bytea *)DatumGetPointer(origentry->key);
  bytea *newval = (bytea *)DatumGetPointer(newentry->key);

  *penalty = (float)hemdist(origval, newval);

  PG_RETURN_POINTER(penalty);
}

/*
 * Picksplit method
 */

#define WISH_F(a, b, c) \
  (double)(-(double)(((a) - (b)) * ((a) - (b)) * ((a) - (b))) * (c))

typedef struct {
  OffsetNumber pos;
  int32 cost;
} SPLITCOST;

static int comparecost(const void *va, const void *vb) {
  SPLITCOST *a = (SPLITCOST *)va;
  SPLITCOST *b = (SPLITCOST *)vb;

  if (a->cost == b->cost) {
    return 0;
  }

  return (a->cost > b->cost) ? 1 : -1;
}

PGDLLEXPORT Datum gmol_picksplit(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gmol_picksplit);
Datum gmol_picksplit(PG_FUNCTION_ARGS) {
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
  int i, signlen = 0;
  SPLITCOST *costvector;

  maxoff = entryvec->n - 1;
  nbytes = (maxoff + 2) * sizeof(OffsetNumber);
  v->spl_left = (OffsetNumber *)palloc(nbytes);
  v->spl_right = (OffsetNumber *)palloc(nbytes);

  for (k = FirstOffsetNumber; k < maxoff; k = OffsetNumberNext(k)) {
    if (signlen == 0) {
      signlen = SIGLEN(GETENTRY(entryvec, k));
    }
    for (j = OffsetNumberNext(k); j <= maxoff; j = OffsetNumberNext(j)) {
      size_waste = hemdist(GETENTRY(entryvec, j), GETENTRY(entryvec, k));
      if (size_waste > waste) {
        waste = size_waste;
        seed_1 = k;
        seed_2 = j;
      }
    }
  }

  if (signlen == 0) {
    signlen = SIGLEN(GETENTRY(entryvec, maxoff));
  }

  left = v->spl_left;
  v->spl_nleft = 0;
  right = v->spl_right;
  v->spl_nright = 0;

  if (signlen == 0 || waste == 0) {
    /* all entries a alltrue  or all the same */

    for (k = FirstOffsetNumber; k <= maxoff; k = OffsetNumberNext(k)) {
      if (k <= (maxoff - FirstOffsetNumber + 1) / 2) {
        v->spl_left[v->spl_nleft] = k;
        v->spl_nleft++;
      } else {
        v->spl_right[v->spl_nright] = k;
        v->spl_nright++;
      }
    }

    signlen = VARSIZE(GETENTRY(entryvec, FirstOffsetNumber));

    datum_l = palloc(signlen);
    memcpy(datum_l, GETENTRY(entryvec, FirstOffsetNumber), signlen);
    v->spl_ldatum = PointerGetDatum(datum_l);
    datum_r = palloc(signlen);
    memcpy(datum_r, GETENTRY(entryvec, FirstOffsetNumber), signlen);
    v->spl_rdatum = PointerGetDatum(datum_r);

    Assert(v->spl_nleft + v->spl_nright == maxoff);
    PG_RETURN_POINTER(v);
  }

  if (seed_1 == 0 || seed_2 == 0) {
    seed_1 = 1;
    seed_2 = 2;
  }

  /* form initial .. */
  if (ISALLTRUE(GETENTRY(entryvec, seed_1))) {
    datum_l = palloc(VARHDRSZ);
    SET_VARSIZE(datum_l, VARHDRSZ);
  } else {
    datum_l = palloc(signlen + VARHDRSZ);
    memcpy(datum_l, GETENTRY(entryvec, seed_1), signlen + VARHDRSZ);
  }

  if (ISALLTRUE(GETENTRY(entryvec, seed_2))) {
    datum_r = palloc(VARHDRSZ);
    SET_VARSIZE(datum_r, VARHDRSZ);
  } else {
    datum_r = palloc(signlen + VARHDRSZ);
    memcpy(datum_r, GETENTRY(entryvec, seed_2), signlen + VARHDRSZ);
  }

  /* sort before ... */
  costvector = (SPLITCOST *)palloc(sizeof(SPLITCOST) * maxoff);
  for (j = FirstOffsetNumber; j <= maxoff; j = OffsetNumberNext(j)) {
    costvector[j - 1].pos = j;
    size_alpha = hemdist(datum_l, GETENTRY(entryvec, j));
    size_beta = hemdist(datum_r, GETENTRY(entryvec, j));
    costvector[j - 1].cost = Abs(size_alpha - size_beta);
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

    size_alpha = hemdist(GETENTRY(entryvec, j), datum_l);
    size_beta = hemdist(GETENTRY(entryvec, j), datum_r);

    if (size_alpha < size_beta + WISH_F(v->spl_nleft, v->spl_nright, 0.1)) {
      if (!ISALLTRUE(datum_l)) {
        if (ISALLTRUE(GETENTRY(entryvec, j))) {
          datum_l = palloc(VARHDRSZ);
          SET_VARSIZE(datum_l, VARHDRSZ);
        } else {
          unsigned char *as = (unsigned char *)VARDATA(datum_l),
                        *bs = (unsigned char *)VARDATA(GETENTRY(entryvec, j));

          for (i = 0; i < signlen; i++) {
            as[i] |= bs[i];
          }
        }
      }
      *left++ = j;
      v->spl_nleft++;
    } else {
      if (!ISALLTRUE(datum_r)) {
        if (ISALLTRUE(GETENTRY(entryvec, j))) {
          datum_r = palloc(VARHDRSZ);
          SET_VARSIZE(datum_r, VARHDRSZ);
        } else {
          unsigned char *as = (unsigned char *)VARDATA(datum_r),
                        *bs = (unsigned char *)VARDATA(GETENTRY(entryvec, j));

          for (i = 0; i < signlen; i++) {
            as[i] |= bs[i];
          }
        }
      }
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

PGDLLEXPORT Datum gmol_consistent(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gmol_consistent);
Datum gmol_consistent(PG_FUNCTION_ARGS) {
  GISTENTRY *entry = (GISTENTRY *)PG_GETARG_POINTER(0);
  StrategyNumber strategy = (StrategyNumber)PG_GETARG_UINT16(2);
  bool *recheck = (bool *)PG_GETARG_POINTER(4);
  bytea *key = (bytea *)DatumGetPointer(entry->key);
  bytea *query;
  bool res = true;

  int siglen = SIGLEN(key);

  fcinfo->flinfo->fn_extra =
      searchMolCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(1), NULL, NULL, &query);

  /*
  ** recheck is required for all strategies
  */
  *recheck = true;

  switch (strategy) {
    case RDKitContains:
      if (!ISALLTRUE(key)) {
        if (siglen != SIGLEN(query)) {
          elog(ERROR, "All fingerprints should be the same length");
        }

        uint8 *k = (uint8 *)VARDATA(key);
        uint8 *q = (uint8 *)VARDATA(query);

        res = bitstringContains(siglen, k, q);
      }
      break;
    case RDKitContained:
      if (!ISALLTRUE(key)) {
        if (siglen != SIGLEN(query)) {
          elog(ERROR, "All fingerprints should be the same length");
        }

        uint8 *k = (uint8 *)VARDATA(key);
        uint8 *q = (uint8 *)VARDATA(query);

        if (GIST_LEAF(entry)) {
          res = bitstringContains(siglen, q, k);
        } else {
          /*
           * Due to superimposed key on inner page we could only check
           * overlapping
           */
          res = bitstringIntersects(siglen, q, k);
        }
      } else if (GIST_LEAF(entry)) {
        /*
         * key is all true, it may be contained in query, iff query is also
         * all true
         */
        res = bitstringAllTrue(siglen, (uint8 *)VARDATA(query));
      }
      break;
    case RDKitEquals:
      if (!ISALLTRUE(key)) {
        /*
        ** verify the necessary condition that key should contain the query
        ** (on leaf nodes, couldn't it also verify that query contains key?)
        */
        if (siglen != SIGLEN(query)) {
          elog(ERROR, "All fingerprints should be the same length");
        }

        uint8 *k = (uint8 *)VARDATA(key);
        uint8 *q = (uint8 *)VARDATA(query);

        res = bitstringContains(siglen, k, q);
      }
      break;
    default:
      elog(ERROR, "Unknown strategy: %d", strategy);
  }

  PG_RETURN_BOOL(res);
}


/*
 * Sortsupport function
 *
 * Returns a comparator function to sort data in a way that preserves locality.
 */
static int
gmol_cmp(Datum x, Datum y, SortSupport ssup)
{
  /* establish order between x and y */
  bytea *a = (bytea*)PG_DETOAST_DATUM(x);
  bytea *b = (bytea*)PG_DETOAST_DATUM(y);

  Assert(!ISALLTRUE(a));
  Assert(!ISALLTRUE(b));

  int siglen = SIGLEN(a);
  Assert(siglen == SIGLEN(b));

  int retval = bitstringGrayCmp(siglen, VARDATA(a), VARDATA(b));
  RDKIT_FREE_IF_COPY_P(a, x);
  RDKIT_FREE_IF_COPY_P(b, y);
  return retval;
}

PGDLLEXPORT Datum  gmol_sortsupport(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gmol_sortsupport);
Datum
gmol_sortsupport(PG_FUNCTION_ARGS)
{
  SortSupport ssup = (SortSupport) PG_GETARG_POINTER(0);
  ssup->comparator = gmol_cmp;
  PG_RETURN_VOID();
}

bool calcConsistency(bool isLeaf, uint16 strategy, double nCommonUp,
                     double nCommonDown, double nKey, double nQuery) {
  bool res = false;

  /*
   * We don't wish to use RDKit's functions to compute similarity
   * for two reasons:
   *  - too expensive
   *  - on inner page it will give wrong result, because inner keys
   *    are not a real fingerprints, they are OR-ed all subsequent
   *    fingerprints
   *
   *  For inner key we should compule maximum possible numerator and
   * minimum possible denominator to find an upper bound of similarity.
   *
   * nCommon and nKey could not become greater on child - so they are
   * an upper bound of number of common bits and number of set bit
   * correspondingly. Lower bound of nKey is nCommon.
   */

  switch (strategy) {
    case RDKitTanimotoStrategy:
      /*
       * Nsame / (Na + Nb - Nsame)
       */
      if (isLeaf) {
        if (nCommonUp / (nKey + nQuery - nCommonUp) >= getTanimotoLimit()) {
          res = true;
        }
      } else {
        if (nCommonUp / nQuery >= getTanimotoLimit()) {
          res = true;
        }
      }
      break;
    case RDKitDiceStrategy:
      /*
       * 2 * Nsame / (Na + Nb)
       */

      if (isLeaf) {
        if (2.0 * nCommonUp / (nKey + nQuery) >= getDiceLimit()) {
          res = true;
        }
      } else {
        if (2.0 * nCommonUp / (nCommonDown + nQuery) >= getDiceLimit()) {
          res = true;
        }
      }
      break;
    default:
      elog(ERROR, "Unknown strategy: %d", strategy);
  }

  return res;
}

static bool rdkit_consistent(GISTENTRY *entry, StrategyNumber strategy,
                             bytea *key, bytea *query) {
  double nCommon, nQuery, nKey = 0.0;

  if (ISALLTRUE(query)) {
    elog(ERROR, "Query malformed");
  }

  /*
   * Counts basic numbers, but don't count nKey on inner
   * page (see comments below)
   */
  int siglen = SIGLEN(query);
  uint8 *q = (uint8 *)VARDATA(query);
  nQuery = (double)bitstringWeight(siglen, q);

  if (ISALLTRUE(key)) {
    if (GIST_LEAF(entry)) {
      nKey = (double)SIGLENBIT(query);
    }
    nCommon = nQuery;
  } else {
    if (siglen != SIGLEN(key)) {
      elog(ERROR, "All fingerprints should be the same length");
    }

    uint8 *k = (uint8 *)VARDATA(key);
    nCommon = bitstringIntersectionWeight(siglen, k, q);

    if (GIST_LEAF(entry)) {
      nKey = (double)bitstringWeight(siglen, k);
    }
  }

  return calcConsistency(GIST_LEAF(entry), strategy, nCommon, nCommon, nKey,
                         nQuery);
}

PGDLLEXPORT Datum gsfp_consistent(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gsfp_consistent);
Datum gsfp_consistent(PG_FUNCTION_ARGS) {
  GISTENTRY *entry = (GISTENTRY *)PG_GETARG_POINTER(0);
  StrategyNumber strategy = (StrategyNumber)PG_GETARG_UINT16(2);
  bool *recheck = (bool *)PG_GETARG_POINTER(4);
  bytea *key = (bytea *)DatumGetPointer(entry->key);
  bytea *query;
  CSfp data;

  fcinfo->flinfo->fn_extra =
      searchSfpCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(1), NULL, &data, &query);

  *recheck = true; /* we use signature, so it's needed to recheck */

  if (ISALLTRUE(key) && !GIST_LEAF(entry)) {
    PG_RETURN_BOOL(true);
  }

  int sum, overlapSum, overlapN;
  countOverlapValues((ISALLTRUE(key)) ? NULL : key, data, NUMBITS, &sum,
                     &overlapSum, &overlapN);

  int nKey = (ISALLTRUE(key))
                 ? NUMBITS
                 : bitstringWeight(SIGLEN(key), (uint8 *)VARDATA(key));

  PG_RETURN_BOOL(calcConsistency(GIST_LEAF(entry), strategy,
                                 overlapSum, /* nCommonUp */
                                 overlapN,   /* nCommonDown */
                                 nKey, sum   /* nQuery */
                                 ));
}

PGDLLEXPORT Datum greaction_compress(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(greaction_compress);
Datum greaction_compress(PG_FUNCTION_ARGS) {
  GISTENTRY *entry = (GISTENTRY *)PG_GETARG_POINTER(0);
  GISTENTRY *retval = entry;

  if (entry->leafkey) {
    CChemicalReaction rxn = constructChemReact(DatumGetMolP(entry->key));

    retval = (GISTENTRY *)palloc(sizeof(GISTENTRY));

    gistentryinit(*retval, PointerGetDatum(makeReactionSign(rxn)), entry->rel,
                  entry->page, entry->offset, false);
    freeChemReaction(rxn);
  } else if (!ISALLTRUE(DatumGetPointer(entry->key))) {
    retval = compressAllTrue(entry);
  }

  PG_RETURN_POINTER(retval);
}

PGDLLEXPORT Datum greaction_consistent(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(greaction_consistent);
Datum greaction_consistent(PG_FUNCTION_ARGS) {
  GISTENTRY *entry = (GISTENTRY *)PG_GETARG_POINTER(0);
  StrategyNumber strategy = (StrategyNumber)PG_GETARG_UINT16(2);
  bool *recheck = (bool *)PG_GETARG_POINTER(4);
  bytea *key = (bytea *)DatumGetPointer(entry->key);
  bytea *query;
  bool res = true;

  fcinfo->flinfo->fn_extra =
      searchReactionCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                          PG_GETARG_DATUM(1), NULL, NULL, &query);
  /*
  ** RDKitContains, RDKitContained, RDKitEquals require a recheck, but
  ** it defaults to false, so that RDkitSmaller and RDKitGreater can reuse
  ** the RDKitContains, RDKitContained implementation.
  */
  *recheck = false;

  switch (strategy) {
    case RDKitContains:
      *recheck = true;
      /* fallthrough */
    case RDKitSmaller:
      if (!ISALLTRUE(key)) {
        int siglen = SIGLEN(key);

        if (siglen != SIGLEN(query)) {
          elog(ERROR, "All fingerprints should be the same length");
        }

        uint8 *k = (uint8 *)VARDATA(key);
        uint8 *q = (uint8 *)VARDATA(query);

        res = bitstringContains(siglen, k, q);
      }
      break;
    case RDKitContained:
      *recheck = true;
      /* fallthrough */
    case RDKitGreater:
      if (!ISALLTRUE(key)) {
        int siglen = SIGLEN(key);

        if (siglen != SIGLEN(query)) {
          elog(ERROR, "All fingerprints should be the same length");
        }

        uint8 *k = (uint8 *)VARDATA(key);
        uint8 *q = (uint8 *)VARDATA(query);

        if (GIST_LEAF(entry)) {
          res = bitstringContains(siglen, q, k);
        } else {
          /*
           * Due to superimposed key on inner page we could only check
           * overlapping
           */
          res = bitstringIntersects(siglen, q, k);
        }
      } else if (GIST_LEAF(entry)) {
        res = bitstringAllTrue(SIGLEN(query), (uint8 *)VARDATA(query));
      }
      break;
    case RDKitEquals:
      *recheck = true;

      if (!ISALLTRUE(key)) {
        int siglen = SIGLEN(key);

        if (siglen != SIGLEN(query)) {
          elog(ERROR, "All fingerprints should be the same length");
        }

        uint8 *k = (uint8 *)VARDATA(key);
        uint8 *q = (uint8 *)VARDATA(query);

        res = bitstringContains(siglen, k, q)
              /*
              ** the original implementation also required the query to
              ** contain the key, but (I think) this is only true on the
              ** leaves (FIXME?)
              */
              && bitstringContains(siglen, q, k);
      }
      break;
    default:
      elog(ERROR, "Unknown strategy: %d", strategy);
  }

  PG_RETURN_BOOL(res);
}
