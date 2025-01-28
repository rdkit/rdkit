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

/*
 * Define the compressed Bfp datum representation (GBfp) to be used
 * as entry in the GiST index
 */

typedef struct {
  char vl_len_[4];
  uint8 flag;                        /* entry data type */
  uint8 data[FLEXIBLE_ARRAY_MEMBER]; /* leaf or inner data */
} GBfp;

typedef struct {
  uint32 weight;
  uint8 fp[FLEXIBLE_ARRAY_MEMBER]; /* indexed bfp */
} GBfpLeafData;

typedef struct {
  uint16 minWeight;
  uint16 maxWeight;
  uint8 fp[FLEXIBLE_ARRAY_MEMBER]; /* union + intersection bfps */
} GBfpInnerData;

#define INNER_KEY 0x01
#define ALL1_UNION 0x02        /* NOT USED YET */
#define ALL0_INTERSECTION 0x04 /* NOT USED YET */

#define IS_INNER_KEY(x) (((x)->flag & INNER_KEY) != 0x00)
#define IS_LEAF_KEY(x) (!IS_INNER_KEY(x))

#define GET_LEAF_DATA(x) ((GBfpLeafData *)((x)->data))
#define GET_INNER_DATA(x) ((GBfpInnerData *)((x)->data))

/* estimate the memory size of the index entries based on the
** signature size
*/
#define GBFP_LEAF_VARSIZE(x) (sizeof(GBfp) + sizeof(GBfpLeafData) + (x))
#define GBFP_INNER_VARSIZE(x) (sizeof(GBfp) + sizeof(GBfpInnerData) + 2 * (x))

#define GBFP_LEAF_SIGLEN(x) (VARSIZE(x) - sizeof(GBfp) - sizeof(GBfpLeafData))
#define GBFP_INNER_SIGLEN(x) \
  ((VARSIZE(x) - sizeof(GBfp) - sizeof(GBfpInnerData)) / 2)
#define GBFP_SIGLEN(x) \
  (IS_INNER_KEY(x) ? GBFP_INNER_SIGLEN(x) : GBFP_LEAF_SIGLEN(x))

#define GETENTRY(vec, pos) ((GBfp *)DatumGetPointer((vec)->vector[(pos)].key))

/* collect 'key' into 'result' */
static void merge_key(GBfp *result, GBfp *key);

/* estimate the "distance"/"difference" between two keys */
static int keys_distance(GBfp *v1, GBfp *v2);

/* clone 'key' into a new inner key */
static GBfp *copy_key(GBfp *key);

/*
 * Compress method
 */

PGDLLEXPORT Datum gbfp_compress(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_compress);
Datum gbfp_compress(PG_FUNCTION_ARGS) {
  GISTENTRY *entry = (GISTENTRY *)PG_GETARG_POINTER(0);
  GISTENTRY *retval;

  /*
  ** On leaf entries, the Bfp is replaced by a GBfp instance, where
  ** the fingerprint is annotated with the precomputed weight (popcount).
  */
  if (entry->leafkey) {
    Bfp *bfp;
    GBfp *gbfp;
    GBfpLeafData *data;
    int size, siglen, weight;

    bfp = DatumGetBfpP(entry->key);

    siglen = BFP_SIGLEN(bfp);
    weight = bitstringWeight(siglen, (uint8 *)VARDATA(bfp));

    retval = (GISTENTRY *)palloc(sizeof(GISTENTRY));

    size = GBFP_LEAF_VARSIZE(siglen);

    gbfp = palloc0(size);
    SET_VARSIZE(gbfp, size);
    data = GET_LEAF_DATA(gbfp);
    data->weight = weight;
    memcpy(data->fp, VARDATA(bfp), siglen);

    gistentryinit(*retval, PointerGetDatum(gbfp), entry->rel, entry->page,
                  entry->offset, false);
  }
  /* no change should be required on inner nodes */
  else {
    retval = entry;
  }

  PG_RETURN_POINTER(retval);
}

/*
 * Decompress method
 */

PGDLLEXPORT Datum gbfp_decompress(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_decompress);
Datum gbfp_decompress(PG_FUNCTION_ARGS) {
  GISTENTRY *entry = (GISTENTRY *)PG_GETARG_POINTER(0);
  GISTENTRY *retval;
  GBfp *key;

  key = (GBfp *)PG_DETOAST_DATUM(entry->key);

  if (key != (GBfp *)DatumGetPointer(entry->key)) {
    retval = (GISTENTRY *)palloc(sizeof(GISTENTRY));
    gistentryinit(*retval, PointerGetDatum(key), entry->rel, entry->page,
                  entry->offset, false);
    PG_RETURN_POINTER(retval);
  }

  PG_RETURN_POINTER(entry);
}

/*
 * Union method
 *
 * summarize the information from a set of keys into a single one
 * the result is always an inner key
 */

PGDLLEXPORT Datum gbfp_union(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_union);
Datum gbfp_union(PG_FUNCTION_ARGS) {
  GistEntryVector *entryvec = (GistEntryVector *)PG_GETARG_POINTER(0);
  int *size = (int *)PG_GETARG_POINTER(1);

  int i;
  GBfp *result, *key;

  key = GETENTRY(entryvec, 0);
  result = copy_key(key);
  *size = VARSIZE(result);

  for (i = 1; i < entryvec->n; ++i) {
    key = GETENTRY(entryvec, i);
    merge_key(result, key);
  }

  PG_RETURN_POINTER(result);
}

/*
 * Same method
 *
 * check if two entries represent the same data
 */

PGDLLEXPORT Datum gbfp_same(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_same);
Datum gbfp_same(PG_FUNCTION_ARGS) {
  GBfp *a = (GBfp *)PG_GETARG_POINTER(0);
  GBfp *b = (GBfp *)PG_GETARG_POINTER(1);
  bool *result = (bool *)PG_GETARG_POINTER(2);

  *result = (VARSIZE(a) == VARSIZE(b)) &&
            (memcmp(VARDATA(a), VARDATA(b), VARSIZE(a) - VARHDRSZ) == 0);

  PG_RETURN_POINTER(result);
}

/*
 * Penalty method
 *
 * estimate the cost of inserting a new key into an existing one
 * this latter is always an inner key
 */

PGDLLEXPORT Datum gbfp_penalty(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_penalty);
Datum gbfp_penalty(PG_FUNCTION_ARGS) {
  GISTENTRY *origentry = (GISTENTRY *)PG_GETARG_POINTER(0);
  GISTENTRY *newentry = (GISTENTRY *)PG_GETARG_POINTER(1);
  float *penalty = (float *)PG_GETARG_POINTER(2);

  GBfp *origval = (GBfp *)DatumGetPointer(origentry->key);
  GBfp *newval = (GBfp *)DatumGetPointer(newentry->key);

  Assert(IS_INNER_KEY(origval));

  *penalty = (float)keys_distance(origval, newval);

  PG_RETURN_POINTER(penalty);
}

/*
 * Picksplit method
 *
 * partition a collection of entries into two clusters
 */

typedef struct {
  OffsetNumber pos;
  int32 cost;
} SPLITCOST;

static int comparecost(const void *va, const void *vb) {
  const SPLITCOST *a = (const SPLITCOST *)va;
  const SPLITCOST *b = (const SPLITCOST *)vb;

  if (a->cost == b->cost) {
    return 0;
  }

  return (a->cost > b->cost) ? 1 : -1;
}

PGDLLEXPORT Datum gbfp_picksplit(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_picksplit);
Datum gbfp_picksplit(PG_FUNCTION_ARGS) {
  GistEntryVector *entryvec = (GistEntryVector *)PG_GETARG_POINTER(0);
  GIST_SPLITVEC *v = (GIST_SPLITVEC *)PG_GETARG_POINTER(1);

  OffsetNumber j, k;
  OffsetNumber *left, *right;
  OffsetNumber maxoff;
  int32 nbytes;

  OffsetNumber seed_1 = 0, seed_2 = 0;

  GBfp *datum_l, *datum_r;
  GBfp *gbfpk, *gbfpj;
  int siglen = 0;

  int32 size_alpha, size_beta;
  int32 size_waste, waste = -1;
  SPLITCOST *costvector;

  maxoff = entryvec->n - 1;
  nbytes = (maxoff + 2) * sizeof(OffsetNumber);

  v->spl_left = (OffsetNumber *)palloc(nbytes);
  left = v->spl_left;
  v->spl_nleft = 0;

  v->spl_right = (OffsetNumber *)palloc(nbytes);
  right = v->spl_right;
  v->spl_nright = 0;

  /*
  ** select the GBfp pair that are most dissimilar.
  **
  ** these will be the seeds of the two subsets.
  */
  for (k = FirstOffsetNumber; k < maxoff; k = OffsetNumberNext(k)) {
    gbfpk = GETENTRY(entryvec, k);
    if (siglen == 0) {
      siglen = GBFP_SIGLEN(gbfpk);
    }
    for (j = OffsetNumberNext(k); j <= maxoff; j = OffsetNumberNext(j)) {
      gbfpj = GETENTRY(entryvec, j);
      size_waste = keys_distance(gbfpk, gbfpj);
      if (size_waste > waste) {
        waste = size_waste;
        seed_1 = k;
        seed_2 = j;
      }
    }
  }

  /*
  ** initialize two empty subsets
  */

  if (seed_1 == 0 || seed_2 == 0) {
    /*
    ** all fps were identical and no waste was measured allowing the seeds to
    ** be assigned, so let's just pick the first two
    */
    seed_1 = 1;
    seed_2 = 2;
  }

  /* form initial .. */
  datum_l = copy_key(GETENTRY(entryvec, seed_1));
  datum_r = copy_key(GETENTRY(entryvec, seed_2));

  /* sort before ... */
  costvector = (SPLITCOST *)palloc(sizeof(SPLITCOST) * maxoff);
  for (j = FirstOffsetNumber; j <= maxoff; j = OffsetNumberNext(j)) {
    costvector[j - 1].pos = j;
    gbfpj = GETENTRY(entryvec, j);
    size_alpha = keys_distance(datum_l, gbfpj);
    size_beta = keys_distance(datum_r, gbfpj);
    costvector[j - 1].cost = abs(size_alpha - size_beta);
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

    gbfpj = GETENTRY(entryvec, j);

    size_alpha = keys_distance(datum_l, gbfpj);
    size_beta = keys_distance(datum_r, gbfpj);

    if ((size_alpha < size_beta) ||
        ((size_alpha == size_beta) && (v->spl_nleft < v->spl_nright))) {
      merge_key(datum_l, gbfpj);

      *left++ = j;
      v->spl_nleft++;
    } else {
      merge_key(datum_r, gbfpj);

      *right++ = j;
      v->spl_nright++;
    }
  }

  v->spl_ldatum = PointerGetDatum(datum_l);
  v->spl_rdatum = PointerGetDatum(datum_r);

  Assert(v->spl_nleft + v->spl_nright == maxoff);

  PG_RETURN_POINTER(v);
}

static bool gbfp_inner_consistent(BfpSignature *query, GBfpInnerData *keyData,
                                  int siglen, StrategyNumber strategy) {
  bool result;
  double t;
  double nCommon, nDelta;
  double nQuery = (double)query->weight;

  switch (strategy) {
    case RDKitTanimotoStrategy:
      /*
       * Nsame / (Na + Nb - Nsame)
       */
      t = getTanimotoLimit();
      /*
      ** The following inequalities hold
      **
      ** Na*t <= Nb <= Na/t
      **
      ** And for the fingerprints in key we have that
      **
      ** minWeight <= Nb <= maxWeight
      **
      ** so if (Na*t > maxWeight) or (Na/t < minWeight) this subtree can be
      ** discarded.
      */
      if ((keyData->maxWeight < t * nQuery) ||
          (nQuery < t * keyData->minWeight)) {
        result = false;
      }
      /* The key in the inner node stores the union of the fingerprints
      ** that populate the child nodes. We use this union to compute an
      ** upper bound to the similarity. If this upper bound is lower than the
      ** threashold value the subtree may be pruned.
      **
      ** T = Ncommon / (Na + Nb - Ncommon) <= Ncommon / Na
      **
      ** The key also stores the intersection of the fingerprints in the child
      ** nodes. The bits in this intersection that are not in the query provide
      ** a lower bound Ndelta to (Nb - Ncommon) and allow to refine the estimate
      ** above to
      **
      ** T <= Ncommon / (Na + Ndelta)
      */
      else {
        nCommon =
            (double)bitstringIntersectionWeight(siglen, keyData->fp, query->fp);
        nDelta = (double)bitstringDifferenceWeight(siglen, query->fp,
                                                   keyData->fp + siglen);
        result = nCommon >= t * (nQuery + nDelta);
      }
      break;
    case RDKitDiceStrategy:
      /*
       * 2 * Nsame / (Na + Nb)
       */
      t = getDiceLimit();
      nCommon =
          (double)bitstringIntersectionWeight(siglen, keyData->fp, query->fp);
      nDelta = (double)bitstringDifferenceWeight(siglen, query->fp,
                                                 keyData->fp + siglen);
      result = 2.0 * nCommon >= t * (nQuery + nCommon + nDelta);
      break;
    default:
      elog(ERROR, "Unknown strategy: %d", strategy);
  }

  PG_RETURN_BOOL(result);
}

static bool gbfp_leaf_consistent(BfpSignature *query, GBfpLeafData *keyData,
                                 int siglen, StrategyNumber strategy) {
  bool result;
  double t;
  double nCommon;
  double nQuery = (double)query->weight;
  double nKey = (double)keyData->weight;

  switch (strategy) {
    case RDKitTanimotoStrategy:
      /*
       * Nsame / (Na + Nb - Nsame)
       */
      t = getTanimotoLimit();
      if ((nKey < t * nQuery) || (nQuery < t * nKey)) {
        result = false;
      } else {
        nCommon =
            (double)bitstringIntersectionWeight(siglen, keyData->fp, query->fp);
        result = nCommon / (nKey + nQuery - nCommon) >= t;
      }
      break;
    case RDKitDiceStrategy:
      /*
       * 2 * Nsame / (Na + Nb)
       */
      t = getDiceLimit();
      nCommon =
          (double)bitstringIntersectionWeight(siglen, keyData->fp, query->fp);
      result = 2.0 * nCommon / (nKey + nQuery) >= t;
      break;
    default:
      elog(ERROR, "Unknown strategy: %d", strategy);
  }

  PG_RETURN_BOOL(result);
}

PGDLLEXPORT Datum gbfp_consistent(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_consistent);
Datum gbfp_consistent(PG_FUNCTION_ARGS) {
  GISTENTRY *entry = (GISTENTRY *)PG_GETARG_POINTER(0);
  StrategyNumber strategy = (StrategyNumber)PG_GETARG_UINT16(2);
  bool *recheck = (bool *)PG_GETARG_POINTER(4);
  bool result;

  GBfp *key;
  BfpSignature *query;
  int siglen;

  *recheck = false;

  fcinfo->flinfo->fn_extra =
      searchBfpCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(1), NULL, NULL, &query);

  siglen = VARSIZE(query) - sizeof(BfpSignature);

  key = (GBfp *)DatumGetPointer(entry->key);

  if (siglen != GBFP_SIGLEN(key)) {
    elog(ERROR, "All fingerprints should be the same length");
  }

  if (GIST_LEAF(entry)) {
    Assert(IS_LEAF_KEY(key));
    result = gbfp_leaf_consistent(query, GET_LEAF_DATA(key), siglen, strategy);
  } else {
    Assert(IS_INNER_KEY(key));
    result =
        gbfp_inner_consistent(query, GET_INNER_DATA(key), siglen, strategy);
  }

  PG_RETURN_BOOL(result);
}

static double gbfp_inner_distance(BfpSignature *query, GBfpInnerData *keyData,
                                  int siglen, StrategyNumber strategy) {
  double nDelta, nQuery, nCommon, similarity;

  nQuery = (double)query->weight;
  nCommon = (double)bitstringIntersectionWeight(siglen, keyData->fp, query->fp);
  nDelta = (double)bitstringDifferenceWeight(siglen, query->fp,
                                             keyData->fp + siglen);

  switch (strategy) {
    case RDKitOrderByTanimotoStrategy:
      /*
       * Nsame / (Na + Nb - Nsame)
       */
      similarity = nCommon / (nQuery + nDelta);
      break;
    case RDKitOrderByDiceStrategy:
      /*
       * 2 * Nsame / (Na + Nb)
       */
      similarity = 2.0 * nCommon / (nQuery + nCommon + nDelta);
      break;
    default:
      elog(ERROR, "Unknown strategy: %d", strategy);
  }

  PG_RETURN_FLOAT8(1.0 - similarity);
}

static double gbfp_leaf_distance(BfpSignature *query, GBfpLeafData *keyData,
                                 int siglen, StrategyNumber strategy) {
  double nKey, nQuery, nCommon, similarity;

  nKey = (double)keyData->weight;
  nQuery = (double)query->weight;
  nCommon = (double)bitstringIntersectionWeight(siglen, keyData->fp, query->fp);

  switch (strategy) {
    case RDKitOrderByTanimotoStrategy:
      /*
       * Nsame / (Na + Nb - Nsame)
       */
      similarity = nCommon / (nKey + nQuery - nCommon);
      break;
    case RDKitOrderByDiceStrategy:
      /*
       * 2 * Nsame / (Na + Nb)
       */
      similarity = 2.0 * nCommon / (nKey + nQuery);
      break;
    default:
      elog(ERROR, "Unknown strategy: %d", strategy);
  }

  return 1.0 - similarity;
}

PGDLLEXPORT Datum gbfp_distance(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_distance);
Datum gbfp_distance(PG_FUNCTION_ARGS) {
  GISTENTRY *entry = (GISTENTRY *)PG_GETARG_POINTER(0);
  StrategyNumber strategy = (StrategyNumber)PG_GETARG_UINT16(2);

  GBfp *key = (GBfp *)DatumGetPointer(entry->key);

  BfpSignature *query;
  int siglen;
  double distance;

  fcinfo->flinfo->fn_extra =
      searchBfpCache(fcinfo->flinfo->fn_extra, fcinfo->flinfo->fn_mcxt,
                     PG_GETARG_DATUM(1), NULL, NULL, &query);

  siglen = VARSIZE(query) - sizeof(BfpSignature);

  if (siglen != GBFP_SIGLEN(key)) {
    elog(ERROR, "All fingerprints should be the same length");
  }

  distance =
      GIST_LEAF(entry)
          ? gbfp_leaf_distance(query, GET_LEAF_DATA(key), siglen, strategy)
          : gbfp_inner_distance(query, GET_INNER_DATA(key), siglen, strategy);

  PG_RETURN_FLOAT8(distance);
}

PGDLLEXPORT Datum gbfp_fetch(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_fetch);
Datum gbfp_fetch(PG_FUNCTION_ARGS) {
  GISTENTRY *entry = (GISTENTRY *)PG_GETARG_POINTER(0);
  GBfp *gbfp = (GBfp *)PG_DETOAST_DATUM(entry->key);

  GBfpLeafData *data;

  int siglen, size;
  Bfp *bfp;
  GISTENTRY *retval;

  Assert(IS_LEAF_KEY(gbfp));

  siglen = GBFP_LEAF_SIGLEN(gbfp);
  data = GET_LEAF_DATA(gbfp);

  size = VARHDRSZ + siglen;

  bfp = palloc(size);
  SET_VARSIZE(bfp, size);
  memcpy(VARDATA(bfp), data->fp, siglen);

  retval = palloc(sizeof(GISTENTRY));

  gistentryinit(*retval, PointerGetDatum(bfp), entry->rel, entry->page,
                entry->offset, false);

  PG_RETURN_POINTER(retval);
}


/*
 * Sortsupport method
 *
 * Returns a comparator function to sort data in a way that preserves locality.
 */
static int
gbfp_cmp(Datum x, Datum y, SortSupport ssup)
{
  /* establish order between x and y */
  GBfp *gbfp1 = (GBfp *)PG_DETOAST_DATUM(x);
  Assert(IS_LEAF_KEY(gbfp1));
  GBfp *gbfp2 = (GBfp *)PG_DETOAST_DATUM(y);
  Assert(IS_LEAF_KEY(gbfp2));

  int siglen = GBFP_LEAF_SIGLEN(gbfp1);
  Assert(siglen == GBFP_LEAF_SIGLEN(gbfp2));

  int retval = bitstringGrayCmp(siglen, GET_LEAF_DATA(gbfp1)->fp, GET_LEAF_DATA(gbfp2)->fp);
  RDKIT_FREE_IF_COPY_P(gbfp1, x);
  RDKIT_FREE_IF_COPY_P(gbfp2, y);
  return retval;
}

PGDLLEXPORT Datum  gbfp_sortsupport(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_sortsupport);
Datum
gbfp_sortsupport(PG_FUNCTION_ARGS)
{
  SortSupport ssup = (SortSupport) PG_GETARG_POINTER(0);
  ssup->comparator = gbfp_cmp;
  PG_RETURN_VOID();
}


/* utility functions */

static void merge_key(GBfp *result, GBfp *key) {
  if (IS_LEAF_KEY(result)) {
    elog(ERROR, "Unexpected leaf key");
  }

  int siglen = GBFP_INNER_SIGLEN(result);
  GBfpInnerData *resultData = GET_INNER_DATA(result);

  if (IS_INNER_KEY(key)) {
    GBfpInnerData *innerData = GET_INNER_DATA(key);

    if (GBFP_INNER_SIGLEN(key) != siglen) {
      elog(ERROR, "All fingerprints should be the same length");
    }

    if (innerData->minWeight < resultData->minWeight) {
      resultData->minWeight = innerData->minWeight;
    }
    if (innerData->maxWeight > resultData->maxWeight) {
      resultData->maxWeight = innerData->maxWeight;
    }
    bitstringUnion(siglen, resultData->fp, innerData->fp);
    bitstringIntersection(siglen, resultData->fp + siglen,
                          innerData->fp + siglen);
  } else {
    GBfpLeafData *leafData = GET_LEAF_DATA(key);

    if (GBFP_LEAF_SIGLEN(key) != siglen) {
      elog(ERROR, "All fingerprints should be the same length");
    }

    if (leafData->weight < resultData->minWeight) {
      resultData->minWeight = leafData->weight;
    }
    if (leafData->weight > resultData->maxWeight) {
      resultData->maxWeight = leafData->weight;
    }
    bitstringUnion(siglen, resultData->fp, leafData->fp);
    bitstringIntersection(siglen, resultData->fp + siglen, leafData->fp);
  }
}

static int keys_distance(GBfp *v1, GBfp *v2) {
  uint8 *u1, *i1, *u2, *i2;
  int32 minw1, maxw1, minw2, maxw2;
  GBfpInnerData *innerData;
  GBfpLeafData *leafData;
  int distance;

  int siglen = GBFP_SIGLEN(v1);

  if (GBFP_SIGLEN(v2) != siglen) {
    elog(ERROR, "All fingerprints should be the same length");
  }

  if (IS_INNER_KEY(v1)) {
    innerData = GET_INNER_DATA(v1);
    u1 = innerData->fp;
    i1 = innerData->fp + siglen;
    minw1 = innerData->minWeight;
    maxw1 = innerData->maxWeight;
  } else {
    leafData = GET_LEAF_DATA(v1);
    u1 = i1 = leafData->fp;
    minw1 = maxw1 = leafData->weight;
  }

  if (IS_INNER_KEY(v2)) {
    innerData = GET_INNER_DATA(v2);
    u2 = innerData->fp;
    i2 = innerData->fp + siglen;
    minw2 = innerData->minWeight;
    maxw2 = innerData->maxWeight;
  } else {
    leafData = GET_LEAF_DATA(v2);
    u2 = i2 = leafData->fp;
    minw2 = maxw2 = leafData->weight;
  }

  distance = abs(minw1 - minw2) + abs(maxw1 - maxw2);
  distance *= siglen;

  distance += bitstringHemDistance(siglen, u1, u2);
  distance += bitstringHemDistance(siglen, i1, i2);

  return distance;
}

static GBfp *copy_inner_key(GBfp *key) {
  int size = VARSIZE(key);
  GBfp *result = palloc(size);
  memcpy(result, key, size);
  return result;
}

static GBfp *copy_leaf_key(GBfp *key) {
  GBfp *result;
  GBfpLeafData *leafData;
  GBfpInnerData *innerData;
  int siglen, size;

  siglen = GBFP_LEAF_SIGLEN(key);
  leafData = GET_LEAF_DATA(key);

  size = GBFP_INNER_VARSIZE(siglen);

  result = palloc0(size);
  SET_VARSIZE(result, size);
  result->flag = INNER_KEY;

  innerData = GET_INNER_DATA(result);
  innerData->minWeight = innerData->maxWeight = leafData->weight;
  memcpy(innerData->fp, leafData->fp, siglen);
  memcpy(innerData->fp + siglen, leafData->fp, siglen);

  return result;
}

static GBfp *copy_key(GBfp *key) {
  return IS_INNER_KEY(key) ? copy_inner_key(key) : copy_leaf_key(key);
}
