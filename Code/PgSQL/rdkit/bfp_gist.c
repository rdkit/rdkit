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
#include <access/tuptoaster.h>
#include <utils/memutils.h>

#include <math.h>

#include "rdkit.h"
#include "guc.h"
#include "cache.h"
#include "bitstring.h"

/*
 * Define the compressed Bfp datum representation (GBfp)  to be used
 * as entry in the GiST index
 */

typedef struct {
  char vl_len_[4];
  uint16 minWeight;
  uint16 maxWeight;
  uint8 fp[FLEXIBLE_ARRAY_MEMBER];
} GBfp;

#define GBFP_SIGLEN(x)  ((VARSIZE(x) - sizeof(GBfp))/2)

#define GETENTRY(vec,pos) ((GBfp *) DatumGetPointer((vec)->vector[(pos)].key))

/*
 * Compress method
 */

PGDLLEXPORT Datum gbfp_compress(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_compress);
Datum
gbfp_compress(PG_FUNCTION_ARGS)
{
  GISTENTRY *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  GISTENTRY *retval;

  Bfp *bfp;
  GBfp *gbfp;
  int size, siglen, weight;
  
  /*
  ** On leaf entries, the Bfp is replaced by a GBfp instance, where
  ** the fingerprint is annotated with the precomputed weight (popcount).
  */
  if (entry->leafkey) {
    bfp = DatumGetBfpP(entry->key);
    
    siglen = BFP_SIGLEN(bfp);
    weight = bitstringWeight(siglen, (uint8 *)VARDATA(bfp));
    
    retval = (GISTENTRY *) palloc(sizeof(GISTENTRY));
    
    size = sizeof(GBfp) + 2*siglen;
    
    gbfp = palloc0(size);
    SET_VARSIZE(gbfp, size);
    gbfp->minWeight = gbfp->maxWeight = weight;
    memcpy(gbfp->fp, VARDATA(bfp), siglen);
    memcpy(gbfp->fp+siglen, VARDATA(bfp), siglen);
    
    gistentryinit(*retval, PointerGetDatum(gbfp),
		  entry->rel, entry->page,
		  entry->offset, FALSE);
  }
  /* no change should be required on internal nodes */
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
Datum
gbfp_decompress(PG_FUNCTION_ARGS)
{
  GISTENTRY *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  GISTENTRY  *retval;
  GBfp *key;

  key = (GBfp *)DatumGetPointer(PG_DETOAST_DATUM(entry->key));

  if (key != (GBfp *)DatumGetPointer(entry->key)) {
    retval = (GISTENTRY *) palloc(sizeof(GISTENTRY));
    gistentryinit(*retval, PointerGetDatum(key),
		  entry->rel, entry->page,
		  entry->offset, FALSE);
    PG_RETURN_POINTER(retval);
  }

  PG_RETURN_POINTER(entry);
}

/*
 * Union method
 */

PGDLLEXPORT Datum gbfp_union(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_union);
Datum
gbfp_union(PG_FUNCTION_ARGS) {
  GistEntryVector *entryvec = (GistEntryVector *) PG_GETARG_POINTER(0);
  int *size = (int *) PG_GETARG_POINTER(1);
  
  int i, siglen;
  GBfp *result, *key;
  
  key = GETENTRY(entryvec, 0);
  *size = VARSIZE(key);
  
  result = palloc(*size);
  memcpy(result, key, *size);

  siglen = GBFP_SIGLEN(key);
  
  for (i = 1; i < entryvec->n; ++i) {
    key = GETENTRY(entryvec, i);
    
    if (GBFP_SIGLEN(key) != siglen) {
      elog(ERROR, "All fingerprints should be the same length");
    }

    if (key->minWeight < result->minWeight) {
      result->minWeight = key->minWeight;
    }
    if (key->maxWeight > result->maxWeight) {
      result->maxWeight = key->maxWeight;
    }
    bitstringUnion(siglen, result->fp, key->fp);
    bitstringIntersection(siglen, result->fp+siglen, key->fp+siglen);
  }

  PG_RETURN_POINTER(result);
}

/*
 * Same method
 */

PGDLLEXPORT Datum gbfp_same(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_same);
Datum
gbfp_same(PG_FUNCTION_ARGS)
{
  GBfp *a = (GBfp *)PG_GETARG_POINTER(0);
  GBfp *b = (GBfp *)PG_GETARG_POINTER(1);
  bool *result = (bool *) PG_GETARG_POINTER(2);

  *result =
    (VARSIZE(a) == VARSIZE(b))
    &&
    (memcmp(VARDATA(a), VARDATA(b), VARSIZE(a) - VARHDRSZ) == 0)
    ;

  PG_RETURN_POINTER(result);
}

/*
 * Penalty method
 */

PGDLLEXPORT Datum gbfp_penalty(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_penalty);
Datum
gbfp_penalty(PG_FUNCTION_ARGS)
{
  GISTENTRY  *origentry = (GISTENTRY *) PG_GETARG_POINTER(0);
  GISTENTRY  *newentry = (GISTENTRY *) PG_GETARG_POINTER(1);
  float      *penalty = (float *) PG_GETARG_POINTER(2);
  
  GBfp *origval = (GBfp *) DatumGetPointer(origentry->key);
  GBfp *newval = (GBfp *) DatumGetPointer(newentry->key);
  
  int siglen = GBFP_SIGLEN(origval);
  
  if (GBFP_SIGLEN(newval) != siglen) {
    elog(ERROR, "All fingerprints should be the same length");
  }

  /* *penalty = 1.0f - (float) bitstringTanimotoSimilarity(siglen,
     origval->fp,
     newval->fp); */
  
  *penalty = (float)
    (bitstringHemDistance(siglen, origval->fp, newval->fp)
     +  bitstringHemDistance(siglen, origval->fp+siglen, newval->fp+siglen))
    ;
  
  PG_RETURN_POINTER(penalty);
}

/*
 * Picksplit method
 */

#define WISH_F(a,b,c) (double)( -(double)(((a)-(b))*((a)-(b))*((a)-(b)))*(c) )

typedef struct {
  OffsetNumber pos;
  int32        cost;
} SPLITCOST;

static int
comparecost(const void *va, const void *vb) {
  const SPLITCOST  *a = (const SPLITCOST *) va;
  const SPLITCOST  *b = (const SPLITCOST *) vb;

  if (a->cost == b->cost) {
    return 0;
  }
  
  return (a->cost > b->cost) ? 1 : -1; 
}

PGDLLEXPORT Datum gbfp_picksplit(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_picksplit);
Datum
gbfp_picksplit(PG_FUNCTION_ARGS)
{
  GistEntryVector *entryvec = (GistEntryVector *) PG_GETARG_POINTER(0);
  GIST_SPLITVEC *v = (GIST_SPLITVEC *) PG_GETARG_POINTER(1);

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
  
  v->spl_left = (OffsetNumber *) palloc(nbytes);
  left = v->spl_left;
  v->spl_nleft = 0;

  v->spl_right = (OffsetNumber *) palloc(nbytes);
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
      size_waste =
	(bitstringHemDistance(siglen, gbfpk->fp, gbfpj->fp)
	 + bitstringHemDistance(siglen, gbfpk->fp+siglen, gbfpj->fp+siglen));
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
    ** be assigned, so just let's pick the first two
    */
    seed_1 = 1;
    seed_2 = 2;
  }

  /* form initial .. */
  datum_l = palloc(2*siglen + sizeof(GBfp));
  memcpy(datum_l, GETENTRY(entryvec, seed_1), 2*siglen + sizeof(GBfp));
  
  datum_r = palloc(2*siglen + sizeof(GBfp));
  memcpy(datum_r, GETENTRY(entryvec, seed_2), 2*siglen + sizeof(GBfp));

  /* sort before ... */
  costvector = (SPLITCOST *) palloc(sizeof(SPLITCOST) * maxoff);
  for (j = FirstOffsetNumber; j <= maxoff; j = OffsetNumberNext(j)) {
    costvector[j - 1].pos = j;
    gbfpj = GETENTRY(entryvec, j);
    size_alpha =
      (bitstringHemDistance(siglen, datum_l->fp, gbfpj->fp)
       + bitstringHemDistance(siglen, datum_l->fp+siglen, gbfpj->fp+siglen));
    size_beta =
      (bitstringHemDistance(siglen, datum_r->fp, gbfpj->fp)
       + bitstringHemDistance(siglen, datum_r->fp+siglen, gbfpj->fp+siglen));
    costvector[j - 1].cost = Abs(size_alpha - size_beta);
  }
  qsort((void *) costvector, maxoff, sizeof(SPLITCOST), comparecost);

  for (k = 0; k < maxoff; k++) {
    j = costvector[k].pos;
    
    if (j == seed_1) {
      *left++ = j;
      v->spl_nleft++;
      continue;
    }
    else if (j == seed_2) {
      *right++ = j;
      v->spl_nright++;
      continue;
    }
    
    gbfpj = GETENTRY(entryvec, j);
      
    size_alpha =
      (bitstringHemDistance(siglen, datum_l->fp, gbfpj->fp)
       + bitstringHemDistance(siglen, datum_l->fp+siglen, gbfpj->fp+siglen));
    size_beta =
      (bitstringHemDistance(siglen, datum_r->fp, gbfpj->fp)
       + bitstringHemDistance(siglen, datum_r->fp+siglen, gbfpj->fp+siglen));
    
    if (size_alpha < size_beta + WISH_F(v->spl_nleft, v->spl_nright, 0.1)) {
      if (gbfpj->minWeight < datum_l->minWeight) {
      	datum_l->minWeight = gbfpj->minWeight;
      }
      if (gbfpj->maxWeight > datum_l->maxWeight) {
      	datum_l->maxWeight = gbfpj->maxWeight;
      }
      bitstringUnion(siglen, datum_l->fp, gbfpj->fp);
      bitstringIntersection(siglen, datum_l->fp+siglen, gbfpj->fp+siglen);
      
      *left++ = j;
      v->spl_nleft++;
    }
    else {
      if (gbfpj->minWeight < datum_r->minWeight) {
      	datum_r->minWeight = gbfpj->minWeight;
      }
      if (gbfpj->maxWeight > datum_r->maxWeight) {
      	datum_r->maxWeight = gbfpj->maxWeight;
      }
      bitstringUnion(siglen, datum_r->fp, gbfpj->fp);
      bitstringIntersection(siglen, datum_r->fp+siglen, gbfpj->fp+siglen);
      
      *right++ = j;
      v->spl_nright++;
    }
  }
  
  v->spl_ldatum = PointerGetDatum(datum_l);
  v->spl_rdatum = PointerGetDatum(datum_r);
  
  Assert( v->spl_nleft + v->spl_nright == maxoff );

  PG_RETURN_POINTER(v);
}

static bool
gbfp_internal_consistent(BfpSignature *query, GBfp *key, int siglen,
			 StrategyNumber strategy)
{
  bool result;
  double t;
  double nCommon, nDelta;
  double nQuery = (double) query->weight;

  switch(strategy) {
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
    if ((key->maxWeight < t*nQuery) || (nQuery < t*key->minWeight)) {
      result = false;
    }
    /* The key in the internal node stores the union of the fingerprints 
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
      nCommon = (double)bitstringIntersectionWeight(siglen, key->fp, query->fp);
      nDelta = (double)bitstringDifferenceWeight(siglen,
						 query->fp, key->fp+siglen);
      result = nCommon >= t*(nQuery + nDelta);
    }
    break;
  case RDKitDiceStrategy:
    /*
     * 2 * Nsame / (Na + Nb)
     */
    t = getDiceLimit();
    nCommon = (double)bitstringIntersectionWeight(siglen, key->fp, query->fp);
    nDelta = (double)bitstringDifferenceWeight(siglen,
					       query->fp, key->fp+siglen);
    result = 2.0 * nCommon >= t*(nQuery + nCommon + nDelta);
    break;
  default:
    elog(ERROR,"Unknown strategy: %d", strategy);
  }
  
  PG_RETURN_BOOL(result);
}

static bool
gbfp_leaf_consistent(BfpSignature *query, GBfp *key, int siglen,
		     StrategyNumber strategy)
{
  bool result;
  double t;
  double nCommon;
  double nQuery = (double) query->weight;
  double nKey = (double)key->minWeight; /* same as maxWeight on a leaf */

  switch(strategy) {
  case RDKitTanimotoStrategy:
    /*
     * Nsame / (Na + Nb - Nsame)
     */
    t = getTanimotoLimit();
    if ((nKey < t*nQuery) || (nQuery < t*nKey)) {
      result = false;
    }
    else {
      nCommon = (double)bitstringIntersectionWeight(siglen, key->fp, query->fp);
      result = nCommon / (nKey + nQuery - nCommon) >= t;
    }
    break;
  case RDKitDiceStrategy:
    /*
     * 2 * Nsame / (Na + Nb)
     */
    t = getDiceLimit();
    nCommon = (double)bitstringIntersectionWeight(siglen, key->fp, query->fp);
    result = 2.0 * nCommon / (nKey + nQuery) >= t;
    break;
  default:
    elog(ERROR,"Unknown strategy: %d", strategy);
  }
  
  PG_RETURN_BOOL(result);
}


PGDLLEXPORT Datum gbfp_consistent(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_consistent);
Datum
gbfp_consistent(PG_FUNCTION_ARGS)
{
  GISTENTRY *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  StrategyNumber strategy = (StrategyNumber) PG_GETARG_UINT16(2);
  bool *recheck = (bool *) PG_GETARG_POINTER(4);
  bool result;

  GBfp *key;
  BfpSignature *query;
  int siglen;
  
  *recheck = false;
  
  fcinfo->flinfo->fn_extra = searchBfpCache(
					    fcinfo->flinfo->fn_extra,
					    fcinfo->flinfo->fn_mcxt,
					    PG_GETARG_DATUM(1), 
					    NULL, NULL,&query);

  siglen = VARSIZE(query) - sizeof(BfpSignature);
  
  key = (GBfp*) DatumGetPointer(entry->key);

  if (siglen != GBFP_SIGLEN(key)) {
    elog(ERROR, "All fingerprints should be the same length");
  }

  if (GIST_LEAF(entry)) {
    result = gbfp_leaf_consistent(query, key, siglen, strategy);
  }
  else {
    result = gbfp_internal_consistent(query, key, siglen, strategy);
  }
  
  PG_RETURN_BOOL(result);
}


PGDLLEXPORT Datum  gbfp_distance(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_distance);
Datum
gbfp_distance(PG_FUNCTION_ARGS)
{
  GISTENTRY *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  StrategyNumber strategy = (StrategyNumber) PG_GETARG_UINT16(2);
  
  bool isLeaf = GIST_LEAF(entry);
  GBfp *key = (GBfp *)DatumGetPointer(entry->key);
  double nKey;

  BfpSignature *query;
  int siglen;
  uint8 *q;
  
  double nCommon, nDelta, nQuery, similarity;
  
  fcinfo->flinfo->fn_extra = searchBfpCache(
					    fcinfo->flinfo->fn_extra,
					    fcinfo->flinfo->fn_mcxt,
					    PG_GETARG_DATUM(1),
					    NULL, NULL,&query);

  siglen = VARSIZE(query) - sizeof(BfpSignature);
  q = query->fp;
  nQuery = (double)query->weight;
  
  if (siglen != GBFP_SIGLEN(key)) {
    elog(ERROR, "All fingerprints should be the same length");
  }

  nCommon = (double)bitstringIntersectionWeight(siglen, key->fp, q);
    
  switch (strategy) {
  case RDKitOrderByTanimotoStrategy:
    /*
     * Nsame / (Na + Nb - Nsame)
     */
    if (isLeaf) {
      nKey = (double)key->minWeight;
      similarity = nCommon / (nKey + nQuery - nCommon);
    }
    else {
      nDelta = (double)bitstringDifferenceWeight(siglen, q, key->fp+siglen);
      similarity = nCommon / (nQuery + nDelta);
    }
    break;
  case RDKitOrderByDiceStrategy:
    /*
     * 2 * Nsame / (Na + Nb)
     */
    if (isLeaf) {
      nKey = (double)key->minWeight;
      similarity = 2.0 * nCommon / (nKey + nQuery);
    }
    else {
      nDelta = (double)bitstringDifferenceWeight(siglen, q, key->fp+siglen);
      similarity =  2.0 * nCommon / (nQuery + nCommon + nDelta);
    }
    break;
  default:
    elog(ERROR,"Unknown strategy: %d", strategy);
  }
  
  PG_RETURN_FLOAT8(1.0 - similarity);
}


PGDLLEXPORT Datum  gbfp_fetch(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_fetch);
Datum
gbfp_fetch(PG_FUNCTION_ARGS)
{
  GISTENTRY  *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  GBfp *gbfp = (GBfp *) DatumGetPointer(entry->key);
  
  int siglen, size;
  Bfp *bfp;
  GISTENTRY  *retval;

  siglen = GBFP_SIGLEN(gbfp);
  size = VARHDRSZ + siglen;
    
  bfp = palloc(size);
  SET_VARSIZE(bfp, size);
  memcpy(VARDATA(bfp), gbfp->fp, siglen);
    
  retval = palloc(sizeof(GISTENTRY));
  
  gistentryinit(*retval, PointerGetDatum(bfp),
		entry->rel, entry->page, entry->offset, FALSE);
  
  PG_RETURN_POINTER(retval);
}

