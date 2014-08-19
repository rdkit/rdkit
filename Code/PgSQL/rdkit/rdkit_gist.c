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
#include "postgres.h"
#include "fmgr.h"
#include "access/gist.h"
#include "access/skey.h"
#include "access/tuptoaster.h"
#include "utils/memutils.h"

#include "rdkit.h"

#define BITBYTE 8
#define SIGLEN(x)       (VARSIZE(x) - VARHDRSZ)
#define SIGLENBIT(x) (SIGLEN(x)*BITBYTE)  /* see makesign */
typedef char *BITVECP;
#define LOOPBYTE                                \
  for(i=0;i<SIGLEN;i++)

#define GETBYTE(x,i) ( *( (BITVECP)(x) + (int)( (i) / BITBYTE ) ) )
#define GETBITBYTE(x,i) ( ((char)(x)) >> i & 0x01 )
#define CLRBIT(x,i)   GETBYTE(x,i) &= ~( 0x01 << ( (i) % BITBYTE ) )
#define SETBIT(x,i)   GETBYTE(x,i) |=  ( 0x01 << ( (i) % BITBYTE ) )
#define GETBIT(x,i) ( (GETBYTE(x,i) >> ( (i) % BITBYTE )) & 0x01 )

#define ISALLTRUE(x)    ( VARSIZE(x) <= VARHDRSZ )

#define GETENTRY(vec,pos) ((bytea *) DatumGetPointer((vec)->vector[(pos)].key))

/* Number of one-bits in an unsigned byte */
static const uint8 number_of_ones[256] = {
  0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
  4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};

static int32
sizebitvec(bytea *b)
{
  int32    size = 0,
    i;
  unsigned char *sign = (unsigned char*)VARDATA(b);

  for(i=0; i<SIGLEN(b); i++)
    size += number_of_ones[sign[i]];

  return size;
}

/*
 * Compress/decompress
 */
static GISTENTRY*
compressAllTrue(GISTENTRY *entry) 
{
  GISTENTRY  *retval = entry;

  bytea   *b = (bytea*)DatumGetPointer(entry->key);
  unsigned char *sign = (unsigned char*)VARDATA(b);
  int i;
                

  for(i=0; i<SIGLEN(b); i++)
    if ( sign[i] != 0xff )
      return retval;

  retval = (GISTENTRY *) palloc(sizeof(GISTENTRY));
  b = palloc(VARHDRSZ);
  SET_VARSIZE(b, VARHDRSZ);

  gistentryinit(*retval, PointerGetDatum(b),
                entry->rel, entry->page,
                entry->offset, FALSE);

  return retval;
}

PG_FUNCTION_INFO_V1(gmol_compress);
Datum gmol_compress(PG_FUNCTION_ARGS);
Datum
gmol_compress(PG_FUNCTION_ARGS)
{
  GISTENTRY  *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  GISTENTRY  *retval = entry;

  if (entry->leafkey) {
    CROMol m = constructROMol(DatumGetMolP(entry->key));

    retval = (GISTENTRY *) palloc(sizeof(GISTENTRY));

    gistentryinit(*retval, PointerGetDatum(makeMolSign(m)),
                  entry->rel, entry->page,
                  entry->offset, FALSE);
    freeCROMol(m);
  }       
  else if ( !ISALLTRUE(DatumGetPointer(entry->key)) )
    {
      retval = compressAllTrue(entry);
    }
                
  PG_RETURN_POINTER(retval);
}

PG_FUNCTION_INFO_V1(gbfp_compress);
Datum gbfp_compress(PG_FUNCTION_ARGS);
Datum
gbfp_compress(PG_FUNCTION_ARGS)
{
  GISTENTRY  *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  GISTENTRY  *retval = entry;

  if (entry->leafkey) {
    MolBitmapFingerPrint fp = constructMolBitmapFingerPrint(DatumGetBitmapFingerPrintP(entry->key));

    retval = (GISTENTRY *) palloc(sizeof(GISTENTRY));

    gistentryinit(*retval, PointerGetDatum(makeSignatureBitmapFingerPrint(fp)),
                  entry->rel, entry->page,
                  entry->offset, FALSE);
    freeMolBitmapFingerPrint(fp);
  }       
  else if ( !ISALLTRUE(DatumGetPointer(entry->key)) )
    {
      retval = compressAllTrue(entry);
    }
                
  PG_RETURN_POINTER(retval);
}

PG_FUNCTION_INFO_V1(gsfp_compress);
Datum gsfp_compress(PG_FUNCTION_ARGS);
Datum
gsfp_compress(PG_FUNCTION_ARGS)
{
  GISTENTRY  *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  GISTENTRY  *retval = entry;

  if (entry->leafkey) {
    MolSparseFingerPrint fp = constructMolSparseFingerPrint(DatumGetSparseFingerPrintP(entry->key));

    retval = (GISTENTRY *) palloc(sizeof(GISTENTRY));

    gistentryinit(*retval, PointerGetDatum(makeSignatureSparseFingerPrint(fp, NUMBITS)),
                  entry->rel, entry->page,
                  entry->offset, FALSE);
    freeMolSparseFingerPrint(fp);
  }       
  else if ( !ISALLTRUE(DatumGetPointer(entry->key)) )
    {
      retval = compressAllTrue(entry);
    }
                
  PG_RETURN_POINTER(retval);
}

PG_FUNCTION_INFO_V1(gmol_decompress);
Datum gmol_decompress(PG_FUNCTION_ARGS);
Datum
gmol_decompress(PG_FUNCTION_ARGS)
{
  GISTENTRY       *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  bytea   *key =  (bytea*)DatumGetPointer(PG_DETOAST_DATUM(entry->key));

  if (key != (bytea *) DatumGetPointer(entry->key))
    {
      GISTENTRY  *retval = (GISTENTRY *) palloc(sizeof(GISTENTRY));

      gistentryinit(*retval, PointerGetDatum(key),
                    entry->rel, entry->page,
                    entry->offset, FALSE);

      PG_RETURN_POINTER(retval);
    }

  PG_RETURN_POINTER(entry);
}

PG_FUNCTION_INFO_V1(gmol_union);
Datum gmol_union(PG_FUNCTION_ARGS);
Datum
gmol_union(PG_FUNCTION_ARGS)
{
  GistEntryVector *entryvec = (GistEntryVector *) PG_GETARG_POINTER(0);
  int        *size = (int *) PG_GETARG_POINTER(1);
  int32            i,j;
  int         signlen;
  bytea      *result, *key;
  unsigned char *s, *k;

  key = GETENTRY(entryvec, 0);
  if (ISALLTRUE(key)) {
    *size = VARHDRSZ;
    result = palloc(VARHDRSZ);
    SET_VARSIZE(result, VARHDRSZ);

    PG_RETURN_POINTER(result);
  }

  signlen = SIGLEN(key);
  *size = VARHDRSZ + signlen;
  result = palloc(VARHDRSZ + signlen);
  SET_VARSIZE(result, VARHDRSZ + signlen);
  memcpy( VARDATA(result), VARDATA(key), signlen );

  s = (unsigned char *)VARDATA(result);
  for (i = 1; i < entryvec->n; i++)
    {
      key = GETENTRY(entryvec, i);
      k = (unsigned char *)VARDATA(key);

      if (ISALLTRUE(key)) {
        *size = VARHDRSZ;
        SET_VARSIZE(result, VARHDRSZ);

        PG_RETURN_POINTER(result);
      }

      if (SIGLEN(key) != signlen)
        elog(ERROR, "All fingerprints should be the same length");

      for(j=0;j<signlen;j++)
        s[j] |= k[j];
    }

  PG_RETURN_POINTER(result);
}

/*
 * Same method
 */

PG_FUNCTION_INFO_V1(gmol_same);
Datum gmol_same(PG_FUNCTION_ARGS);
Datum
gmol_same(PG_FUNCTION_ARGS)
{
  bytea   *a = (bytea*)PG_GETARG_POINTER(0);
  bytea   *b = (bytea*)PG_GETARG_POINTER(1);
  bool    *result = (bool *) PG_GETARG_POINTER(2);

  if (VARSIZE(a) != VARSIZE(b))
    {
      /* alltrue ot not */
      *result = false;
    }
  else
    {
      *result = (memcmp(VARDATA(a), VARDATA(b), VARSIZE(a) - VARHDRSZ) == 0) ? true : false;
    }

  PG_RETURN_POINTER(result);
}

/*
 * Penalty method
 */
static int
hemdistsign(bytea *a, bytea *b)
{
  unsigned i,
    dist = 0;
  unsigned char   *as = (unsigned char *)VARDATA(a),
    *bs = (unsigned char *)VARDATA(b);

  if (SIGLEN(a) != SIGLEN(b))
    elog(ERROR, "All fingerprints should be the same length");
#ifndef USE_BUILTIN_POPCOUNT
  for(i=0;i<SIGLEN(a);i++)
    {
      int diff = as[i] ^ bs[i];
      dist += number_of_ones[diff];
    }
#else
  unsigned eidx=SIGLEN(a)/sizeof(unsigned int);
  for(i=0;i<eidx;++i){
    dist += __builtin_popcount(((unsigned int *)as)[i] ^ ((unsigned int *)bs)[i]);
  }
  for(i=eidx*sizeof(unsigned);i<SIGLEN(a);++i){
    int diff = as[i] ^ bs[i];
    dist += number_of_ones[diff];
  }
#endif
  return dist;
}
                                                                                                                                         
static int
soergeldistsign(bytea *a, bytea *b) {
  if (SIGLEN(a) != SIGLEN(b))
    elog(ERROR, "All fingerprints should be the same length");
  unsigned int union_popcount=0,intersect_popcount=0;
  unsigned int i;
#ifndef USE_BUILTIN_POPCOUNT
  unsigned char   *as = (unsigned char *)VARDATA(a);
  unsigned char   *bs = (unsigned char *)VARDATA(b);
  for (i=0; i<SIGLEN(a); i++) {
    union_popcount += number_of_ones[as[i] | bs[i]];
    intersect_popcount += number_of_ones[as[i] & bs[i]];
  }
#else
  unsigned *as = (unsigned *)VARDATA(a);
  unsigned *bs = (unsigned *)VARDATA(b);
  unsigned eidx=SIGLEN(a)/sizeof(unsigned);
  for(i=0;i<eidx;++i){
    union_popcount += __builtin_popcount(as[i] | bs[i]);
    intersect_popcount += __builtin_popcount(as[i] & bs[i]);
  }
  for(i=eidx*sizeof(unsigned);i<SIGLEN(a);++i){
    union_popcount += number_of_ones[as[i] | bs[i]];
    intersect_popcount += number_of_ones[as[i] & bs[i]];
  }
#endif
  if (union_popcount == 0) {
    return 1;
  }
  return (int)floor(10000*(1.0-intersect_popcount / union_popcount));
}

static int
hemdist(bytea *a, bytea *b)
{
  if (ISALLTRUE(a))
    {
      if (ISALLTRUE(b))
        return 0;
      else
        return SIGLENBIT(b) - sizebitvec(b);
    }
  else if (ISALLTRUE(b))
    return SIGLENBIT(a) - sizebitvec(a);

  return hemdistsign(a, b);
}

static int
soergeldist(bytea *a, bytea *b) {
  if (ISALLTRUE(a)) {
    if (ISALLTRUE(b))
      return 0;
    else
      return SIGLENBIT(b) - sizebitvec(b);
  }
  else if (ISALLTRUE(b)) {
    return SIGLENBIT(a) - sizebitvec(a);
  }
  return soergeldistsign(a, b);
}

PG_FUNCTION_INFO_V1(gmol_penalty);
Datum gmol_penalty(PG_FUNCTION_ARGS);
Datum
gmol_penalty(PG_FUNCTION_ARGS)
{
  GISTENTRY  *origentry = (GISTENTRY *) PG_GETARG_POINTER(0); /* always ISSIGNKEY */
  GISTENTRY  *newentry = (GISTENTRY *) PG_GETARG_POINTER(1);
  float      *penalty = (float *) PG_GETARG_POINTER(2);
  bytea *origval = (bytea *) DatumGetPointer(origentry->key);
  bytea *newval = (bytea *) DatumGetPointer(newentry->key);

  *penalty = hemdist(origval, newval);

  PG_RETURN_POINTER(penalty);
}

PG_FUNCTION_INFO_V1(gbfp_penalty);
Datum gbfp_penalty(PG_FUNCTION_ARGS);
Datum
gbfp_penalty(PG_FUNCTION_ARGS)
{
  GISTENTRY  *origentry = (GISTENTRY *) PG_GETARG_POINTER(0); /* always ISSIGNKEY */
  GISTENTRY  *newentry = (GISTENTRY *) PG_GETARG_POINTER(1);
  float      *penalty = (float *) PG_GETARG_POINTER(2);
  bytea *origval = (bytea *) DatumGetPointer(origentry->key);
  bytea *newval = (bytea *) DatumGetPointer(newentry->key);

  *penalty = soergeldist(origval, newval);

  PG_RETURN_POINTER(penalty);
}

/*
 * Picksplit method
 */

#define WISH_F(a,b,c) (double)( -(double)(((a)-(b))*((a)-(b))*((a)-(b)))*(c) )

typedef struct
{
  OffsetNumber pos;
  int32        cost;
} SPLITCOST;

static int
comparecost(const void *va, const void *vb)
{
  SPLITCOST  *a = (SPLITCOST *) va;
  SPLITCOST  *b = (SPLITCOST *) vb;

  if (a->cost == b->cost)
    return 0;
  else
    return (a->cost > b->cost) ? 1 : -1;
}

PG_FUNCTION_INFO_V1(gmol_picksplit);
Datum gmol_picksplit(PG_FUNCTION_ARGS);
Datum
gmol_picksplit(PG_FUNCTION_ARGS)
{
  GistEntryVector *entryvec = (GistEntryVector *) PG_GETARG_POINTER(0);
  GIST_SPLITVEC *v = (GIST_SPLITVEC *) PG_GETARG_POINTER(1);
  OffsetNumber k,
    j;
  bytea                   *datum_l,
    *datum_r;
  int32        size_alpha,
    size_beta;
  int32        size_waste,
    waste = -1;
  int32        nbytes;
  OffsetNumber seed_1 = 0,
    seed_2 = 0;
  OffsetNumber *left,
    *right;
  OffsetNumber maxoff;
  int         i, signlen = 0;
  SPLITCOST  *costvector;

  maxoff = entryvec->n - 1;
  nbytes = (maxoff + 2) * sizeof(OffsetNumber);
  v->spl_left = (OffsetNumber *) palloc(nbytes);
  v->spl_right = (OffsetNumber *) palloc(nbytes);

  for (k = FirstOffsetNumber; k < maxoff; k = OffsetNumberNext(k))
    {
      if (signlen == 0)
        signlen = SIGLEN(GETENTRY(entryvec, k));
      for (j = OffsetNumberNext(k); j <= maxoff; j = OffsetNumberNext(j))
        {
          size_waste = hemdist(GETENTRY(entryvec, j), GETENTRY(entryvec, k));
          if (size_waste > waste)
            {
              waste = size_waste;
              seed_1 = k;
              seed_2 = j;
            }
        }
    }

  if (signlen == 0)
    signlen = SIGLEN(GETENTRY(entryvec, maxoff));


  left = v->spl_left;
  v->spl_nleft = 0;
  right = v->spl_right;
  v->spl_nright = 0;

  if (signlen == 0 || waste == 0)
    {
      /* all entries a alltrue  or all the same */

      for (k = FirstOffsetNumber; k <= maxoff; k = OffsetNumberNext(k))
        {
          if (k <= (maxoff - FirstOffsetNumber + 1) / 2)
            {
              v->spl_left[v->spl_nleft] = k;
              v->spl_nleft++;
            }
          else
            {
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

      Assert( v->spl_nleft + v->spl_nright == maxoff );
      PG_RETURN_POINTER(v);
    }

  if (seed_1 == 0 || seed_2 == 0)
    {
      seed_1 = 1;
      seed_2 = 2;
    }

  /* form initial .. */
  if (ISALLTRUE(GETENTRY(entryvec, seed_1)))
    {
      datum_l = palloc(VARHDRSZ);
      SET_VARSIZE(datum_l, VARHDRSZ);
    }
  else
    {
      datum_l = palloc(signlen + VARHDRSZ);
      memcpy(datum_l , GETENTRY(entryvec, seed_1) , signlen + VARHDRSZ);
    }
  if (ISALLTRUE(GETENTRY(entryvec, seed_2)))
    {
      datum_r = palloc(VARHDRSZ);
      SET_VARSIZE(datum_r, VARHDRSZ);
    }
  else
    {
      datum_r = palloc(signlen + VARHDRSZ);
      memcpy(datum_r , GETENTRY(entryvec, seed_2) , signlen + VARHDRSZ);
    }

  /* sort before ... */
  costvector = (SPLITCOST *) palloc(sizeof(SPLITCOST) * maxoff);
  for (j = FirstOffsetNumber; j <= maxoff; j = OffsetNumberNext(j))
    {
      costvector[j - 1].pos = j;
      size_alpha = hemdist(datum_l, GETENTRY(entryvec, j));
      size_beta  = hemdist(datum_r, GETENTRY(entryvec, j));
      costvector[j - 1].cost = Abs(size_alpha - size_beta);
    }
  qsort((void *) costvector, maxoff, sizeof(SPLITCOST), comparecost);

  for (k = 0; k < maxoff; k++)
    {
      j = costvector[k].pos;
      if (j == seed_1)
        {
          *left++ = j;
          v->spl_nleft++;
          continue;
        }
      else if (j == seed_2)
        {
          *right++ = j;
          v->spl_nright++;
          continue;
        }

      size_alpha = hemdist(GETENTRY(entryvec, j), datum_l);
      size_beta =  hemdist(GETENTRY(entryvec, j), datum_r);

      if (size_alpha < size_beta + WISH_F(v->spl_nleft, v->spl_nright, 0.1))
        {
          if (!ISALLTRUE(datum_l)) 
            {
              if (ISALLTRUE(GETENTRY(entryvec, j)))
                {
                  datum_l = palloc(VARHDRSZ);
                  SET_VARSIZE(datum_l, VARHDRSZ);
                }
              else
                {
                  unsigned char   *as = (unsigned char *)VARDATA(datum_l),
                    *bs = (unsigned char *)VARDATA(GETENTRY(entryvec, j));

                  for(i=0;i<signlen;i++)
                    as[i] |= bs[i];
                }
            }
          *left++ = j;
          v->spl_nleft++;
        }
      else
        {
          if (!ISALLTRUE(datum_r)) 
            {
              if (ISALLTRUE(GETENTRY(entryvec, j)))
                {
                  datum_r = palloc(VARHDRSZ);
                  SET_VARSIZE(datum_r, VARHDRSZ);
                }
              else
                {
                  unsigned char   *as = (unsigned char *)VARDATA(datum_r),
                    *bs = (unsigned char *)VARDATA(GETENTRY(entryvec, j));

                  for(i=0;i<signlen;i++)
                    as[i] |= bs[i];
                }
            }
          *right++ = j;
          v->spl_nright++;
        }
    }
  *right = *left = FirstOffsetNumber;
  v->spl_ldatum = PointerGetDatum(datum_l);
  v->spl_rdatum = PointerGetDatum(datum_r);

  Assert( v->spl_nleft + v->spl_nright == maxoff );

  PG_RETURN_POINTER(v);
}

/*
 * Consistent function
 */

PG_FUNCTION_INFO_V1(gmol_consistent);
Datum gmol_consistent(PG_FUNCTION_ARGS);
Datum
gmol_consistent(PG_FUNCTION_ARGS)
{
  GISTENTRY               *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  StrategyNumber  strategy = (StrategyNumber) PG_GETARG_UINT16(2);
  bool                    *recheck = (bool *) PG_GETARG_POINTER(4);
  bytea                   *key = (bytea*)DatumGetPointer(entry->key);
  bytea                   *query;
  bool                    res = true;

  fcinfo->flinfo->fn_extra = SearchMolCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(1), 
                                            NULL, NULL,&query);

  switch(strategy)
    {
    case RDKitContains:
      *recheck = true;

      if (!ISALLTRUE(key))
        {
          int i;
          unsigned char   *k = (unsigned char*)VARDATA(key),
            *q = (unsigned char*)VARDATA(query);

          if (SIGLEN(key) != SIGLEN(query))
            elog(ERROR, "All fingerprints should be the same length");

          for(i=0; res && i<SIGLEN(key); i++)
            if ( (k[i] & q[i]) != q[i])
              res = false;
        }
      break;
    case RDKitContained:
      *recheck = true;

      if (!ISALLTRUE(key))
        {
          int i;
          unsigned char   *k = (unsigned char*)VARDATA(key),
            *q = (unsigned char*)VARDATA(query);

          if (SIGLEN(key) != SIGLEN(query))
            elog(ERROR, "All fingerprints should be the same length");

          if ( GIST_LEAF(entry) )
            {
              for(i=0; res && i<SIGLEN(key); i++)
                if ( (k[i] & q[i]) != k[i])
                  res = false;
            }
          else
            {
              /*
               * Due to superimposed key on inner page we could only check
               * overlapping
               */
              res = false;
              for(i=0; res == false && i<SIGLEN(key); i++)
                if ( k[i] & q[i] )
                  res = true;
            }
        } 
      else if (GIST_LEAF(entry))
        {
          int i;
          unsigned char *q = (unsigned char*)VARDATA(query);

          res = true;
          for(i=0; res && i<SIGLEN(query); i++)
            if ( q[i] != 0xff )
              res = false;
        }
      break;
    case RDKitEquals:
      *recheck = true;

      if (!ISALLTRUE(key))
        {
          int i;
          unsigned char   *k = (unsigned char*)VARDATA(key),
            *q = (unsigned char*)VARDATA(query);

          if (SIGLEN(key) != SIGLEN(query))
            elog(ERROR, "All fingerprints should be the same length");

          for(i=0; res && i<SIGLEN(key); i++){
        	unsigned char temp = k[i] & q[i];
            if ( temp != q[i] || temp != k[i])
              res = false;
          }
        }
      break;
    default:
      elog(ERROR,"Unknown strategy: %d", strategy);
    }

  PG_RETURN_BOOL(res);
}

bool
calcConsistency(bool isLeaf, uint16 strategy, 
                double nCommonUp, double nCommonDown, double nKey, double nQuery)
{
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

  switch(strategy) {
  case RDKitTanimotoStrategy:
    /*
     * Nsame / (Na + Nb - Nsame)
     */
    if ( isLeaf )
      {
        if ( nCommonUp / (nKey + nQuery - nCommonUp) >= getTanimotoLimit() )
          res = true;
      }
    else
      {
        if ( nCommonUp / nQuery >= getTanimotoLimit() )
          res = true;
      }

    break;
  case RDKitDiceStrategy:
    /*
     * 2 * Nsame / (Na + Nb)
     */

    if ( isLeaf )
      {
        if ( 2.0 * nCommonUp / (nKey + nQuery) >= getDiceLimit() )
          res = true;
      }
    else
      {
        if ( 2.0 * nCommonUp / (nCommonDown + nQuery) >= getDiceLimit() )
          res = true;
      }

    break;
  default:
    elog(ERROR,"Unknown strategy: %d", strategy);
  }

  return res;
}

static bool
rdkit_consistent(GISTENTRY *entry, StrategyNumber strategy, bytea *key, bytea *query)
{
  double nCommon, nQuery, nKey = 0.0;

  if (ISALLTRUE(query))
    elog(ERROR, "Query malformed");

  /* 
   * Counts basic numbers, but don't count nKey on inner
   * page (see comments below)  
   */
  nQuery = (double)sizebitvec(query);
  if (ISALLTRUE(key)) {
    if (GIST_LEAF(entry))
      nKey = (double)SIGLENBIT(query);
    nCommon = nQuery;
  } else {
    int i, cnt = 0;
    unsigned char *pk = (unsigned char*)VARDATA(key);
    unsigned char *pq = (unsigned char*)VARDATA(query);

    if (SIGLEN(key) != SIGLEN(query))
      elog(ERROR, "All fingerprints should be the same length");

#ifndef USE_BUILTIN_POPCOUNT
    for(i=0;i<SIGLEN(key);i++)
      cnt += number_of_ones[ pk[i] & pq[i] ];
#else
    unsigned eidx=SIGLEN(key)/sizeof(unsigned int);
    for(i=0;i<eidx;++i){
      cnt += __builtin_popcount(((unsigned int *)pk)[i] & ((unsigned int *)pq)[i]);
    }
    for(i=eidx*sizeof(unsigned);i<SIGLEN(key);++i){
      cnt += number_of_ones[ pk[i] & pq[i] ];
    }
#endif      

    nCommon = (double)cnt;
    if (GIST_LEAF(entry))
      nKey = (double)sizebitvec(key);
  }

  return calcConsistency(GIST_LEAF(entry), strategy, nCommon, nCommon, nKey, nQuery);
}

#if PG_VERSION_NUM >= 90100
Datum  gbfp_distance(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gbfp_distance);

Datum
gbfp_distance(PG_FUNCTION_ARGS)
{
    GISTENTRY      *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
    // bytea          *query = PG_GETARG_DATA_TYPE_P(1);
    StrategyNumber  strategy = (StrategyNumber) PG_GETARG_UINT16(2);
    bytea          *key = (bytea*)DatumGetPointer(entry->key);

    bytea          *query;
    double          nCommon, nCommonUp, nCommonDown, nQuery, distance;
    double          nKey = 0.0;

    fcinfo->flinfo->fn_extra = SearchBitmapFPCache(
                                                   fcinfo->flinfo->fn_extra,
                                                   fcinfo->flinfo->fn_mcxt,
                                                   PG_GETARG_DATUM(1),
                                                   NULL, NULL,&query);

    if (ISALLTRUE(query))
        elog(ERROR, "Query malformed");

    /*
    * Counts basic numbers, but don't count nKey on inner
    * page (see comments below)
    */
    nQuery = (double)sizebitvec(query);
    if (ISALLTRUE(key))
        {

        if (GIST_LEAF(entry)) nKey = (double)SIGLENBIT(query);

        nCommon = nQuery;
        }
    else
        {
        int i, cnt = 0;
        unsigned char *pk = (unsigned char*)VARDATA(key),
            *pq = (unsigned char*)VARDATA(query);

        if (SIGLEN(key) != SIGLEN(query))
            elog(ERROR, "All fingerprints should be the same length");

#ifndef USE_BUILTIN_POPCOUNT
        for(i=0;i<SIGLEN(key);i++)
            cnt += number_of_ones[ pk[i] & pq[i] ];
#else
        unsigned eidx=SIGLEN(key)/sizeof(unsigned int);
        for(i=0;i<SIGLEN(key)/sizeof(unsigned int);++i){
          cnt += __builtin_popcount(((unsigned int *)pk)[i] & ((unsigned int *)pq)[i]);
        }
        for(i=eidx*sizeof(unsigned);i<SIGLEN(key);++i){
          cnt += number_of_ones[ pk[i] & pq[i] ];
        }
#endif        

        nCommon = (double)cnt;
        if (GIST_LEAF(entry))
            nKey = (double)sizebitvec(key);
        }

    nCommonUp = nCommon;
    nCommonDown = nCommon;

    switch(strategy)
    {
        case RDKitOrderByTanimotoStrategy:
        /*
        * Nsame / (Na + Nb - Nsame)
        */
        if (GIST_LEAF(entry))
        {
            distance = nCommonUp / (nKey + nQuery - nCommonUp);
        }

        else
        {
            distance = nCommonUp / nQuery;
        }

        break;

        case RDKitOrderByDiceStrategy:
        /*
        * 2 * Nsame / (Na + Nb)
        */
        if (GIST_LEAF(entry))
        {
            distance = 2.0 * nCommonUp / (nKey + nQuery);
        }

        else
        {
            distance =  2.0 * nCommonUp / (nCommonDown + nQuery);
        }

        break;

        default:
        elog(ERROR,"Unknown strategy: %d", strategy);
    }

    PG_RETURN_FLOAT8(1.0 - distance);
}
#endif

PG_FUNCTION_INFO_V1(gbfp_consistent);
Datum gbfp_consistent(PG_FUNCTION_ARGS);
Datum
gbfp_consistent(PG_FUNCTION_ARGS)
{
  GISTENTRY               *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  StrategyNumber  strategy = (StrategyNumber) PG_GETARG_UINT16(2);
  bool                    *recheck = (bool *) PG_GETARG_POINTER(4);
  bytea                   *key = (bytea*)DatumGetPointer(entry->key);
  bytea                   *query;

  fcinfo->flinfo->fn_extra = SearchBitmapFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(1), 
                                                 NULL, NULL,&query);

  *recheck = false;

  PG_RETURN_BOOL(rdkit_consistent(entry, strategy, key, query));
}

PG_FUNCTION_INFO_V1(gsfp_consistent);
Datum gsfp_consistent(PG_FUNCTION_ARGS);
Datum
gsfp_consistent(PG_FUNCTION_ARGS)
{
  GISTENTRY               *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  StrategyNumber  strategy = (StrategyNumber) PG_GETARG_UINT16(2);
  bool                    *recheck = (bool *) PG_GETARG_POINTER(4);
  bytea                   *key = (bytea*)DatumGetPointer(entry->key);
  bytea                   *query;
  MolSparseFingerPrint data;

  fcinfo->flinfo->fn_extra = SearchSparseFPCache(
                                                 fcinfo->flinfo->fn_extra,
                                                 fcinfo->flinfo->fn_mcxt,
                                                 PG_GETARG_DATUM(1), 
                                                 NULL, &data, &query);

  *recheck = true; /* we use signature, so it's needed to recheck */

  if (ISALLTRUE(key) && !GIST_LEAF(entry))
    {
      PG_RETURN_BOOL(true);
    }
  else
    {
      int sum,
        overlapSum,
        overlapN;

      countOverlapValues(
                         (ISALLTRUE(key)) ? NULL : key, data, NUMBITS,
                         &sum, &overlapSum, &overlapN
                         );
      PG_RETURN_BOOL(
                     calcConsistency(
                                     GIST_LEAF(entry), strategy,
                                     overlapSum, /* nCommonUp */
                                     overlapN, /* nCommonDown */
                                     (ISALLTRUE(key)) ? NUMBITS : sizebitvec(key),  /* nKey */
                                     sum                                                             /* nQuery */
                                     )
                     );
    }
}

PG_FUNCTION_INFO_V1(greaction_compress);
Datum greaction_compress(PG_FUNCTION_ARGS);
Datum
greaction_compress(PG_FUNCTION_ARGS)
{
  GISTENTRY  *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  GISTENTRY  *retval = entry;

  if (entry->leafkey) {
    CChemicalReaction rxn = constructChemReact(DatumGetMolP(entry->key));

    retval = (GISTENTRY *) palloc(sizeof(GISTENTRY));

    gistentryinit(*retval, PointerGetDatum(makeReactionSign(rxn)),
                  entry->rel, entry->page,
                  entry->offset, FALSE);
    freeChemReaction(rxn);
  }
  else if ( !ISALLTRUE(DatumGetPointer(entry->key)) )
    {
      retval = compressAllTrue(entry);
    }

  PG_RETURN_POINTER(retval);
}

PG_FUNCTION_INFO_V1(greaction_consistent);
Datum greaction_consistent(PG_FUNCTION_ARGS);
Datum
greaction_consistent(PG_FUNCTION_ARGS)
{
  GISTENTRY               *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
  StrategyNumber  strategy = (StrategyNumber) PG_GETARG_UINT16(2);
  bool                    *recheck = (bool *) PG_GETARG_POINTER(4);
  bytea                   *key = (bytea*)DatumGetPointer(entry->key);
  bytea                   *query;
  bool                    res = true;

  fcinfo->flinfo->fn_extra = SearchChemReactionCache(
                                            fcinfo->flinfo->fn_extra,
                                            fcinfo->flinfo->fn_mcxt,
                                            PG_GETARG_DATUM(1),
                                            NULL, NULL,&query);

  switch(strategy)
    {
    case RDKitContains:
      *recheck = true;

      if (!ISALLTRUE(key))
        {
          int i;
          unsigned char   *k = (unsigned char*)VARDATA(key),
            *q = (unsigned char*)VARDATA(query);

          if (SIGLEN(key) != SIGLEN(query))
            elog(ERROR, "All fingerprints should be the same length");

          for(i=0; res && i<SIGLEN(key); i++)
            if ( (k[i] & q[i]) != q[i])
              res = false;
        }
      break;
    case RDKitContained:
      *recheck = true;

      if (!ISALLTRUE(key))
        {
          int i;
          unsigned char   *k = (unsigned char*)VARDATA(key),
            *q = (unsigned char*)VARDATA(query);

          if (SIGLEN(key) != SIGLEN(query))
            elog(ERROR, "All fingerprints should be the same length");

          if ( GIST_LEAF(entry) )
            {
              for(i=0; res && i<SIGLEN(key); i++)
                if ( (k[i] & q[i]) != k[i])
                  res = false;
            }
          else
            {
              /*
               * Due to superimposed key on inner page we could only check
               * overlapping
               */
              res = false;
              for(i=0; res == false && i<SIGLEN(key); i++)
                if ( k[i] & q[i] )
                  res = true;
            }
        }
      else if (GIST_LEAF(entry))
        {
          int i;
          unsigned char *q = (unsigned char*)VARDATA(query);

          res = true;
          for(i=0; res && i<SIGLEN(query); i++)
            if ( q[i] != 0xff )
              res = false;
        }
      break;
    case RDKitEquals:
      *recheck = true;

      if (!ISALLTRUE(key))
        {
          int i;
          unsigned char   *k = (unsigned char*)VARDATA(key),
            *q = (unsigned char*)VARDATA(query);

          if (SIGLEN(key) != SIGLEN(query))
            elog(ERROR, "All fingerprints should be the same length");

          for(i=0; res && i<SIGLEN(key); i++){
          	unsigned char temp = k[i] & q[i];
              if ( temp != q[i] || temp != k[i])
                res = false;
            }
        }
      break;
    case RDKitSmaller:
      *recheck = false;

      if (!ISALLTRUE(key))
        {
          int i;
          unsigned char   *k = (unsigned char*)VARDATA(key),
            *q = (unsigned char*)VARDATA(query);

          if (SIGLEN(key) != SIGLEN(query))
            elog(ERROR, "All fingerprints should be the same length");

          for(i=0; res && i<SIGLEN(key); i++)
            if ( (k[i] & q[i]) != q[i])
              res = false;
        }
      break;
    case RDKitGreater:
      *recheck = false;

      if (!ISALLTRUE(key))
        {
          int i;
          unsigned char   *k = (unsigned char*)VARDATA(key),
            *q = (unsigned char*)VARDATA(query);

          if (SIGLEN(key) != SIGLEN(query))
            elog(ERROR, "All fingerprints should be the same length");

          if ( GIST_LEAF(entry) )
            {
              for(i=0; res && i<SIGLEN(key); i++)
                if ( (k[i] & q[i]) != k[i])
                  res = false;
            }
          else
            {
              /*
               * Due to superimposed key on inner page we could only check
               * overlapping
               */
              res = false;
              for(i=0; res == false && i<SIGLEN(key); i++)
                if ( k[i] & q[i] )
                  res = true;
            }
        }
      else if (GIST_LEAF(entry))
        {
          int i;
          unsigned char *q = (unsigned char*)VARDATA(query);

          res = true;
          for(i=0; res && i<SIGLEN(query); i++)
            if ( q[i] != 0xff )
              res = false;
        }
      break;
    default:
      elog(ERROR,"Unknown strategy: %d", strategy);
    }

  PG_RETURN_BOOL(res);
}



