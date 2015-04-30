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
#include "utils/memutils.h"

#define MAGICKNUMBER  0xBEEC0DED
/*
 * Deconstructed values cache
 */

typedef enum EntryKind {
  MolKind,
  BitmapFpKind,
  SparseFpKind,
  ChemReactionKind
} EntryKind;

typedef struct ValueCacheEntry {
  Datum                   toastedValue;
  EntryKind               kind;

  union {
    struct {
      Mol                             *value;
      CROMol                  mol;
    } mol;
    struct {
      BitmapFingerPrint               *value;
      MolBitmapFingerPrint    fp;
    } bitmap;
    struct {
      SparseFingerPrint               *value;
      MolSparseFingerPrint    fp;
    } sparse;
    struct {
      ChemReactionBA                  *value;
      CChemicalReaction       rxn;
    } chemreact;
  }                               detoasted;
  bytea                   *sign;
        
  struct ValueCacheEntry  *prev;
  struct ValueCacheEntry  *next;
} ValueCacheEntry;

#define NENTRIES        (16)

typedef struct ValueCache {
  uint32                          magickNumber;
  MemoryContext           ctx;
  int32                           nentries;
  ValueCacheEntry         *head;
  ValueCacheEntry         *tail;
  ValueCacheEntry*        entries[ NENTRIES ];

  void                            (*resetOrig) (MemoryContext context);
  void                                (*deleteOrig) (MemoryContext context);
} ValueCache;

/*********** Managing LRU **********/
static void 
moveFirst(ValueCache *ac, ValueCacheEntry *entry)
{
  /*
   * delete entry form a list
   */
  Assert( entry != ac->head );
  if ( entry == ac->tail )
    {
      Assert( entry->next == NULL );
      ac->tail = entry->prev;
      if ( ac->tail )
        ac->tail->next = NULL;
      else
        ac->head = NULL;
    }
  else
    {
      entry->prev->next = entry->next;
      entry->next->prev = entry->prev;
    }

  /*
   * Install into head 
   */

  Assert( ac->head != NULL );
  Assert( ac->tail != NULL );

  entry->next = ac->head;
  entry->prev = NULL;
  ac->head->prev = entry;
  ac->head = entry;
}

#define DATUMSIZE(d)    VARSIZE_ANY(DatumGetPointer(d)) 
static int
cmpDatum(Datum a, Datum b)
{
  int32   la = DATUMSIZE(a);
  int32   lb = DATUMSIZE(b);

  if ( la == lb )
    return memcmp( DatumGetPointer(a), DatumGetPointer(b), la );

  return (la > lb) ? 1 : -1;
}

static void
cleanupData(ValueCacheEntry *entry)
{
  pfree( DatumGetPointer(entry->toastedValue) );
  switch(entry->kind) 
    {
    case MolKind:
      if (entry->detoasted.mol.value) pfree( entry->detoasted.mol.value );
      if (entry->detoasted.mol.mol) freeCROMol( entry->detoasted.mol.mol );
      break;
    case BitmapFpKind:
      if (entry->detoasted.bitmap.value) pfree( entry->detoasted.bitmap.value );
      if (entry->detoasted.bitmap.fp) freeMolBitmapFingerPrint( entry->detoasted.bitmap.fp );
      break;
    case SparseFpKind:
      if (entry->detoasted.sparse.value) pfree( entry->detoasted.sparse.value );
      if (entry->detoasted.sparse.fp) freeMolSparseFingerPrint( entry->detoasted.sparse.fp );
      break;
    case ChemReactionKind:
      if (entry->detoasted.chemreact.value) pfree( entry->detoasted.chemreact.value );
      if (entry->detoasted.chemreact.rxn) freeChemReaction( entry->detoasted.chemreact.rxn );
      break;
    default:
      elog(ERROR, "Unknown kind: %d", entry->kind);
    }
  if (entry->sign) pfree( entry->sign );

  memset(entry, 0, offsetof(ValueCacheEntry, prev));
}

static void
makeEntry(ValueCache *ac, ValueCacheEntry *entry, Datum value, EntryKind kind)
{
  entry->toastedValue = (Datum)MemoryContextAlloc( ac->ctx, DATUMSIZE(value) );
  entry->kind = kind;
  memcpy( DatumGetPointer(entry->toastedValue), DatumGetPointer(value), DATUMSIZE(value) );
}

static int
cmpEntry(const void *a, const void *b)
{
  return cmpDatum( (*(ValueCacheEntry**)a)->toastedValue, (*(ValueCacheEntry**)b)->toastedValue ); 
}

/*********** Managing cache structure **********/

typedef struct CacheHolder {
  MemoryContext                   ctx;
  ValueCache                              *cache;
  struct CacheHolder *next;
} CacheHolder;

static  CacheHolder *holder = NULL;

static void
cleanupRDKitCache(MemoryContext context)
{
  CacheHolder *h = holder,
    *p = NULL;

  /*
   * Find holder and clean non-postgres values.
   * Note, one context could contains several caches
   */
  while( h ) {
    if (h->ctx == context)
      {
        if (h->cache->ctx != context || h->cache->magickNumber != MAGICKNUMBER)
          {
            elog(WARNING, "Something wrong in cleanupRDKitCache");
          }
        else
          {
            int i;

            for(i=0;i<h->cache->nentries;i++)
              cleanupData(h->cache->entries[i]);
            h->cache->nentries = 0;
          }

        /* remove current holder from list */
        if (p==NULL)
          {
            holder = h->next;
            free(h);
            h = holder;
          }
        else
          {
            p->next = h->next;
            free(h);
            h = p->next;
          }
        continue;
      }

    p = h;
    h = h->next;
  }
}

MemoryContextMethods    *methodsOrig = NULL;
MemoryContextMethods    methodsCache;

static void
resetCacheContext(MemoryContext context)
{
  cleanupRDKitCache(context);
  methodsOrig->reset(context);
}

static void
deleteCacheContext(MemoryContext context)
{
  cleanupRDKitCache(context);
#if PG_VERSION_NUM>=90000
  methodsOrig->delete_context(context);
#else
  methodsOrig->delete(context);
#endif
}

static ValueCache*
createCache(void *cache, struct MemoryContextData * ctx) 
{
  ValueCache              *ac;
  CacheHolder     *newholder;

  if (cache != NULL)
    {
      ac = (ValueCache*)cache;
      if (ac->ctx != ctx)
        elog(ERROR, "We can't use our approache with cache :(");
    }

  /*
   * We need to add cleanup data to delete and reset of out memory context,
   * for that we install new handlers, but we need to store old ones.
   * HACK!: we believe that there is only single memory context type 
   * in postgres 
   */

  /* define new methods */
  if ( !methodsOrig ) {
    methodsOrig = ctx->methods;
    methodsCache = *methodsOrig;
    methodsCache.reset = resetCacheContext;
#if PG_VERSION_NUM>=90000
    methodsCache.delete_context = deleteCacheContext;
#else
    methodsCache.delete = deleteCacheContext;
#endif
  }

  /*
   * Try to connect to existing cache!
   */
  newholder = holder;
  while(newholder) {
    if (newholder->ctx == ctx)
      {
        cache = (void*)newholder->cache;
        break;
      }
    newholder = newholder->next;
  }

  if (!cache)
    {
      /*
       * did not found a cache in current context, make new one
       */
      cache = MemoryContextAllocZero(ctx, sizeof(ValueCache));
      ac = (ValueCache*) cache;
      ac->magickNumber = MAGICKNUMBER;
      ac->ctx = ctx;

      newholder = malloc(sizeof(*newholder));
      if (!newholder)
        elog(ERROR, "Could not allocate %ld bytes", sizeof(*newholder)); 

      /* init holder */
      newholder->ctx = ctx;
      newholder->cache = ac;

      if ( !(ctx->methods == methodsOrig || ctx->methods == &methodsCache /* already used for another cache */) )
        elog(ERROR, "We can't use our approache with cache :((");       
                        
      ctx->methods = &methodsCache;

      /* store holder */
      newholder->next = holder;
      holder = newholder;
    }

  return (ValueCache*)cache;
}

/*********** SEARCHING **********/

/*
 * entry->kind == MolKind
 *      m, mol, fp, val
 * entry->kind == FpKind
 *      f, fp, val
 */

static void
fetchData(ValueCache *ac, ValueCacheEntry *entry,
          void **detoasted, void **internal, void **sign)
{
  MemoryContext   old;
  void * _tmp;
        
  switch(entry->kind) 
    {
    case MolKind:
      if (detoasted) 
        {
          if (entry->detoasted.mol.value == NULL)
            {
              Mol *detoastedMol;

              detoastedMol = DatumGetMolP(entry->toastedValue);
              entry->detoasted.mol.value = MemoryContextAlloc( ac->ctx, VARSIZE(detoastedMol));
              memcpy( entry->detoasted.mol.value, detoastedMol, VARSIZE(detoastedMol));
            }
          *detoasted = entry->detoasted.mol.value;
        }

      if (internal)
        {
          if (entry->detoasted.mol.mol == NULL)
            {
              fetchData(ac, entry, &_tmp, NULL, NULL);
              entry->detoasted.mol.mol = constructROMol(entry->detoasted.mol.value);
            }
          *internal = entry->detoasted.mol.mol;
        }

      if (sign)
        {
          if (entry->sign == NULL)
            {
              fetchData(ac, entry, NULL, &_tmp, NULL);
              old = MemoryContextSwitchTo( ac->ctx );
              entry->sign = makeMolSign(entry->detoasted.mol.mol);
              MemoryContextSwitchTo(old);
            }
          *sign = entry->sign;
        }
      break;
    case BitmapFpKind:
      if (detoasted) 
        {
          if (entry->detoasted.bitmap.value == NULL)
            {
              BitmapFingerPrint *detoastedFP;

              detoastedFP = DatumGetBitmapFingerPrintP(entry->toastedValue);
              entry->detoasted.bitmap.value = MemoryContextAlloc( ac->ctx, VARSIZE(detoastedFP));
              memcpy( entry->detoasted.bitmap.value, detoastedFP, VARSIZE(detoastedFP));
            }
          *detoasted = entry->detoasted.bitmap.value;
        }

      if (internal)
        {
          if (entry->detoasted.bitmap.fp == NULL)
            {
              fetchData(ac, entry, &_tmp, NULL, NULL);
              entry->detoasted.bitmap.fp = constructMolBitmapFingerPrint(entry->detoasted.bitmap.value);
            }
          *internal = entry->detoasted.bitmap.fp;
        }

      if (sign)
        {
          if (entry->sign == NULL)
            {
              fetchData(ac, entry, NULL, &_tmp, NULL);
              old = MemoryContextSwitchTo( ac->ctx );
              entry->sign = makeSignatureBitmapFingerPrint(entry->detoasted.bitmap.fp);
              MemoryContextSwitchTo(old);
            }
          *sign = entry->sign;
        }
      break;
      break;
    case SparseFpKind:
      if (detoasted) 
        {
          if (entry->detoasted.sparse.value == NULL)
            {
              SparseFingerPrint *detoastedFP;

              detoastedFP = DatumGetSparseFingerPrintP(entry->toastedValue);
              entry->detoasted.sparse.value = MemoryContextAlloc( ac->ctx, VARSIZE(detoastedFP));
              memcpy( entry->detoasted.sparse.value, detoastedFP, VARSIZE(detoastedFP));
            }
          *detoasted = entry->detoasted.sparse.value;
        }

      if (internal)
        {
          if (entry->detoasted.sparse.fp == NULL)
            {
              fetchData(ac, entry, &_tmp, NULL, NULL);
              entry->detoasted.sparse.fp = constructMolSparseFingerPrint(entry->detoasted.sparse.value);
            }
          *internal = entry->detoasted.sparse.fp;
        }

      if (sign)
        {
          if (entry->sign == NULL)
            {
              fetchData(ac, entry, NULL, &_tmp, NULL);
              old = MemoryContextSwitchTo( ac->ctx );
              entry->sign = makeSignatureSparseFingerPrint(entry->detoasted.sparse.fp, NUMBITS);
              MemoryContextSwitchTo(old);
            }
          *sign = entry->sign;
        }
      break;
      break;
    case ChemReactionKind:
      if (detoasted) 
        {
          if (entry->detoasted.chemreact.value == NULL)
            {
              ChemReactionBA *detoastedRxn;

              detoastedRxn = DatumGetChemReactionP(entry->toastedValue);
              entry->detoasted.chemreact.value = MemoryContextAlloc( ac->ctx, VARSIZE(detoastedRxn));
              memcpy( entry->detoasted.chemreact.value, detoastedRxn, VARSIZE(detoastedRxn));
            }
          *detoasted = entry->detoasted.chemreact.value;
        }

      if (internal)
        {
          if (entry->detoasted.chemreact.rxn == NULL)
            {
              fetchData(ac, entry, &_tmp, NULL, NULL);
              entry->detoasted.chemreact.rxn = constructChemReact(entry->detoasted.chemreact.value);
            }
          *internal = entry->detoasted.chemreact.rxn;
        }

      if (sign)
      {
    	  if (entry->sign == NULL)
    	  {
    		  fetchData(ac, entry, NULL, &_tmp, NULL);
    		  old = MemoryContextSwitchTo( ac->ctx );
    		  entry->sign = makeReactionSign(entry->detoasted.chemreact.rxn);
    		  MemoryContextSwitchTo(old);
    	  }
    	  *sign = entry->sign;
      }
      break;
    default:
      elog(ERROR, "Unknown kind: %d", entry->kind);
    }
}

static void*
SearchValueCache( void *cache, struct MemoryContextData * ctx, 
                  /*  input: */ Datum a, EntryKind kind, 
                  /* output: */ void **detoasted, void **internal, void **sign)
{
  ValueCache              *ac = createCache(cache, ctx);
  ValueCacheEntry *entry;

  /*
   * Fast check of recent used value 
   */
  if ( ac->head && cmpDatum(ac->head->toastedValue, a) == 0 )
    {
      Assert( ac->head->kind == kind );
      fetchData(ac, ac->head, detoasted, internal, sign);
      return cache;
    }

  if ( ac->head == NULL  )
    {
      ac->entries[0] = ac->head = ac->tail = MemoryContextAllocZero(ctx, sizeof(ValueCacheEntry));
      ac->nentries = 1;
      makeEntry(ac, ac->head, a, kind);
      fetchData(ac, ac->head, detoasted, internal, sign);
      return cache;
    }

  do {
    ValueCacheEntry **StopLow = ac->entries;
    ValueCacheEntry **StopHigh = ac->entries + ac->nentries;
    ValueCacheEntry **StopMiddle;
    int cmp;

    while (StopLow < StopHigh) {
      StopMiddle = StopLow + ((StopHigh - StopLow) >> 1);
      entry = *StopMiddle;
      cmp = cmpDatum(entry->toastedValue, a); 

      if ( cmp == 0 )
        {
          moveFirst(ac, entry);
          Assert( ac->head->kind == kind );
          fetchData(ac, ac->head, detoasted, internal, sign);
          return cache;
        }
      else if ( cmp < 0 )
        StopLow = StopMiddle + 1;
      else
        StopHigh = StopMiddle;
    }
  } while(0);

  /*
   * Not found 
   */

  if ( ac->nentries < NENTRIES )
    {
      entry = ac->entries[ ac->nentries ] = MemoryContextAllocZero(ctx, sizeof(ValueCacheEntry));

      /* install first */
      entry->next = ac->head;
      entry->prev = NULL;
      ac->head->prev = entry;
      ac->head = entry;

      ac->nentries ++;

      makeEntry(ac, ac->head, a, kind);
      fetchData(ac, ac->head, detoasted, internal, sign);
    } 
  else
    {
      cleanupData( ac->tail );
      moveFirst(ac, ac->tail );
      makeEntry(ac, ac->head, a, kind);
      fetchData(ac, ac->head, detoasted, internal, sign);
    }

  qsort(ac->entries, ac->nentries, sizeof(ValueCacheEntry*), cmpEntry);   
  return cache;
}

void* 
SearchMolCache(void *cache, struct MemoryContextData * ctx, Datum a, 
               Mol **m, CROMol *mol, bytea ** val)
{
  return  SearchValueCache( 
                           cache, ctx, 
                           /*  input: */ a, MolKind, 
                           /* output: */ (void**)m, (void**)mol, (void**)val);
}

void* 
SearchBitmapFPCache(void *cache, struct MemoryContextData * ctx, Datum a, 
                    BitmapFingerPrint **f, MolBitmapFingerPrint *fp, bytea **val)
{
  return  SearchValueCache( 
                           cache, ctx, 
                           /*  input: */ a, BitmapFpKind, 
                           /* output: */ (void**)f, (void**)fp, (void**)val);
}

void* 
SearchSparseFPCache(void *cache, struct MemoryContextData * ctx, Datum a, 
                    SparseFingerPrint **f, MolSparseFingerPrint *fp, bytea **val)
{
  return  SearchValueCache( 
                           cache, ctx, 
                           /*  input: */ a, SparseFpKind, 
                           /* output: */ (void**)f, (void**)fp, (void**)val);
}

void* 
SearchChemReactionCache(void *cache, struct MemoryContextData * ctx, Datum a, 
               ChemReactionBA **rxnBA, CChemicalReaction *rxn, bytea ** val)
{
  return  SearchValueCache( 
                           cache, ctx, 
                           /*  input: */ a, ChemReactionKind, 
                           /* output: */ (void**)rxnBA, (void**)rxn, (void**)val);
}
