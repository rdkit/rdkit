// $Id$
//
//  Copyright (c) 2010-2015, Novartis Institutes for BioMedical Research Inc.
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
#include <RDGeneral/export.h>
#ifndef _RDKIT_CACHE_H_
#define _RDKIT_CACHE_H_

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Cache subsystem. Molecules and fingerprints I/O is extremely expensive.
 *
 * This module provides support for caching the toasted data values that
 * functions receive as argument, together with their detoasted,
 * deserialized, or signature representations. This may help saving the time
 * that would be otherwise required to some queries to process the same
 * argument values multiple times (for example repeatedly unpickling the
 * same RDKit molecule).
 *
 * This cache mechanism also provides the functions code with a high level
 * API that encapsulates the details of the conversion from the toasted
 * PostgreSQL Datum and the different representations of the supported
 * value types.
 *
 * A cache instance is associated with a memory context where the
 * cached values are allocated, and it's therefore automatically cleared
 * when the same context is reset or destroyed. For this reason the client
 * code doesn't usually create/use this cache with the "current" memory
 * context, but rather with the longer-living context that is available
 * to the function as fcinfo->flinfo->fn_mcxt.
 */

struct MemoryContextData; /* forward declaration to prevent conflicts with C++
                           */

void *searchMolCache(void *cache, struct MemoryContextData *ctx, Datum a,
                     Mol **m, CROMol *mol, bytea **sign);

void *searchXQMolCache(void *cache, struct MemoryContextData *ctx, Datum a,
                     XQMol **m, CXQMol *mol, bytea **sign);

void *searchBfpCache(void *cache, struct MemoryContextData *ctx, Datum a,
                     Bfp **f, CBfp *fp, BfpSignature **sign);

void *searchSfpCache(void *cache, struct MemoryContextData *ctx, Datum a,
                     Sfp **f, CSfp *fp, bytea **sign);

void *searchReactionCache(void *cache, struct MemoryContextData *ctx, Datum a,
                          Reaction **r, CChemicalReaction *rxn, bytea **sign);

#ifdef __cplusplus
}
#endif
#endif
