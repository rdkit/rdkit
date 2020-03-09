//
//  Copyright (c) 2016, Riccardo Vianello
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
//     * Neither the name of the authors nor the names of their contributors
//       may be used to endorse or promote products derived from this software
//       without specific prior written permission.
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

#include <access/gin.h>
#if PG_VERSION_NUM >= 90500
#include <access/stratnum.h>
#else
#include <access/skey.h>
#endif
#include <fmgr.h>

#include "rdkit.h"
#include "bitstring.h"
#include "guc.h"

static Datum *gin_bfp_extract(Bfp *bfp, int32 *nkeys) {
  Datum *keys = NULL;

  int32 weight, siglen = BFP_SIGLEN(bfp);
  uint8 *fp = (uint8 *)VARDATA(bfp);

  *nkeys = weight = bitstringWeight(siglen, fp);

  if (weight != 0) {
    int32 i, j, keycount;

    keys = palloc(sizeof(Datum) * weight);

    for (keycount = 0, i = 0; i < siglen; ++i) {
      uint8 byte = fp[i];
      for (j = 0; j < 8; ++j) {
        if (byte & 0x01) {
          int32 key = 8 * i + j;
          keys[keycount++] = Int32GetDatum(key);
        }
        byte >>= 1;
      }
    }
  }
  return keys;
}

PGDLLEXPORT Datum gin_bfp_extract_value(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gin_bfp_extract_value);
Datum gin_bfp_extract_value(PG_FUNCTION_ARGS) {
  Bfp *bfp = PG_GETARG_BFP_P(0);
  int32 *nkeys = (int32 *)PG_GETARG_POINTER(1);

  PG_RETURN_POINTER(gin_bfp_extract(bfp, nkeys));
}

PGDLLEXPORT Datum gin_bfp_extract_query(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gin_bfp_extract_query);
Datum gin_bfp_extract_query(PG_FUNCTION_ARGS) {
  Bfp *bfp = PG_GETARG_BFP_P(0);
  int32 *nkeys = (int32 *)PG_GETARG_POINTER(1);
  /* StrategyNumber strategy = PG_GETARG_UINT16(2); */
  /* bool **pmatch = (bool **) PG_GETARG_POINTER(3); */
  /* Pointer **extra_data = (Pointer **) PG_GETARG_POINTER(4); */
  /* bool **nullFlags = (bool **) PG_GETARG_POINTER(5); */
  int32 *searchMode = (int32 *)PG_GETARG_POINTER(6);

  Datum *keys = gin_bfp_extract(bfp, nkeys);

  if (*nkeys == 0) {
    *searchMode = GIN_SEARCH_MODE_ALL;
  }

  PG_RETURN_POINTER(keys);
}

PGDLLEXPORT Datum gin_bfp_consistent(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gin_bfp_consistent);
Datum gin_bfp_consistent(PG_FUNCTION_ARGS) {
  bool *check = (bool *)PG_GETARG_POINTER(0);
  StrategyNumber strategy = PG_GETARG_UINT16(1);
  /* Bfp *query = PG_GETARG_BFP_P(2); */
  int32 nkeys = PG_GETARG_INT32(3);
  /* Pointer *extra_data = (Pointer *) PG_GETARG_POINTER(4); */
  bool *recheck = (bool *)PG_GETARG_POINTER(5);
  /* Datum * queryKeys = PG_GETARG_POINTER(6); */
  /* bool *nullFlags = (bool *) PG_GETARG_POINTER(7); */

  double threshold;
  bool result;

  int32 i, nCommon = 0;
  for (i = 0; i < nkeys; ++i) {
    if (check[i] == true) {
      ++nCommon;
    }
  }

  switch (strategy) {
    case RDKitTanimotoStrategy:
      /*
       * Nsame / (Na + Nb - Nsame)
       */
      threshold = getTanimotoLimit();
      result = nCommon >= threshold * nkeys;
      break;
    case RDKitDiceStrategy:
      /*
       * 2 * Nsame / (Na + Nb)
       */
      threshold = getDiceLimit();
      result = 2.0 * nCommon >= threshold * (nCommon + nkeys);
      break;
    default:
      elog(ERROR, "Unknown strategy: %d", strategy);
  }

  *recheck = result;

  PG_RETURN_BOOL(result);
}

PGDLLEXPORT Datum gin_bfp_triconsistent(PG_FUNCTION_ARGS);
PG_FUNCTION_INFO_V1(gin_bfp_triconsistent);
Datum gin_bfp_triconsistent(PG_FUNCTION_ARGS) {
#if PG_VERSION_NUM >= 90300
  /*

   */
  GinTernaryValue *check = (GinTernaryValue *)PG_GETARG_POINTER(0);
  StrategyNumber strategy = PG_GETARG_UINT16(1);
  /* Bfp *query = PG_GETARG_BFP_P(2); */
  int32 nkeys = PG_GETARG_INT32(3);
  /* Pointer *extra_data = (Pointer *) PG_GETARG_POINTER(4); */
  /* Datum * queryKeys = PG_GETARG_POINTER(5); */
  /* bool *nullFlags = (bool *) PG_GETARG_POINTER(6); */

  double threshold;
  GinTernaryValue result = GIN_MAYBE;

  int32 i, nCommon = 0, nCommonMaybe = 0;
  for (i = 0; i < nkeys; ++i) {
    if (check[i] == GIN_TRUE) {
      ++nCommon;
      ++nCommonMaybe;
    } else if (check[i] == GIN_MAYBE) {
      ++nCommonMaybe;
    }
  }

  switch (strategy) {
    case RDKitTanimotoStrategy:
      /*
       * Nsame / (Na + Nb - Nsame)
       */
      threshold = getTanimotoLimit();
      if (nCommonMaybe < threshold * nkeys) {
        result = GIN_FALSE;
      }
      break;
    case RDKitDiceStrategy:
      /*
       * 2 * Nsame / (Na + Nb)
       */
      threshold = getDiceLimit();
      if (2.0 * nCommonMaybe < threshold * (nCommonMaybe + nkeys)) {
        result = GIN_FALSE;
      }
      break;
    default:
      elog(ERROR, "Unknown strategy: %d", strategy);
  }

  PG_RETURN_GIN_TERNARY_VALUE(result);
#endif
}
