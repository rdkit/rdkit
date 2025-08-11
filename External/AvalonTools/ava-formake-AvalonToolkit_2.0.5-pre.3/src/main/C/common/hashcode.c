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

/**********************************************************************/
/*                                                                    */
/*      File:    hashcode.c                                           */
/*                                                                    */
/*      Purpose: implements a one-at-a-time hashfunction to be        */
/*               used for fingerprint hashing.                        */
/*               The algorithm is adapted from one_at_a_time()        */
/*               http://en.wikipedia.org/wiki/Jenkins_hash_function   */
/*               by Bob Jenkins. It now used 64-bit intermediates.    */
/*                                                                    */
//      History: creation:       14-May-2011                          */
//                                                                    */
/**********************************************************************/

#include "hashcode.h"

// use old hash function compatibility
// this will break with revision 36889+
#ifdef __OLDFP
uint64_t next_hash(uint64_t hash, uint64_t data)
{
    hash = (uint64_t)((uint64_t)(((hash)&(uint64_t)0x000FFFFL)*(uint64_t)(data))+(uint64_t)(((hash)&(uint64_t)0xFFF0000L)>>(uint64_t)16)+(uint64_t)7);
    return hash;
}

uint64_t hash_position(uint64_t hash, int nslots)
{
    return (hash%(uint64_t)nslots);
}
#else
uint64_t next_hash(uint64_t hash, uint64_t data)
{
    hash += data;
    hash += (uint64_t)(hash << (uint64_t)10);
    hash ^= (uint64_t)(hash >> (uint64_t)6);
    return hash;
}

uint64_t hash_position(uint64_t hash, int nslots)
{
    hash += (uint64_t)(hash << (uint64_t)3);
    hash ^= (uint64_t)(hash >> (uint64_t)11);
    hash += (uint64_t)(hash << (uint64_t)15);
    return (hash%(uint64_t)nslots);
}
#endif

/**
 * Returns a 64 bit hashcode derived from str[].
 */
uint64_t hash_string(char *str)
{
   uint64_t hash = 1001;
   for (; (*str); str++)
      hash = next_hash(hash, *str);
   return hash;
}
