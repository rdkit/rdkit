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

/************************************************************************/
/*                                                                      */
/*    File:           set.h                                             */
/*                                                                      */
/*    Purpose:        Defines the data types and macros to implements   */
/*                    set operations. Sets are implemented as arrays    */
/*                    of unsigned integers wherein the individual       */
/*                    bits represent the members of a set. Member       */
/*                    number start with 0!!. The functions on sets      */
/*                    work like the C "string"-library, i.e. the        */
/*                    first pointer parameter is returned as the        */
/*                    function result. This means the first set is      */
/*                    overwritten!                                      */
/*                                                                      */
/************************************************************************/

#ifndef _SET_H_
#define _SET_H_

/* Only 16 bits are used per unsigned int even for 64-bit architectures. */
#define BASE_BITS      16
#define BASE_BIT_SHIFT 4
#define BASE_BIT_MASK  0xF

typedef unsigned int set_base_t;

typedef struct BIT_SET_T
   {
      unsigned    max_member;
      set_base_t *bit_array;
   } bit_set_t;


/**********************   function prototypes   *************************/

extern
bit_set_t *NewSet(unsigned int max_member);
/*
 * Allocate a new set. The set is initialized to cover values from
 * 0 to max_member.
 */

extern
void DisposeSet(bit_set_t *set);
/*
 * Deallocate the bit set pointed to by set.
 */


extern
unsigned int MaxMember(bit_set_t *set);
/*
 * Returns the allocated size of *set
 */

extern
int IsMember(bit_set_t *set, unsigned int member);
/*
 * TRUE if member is in *set
 */

extern
int NextMember(bit_set_t *set, unsigned int start_member);
/*
 * Returns the next member of this bit set starting search with start_member.
 * (-1) is returned if no further member is found.
 */

extern
bit_set_t *PutMember(bit_set_t *set,
                 unsigned int member);
/*
 * Adds member to *set.
 */

extern
bit_set_t *ClearSet(bit_set_t *set);
/*
 * Clears all bits in *set.
 */

extern
bit_set_t *CopySet(bit_set_t *dest, bit_set_t *src);
/*
 * Copies *src onto *dest. However, it does _not_ allocate the
 * set, it only checks if destination is large enough.
 */

extern
bit_set_t *SetUnion(bit_set_t *set1, bit_set_t *set2);
/*
 * Computes the set union of *set1 and *set2. The result is a pointer
 * to the changed *set1. The result is not allocated.
 */

extern
bit_set_t *SetExclusiveUnion(bit_set_t *set1, bit_set_t *set2);
/*
 * Computes the exclusive set union of *set1 and *set2. The result is
 * a pointer to the changed *set1. The result is not allocated.
 */

extern
int SetIsEmpty(bit_set_t *set);
/*
 * Returns TRUE if there are no bits in set.
 */

extern
int IntersectionIsEmpty(bit_set_t *set1, bit_set_t *set2);
/*
 * Tests if set1 and set2 don't share any bits.
 */

extern
bit_set_t *SetIntersection(bit_set_t *set1, bit_set_t *set2);
/*
 * Computes the set intersection of *set1 and *set2. The result is a
 * pointer to the changed *set1. The result is not allocated.
 */

extern
bit_set_t *SetDifference(bit_set_t *from, bit_set_t *subtract);
/*
 * Computes the set difference of *set1 and *set2. The result is a
 * pointer to the changed *set1. The result is not allocated.
 */

extern
int Cardinality(bit_set_t *set);
/*
 * Returns the number of elements in *set.
 */

extern
int CompareSets(bit_set_t *set1, bit_set_t *set2);
/*
 * Compares the two sets and returns 0 if the two are equal.
 * It returns an integer that is negative if set1 is considered lower
 * and an number greater than zero if set1 is considered higher than set2.
 */

/**************************   macros   *****************************/

#define ISMEMBER(set, member) \
   ((set)->bit_array[(member)/BASE_BITS]&(((set_base_t)1)<<(member)%BASE_BITS))

#endif

