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
/*    File:           set.c                                             */
/*                                                                      */
/*    Purpose:        Implements the functions to handle set            */
/*                    allocation and set operations.                    */
/*                                                                      */
//    Histroy:        04-Jan-93       Changed to use ANSI function      */
//                                    prototypes.                       */
//                                                                      */
/************************************************************************/

#include "set.h"

#include <stdio.h>

#include "local.h"

#define FALSE 0
#define TRUE  1

bit_set_t *NewSet(unsigned int max_member)
/*
 * Allocate a new set. The set is initialized to cover values from
 * 0 to max_member.
 */
{
   bit_set_t *result;

   result = TypeAlloc(1, bit_set_t);
   result->max_member = max_member;
   result->bit_array = TypeAlloc((max_member/BASE_BITS)+1, set_base_t);
   return(result);
}

void DisposeSet(bit_set_t *set)
/*
 * Deallocate the bit set pointed to by set.
 */
{
   free((char *)set->bit_array);
   free((char *)set);
}

unsigned int MaxMember(bit_set_t *set)
/*
 * Returns the allocated size of *set
 */
{
   return(set->max_member);
}

int IsMember(bit_set_t *set, unsigned member)
/*
 * TRUE if member is in *set
 */
{
   int result;

   if (set == (bit_set_t *)NULL) {
      ShowMessage("globbered set pointer", "IsMember");
      return (FALSE);
   }
   else if (set->max_member < member  ||  set->bit_array[member/BASE_BITS] == 0)
      return (FALSE);
   else
      return (0 != (set->bit_array[member/BASE_BITS] & ((set_base_t)1<<(member%BASE_BITS))));
}

int NextMember(bit_set_t *set, unsigned start_member)
/*
 * Returns the next member of this bit set starting search with start_member.
 * (-1) is returned if no further member is found.
 */
{
   if (set == (bit_set_t *)NULL)
   {
      ShowMessage("globbered set pointer", "IsMember");
      return (-1);
   }
   for (; start_member <= set->max_member; start_member++)
      if (0 != set->bit_array[start_member/BASE_BITS]  &&
          0 != (set->bit_array[start_member/BASE_BITS] & ((set_base_t)1<<(start_member%BASE_BITS))))
            return start_member;
   return (-1);
}

bit_set_t *PutMember(bit_set_t *set, unsigned int member)
/*
 * Adds member to *set.
 */
{
   bit_set_t *sp;

   if (set == (bit_set_t *)NULL || set->bit_array == (set_base_t *)NULL)
      ShowMessage("globbered set pointer","PutMember");
   else if (set->max_member < member)
   {
      sp = NewSet(member);
      sp = CopySet(sp, set);
      DisposeSet(set);
      set = sp;
      set->bit_array[member/BASE_BITS] |= ((set_base_t)1<<(member%BASE_BITS));
   }
   else
      set->bit_array[member/BASE_BITS] |= ((set_base_t)1<<(member%BASE_BITS));
   return(set);
}

bit_set_t *ClearSet(bit_set_t *set)
/*
 * Clears all bits in *set.
 */
{
   register int i;
   register set_base_t *ip;

   if (set == (bit_set_t *)NULL)
      ShowMessage("globbered set pointer","ClearSet");
   else
      for (i=set->max_member/BASE_BITS, ip=set->bit_array;
           i>=0;
           i--, ip++)
         *ip = 0;
   return(set);
}

bit_set_t *CopySet(bit_set_t *dest, bit_set_t *src)
/*
 * Copies *src onto *dest. However, it does _not_ allocate the
 * set, it only checks if destination is large enough.
 */
{
   register int i;
   register set_base_t *ip1, *ip2;

   if (dest == (bit_set_t *)NULL || src == (bit_set_t *)NULL)
      ShowMessage("globbered set pointer","CopySet");
   else if (dest->max_member < src->max_member)
      ShowMessage("destination set size < source set size","CopySet");
   else
      for (i=src->max_member/BASE_BITS,
             ip1=dest->bit_array, ip2=src->bit_array;
           i>=0;
           i--, ip1++, ip2++)
         *ip1 = *ip2;
   return(dest);
}

int CompareSets(bit_set_t *set1, bit_set_t *set2)
/*
 * Compares the two sets and returns 0 if the two are equal.
 * It returns an integer that is negative if set1 is considered lower
 * and an number greater than zero if set1 is considered higher than set2.
 */
{
   register int i;
   int c1, c2;

   if (set1 == (bit_set_t *)NULL || set2 == (bit_set_t *)NULL)
      ShowMessage("globbered set pointer","CopySet");
   else if (set1->max_member < set2->max_member)
      ShowMessage("destination set size < source set size","CopySet");
   else
   {
      c1 = Cardinality(set1);
      c2 = Cardinality(set2);
      if (c1 != c2) return (c1-c2);
      for (i=0; i<set1->max_member/BASE_BITS; i++)
         if (set1->bit_array[i] != set2->bit_array[i])
            return (((long)set1->bit_array[i])-((long)set2->bit_array[i]));
   }
   return(0);
}

bit_set_t *SetUnion(bit_set_t *set1, bit_set_t *set2)
/*
 * Computes the set union of *set1 and *set2. The result is a pointer
 * to the changed *set1. The result is not allocated.
 */
{
   register int i;

   if (set1 == (bit_set_t *)NULL || set2 == (bit_set_t *)NULL)
      ShowMessage("globbered set pointer","SetUnion");
   else if (set1->max_member < set2->max_member)
      ShowMessage("destination set size < source set size","SetUnion");
   else
      for (i=0; i<(set2->max_member/BASE_BITS)+1; i++)
         set1->bit_array[i] |= set2->bit_array[i];
   return(set1);
}

bit_set_t *SetExclusiveUnion(bit_set_t *set1, bit_set_t *set2)
/*
 * Computes the exclusive set union of *set1 and *set2. The result is
 * a pointer to the changed *set1. The result is not allocated.
 */
{
   register int i;

   if (set1 == (bit_set_t *)NULL || set2 == (bit_set_t *)NULL)
      ShowMessage("globbered set pointer","SetExclusiveUnion");
   else if (set1->max_member < set2->max_member)
      ShowMessage("destination set size < source set size",
                  "SetExclusiveUnion");
   else
      for (i=0; i<(set2->max_member/BASE_BITS)+1; i++)
         set1->bit_array[i] ^= set2->bit_array[i];
   return(set1);
}

int SetIsEmpty(bit_set_t *set)
/*
 * Returns TRUE if there are no bits in set.
 */
{
   int i;

   for (i=0; i<(set->max_member/BASE_BITS)+1; i++)
      if (set->bit_array[i] != 0) return FALSE;
   return TRUE;
}

int IntersectionIsEmpty(bit_set_t *set1, bit_set_t *set2)
/*
 * Tests if set1 and set2 don't share any bits.
 */
{
   register int i, n;

   if (set1 == (bit_set_t *)NULL || set2 == (bit_set_t *)NULL)
      ShowMessage("globbered set pointer", "IntersectionIsEmpty");
   else
   {
      n = set1->max_member/BASE_BITS;
      if ((set2->max_member/BASE_BITS)>n) n = set2->max_member/BASE_BITS;
      n++;

      for (i=0; i<n; i++)
          if (set1->bit_array[i] & set2->bit_array[i]) return (FALSE);
   }
   return(TRUE);
}

bit_set_t *SetIntersection(bit_set_t *set1, bit_set_t *set2)
/*
 * Computes the set intersection of *set1 and *set2. The result is a
 * pointer to the changed *set1. The result is not allocated.
 */
{
   register int i;

   if (set1 == (bit_set_t *)NULL || set2 == (bit_set_t *)NULL)
      ShowMessage("globbered set pointer", "SetIntersection");
   else if (set1->max_member < set2->max_member)
      ShowMessage("destination set size < source set size","SetIntersection");
   else
   {
      for (i=0; i<(set2->max_member/BASE_BITS)+1; i++)
         set1->bit_array[i] &= set2->bit_array[i];
      for (i=(set2->max_member/BASE_BITS)+1;
           i<(set1->max_member/BASE_BITS)+1; i++)
         set1->bit_array[i] = 0;
   }
   return(set1);
}

bit_set_t *SetDifference(bit_set_t *set1, bit_set_t *set2)
/*
 * Computes the set difference of *set1 and *set2. The result is a
 * pointer to the changed *set1. The result is not allocated.
 */
{
   register int i;

   if (set1 == (bit_set_t *)NULL || set2 == (bit_set_t *)NULL)
      ShowMessage("globbered set pointer", "SetDifference");
   else if (set1->max_member < set2->max_member)
      ShowMessage("destination set size < source set size","SetDifference");
   else
      for (i=0; i<(set2->max_member/BASE_BITS)+1; i++)
         set1->bit_array[i] &= ~set2->bit_array[i];
   return(set1);
}

int Cardinality(bit_set_t *set)
/*
 * Returns the number of elements in *set.
 */
{
   register unsigned int i, n;

   for (i=0, n=0; i<=set->max_member; i++)
      if (ISMEMBER(set, i)) n++;

   return(n);
}

void FPrintSet(FILE *fp, bit_set_t *set)
/*
 * Prints the contents of *set in ASCII format.
 */
{
   register unsigned int i;

   for (i=0; i<=set->max_member; i++)
      if (ISMEMBER(set, i))
         fprintf(fp, " %d", i);
}

#ifdef SETTEST
define MAX(a,b) ((a)>(b)?(a):(b))

#ifdef __TURBOC__
#include <process.h>
unsigned _stklen = 0xFEEE;
#endif

main(int argc, char *argv[])         /* Test driver for set routines */
{
   char str[80];
   int l[10];
   int i, n;

   bit_set_t *set1, *set2, *result;

   do
   {
      gets(str);
      n = sscanf(str,"%d%d%d%d%d%d%d%d%d%d",
                     l+0,l+1,l+2,l+3,l+4,l+5,l+6,l+7,l+8,l+9);
      if (n<=0) break;
      else printf("n=%d\n",n);
      set1 = NewSet(l[n-1]);
      for (i=0; i<n; i++)
         set1 = PutMember(set1,l[i]);

      gets(str);
      n = sscanf(str,"%d%d%d%d%d%d%d%d%d%d",
                     l+0,l+1,l+2,l+3,l+4,l+5,l+6,l+7,l+8,l+9);
      if (n<=0) break;
      else printf("n=%d\n",n);
      set2 = NewSet(l[n-1]);
      for (i=0; i<n; i++)
         set2 = PutMember(set2,l[i]);

      result = NewSet(MAX(MaxMember(set1),MaxMember(set2)));
      printf("set1 : "); FPrintSet(stdout,set1); printf("\n");
      printf("set2 : "); FPrintSet(stdout,set2); printf("\n");

      printf("SetUnion(set1,set2) : ");
      FPrintSet(stdout,SetUnion(CopySet(result,set1),set2));
      printf("\n");

      printf("SetIntersection(set1,set2) : ");
      FPrintSet(stdout,SetIntersection(CopySet(result,set1),set2));
      printf("\n");

      printf("SetDifference(set1,set2) : ");
      FPrintSet(stdout,SetDifference(CopySet(result,set1),set2));
      printf("\n");
   } while(1);
}
#endif
