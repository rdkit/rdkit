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
/*    File:           mymacros.h                                        */
/*                                                                      */
/*    Author:         B. Rohde                                          */
//                                                                      */
//    History:        02-AUG-89       Creation                          */
/*                                                                      */
/*    Purpose:        Defines some useful macros.                       */
/*                    It is not meant to be a standard set, but         */
/*                    only to save typeing and improve readability.     */
/*                                                                      */
/************************************************************************/

#ifndef _MYMACROS_H_
#define _MYMACROS_H_ 1

#define IsNULL(p) ((char *)(p) == (char *)NULL)

#define MIN(x,y)        ((x)<(y) ? (x) : (y))
#define MAX(x,y)        ((x)>(y) ? (x) : (y))
#define ABS(x)           (((x)<0) ? ((-1)*(x)) : (x))


#define TypeAlloc(n, type) ((type *)MyCalloc(n, sizeof(type)))
#define TypeRealloc(ptr, nold, nnew, type) ((type *)MyRealloc((char *)ptr, (unsigned int)nold, (unsigned int)nnew, sizeof(type)))

#define STRING_BEGINS(str, begin) (0==strncmp(str, begin, strlen(begin)))

#define Error(message)       (fprintf(stderr, "%s\n", message), exit(1))

#define Warning(message) fprintf(stderr, "%s\n", message)

#define SAFEGUARD(strp) ((strp)?(strp):("(NULL)"))

#define COUNT_LIST(counter, list, _hp) \
   for (_hp=list, counter=0; \
        !IsNULL(_hp); \
     _hp=_hp->next, counter++)
#endif
