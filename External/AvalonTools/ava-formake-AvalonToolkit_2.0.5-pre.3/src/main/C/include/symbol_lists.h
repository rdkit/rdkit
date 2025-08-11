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
/*    File:           symbol_lists.h                                    */
/*                                                                      */
/*    Author:         B. Rohde                                          */
/*                                                                      */
/*    Purpose:        Contains the declarations of types used to        */
/*                    handle lists of chemical symbols for atoms in     */
/*                    MDL CT-files.                                     */
/*                                                                      */
/*    History:        12-Apr-1991     Creation                          */
/*                    12-May-1991     Renamed from symlist.h to         */
/*                                    symbol_lists.h                    */
/*                                                                      */
/************************************************************************/

#ifndef _SYMBOL_LIST_H_
#define _SYMBOL_LIST_H_      1

#define MAXSYMBOL       80
#define EXCLUSIVE        0
#define INCLUSIVE        1

struct symbol_list_t
   {
      int  atom;                   /* number of atom to which the list applies    */
      int  logic;                  /* logic of atom list (inclusive or exclusive) */
      char string[MAXSYMBOL+1];    /* atom symbol string (comma as separator)     */
      struct symbol_list_t *next;  /* pointer to next symbol list                 */
   };

#include "forio.h"

extern
void FreeSymbolLists(struct symbol_list_t *symbol_lists);
/*
 * Frees the storage occupied by the symbol lists *symbol_lists.
 */

extern
struct symbol_list_t *CopySymbolLists(struct symbol_list_t *symbol_lists);
/*
 * Copies the symbol lists *symbol_lists and returns a pointer to the
 * result.
 */

extern
struct symbol_list_t *ReadSymbolLists(Fortran_FILE *fp, int nlists);
/*
 * Reads nlists symbol list descriptions off the file *fp and returns
 * a pointer to the linkend descriptions.
 */

extern
void PrintSymbolLists(FILE *fp, struct symbol_list_t *symbol_lists);
/*
 * Prints the symbol lists pointed to by symbol_lists to the file *fp.
 */

#endif
