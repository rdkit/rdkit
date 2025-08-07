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
/*    file:           symboltable.h                                     */
/*                                                                      */
/*    author:         B. Rohde                                          */
/*                                                                      */
/*    history:        02-AUG-89       Creation                          */
/*                                                                      */
/*    purpose:        Provides the header for the services dealing      */
/*                    with symbol tables. Defines the symbol table      */
/*                    entry and function prototypes.                    */
/*                                                                      */
/************************************************************************/

#ifndef _SYMBOLTABLE_H_
#define _SYMBOLTABLE_H_      1

typedef struct SYMBOL_ENTRY_T
   {
      int   symbol_id;
      char *symbol_string;
   } symbol_entry_t;

extern
int StringToId(symbol_entry_t table[], char *string);
/*
 * Searches a symbol table for entries with the character string of the
 * second parameter and returns the corresponding id-number.
 * If the string is not found, the function returns the id of the NULL
 * symbol terminating the table.
 */

extern
char * IdToString(symbol_entry_t table[], int id);
/*
 * Searches a symbol table for entries with the id-number of the
 * second parameter and returns the corresponding symbol string.
 * If the id-number is not found, the function returns (char *)NULL.
 */

#endif
