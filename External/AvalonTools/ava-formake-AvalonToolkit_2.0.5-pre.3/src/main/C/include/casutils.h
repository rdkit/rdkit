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
/*  File:   casutils.h                                                  */
/*                                                                      */
/*  Author: B. Rohde                                                    */
/*                                                                      */
/*  Purpose:  Declares utlity functions, data structures, and macros    */
/*                needed to deal with CAS connection tables.            */
/*                                                                      */
//  History:        26-Apr-94       Start of development.               */
//                                                                      */
/************************************************************************/

struct cas_node_entry_t
   {
      char *symbol;
      short valence;
      char *chem_name;
      short hex_value;
   };

#define HCO_MIN(d) (0xF & (d))
#define HCO_MAX(d)  (0xF & ((d)>>4))

/* Labels used to indicate that a charge tautomerizes in stead of a hydrogen */
#define PLUS_TAUTOMER   0x100
#define MINUS_TAUTOMER  0x200

extern
struct cas_node_entry_t *ToNodeEntry(unsigned char hex_value);
/*
 * Translates the numeric CAS node hex_value to a pointer to the
 * corresponding entry in cas_node_table. If the node hex_value is not
 * found in the table, a pointer to the last element is returned.
 */

extern
int CASValence(char *symbol);
/*
 * Returns the preferred valence of the element *symbol in the CAS
 * convention system.
 *
 * The function returns 0 for unknown elements;
 */

extern
unsigned char etoa(register unsigned int chr);
/*
 * Converts the EBCDIC character to its ASCII equivalent. It returns
 * ('?' | 0x80) if vharacter is not convertable.
 */

