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
#include <windows.h>

#ifdef __cplusplus
extern "C" {
#endif

/* This set of macros is necessary to make a source file that exports */
/* the DLL entry points in uppercase on both VC++ and Borland C++     */
#define HexFingerprintFromSmiles    HEXFINGERPRINTFROMSMILES
#define HexFingerprintFromMolString HEXFINGERPRINTFROMMOLSTRING
#define CanonicalizeSmiles          CANONICALIZESMILES
#define TransformSmiles             TRANSFORMSMILES
#define MOLStringToSmiles           MOLSTRINGTOSMILES
#define SmilesMatchesQueryCT        SMILESMATCHESQUERYCT
#define SmilesSimilarity            SMILESSIMILARITY
#define HeavyAtomCount              HEAVYATOMCOUNT

#define COMPONENTS              1
#define MAIN_PARTS              2
#define NON_STEREO_MAIN_PARTS   4

int _export FAR PASCAL CanonicalizeSmiles(LPSTR smiles,
                                          LPSTR cansmi,
                                          int nsmi,
                                          unsigned int level);
/*
 * Computes the canonical form of smiles and places it in cansmi[0..nsmi-1].
 * It returns the length of the canonical smiles if it worked, the negative
 * length if the size of the buffer was too small, and 0 in case of error.
 *
 * grouping can be NONE, COMPONENTS, MAIN_PARTS, or NON_STEREO_MAIN_PARTS
 * to indicate what the level of precision for identity is.
 */

int _export FAR PASCAL HexFingerprintFromSmiles(LPSTR smiles,
                                                LPSTR hexfp,
                                                int nchar);
/*
 * Computes a hexadecimal representation of the Avalon fingerprint of smiles.
 */

#define APPLY_SHORTCUTS (1<<4)

int _export FAR PASCAL TransformSmiles(LPSTR insmiles,
                                       LPSTR outsmiles,
                                       int nsmi,
                                       unsigned int flags);
/*
 * Transforms the structure denoted by insmiles into an outsmiles applying the methods indicated by flags.
 *
 * The flag APPLY_SHORTCUTS abbreviates the common shortcuts and returns an extended SMILES format with atom texts
 * enclised in {} pairs.
 */

int _export FAR PASCAL HexFingerprintFromMOLString(LPSTR smiles,
                                                   LPSTR hexfp,
                                                   int nchar,
                                                   int asQuery);
/*
 * Computes a hexadecimal representation of the Avalon fingerprint of CT.
 */

int _export FAR PASCAL HeavyAtomCount(LPSTR smiles);
/*
 * Computes the number of heavy atoms in the molecule represented by smiles.
 */

int _export FAR PASCAL SmilesSimilarity(LPSTR smiles1,
                                        LPSTR smiles2);
/*
 * Computes the fingerprint-based similarity of the structures represented by
 * the two SMILES smiles1 and smiles2 and returns an integer percentage.
 */

int _export FAR PASCAL SmilesMatchesQueryCT(LPSTR smiles,
                                            LPSTR mol_string);
/*
 * Returns TRUE if the input smiles matches the query mol_string and false
 * otherwise.
 */

/*      Not Yet Implemented !!!
int _export FAR PASCAL MOLStringToSmiles(LPSTR mol_string,
                                         LPSTR smiles,
                                         int nsmi,
                                         int with_coordinates);
 * Converts the MDL connection table contained in mol_string to a SMILES.
 * It returns the length of the smiles if it worked, the negative
 * length if the size of the buffer was too small, and 0 in case of error.
 *
 * with_coordinates indicates if the 2D coordinates are to be appended
 * to the SMILES on output.
 */

#ifdef __cplusplus
}
#endif
