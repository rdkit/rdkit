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
/*    File:           mdl.h                                             */
/*                                                                      */
/*    Author:         B. Rohde                                          */
/*                                                                      */
/*    Purpose:        Contains definitions valid for all work with      */
/*                    MDL data structures.                              */
/*                                                                      */
/*    History:        11-Apr-1991     Creation                          */
/*                                                                      */
/************************************************************************/

#ifndef _MDL_H_
#define _MDL_H_              1

#include "local.h"

#define MDL_VERSION    1      /* version number of MOL file format    */
#define MDL_MAXLINE   80      /* maximum line length in MDL files     */
#define NONE           0      /* general default value for most data  */
                              /* fields in MDL data structures        */

/* constants for atom definitions */
/***********      charges    ********/
#define  RADICAL        4
#define  ANY_CHARGE     8

#define  PLUS3          1
#define  PLUS2          2
#define  PLUS1          3
#define  MINUS1         5
#define  MINUS2         6
#define  MINUS3         7
#define      xCHARGE(cr)    ((cr) <= 0  || (cr) >= ANY_CHARGE ?\
                     (0) : RADICAL - (cr))

/***********      radicals           ********/
#define  SINGLET       1
#define  DOUBLET       2
#define  TRIPLET       3
#define  ANY_RADICAL   8

/*********** stereo parity ********/
#define  ODD            1
#define  EVEN           2
#define  UNMARKED       3

/*********** hydrogen count *******/
#define  ZERO_COUNT     1
#define  AT_LEAST_1     2
#define  AT_LEAST_2     3
#define  AT_LEAST_3     4

/* constants for bond definitions */
#define  SINGLE            1
#define  DOUBLE            2
#define  TRIPLE            3
#define  AROMATIC          4
#define  SINGLE_DOUBLE     5
#define  SINGLE_AROMATIC   6
#define  DOUBLE_AROMATIC   7
#define  ANY_BOND          8
#define  ALL_BOND_TYPES    0xF

#define  UP                0x01
#define  DOWN              0x06
#define  EITHER            0x04
#define  CIS_TRANS_EITHER  0x03
#define  CIS_TRANS_SWAPPED 0x08

#define  TRANS_MASK        0x10
#define  CIS_MASK   0x20

#define  RING              1
#define  CHAIN             2

#endif
