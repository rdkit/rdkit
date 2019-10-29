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

#include <RDGeneral/export.h>
#ifndef BINARY_STRING_OPERATIONS_INCLUDED
#define BINARY_STRING_OPERATIONS_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

void bitstringUnion(int length, uint8 *bstr1, uint8 *bstr2);
void bitstringIntersection(int length, uint8 *bstr1, uint8 *bstr2);

int bitstringWeight(int length, uint8 *bstr);
int bitstringIntersectionWeight(int length, uint8 *bstr1, uint8 *bstr2);
int bitstringDifferenceWeight(int length, uint8 *bstr1, uint8 *bstr2);

int bitstringHemDistance(int length, uint8 *bstr1, uint8 *bstr2);
double bitstringTanimotoSimilarity(int length, uint8 *bstr1, uint8 *bstr2);
double bitstringTanimotoDistance(int length, uint8 *bstr1, uint8 *bstr2);

bool bitstringContains(int length, uint8 *bstr1, uint8 *bstr2);
bool bitstringIntersects(int length, uint8 *bstr1, uint8 *bstr2);
bool bitstringAllTrue(int length, uint8 *bstr);

void bitstringSimpleSubset(int length, uint8 *bstr, int sub_weight,
                           uint8 *sub_bstr);
void bitstringRandomSubset(int length, int weight, uint8 *bstr, int sub_weight,
                           uint8 *sub_bstr);

#ifdef __cplusplus
} /* extern "C" { */
#endif

#endif
