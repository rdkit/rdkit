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

#include <postgres.h>

#ifdef RDK_OPTIMIZE_POPCNT

#ifdef _MSC_VER
#include <intrin.h>
#define POPCNT __popcnt
typedef unsigned int POPCNT_TYPE;
#else
#define POPCNT __builtin_popcountll
typedef unsigned long long POPCNT_TYPE;
#endif

#endif

#include "bitstring.h"

/* Number of one-bits in an unsigned byte */
static const uint8 number_of_ones[256] = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};

void bitstringUnion(int length, uint8 *bstr1, uint8 *bstr2) {
  int i;
  uint8 *bstr1_end = bstr1 + length;

#ifdef RDK_OPTIMIZE_POPCNT
  int ilength = length / sizeof(POPCNT_TYPE);

  POPCNT_TYPE *ibstr1 = (POPCNT_TYPE *)bstr1;
  POPCNT_TYPE *ibstr2 = (POPCNT_TYPE *)bstr2;
  POPCNT_TYPE *ibstr1_end = ibstr1 + ilength;
  POPCNT_TYPE *ibstr1_end4 = ibstr1_end - (ilength % 4);

  while (ibstr1 < ibstr1_end4) {
    *ibstr1++ |= *ibstr2++;
    *ibstr1++ |= *ibstr2++;
    *ibstr1++ |= *ibstr2++;
    *ibstr1++ |= *ibstr2++;
  }

  while (ibstr1 < ibstr1_end) {
    *ibstr1++ |= *ibstr2++;
  }

  bstr1 = (uint8 *)ibstr1;
  bstr2 = (uint8 *)ibstr2;
#endif

  while (bstr1 < bstr1_end) {
    *bstr1++ |= *bstr2++;
  }
}

void bitstringIntersection(int length, uint8 *bstr1, uint8 *bstr2) {
  int i;
  uint8 *bstr1_end = bstr1 + length;

#ifdef RDK_OPTIMIZE_POPCNT
  int ilength = length / sizeof(POPCNT_TYPE);

  POPCNT_TYPE *ibstr1 = (POPCNT_TYPE *)bstr1;
  POPCNT_TYPE *ibstr2 = (POPCNT_TYPE *)bstr2;
  POPCNT_TYPE *ibstr1_end = ibstr1 + ilength;
  POPCNT_TYPE *ibstr1_end4 = ibstr1_end - (ilength % 4);

  while (ibstr1 < ibstr1_end4) {
    *ibstr1++ &= *ibstr2++;
    *ibstr1++ &= *ibstr2++;
    *ibstr1++ &= *ibstr2++;
    *ibstr1++ &= *ibstr2++;
  }

  while (ibstr1 < ibstr1_end) {
    *ibstr1++ &= *ibstr2++;
  }

  bstr1 = (uint8 *)ibstr1;
  bstr2 = (uint8 *)ibstr2;
#endif

  while (bstr1 < bstr1_end) {
    *bstr1++ &= *bstr2++;
  }
}

int bitstringWeight(int length, uint8 *bstr) {
  int total_popcount = 0;
  uint8 *bstr_end = bstr + length;

#ifdef RDK_OPTIMIZE_POPCNT
  int ilength = length / sizeof(POPCNT_TYPE);

  POPCNT_TYPE *ibstr = (POPCNT_TYPE *)bstr;
  POPCNT_TYPE *ibstr_end = ibstr + ilength;
  POPCNT_TYPE *ibstr_end4 = ibstr_end - (ilength % 4);

  while (ibstr < ibstr_end4) {
    total_popcount += POPCNT(*ibstr++);
    total_popcount += POPCNT(*ibstr++);
    total_popcount += POPCNT(*ibstr++);
    total_popcount += POPCNT(*ibstr++);
  }

  while (ibstr < ibstr_end) {
    total_popcount += POPCNT(*ibstr++);
  }

  bstr = (uint8 *)ibstr;
#endif

  while (bstr < bstr_end) {
    total_popcount += number_of_ones[*bstr++];
  }

  return total_popcount;
}

int bitstringIntersectionWeight(int length, uint8 *bstr1, uint8 *bstr2) {
  int intersect_popcount = 0;
  uint8 *bstr1_end = bstr1 + length;

#ifdef RDK_OPTIMIZE_POPCNT
  int ilength = length / sizeof(POPCNT_TYPE);

  POPCNT_TYPE *ibstr1 = (POPCNT_TYPE *)bstr1;
  POPCNT_TYPE *ibstr2 = (POPCNT_TYPE *)bstr2;
  POPCNT_TYPE *ibstr1_end = ibstr1 + ilength;
  POPCNT_TYPE *ibstr1_end4 = ibstr1_end - (ilength % 4);

  while (ibstr1 < ibstr1_end4) {
    intersect_popcount += POPCNT(*ibstr1++ & *ibstr2++);
    intersect_popcount += POPCNT(*ibstr1++ & *ibstr2++);
    intersect_popcount += POPCNT(*ibstr1++ & *ibstr2++);
    intersect_popcount += POPCNT(*ibstr1++ & *ibstr2++);
  }

  while (ibstr1 < ibstr1_end) {
    intersect_popcount += POPCNT(*ibstr1++ & *ibstr2++);
  }

  bstr1 = (uint8 *)ibstr1;
  bstr2 = (uint8 *)ibstr2;
#endif

  while (bstr1 < bstr1_end) {
    intersect_popcount += number_of_ones[*bstr1++ & *bstr2++];
  }

  return intersect_popcount;
}

int bitstringDifferenceWeight(int length, uint8 *bstr1, uint8 *bstr2) {
  int difference = 0;
  uint8 *bstr1_end = bstr1 + length;

#ifdef RDK_OPTIMIZE_POPCNT
  int ilength = length / sizeof(POPCNT_TYPE);

  POPCNT_TYPE ib1;
  POPCNT_TYPE *ibstr1 = (POPCNT_TYPE *)bstr1;
  POPCNT_TYPE *ibstr2 = (POPCNT_TYPE *)bstr2;
  POPCNT_TYPE *ibstr1_end = ibstr1 + ilength;
  POPCNT_TYPE *ibstr1_end4 = ibstr1_end - (ilength % 4);

  while (ibstr1 < ibstr1_end4) {
    ib1 = *ibstr1++;
    difference += POPCNT(ib1 ^ (ib1 | *ibstr2++));
    ib1 = *ibstr1++;
    difference += POPCNT(ib1 ^ (ib1 | *ibstr2++));
    ib1 = *ibstr1++;
    difference += POPCNT(ib1 ^ (ib1 | *ibstr2++));
    ib1 = *ibstr1++;
    difference += POPCNT(ib1 ^ (ib1 | *ibstr2++));
  }

  while (ibstr1 < ibstr1_end) {
    ib1 = *ibstr1++;
    difference += POPCNT(ib1 ^ (ib1 | *ibstr2++));
  }

  bstr1 = (uint8 *)ibstr1;
  bstr2 = (uint8 *)ibstr2;
#endif

  while (bstr1 < bstr1_end) {
    uint8 b1 = *bstr1++;
    uint8 b2 = *bstr2++;
    difference += number_of_ones[b1 ^ (b1 | b2)];
  }

  return difference;
}

int bitstringHemDistance(int length, uint8 *bstr1, uint8 *bstr2) {
  int difference = 0;
  uint8 *bstr1_end = bstr1 + length;

#ifdef RDK_OPTIMIZE_POPCNT
  int ilength = length / sizeof(POPCNT_TYPE);

  POPCNT_TYPE *ibstr1 = (POPCNT_TYPE *)bstr1;
  POPCNT_TYPE *ibstr2 = (POPCNT_TYPE *)bstr2;
  POPCNT_TYPE *ibstr1_end = ibstr1 + ilength;
  POPCNT_TYPE *ibstr1_end4 = ibstr1_end - (ilength % 4);

  while (ibstr1 < ibstr1_end4) {
    difference += POPCNT(*ibstr1++ ^ *ibstr2++);
    difference += POPCNT(*ibstr1++ ^ *ibstr2++);
    difference += POPCNT(*ibstr1++ ^ *ibstr2++);
    difference += POPCNT(*ibstr1++ ^ *ibstr2++);
  }

  while (ibstr1 < ibstr1_end) {
    difference += POPCNT(*ibstr1++ ^ *ibstr2++);
  }

  bstr1 = (uint8 *)ibstr1;
  bstr2 = (uint8 *)ibstr2;
#endif

  while (bstr1 < bstr1_end) {
    difference += number_of_ones[*bstr1++ ^ *bstr2++];
  }

  return difference;
}

double bitstringTanimotoSimilarity(int length, uint8 *bstr1, uint8 *bstr2) {
  double sim;

  int union_popcount = 0;
  int intersect_popcount = 0;

  uint8 *bstr1_end = bstr1 + length;

#ifdef RDK_OPTIMIZE_POPCNT
  int ilength = length / sizeof(POPCNT_TYPE);

  POPCNT_TYPE *ibstr1 = (POPCNT_TYPE *)bstr1;
  POPCNT_TYPE *ibstr2 = (POPCNT_TYPE *)bstr2;
  POPCNT_TYPE *ibstr1_end = ibstr1 + ilength;
  POPCNT_TYPE *ibstr1_end4 = ibstr1_end - (ilength % 4);

  POPCNT_TYPE ib1, ib2;
  while (ibstr1 < ibstr1_end4) {
    ib1 = *ibstr1++;
    ib2 = *ibstr2++;
    union_popcount += POPCNT(ib1 | ib2);
    intersect_popcount += POPCNT(ib1 & ib2);
    ib1 = *ibstr1++;
    ib2 = *ibstr2++;
    union_popcount += POPCNT(ib1 | ib2);
    intersect_popcount += POPCNT(ib1 & ib2);
    ib1 = *ibstr1++;
    ib2 = *ibstr2++;
    union_popcount += POPCNT(ib1 | ib2);
    intersect_popcount += POPCNT(ib1 & ib2);
    ib1 = *ibstr1++;
    ib2 = *ibstr2++;
    union_popcount += POPCNT(ib1 | ib2);
    intersect_popcount += POPCNT(ib1 & ib2);
  }

  while (ibstr1 < ibstr1_end) {
    ib1 = *ibstr1++;
    ib2 = *ibstr2++;
    union_popcount += POPCNT(ib1 | ib2);
    intersect_popcount += POPCNT(ib1 & ib2);
  }

  bstr1 = (uint8 *)ibstr1;
  bstr2 = (uint8 *)ibstr2;
#endif

  while (bstr1 < bstr1_end) {
    uint8 b1 = *bstr1++;
    uint8 b2 = *bstr2++;
    union_popcount += number_of_ones[b1 | b2];
    intersect_popcount += number_of_ones[b1 & b2];
  }

  if (union_popcount != 0) {
    sim = ((double)intersect_popcount) / union_popcount;
  } else {
    sim = 1.0;
  }

  return sim;
}

double bitstringTanimotoDistance(int length, uint8 *bstr1, uint8 *bstr2) {
  return 1. - bitstringTanimotoSimilarity(length, bstr1, bstr2);
}

bool bitstringContains(int length, uint8 *bstr1, uint8 *bstr2) {
  bool contains = true;
  uint8 *bstr1_end = bstr1 + length;

#ifdef RDK_OPTIMIZE_POPCNT
  int ilength = length / sizeof(POPCNT_TYPE);

  POPCNT_TYPE *ibstr1 = (POPCNT_TYPE *)bstr1;
  POPCNT_TYPE *ibstr2 = (POPCNT_TYPE *)bstr2;
  POPCNT_TYPE *ibstr1_end = ibstr1 + ilength;

  while (contains && ibstr1 < ibstr1_end) {
    POPCNT_TYPE i1 = *ibstr1++;
    POPCNT_TYPE i2 = *ibstr2++;
    contains = i1 == (i1 | i2);
  }

  bstr1 = (uint8 *)ibstr1;
  bstr2 = (uint8 *)ibstr2;
#endif

  while (contains && bstr1 < bstr1_end) {
    uint8 b1 = *bstr1++;
    uint8 b2 = *bstr2++;
    contains = b1 == (b1 | b2);
  }

  return contains;
}

bool bitstringIntersects(int length, uint8 *bstr1, uint8 *bstr2) {
  bool intersects = false;
  uint8 *bstr1_end = bstr1 + length;

#ifdef RDK_OPTIMIZE_POPCNT
  int ilength = length / sizeof(POPCNT_TYPE);

  POPCNT_TYPE *ibstr1 = (POPCNT_TYPE *)bstr1;
  POPCNT_TYPE *ibstr2 = (POPCNT_TYPE *)bstr2;
  POPCNT_TYPE *ibstr1_end = ibstr1 + ilength;

  while (!intersects && (ibstr1 < ibstr1_end)) {
    POPCNT_TYPE i1 = *ibstr1++;
    POPCNT_TYPE i2 = *ibstr2++;
    intersects = (i1 & i2) != 0;
  }

  bstr1 = (uint8 *)ibstr1;
  bstr2 = (uint8 *)ibstr2;
#endif

  while (!intersects && (bstr1 < bstr1_end)) {
    uint8 b1 = *bstr1++;
    uint8 b2 = *bstr2++;
    intersects = (b1 & b2) != 0;
  }

  return intersects;
}

bool bitstringAllTrue(int length, uint8 *bstr) {
  bool allTrue = true;
  uint8 *bstr_end = bstr + length;

#ifdef RDK_OPTIMIZE_POPCNT
  /* TODO consider optimizing */
#endif

  while (allTrue && bstr < bstr_end) {
    allTrue = *bstr++ == 0xFF;
  }

  return allTrue;
}

void bitstringSimpleSubset(int length, uint8 *bstr, int sub_weight,
                           uint8 *sub_bstr) {
  /*
  ** a simple implementation. just pick the first sub_weight 1s from bstr
  ** and set them in sub_bstr
  */

  uint8 *bstr_end = bstr + length;

  int i;
  int bitcount = 0;
  uint8 bit, byte;

  while (bitcount < sub_weight && bstr < bstr_end) {
    byte = *bstr++;
    for (i = 0, bit = 0x01; bitcount < sub_weight && i < 8; ++i, bit <<= 1) {
      if (byte & bit) {
        ++bitcount;
        *sub_bstr |= bit;
      }
    }
    ++sub_bstr;
  }
}

void bitstringRandomSubset(int length, int weight, uint8 *bstr, int sub_weight,
                           uint8 *sub_bstr) {
  int i, j;
  int bitindex, bitcount, *bits;
  uint8 byte;

  Assert(sub_weight <= weight);

  /* fill an array with the indices of the 1s in bstr */
  bits = palloc(weight * sizeof(int));

  for (bitcount = 0, i = 0; i < length; ++i) {
    byte = bstr[i];
    for (j = 0; j < 8; ++j) {
      if (byte & 0x01) {
        Assert(bitcount < weight);
        bitindex = 8 * i + j;
        bits[bitcount++] = bitindex;
      }
      byte >>= 1;
    }
  }

  /* move a random subset of bit indices to the front */
  for (bitcount = 0; bitcount < sub_weight; ++bitcount) {
    /* pick a random index in [bitcount, weight - 1] */
    double r = ((double)rand()) / RAND_MAX * (weight - 1 - bitcount);
    i = bitcount + (int)(r + .5);
    /* swap the values (move the selected entry from i to bitcount) */
    bitindex = bits[i];
    bits[i] = bits[bitcount];
    bits[bitcount] = bitindex;
  }

  /* set the selected bits in the output sub_bstr */
  for (bitcount = 0; bitcount < sub_weight; ++bitcount) {
    bitindex = bits[bitcount];
    sub_bstr[bitindex / 8] |= (0x01 << bitindex % 8);
  }

  /* deallocate the array of bit indices */
  pfree(bits);
}

int bitstringGrayCmp(int length, uint8 *bstr1, uint8 *bstr2) {
  /*
   * Compare two bistrings according to their ordering in a
   * Binary Reflected Gray Code
   *
   * Which digit represents a higher value, 0 or 1 ?
   *
   * 000
   * 001
   * 011
   * 010
   * 110
   * 111
   * 101
   * 100
   *
   * In contrast to the usual representation of binary numbers,
   * the code is not just positional.
   * 
   * When we start from the leftmost position, 1 > 0. In moving
   * to right, every time a pair of 1s is found, it means the
   * remaining parts of the codes resulted from a reflection, and
   * the digit with the highest value changes from 1 to 0 and
   * viceversa.
   */

  uint8 higher = 1;

  const uint8 *bstr1_end = bstr1 + length;

  while (bstr1 < bstr1_end) {
    const uint8 bytea = *bstr1++;
    const uint8 byteb = *bstr2++;
    if (bytea == byteb) {
      /* if the number of 1s is odd, higher is flipped */
      higher ^= (1 & number_of_ones[bytea]);
    }
    else {
      uint8 mask = 0x80;
      while (mask) {
        uint8 bita = (bytea & mask) ? 1 : 0;
        uint8 bitb = (byteb & mask) ? 1 : 0;
        if (bita != bitb) {
          return (bita == higher) ? 1 : -1;
        }
        else {
          /* flip higher if bita is 1 */
          higher ^= bita;
        }
        mask >>= 1;
      }
      Assert(!"It should never get here if bytea != byteb");
    }
  }

  /* same bfp value */
  return 0;
}
