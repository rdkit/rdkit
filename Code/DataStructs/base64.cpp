// $Id$
//
// Copyright (c) 2002-2008 greg landrum and rational discovery llc
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <iostream>
#include <cstring>
#include "base64.h"
// Encoding table from RFC 2045
//
//     Value Encoding  Value Encoding  Value Encoding  Value Encoding
//         0 A            17 R            34 i            51 z
//         1 B            18 S            35 j            52 0
//         2 C            19 T            36 k            53 1
//         3 D            20 U            37 l            54 2
//         4 E            21 V            38 m            55 3
//         5 F            22 W            39 n            56 4
//         6 G            23 X            40 o            57 5
//         7 H            24 Y            41 p            58 6
//         8 I            25 Z            42 q            59 7
//         9 J            26 a            43 r            60 8
//        10 K            27 b            44 s            61 9
//        11 L            28 c            45 t            62 +
//        12 M            29 d            46 u            63 /
//        13 N            30 e            47 v
//        14 O            31 f            48 w         (pad) =
//        15 P            32 g            49 x
//        16 Q            33 h            50 y

char *Base64Encode(const char *inText, const unsigned int inLen) {
  return Base64Encode((const unsigned char *)inText, inLen);
}

char *Base64Encode(const unsigned char *inText, const unsigned int inLen) {
  // Notes:
  //   - whoever calls us is responsible for free'ing the result we return
  //   - we cheat and don't worry about breaking lines
  static unsigned char transTable[64] = {
      'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
      'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
      'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
      'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
      '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '+', '/'};

  int resSize = (4 * inLen) / 3;
  while (resSize % 4) {
    resSize++;
  }
  char *res = new char[resSize + 1];
  unsigned int i = 0;
  int pos = 0;
  while (i < inLen) {
    res[pos++] = transTable[inText[i] >> 2];
    if (i + 1 < inLen) {
      res[pos++] = transTable[((inText[i] & 3) << 4) | (inText[i + 1] >> 4)];
      if (i + 2 < inLen) {
        res[pos++] =
            transTable[((inText[i + 1] & 0xF) << 2) | (inText[i + 2] >> 6)];
        res[pos++] = transTable[inText[i + 2] & 0x3F];
      } else {
        // single padding
        res[pos++] = transTable[((inText[i + 1] & 0xF) << 2)];
        res[pos++] = '=';
      }
    } else {
      // double padding
      res[pos++] = transTable[((inText[i] & 3) << 4)];
      res[pos++] = '=';
      res[pos++] = '=';
    }
    i += 3;
  }
  res[resSize] = 0;
  return res;
}

char *Base64Decode(const char *inText, unsigned int *size) {
  // Notes:
  //   - whoever calls us is responsible for free'ing the result we return

  unsigned char transTable[256];
  size_t inLen = strlen(inText);

  size_t i;
  // FIX: we don't really need to build this table here
  for (i = 0; i < 255; i++) {
    transTable[i] = 0x80;
  }
  for (i = 'A'; i <= 'Z'; i++) {
    transTable[i] = (unsigned char)i - 'A';
  }
  for (i = 'a'; i <= 'z'; i++) {
    transTable[i] = (unsigned char)i - 'a' + 26;
  }
  for (i = '0'; i <= '9'; i++) {
    transTable[i] = (unsigned char)i - '0' + 52;
  }
  transTable[static_cast<int>('+')] = 62;
  transTable[static_cast<int>('/')] = 63;

  size_t outLen = 3 * inLen / 4;
  auto *res = new char[outLen];
  res[outLen - 1] = 0;
  size_t pos = 0;
  i = 0;
  // decode 4 bytes at a time
  unsigned char block[4];
  int nInBlock = 0;
  while (i < inLen) {
    unsigned char c = inText[i];

    // above we set 0x80 as the junk marker in the translation table
    if (!(transTable[c] & 0x80)) {
      block[nInBlock++] = transTable[c];
      if (nInBlock == 4) {
        // finished a block
        res[pos++] = (block[0] << 2) | (block[1] >> 4);
        res[pos++] = (block[1] << 4) | (block[2] >> 2);
        res[pos++] = (block[2] << 6) | block[3];
        nInBlock = 0;
      }
    }
    i++;
  }

  // okay, now there can be 2 or 3 chars remaining to be processed
  //  (before the padding)
  if (nInBlock > 1) {
    res[pos++] = (block[0] << 2) | (block[1] >> 4);
    if (nInBlock > 2) {
      res[pos++] = (block[1] << 4) | (block[2] >> 2);
      res[pos] = (block[2] << 6);
    }
  }
  *size = pos;
  return res;
}
