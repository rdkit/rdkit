// $Id$
//
// Copyright (c) 2002-20`0  greg Landrum, Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "BitVects.h"
#include "BitVectUtils.h"
#include <RDGeneral/Invariant.h>
#include <iostream>

//! Convert a SparseBitVector to an ExplicitBitVector
ExplicitBitVect *convertToExplicit(const SparseBitVect *sbv) {
  unsigned int sl = sbv->getNumBits();
  auto *ebv = new ExplicitBitVect(sl);
  const IntSet *bset = sbv->getBitSet();
  for (int it : *bset) {
    ebv->setBit(it);
  }
  return ebv;
}

void a2b(const char *, char *);

//! \brief Construct a BitVect from the ASCII representation of a
//! Daylight fingerprint string
template <typename T>
void FromDaylightString(T &sbv, const std::string &s) {
  sbv.clearBits();
  size_t length = s.length();
  size_t nBits;

  if (s[length - 1] == '\n') {
    length -= 1;
  }

  // 4 bytes in the ascii correspond to 3 bytes in the binary
  //  plus there's one extra ascii byte for the pad marker
  length -= 1;
  nBits = (3 * length / 4) * 8;

  switch (s[length]) {
    case '1':
      nBits -= 16;
      break;
    case '2':
      nBits -= 8;
      break;
    case '3':
      break;
    default:
      throw "ValueError bad daylight fingerprint string";
  }
  size_t i = 0, nBitsDone = 0;
  while (i < length) {
    char bytes[3];
    a2b(s.c_str() + i, bytes);
    for (size_t j = 0; j < 3 && nBitsDone < nBits; j++) {
      unsigned char query = 0x80;
      for (size_t k = 0; k < 8; k++) {
        if (bytes[j] & query) {
          sbv.setBit(nBitsDone);
        }
        query >>= 1;
        nBitsDone++;
      }
    }
    i += 4;
  }
}

template RDKIT_DATASTRUCTS_EXPORT void FromDaylightString(SparseBitVect &sbv,
                                                          const std::string &s);
template RDKIT_DATASTRUCTS_EXPORT void FromDaylightString(ExplicitBitVect &sbv,
                                                          const std::string &s);

//! \brief Construct a BitVect from the ASCII representation of a
//! BitString
template <typename T>
void FromBitString(T &sbv, const std::string &s) {
  PRECONDITION(s.length() <= sbv.getNumBits(), "bad bitvect length");
  sbv.clearBits();
  for (unsigned int i = 0; i < sbv.getNumBits(); ++i) {
    if (s[i] == '1') {
      sbv.setBit(i);
    }
  }
}

template RDKIT_DATASTRUCTS_EXPORT void FromBitString(SparseBitVect &sbv,
                                                     const std::string &s);
template RDKIT_DATASTRUCTS_EXPORT void FromBitString(ExplicitBitVect &sbv,
                                                     const std::string &s);

//! converts 4 ascii bytes at a4 to 3 binary bytes
/*!
 THE FOLLOWING IS TAKEN FROM THE DAYLIGHT CONTRIB PROGRAM
   ascii2bits.c
*********************************************************************
*** a2b - converts 4 ascii bytes at a4 to 3 binary
***       bytes at b3.
***
***  ASCII:    |=======+=======+=======+=======| etc.
***                                            ^
***    becomes...                      3  <->  4
***                                    v
***  BINARY:   |=====+=====+=====+=====| etc.
********************************************************************
*/
void a2b(const char *a4, char *b3) {
  int i;
  char byte = 0x00, b = 0x00;

  /*********************************************
  *** Use the Daylight mapping to convert each
  *** ascii char to its 6-bit code.
  ***
  *** a4: |xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx (printable)
  ***     |=======+=======+=======+=======|
  ***   becomes...
  *** a4: |00xxxxxx00xxxxxx00xxxxxx00xxxxxx
  ***     |=======+=======+=======+=======|
  *********************************************/
  for (i = 0; i < 4; ++i) {
    switch (a4[i]) {
      case '.':
        byte = 0x00;
        break; /* 00 = __000000 */
      case '+':
        byte = 0x01;
        break; /* 01 = __000001 */
      case '0':
        byte = 0x02;
        break; /* 02 = __000010 */
      case '1':
        byte = 0x03;
        break; /* 03 = __000011 */
      case '2':
        byte = 0x04;
        break; /* 04 = __000100 */
      case '3':
        byte = 0x05;
        break; /* 05 = __000101 */
      case '4':
        byte = 0x06;
        break; /* 06 = __000110 */
      case '5':
        byte = 0x07;
        break; /* 07 = __000111 */
      case '6':
        byte = 0x08;
        break; /* 08 = __001000 */
      case '7':
        byte = 0x09;
        break; /* 09 = __001001 */
      case '8':
        byte = 0x0a;
        break; /* 10 = __001010 */
      case '9':
        byte = 0x0b;
        break; /* 11 = __001011 */
      case 'A':
        byte = 0x0c;
        break; /* 12 = __001100 */
      case 'B':
        byte = 0x0d;
        break; /* 13 = __001101 */
      case 'C':
        byte = 0x0e;
        break; /* 14 = __001110 */
      case 'D':
        byte = 0x0f;
        break; /* 15 = __001111 */
      case 'E':
        byte = 0x10;
        break; /* 16 = __010000 */
      case 'F':
        byte = 0x11;
        break; /* 17 = __010001 */
      case 'G':
        byte = 0x12;
        break; /* 18 = __010010 */
      case 'H':
        byte = 0x13;
        break; /* 19 = __010011 */
      case 'I':
        byte = 0x14;
        break; /* 20 = __010100 */
      case 'J':
        byte = 0x15;
        break; /* 21 = __010101 */
      case 'K':
        byte = 0x16;
        break; /* 22 = __010110 */
      case 'L':
        byte = 0x17;
        break; /* 23 = __010111 */
      case 'M':
        byte = 0x18;
        break; /* 24 = __011000 */
      case 'N':
        byte = 0x19;
        break; /* 25 = __011001 */
      case 'O':
        byte = 0x1a;
        break; /* 26 = __011010 */
      case 'P':
        byte = 0x1b;
        break; /* 27 = __011011 */
      case 'Q':
        byte = 0x1c;
        break; /* 28 = __011100 */
      case 'R':
        byte = 0x1d;
        break; /* 29 = __011101 */
      case 'S':
        byte = 0x1e;
        break; /* 30 = __011110 */
      case 'T':
        byte = 0x1f;
        break; /* 31 = __011111 */
      case 'U':
        byte = 0x20;
        break; /* 32 = __100000 */
      case 'V':
        byte = 0x21;
        break; /* 33 = __100001 */
      case 'W':
        byte = 0x22;
        break; /* 34 = __100010 */
      case 'X':
        byte = 0x23;
        break; /* 35 = __100011 */
      case 'Y':
        byte = 0x24;
        break; /* 36 = __100100 */
      case 'Z':
        byte = 0x25;
        break; /* 37 = __100101 */
      case 'a':
        byte = 0x26;
        break; /* 38 = __100110 */
      case 'b':
        byte = 0x27;
        break; /* 39 = __100111 */
      case 'c':
        byte = 0x28;
        break; /* 40 = __101000 */
      case 'd':
        byte = 0x29;
        break; /* 41 = __101001 */
      case 'e':
        byte = 0x2a;
        break; /* 42 = __101010 */
      case 'f':
        byte = 0x2b;
        break; /* 43 = __101011 */
      case 'g':
        byte = 0x2c;
        break; /* 44 = __101100 */
      case 'h':
        byte = 0x2d;
        break; /* 45 = __101101 */
      case 'i':
        byte = 0x2e;
        break; /* 46 = __101110 */
      case 'j':
        byte = 0x2f;
        break; /* 47 = __101111 */
      case 'k':
        byte = 0x30;
        break; /* 48 = __110000 */
      case 'l':
        byte = 0x31;
        break; /* 49 = __110001 */
      case 'm':
        byte = 0x32;
        break; /* 50 = __110010 */
      case 'n':
        byte = 0x33;
        break; /* 51 = __110011 */
      case 'o':
        byte = 0x34;
        break; /* 52 = __110100 */
      case 'p':
        byte = 0x35;
        break; /* 53 = __110101 */
      case 'q':
        byte = 0x36;
        break; /* 54 = __110110 */
      case 'r':
        byte = 0x37;
        break; /* 55 = __110111 */
      case 's':
        byte = 0x38;
        break; /* 56 = __111000 */
      case 't':
        byte = 0x39;
        break; /* 57 = __111001 */
      case 'u':
        byte = 0x3a;
        break; /* 58 = __111010 */
      case 'v':
        byte = 0x3b;
        break; /* 59 = __111011 */
      case 'w':
        byte = 0x3c;
        break; /* 60 = __111100 */
      case 'x':
        byte = 0x3d;
        break; /* 61 = __111101 */
      case 'y':
        byte = 0x3e;
        break; /* 62 = __111110 */
      case 'z':
        byte = 0x3f;
        break; /* 63 = __111111 */
    }

    /*********************************************
    *** Now copy the 4x6=24 bits from a4 to b3.
    ***
    *** a4: |--000000--111111--222222--333333
    ***     |=======+=======+=======+=======|
    ***
    *** b3: |000000111111222222333333
    ***     |=====+=====+=====+=====|
    *********************************************/
    switch (i) {
      case 0:
        b3[0] = (byte << 2); /*** 6 bits into 1st byte ***/
        break;
      case 1:
        b3[0] |= ((b = byte) >> 4); /*** 2 bits into 1st byte ***/
        b3[1] = ((b = byte) << 4);  /*** 4 bits into 2nd byte ***/
        break;
      case 2:
        b3[1] |= ((b = byte) >> 2); /*** 4 bits into 2nd byte ***/
        b3[2] = ((b = byte) << 6);  /*** 2 bits into 3rd byte ***/
        break;
      case 3:
        b3[2] |= byte; /*** 6 bits into 3rd byte ***/
        break;
    }
  }
  return;
}

// Demo Data:
// 256 bits:
//.b7HEa..ccc+gWEIr89.8lV8gOF3aXFFR.+Ps.mZ6lg.2
//
// 00000010 01110010 01010011 01000010 01100000
// 00000000 10100010 10001010 00000001 10110010
// 00100100 00010100 11011100 10100010 11000000
// 00101011 00011000 01001010 10110001 10100100
// 01000101 10011010 00110100 01010001 01110100
// 00000000 01011011 11100000 00001100 10100101
// 00100011 00011011
