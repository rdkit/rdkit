//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once
namespace RDKit {
    namespace FMCS {
        typedef unsigned long long BitSet;
        class Composition2N { // generator of 2^N-1 possible bit combinations
            BitSet  Bits, InverseBits;
            BitSet  MaxValue, ValueMask; // need for inverse bitset must be 2^N-1
        public:
            Composition2N(BitSet maxValue, BitSet valueMask) : Bits(0), InverseBits(0), MaxValue(maxValue), ValueMask(valueMask) {}

            static void compute2N(unsigned power, BitSet& value) {
                value = 1uLL << power;
            }

            BitSet getBitSet()const {
                return InverseBits;   // inverse to generate biggest seed first and then decrease number of external bonds
            }

            bool generateNext() {
                if((++Bits) <= MaxValue) {
                    InverseBits = (~Bits+1) & ValueMask;
                    return true;
                } else
                    return false;
            }
            bool is2Power()const { // one bit is set only
                BitSet  bits = getBitSet();
                unsigned n = 0;
                while(0==(bits & 1uLL) && ++n < sizeof(bits)*8)    //find lowest bitwise 1
                    bits >>= 1u;    //shift all zero lower bits
                if(0!=(bits & 1uLL))
                    bits >>= 1u;    //shift first set bit too
                return 0==bits;     //remained bits except lowest 1
            }
            //unused:        bool nonZero() {return 0!=getBitSet();}
            bool isSet(unsigned bit)const {
                return 0 != (getBitSet() & (1uLL << bit));
            }
        };
    }
}
