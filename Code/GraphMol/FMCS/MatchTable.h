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
#include <vector>
#include <stdexcept>

namespace RDKit {
    namespace FMCS {
        template<typename T>
        class TArray2D { // for scalar value types ! including bool with special STL implementation (no reference to item - bitset used)
            size_t          XSize;
            size_t          YSize;
            std::vector<T>  Data;
        public:
            inline TArray2D(size_t cy=0, size_t cx=0) : XSize(cx), YSize(cy), Data(cx*cy) {}
            inline size_t getXSize()const {
                return XSize;
            }
            inline size_t getYSize()const {
                return YSize;
            }
            inline bool   empty   ()const {
                return Data.empty();
            }
            inline void   clear   () {
                Data.clear();
                XSize = 0;
                YSize = 0;
            }
            inline void   resize(size_t cy, size_t cx) {
                Data.resize(cx*cy);
                XSize = cx;
                YSize = cy;
            }
            inline void set(size_t row, size_t col, T val) {
                Data[row*XSize + col] = val;
            }
            inline T at(size_t row, size_t col)     {
                return Data[row*XSize + col];
            }
            inline T at(size_t row, size_t col)const {
                return Data[row*XSize + col];
            }
        };

        typedef TArray2D<bool>    MatchTable;   // row is index in QueryMolecule

    }
}
