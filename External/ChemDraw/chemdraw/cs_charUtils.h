// BSD 3-Clause License
// 
// Copyright (c) 2025, Glysade Inc
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#ifndef CS_CHAR_UTILS_H
#define CS_CHAR_UTILS_H

#include "CoreChemistryAPI.h"

namespace cs {
  inline bool IsChargeSign(char s) {
    return s == '-'; // XXX need a lot more logic here....
  }

// Modern symbols - these use UTF-8 where possible to fully encode characters
CORE_CHEMISTRY_API extern const char* kSymbolDegree;
CORE_CHEMISTRY_API extern const char* kSymbolEllipsis;
CORE_CHEMISTRY_API extern const char* kSymbolEnDash;
CORE_CHEMISTRY_API extern const char* kSymbolEmDash;
CORE_CHEMISTRY_API extern const char* kSymbolPlusMinus;
CORE_CHEMISTRY_API extern const char* kSymbolBullet;
CORE_CHEMISTRY_API extern const char* kSymbolCenterDot;
CORE_CHEMISTRY_API extern const char* kSymbolReg;
CORE_CHEMISTRY_API extern const char* kSymbolCopyright;
CORE_CHEMISTRY_API extern const char* kSymbolAngstrom;
CORE_CHEMISTRY_API extern const char* kSymbolMicro;
CORE_CHEMISTRY_API extern const char* kSymbolCent;
CORE_CHEMISTRY_API extern const char* kSymbolPound;
}
#endif
