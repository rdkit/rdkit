// BSD 3-Clause License
// 
// Copyright (c) 1986-2025, CambridgeSoft Corp, Revvity Inc and others.
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
#pragma once

bool  csIsNetVariant (const char *s, long len);
bool  csIsValidBase64 (const char *s, long len);
long  csBase64Encode_Ext (const char	*input,
			  long			inputLen,
			  char			*&output,
			  bool			bNetVariant,
			  unsigned int	lineBreakFrequency = 0);
long  csBase64Decode(const char *input, long inputLen, char *&output, bool bNetVariantOnly);
void csSetBase64DecodeAuthorization (bool bAuthorized);

// Added by Glysade
#include <string>

inline std::string csBase64ToBinary(const std::string &input) {
  char *output;
  auto outputlen = csBase64Decode(input.c_str(), input.size(), output, false);
  auto res = std::string(output, outputlen);
  delete output;
  return res;
}

inline std::string csBinaryToBase64(const std::string &input) {
  char *output;
  auto outputlen = csBase64Encode_Ext(input.c_str(), input.size(), output, false);
  auto res = std::string(output, outputlen);
  delete output;
  return res;
}
