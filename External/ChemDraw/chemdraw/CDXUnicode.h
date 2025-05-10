// CommonCS/LibCommon/Hdr/CDXUnicode.h
// Contains: Program-independent class library for managing CDX objects
// Copyright © 1986-2004, CambridgeSoft Corp., All Rights Reserved

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

#include <string>
#include <RDGeneral/export.h>
#include "CoreChemistryAPI.h"
#include "CDXConstants.h"

CORE_CHEMISTRY_API std::string UnicodeToString(const std::string &fontName, CDXCharSet *charSet, const std::string &theText);
CORE_CHEMISTRY_API std::string UnicodeToString(const std::string &theText);
CORE_CHEMISTRY_API std::string StringToUnicode(const std::string &fontName, const std::string &theText, CDXCharSet charSet = kCDXCharSetUnknown);

#ifdef __linux
typedef std::basic_string<wchar_t> UTF16_string;
#else
typedef std::basic_string<UINT16> UTF16_string;
#endif

CORE_CHEMISTRY_API UTF16_string ConvertUTF8to16(const std::string &utf8);
std::string ConvertUTF16to8(const UTF16_string &utf16);

// Transcode a string from one CDXCharSet to another. Pass kCDXCharSetUTF8 to convert to and from UTF-8.
#if UTF8_STD_STRING
CORE_CHEMISTRY_API std::string TranscodeString(const std::string& inString, CDXCharSet inSourceEncoding, CDXCharSet inDestinationEncoding);
#endif
