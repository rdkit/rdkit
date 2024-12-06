// CommonCS/LibCommon/Hdr/CDXUnicode.h
// Contains: Program-independent class library for managing CDX objects
// Copyright © 1986-2004, CambridgeSoft Corp., All Rights Reserved

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
