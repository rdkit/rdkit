
// CommonCS/LibCommon/Src/CDXUnicode.cpp
// Contains: Program-independent class library for managing CDX objects
// Copyright © 1986-2008, CambridgeSoft Corp., All Rights Reserved

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

#include "CDXUnicode.h"
#include "cs_auto_buffer.h"
#include "cs_stringUtils.h"
#include "cs_specialchardefs.h"
#include "cs_charUtils.h"
#include "UTF8Iterator.h"

#if TARGET_OS_MAC
    // Nothing special
#elif TARGET_WEB
    #include <emscripten.h>
    #include <string.h>
#elif TARGET_OS_WIN32
    #include <windows.h>
#elif defined __linux
    #include <vector>
    #include <wchar.h>
    #include <locale.h>
    #include <boost/locale.hpp>
#else
    #error Unsupported platform
#endif

#include <algorithm>
#include <stdexcept>

using namespace std;

namespace
{
    #if TARGET_OS_WIN32  
    /**
     *  Convert the string to UTF-8 
     *
     *  @param theText The string to be converted to UTF-8
     *  @param charSet The Code page to use in performing the conversion
     *
     *  @return the UTF-8 string
     */
    string StringToUTF8(const string &theText, CDXCharSet charSet)
    {
        // First convert to unicode using UTF-16 encoding
        int nWideChars = 0;
        UINT32 codePage;

        // If we were given a specific encoding, try that first (a CDX CharSet == a Windows Code Page)
        if (charSet != kCDXCharSetUnknown)
        {
            codePage = charSet;
            nWideChars = MultiByteToWideChar(codePage, MB_ERR_INVALID_CHARS, theText.data(), (int)theText.size(), 0, 0);
        }

        // If we weren't given a specific encoding or its translation failed, try using the system default codepage
        if (nWideChars == 0)
        {
            codePage = CP_ACP;
            nWideChars = MultiByteToWideChar(codePage, MB_ERR_INVALID_CHARS, theText.data(), (int)theText.size(), 0, 0);
        }

        // If we can't convert using the system default codepage, try hardcoding to ascii
        if (nWideChars == 0)
        {
            codePage = kCDXCharSetWin31Latin1;
            nWideChars = MultiByteToWideChar(codePage, MB_ERR_INVALID_CHARS, theText.data(), (int)theText.size(), 0, 0);
        }

        // We really had better have gotten something by this point
        ASSERT(nWideChars != 0);
        if (nWideChars == 0)
        {
            throw runtime_error("Text cannot be converted to Unicode");
        }

        auto_buffer<wchar_t> wide_chars(nWideChars);
        nWideChars = MultiByteToWideChar(codePage, MB_ERR_INVALID_CHARS, theText.data(), (int)theText.size(), wide_chars.ptr(), nWideChars);

        // Now convert to UTF-8
        string utf8 = ConvertUTF16to8(UTF16_string((UINT16*)wide_chars.ptr(), nWideChars));

        return utf8;
    }
    #endif // TARGET_OS_WIN32
}

// This variant only works if the input text is always in the standard character set
std::string UnicodeToString(const std::string &theText)
{
	CDXCharSet actualCharSet = kCDXCharSetUnknown;
	// XXX Need implementation return UnicodeToString("", &actualCharSet, theText);
  return theText;
}

UTF16_string ConvertUTF8to16(const std::string &utf8)
{
    UTF16_string utf16;
    utf16.reserve(utf8.length());

    for (UTF8Iterator i(utf8); !i.AtEnd(); i++)
    {
        auto character = i.GetCharacter();

        if (character < 0x10000)
        {
            // Requires one UTF-16 codeunits
            utf16 += (UINT16)character;
        }
        else
        {
            // Requires two UTF-16 codeunits
            uint32_t t = character - 0x10000;
            utf16 += (UINT16)(((t << 12) >> 22) + 0xD800);
            utf16 += (UINT16)(((t << 22) >> 22) + 0xDC00);
        }
    }

	return utf16;
}

std::string ConvertUTF16to8(const UTF16_string &utf16)
{
	string utf8;
    utf8.reserve(utf16.length());

    for (UTF16_string::const_iterator p = utf16.begin();  p != utf16.end();  ++p)
	{
		if (*p <= 0x7F) // low-ascii converts to single byte
        {
			utf8 += (char)(*p);
        }
        else if (*p <= 0x7FF) // 00000xxx xxyyyyyy converts to two bytes: 110xxxxx 10yyyyyy
		{
			utf8 += 0xC0 | (*p >> 6);
			utf8 += 0x80 | (*p & 0x3F);
		}
		else if (*p < (UINT16)0xD800 || *p >= (UINT16)0xE000) 
		{
			// xxxxyyyy yyzzzzzz converts to three bytes: 1110xxxx 10yyyyyy 10zzzzzz
			// but if the utf16 word is in the range 0xD800 through 0xDFFF then it's part of a composite
			utf8 += 0xE0 | (*p >> 12);
			utf8 += 0x80 | ((*p >> 6) & 0x3F);
			utf8 += 0x80 | (*p & 0x3F);
		}
		else
		{
			if (*p >= (UINT16)0xDC00 || (p+1) == utf16.end())
            {
				throw runtime_error("Invalid UTF-16 text");
            }
            
            UINT32 c = ((UINT32)*p - (UINT32)0xD800) << 10;
			++p;
			c += ((UINT32)*p - (UINT32)0xDC00);
            c += 0x10000;
            
			utf8 += 0xF0 | (c >> 18);
			utf8 += 0x80 | ((c >> 12) & 0x3F);
			utf8 += 0x80 | ((c >> 6) & 0x3F);
			utf8 += 0x80 | (c & 0x3F);
		}
	}
	return utf8;
}

/**
 *  Convert a UTF-8 string to a CDX string using font information. The font information can be dropped when
 *  we are using UTF-8 internally on all platforms.
 */
std::string UnicodeToString(const std::string &fontName, CDXCharSet *charSet, const std::string &theText)
{
	if (theText.empty())
    {
		return "";
    }
    
#if UTF8_STD_STRING
    // Make sure that nothing else tries to convert us
    if (charSet)
    {
        *charSet = kCDXCharSetUTF8;
    }

    // Return our UTF-8 (ish) internal string
    return theText;
#endif

#if (TARGET_OS_MAC && !TARGET_OS_IPHONE) || TARGET_WEB
    // Mac and Web use UTF-8
    string result;
    ASSERT(0);
#elif defined __linux || TARGET_OS_IPHONE
    string result;
	string::size_type i,n = theText.size();
	for (i=0; i<n;) {
		unsigned char c = theText[i++];
		if (c == '&' && theText[i] == '#') {
			i++;
			c = atol(theText.c_str()+i);
			i = theText.find(';', i);
			if (i == string::npos)
				i = n;
			else
				i++;
		}
		result += c;
	}
#else
	// First, we have to decode UTF-8 manually
	UTF16_string utf16 = ConvertUTF8to16(theText);

	// Now convert from UTF-16 to mbcs
	auto_buffer<char> outbuf (2 * (int)utf16.length());
	size_t bestDefault = 0;
	static unsigned int codePages[] = {
		CP_ACP,
		kCDXCharSetHebrew,
		kCDXCharSetArabic,
		kCDXCharSetWin31EasternEuropean,
		kCDXCharSetWin31Cyrillic,
		kCDXCharSetJapanese,
		kCDXCharSetChineseTraditional,
		kCDXCharSetChineseSimplified,
		kCDXCharSetKorean
	};
	const int numCodePages = sizeof codePages/sizeof codePages[0];

	
	// Go through each of the code pages looking for one which can translate all the characters.
	BOOL usedDefault;
	char defaultChar[] = { 127, 0 };
	int oOutputLen = WideCharToMultiByte(*charSet, 0, (LPCWSTR)utf16.data(), (int)utf16.length(), NULL, 0, defaultChar, &usedDefault);
	if (oOutputLen == 0 || usedDefault)
	{
		for (UINT *cpp = codePages;  cpp != codePages + numCodePages;  ++cpp)
		{
			oOutputLen = WideCharToMultiByte(*cpp, 0, (LPCWSTR)utf16.data(), (int)utf16.length(), outbuf.ptr(), (int)outbuf.length(), defaultChar, &usedDefault);
			if (!usedDefault)
			{
				*charSet = (CDXCharSet) *cpp;
				break;
			}
			size_t numDefault = std::count(outbuf.ptr(), outbuf.ptr() + outbuf.length(), 127);
			if (cpp == codePages || numDefault < bestDefault)
			{
				bestDefault = numDefault;
				*charSet = (CDXCharSet) *cpp;
			}
		}
	}

	oOutputLen = WideCharToMultiByte(*charSet, 0, (LPCWSTR)utf16.data(), (int)utf16.length(), outbuf.ptr(), (int)outbuf.length(), 0, 0);

	// The converted string might include a null at the end
	if (oOutputLen > 0 && outbuf[oOutputLen-1] == '\0')
		--oOutputLen;
	string result(outbuf.ptr(), oOutputLen);
#endif
	return result;
}

/**
 *  Convert a CDX string to a UTF-8 string using font information. This method can be removed
 *  when all platforms are using UTF-8 std::string as no conversion will be necessary,
 */
std::string StringToUnicode(const std::string &fontName, const std::string &theText, CDXCharSet charSet)
{
	if (theText.empty())
    {
		return theText;
    }
    
#if UTF8_STD_STRING
    return theText;
#endif

	std::string utf8;
#if (TARGET_OS_MAC && !TARGET_OS_IPHONE) || TARGET_WEB
    // Mac and Web use UTF-8
    ASSERT(0);
#elif TARGET_OS_IPHONE
    return theText; // always assume string is utf8 on iOS
#elif defined __linux
#ifndef __hpux
	switch (charSet) {
	case kCDXCharSetJapanese:
		setlocale(LC_ALL, "japanese");
		break;
	case kCDXCharSetChineseSimplified:
		setlocale(LC_ALL, "chinese-simplified");
		break;
	case kCDXCharSetChineseTraditional:
		setlocale(LC_ALL, "chinese-traditional");
		break;
	case kCDXCharSetKorean:
		setlocale(LC_ALL, "korean");
		break;
    default:
        break;  // these cases should be handled also
	}

	UTF16_string utf16;

	const char *p = theText.data();
	long x,i,n = theText.size();
	for (i=0; i<n; i+=x) {
		wchar_t w;
		x = mbrtowc(&w, p+i, n-i, 0);
		if (x<0)
			break;
		utf16 += w;
	}
		
	utf8 = ConvertUTF16to8(utf16);
#else	
/* this is the correct solution, but did not work on Solaris
	size_t n = 0;
	const char *p;
	
	p = theText.data();
	n = mbsrtowcs(0, &p, theText.size(), 0);
	if (n == -1)	// error
		utf8 = theText;
	else {
		vector<wchar_t> w(n);
		p = theText.data();
		mbsrtowcs(&w[0], &p, theText.size(), 0);
		UTF16_string utf16((UINT16 *) &w[0], n);
		utf8 = ConvertUTF16to8(utf16);
	}
*/
// the one below works for high ascii characters, but does not work for chinese and japanese character sets

	int i,n = theText.size();
	int nc = 0;
	
	for (i=0; i<n; i++) {
		unsigned char c = theText[i];
		nc += c < 0x80 ? 1 : 2;
	}
	
	utf8.reserve(nc);

	for (i=0; i<n; i++) {
		unsigned char c = theText[i];
		if (c < 0x80)
			utf8 += c;
		else {
			unsigned char c1 = (c & 0xC0) >> 6;
			unsigned char c2 = (c & 0x3F);
			utf8 += 0xC0 | c1;
			utf8 += 0x80 | c2;
		}
	}
#endif // __hpux

#else
    utf8 = StringToUTF8(theText, charSet);
#endif
	return utf8;
}

/**
 *  Mac implementation of unicode conversions
 */
#if TARGET_OS_MAC

/**
 *  Return our encoding map, creating it if required
 */
static const map<CDXCharSet,CFStringEncoding>& CDXToCFEncodingMap()
{
    static map<CDXCharSet,CFStringEncoding> sEncodingMap;
    
    if (sEncodingMap.size() == 0)
    {
        sEncodingMap[kCDXCharSetEBCDICOEM] = kCFStringEncodingEBCDIC_CP037;
        sEncodingMap[kCDXCharSetMSDOSUS] = kCFStringEncodingDOSLatinUS;
        sEncodingMap[kCDXCharSetEBCDIC500V1] = kCFStringEncodingEBCDIC_US;
        sEncodingMap[kCDXCharSetArabicASMO708] = kCFStringEncodingISOLatinArabic;
        sEncodingMap[kCDXCharSetArabicASMO449P] = kCFStringEncodingISOLatinArabic;
        sEncodingMap[kCDXCharSetArabicTransparent] = kCFStringEncodingISOLatinArabic;
        sEncodingMap[kCDXCharSetArabicTransparentASMO] = kCFStringEncodingISOLatinArabic;
        sEncodingMap[kCDXCharSetGreek437G] = kCFStringEncodingDOSGreek;
        sEncodingMap[kCDXCharSetBalticOEM] = kCFStringEncodingDOSBalticRim;
        sEncodingMap[kCDXCharSetMSDOSLatin1] = kCFStringEncodingDOSLatin1;
        sEncodingMap[kCDXCharSetMSDOSLatin2] = kCFStringEncodingDOSLatin2;
        sEncodingMap[kCDXCharSetIBMCyrillic] = kCFStringEncodingDOSCyrillic;
        sEncodingMap[kCDXCharSetIBMTurkish] = kCFStringEncodingDOSTurkish;
        sEncodingMap[kCDXCharSetMSDOSPortuguese] = kCFStringEncodingDOSPortuguese;
        sEncodingMap[kCDXCharSetMSDOSIcelandic] = kCFStringEncodingDOSIcelandic;
        sEncodingMap[kCDXCharSetHebrewOEM] = kCFStringEncodingDOSHebrew;
        sEncodingMap[kCDXCharSetMSDOSCanadianFrench] = kCFStringEncodingDOSCanadianFrench;
        sEncodingMap[kCDXCharSetArabicOEM] = kCFStringEncodingDOSArabic;
        sEncodingMap[kCDXCharSetMSDOSNordic] = kCFStringEncodingDOSNordic;
        sEncodingMap[kCDXCharSetMSDOSRussian] = kCFStringEncodingDOSRussian;
        sEncodingMap[kCDXCharSetIBMModernGreek] = kCFStringEncodingDOSGreek2;
        sEncodingMap[kCDXCharSetThai] = kCFStringEncodingDOSThai;
        sEncodingMap[kCDXCharSetEBCDIC] = kCFStringEncodingEBCDIC_US;
        sEncodingMap[kCDXCharSetJapanese] = kCFStringEncodingDOSJapanese;
        sEncodingMap[kCDXCharSetChineseSimplified] = kCFStringEncodingDOSChineseSimplif;
        sEncodingMap[kCDXCharSetKorean] = kCFStringEncodingDOSKorean;
        sEncodingMap[kCDXCharSetChineseTraditional] = kCFStringEncodingDOSChineseTrad;
        sEncodingMap[kCDXCharSetUnicodeISO10646] = kCFStringEncodingUTF16;
        sEncodingMap[kCDXCharSetWin31EasternEuropean] = kCFStringEncodingWindowsLatin2;
        sEncodingMap[kCDXCharSetWin31Cyrillic] = kCFStringEncodingWindowsCyrillic;
        sEncodingMap[kCDXCharSetWin31Latin1] = kCFStringEncodingWindowsLatin1;
        sEncodingMap[kCDXCharSetWin31Greek] = kCFStringEncodingWindowsGreek;
        sEncodingMap[kCDXCharSetWin31Turkish] = kCFStringEncodingWindowsLatin5;
        sEncodingMap[kCDXCharSetHebrew] = kCFStringEncodingWindowsHebrew;
        sEncodingMap[kCDXCharSetArabic] = kCFStringEncodingWindowsArabic;
        sEncodingMap[kCDXCharSetBaltic] = kCFStringEncodingWindowsBalticRim;
        sEncodingMap[kCDXCharSetVietnamese] = kCFStringEncodingWindowsVietnamese;
        sEncodingMap[kCDXCharSetKoreanJohab] = kCFStringEncodingWindowsKoreanJohab;
        sEncodingMap[kCDXCharSetMacRoman] = kCFStringEncodingMacRoman;
        sEncodingMap[kCDXCharSetMacJapanese] = kCFStringEncodingMacJapanese;
        sEncodingMap[kCDXCharSetMacTradChinese] = kCFStringEncodingMacChineseTrad;
        sEncodingMap[kCDXCharSetMacKorean] = kCFStringEncodingMacKorean;
        sEncodingMap[kCDXCharSetMacArabic] = kCFStringEncodingMacArabic;
        sEncodingMap[kCDXCharSetMacHebrew] = kCFStringEncodingMacHebrew;
        sEncodingMap[kCDXCharSetMacGreek] = kCFStringEncodingMacGreek;
        sEncodingMap[kCDXCharSetMacCyrillic] = kCFStringEncodingMacCyrillic;
        sEncodingMap[kCDXCharSetMacDevanagari] = kCFStringEncodingMacDevanagari;
        sEncodingMap[kCDXCharSetMacGurmukhi] = kCFStringEncodingMacGurmukhi;
        sEncodingMap[kCDXCharSetMacGujarati] = kCFStringEncodingMacGujarati;
        sEncodingMap[kCDXCharSetMacOriya] = kCFStringEncodingMacOriya;
        sEncodingMap[kCDXCharSetMacBengali] = kCFStringEncodingMacBengali;
        sEncodingMap[kCDXCharSetMacTamil] = kCFStringEncodingMacTamil;
        sEncodingMap[kCDXCharSetMacTelugu] = kCFStringEncodingMacTelugu;
        sEncodingMap[kCDXCharSetMacKannada] = kCFStringEncodingMacKannada;
        sEncodingMap[kCDXCharSetMacMalayalam] = kCFStringEncodingMacMalayalam;
        sEncodingMap[kCDXCharSetMacSinhalese] = kCFStringEncodingMacSinhalese;
        sEncodingMap[kCDXCharSetMacBurmese] = kCFStringEncodingMacBurmese;
        sEncodingMap[kCDXCharSetMacKhmer] = kCFStringEncodingMacKhmer;
        sEncodingMap[kCDXCharSetMacThai] = kCFStringEncodingMacThai;
        sEncodingMap[kCDXCharSetMacLao] = kCFStringEncodingMacLaotian;
        sEncodingMap[kCDXCharSetMacGeorgian] = kCFStringEncodingMacGeorgian;
        sEncodingMap[kCDXCharSetMacArmenian] = kCFStringEncodingMacArmenian;
        sEncodingMap[kCDXCharSetMacSimpChinese] = kCFStringEncodingMacChineseSimp;
        sEncodingMap[kCDXCharSetMacTibetan] = kCFStringEncodingMacTibetan;
        sEncodingMap[kCDXCharSetMacMongolian] = kCFStringEncodingMacMongolian;
        sEncodingMap[kCDXCharSetMacEthiopic] = kCFStringEncodingMacEthiopic;
        sEncodingMap[kCDXCharSetMacCentralEuroRoman] = kCFStringEncodingMacCentralEurRoman;
        sEncodingMap[kCDXCharSetMacVietnamese] = kCFStringEncodingMacVietnamese;
        sEncodingMap[kCDXCharSetMacExtArabic] = kCFStringEncodingMacExtArabic;
        sEncodingMap[kCDXCharSetMacSymbol] = kCFStringEncodingMacSymbol;
        sEncodingMap[kCDXCharSetMacDingbats] = kCFStringEncodingMacDingbats;
        sEncodingMap[kCDXCharSetMacCroatian] = kCFStringEncodingMacCroatian;
        sEncodingMap[kCDXCharSetMacRomanian] = kCFStringEncodingMacRomanian;
        sEncodingMap[kCDXCharSetMacCeltic] = kCFStringEncodingMacCeltic;
        sEncodingMap[kCDXCharSetMacGaelic] = kCFStringEncodingMacGaelic;
        sEncodingMap[kCDXCharSetMacIcelandic] = kCFStringEncodingMacIcelandic;
        sEncodingMap[kCDXCharSetMacTurkish] = kCFStringEncodingMacTurkish;
        sEncodingMap[kCDXCharSetUTF8] = kCFStringEncodingUTF8;
    }
    
    return sEncodingMap;
}

/**
 *  Get a CFEncoding for a CDX encoding
 */
static CFStringEncoding CDXToCFEncoding(CDXCharSet inCDXEncoding)
{
    const map<CDXCharSet,CFStringEncoding>::const_iterator i = CDXToCFEncodingMap().find(inCDXEncoding);
    if (i == CDXToCFEncodingMap().end())
    {
        return CFStringGetSystemEncoding();
    }
    else
    {
        return i->second;
    }
}

/**
 *  Convert inString from one encoding to another
 *
 *  @param inString              String to convert in inSourceEncoding
 *  @param inSourceEncoding      Source encoding
 *  @param inDestinationEncoding Destination encoding
 *
 *  @return inString encoded in the destination encoding
 */
#if UTF8_STD_STRING
string TranscodeString(const std::string& inString, CDXCharSet inSourceEncoding, CDXCharSet inDestinationEncoding)
{
    if (inSourceEncoding == inDestinationEncoding)
    {
        // Nothing to do
        return inString;
    }

    CFStringEncoding cfSourceEncoding = CDXToCFEncoding(inSourceEncoding);
    CFStringEncoding cfDestinationEncoding = CDXToCFEncoding(inDestinationEncoding);
    
    if ((cfSourceEncoding == kCFStringEncodingInvalidId) || (cfSourceEncoding == kCFStringEncodingInvalidId))
    {
        // Nothing we can do
        return inString;
    }

    // Create a cfString from the source string and encoding
    CFStringRef cfString = CFStringCreateWithCStringNoCopy(kCFAllocatorDefault, inString.c_str(), cfSourceEncoding, kCFAllocatorNull);
    if (cfString == NULL)
    {
        // Nothing we can do
        return inString;
    }
    
    // Convert this CFString to the new encoding
    CFIndex cfStringLength = CFStringGetLength(cfString);
    CFIndex destinationMaxLength = CFStringGetMaximumSizeForEncoding(cfStringLength, cfDestinationEncoding);
    std::string destinationString;
    destinationString.resize(destinationMaxLength);
    CFIndex destinationLength = 0;
    CFStringGetBytes(cfString, CFRangeMake(0,cfStringLength), cfDestinationEncoding, 0, false, (UInt8*)destinationString.c_str(), destinationMaxLength, &destinationLength);
    destinationString.resize(destinationLength);
    
    CFRelease(cfString);
    
    return destinationString;
}
#endif
#elif TARGET_OS_WIN32 // TARGET_OS_MAC

/**
 *  Convert inString from one encoding to another
 *
 *  @param inString              String to convert in inSourceEncoding
 *  @param inSourceEncoding      Source encoding
 *  @param inDestinationEncoding Destination encoding
 *
 *  @return inString encoded in the destination encoding
 */
#if UTF8_STD_STRING
string TranscodeString(const string& inString, CDXCharSet inSourceEncoding, CDXCharSet inDestinationEncoding)
{
    if ((inSourceEncoding == inDestinationEncoding) || inString.empty())
    {
        // Nothing to do
        return inString;
    }
   
    // 1. Convert the source string to UTF-16 using given source encoding
    int wideCharsLength = ::MultiByteToWideChar(inSourceEncoding, MB_ERR_INVALID_CHARS, inString.data(), (int)inString.size(), 0, 0);
    ASSERT(wideCharsLength != 0);
    if (wideCharsLength == 0)
    {
        // Nothing we can do
        return inString;
    }

    vector<wchar_t> wideCharBuffer(wideCharsLength + 1);
    wideCharsLength = ::MultiByteToWideChar(inSourceEncoding, MB_ERR_INVALID_CHARS, inString.data(), (int)inString.size(), wideCharBuffer.data(), wideCharsLength);
    wideCharBuffer[wideCharsLength] = '\0';
    
    // 2. Convert this UTF-16 string to the destination encoding
    if (inDestinationEncoding == kCDXCharSetUTF8)
    {
        // This typically used when open CDX file and convert the legacy string to UTF-8
        return ConvertUTF16to8(UTF16_string((UINT16*)wideCharBuffer.data(), wideCharsLength));
    }
    else
    {
        // This typically used when to save CDX file and keep the strings in the legacy format for the backward compatibility
        int length = ::WideCharToMultiByte(inDestinationEncoding, 0, wideCharBuffer.data(), -1, NULL, 0, NULL, NULL);
        ASSERT(length != 0);
        if (length == 0)
        {
            return inString;
        }

        vector<char> multiBytecharBuffer(length);
        ::WideCharToMultiByte(inDestinationEncoding, 0, wideCharBuffer.data(), -1, multiBytecharBuffer.data(), length, NULL, NULL);
        string strBuffer(multiBytecharBuffer.data());

        return strBuffer;
    }
}
#endif
#elif TARGET_WEB // TARGET_OS_WIN32

namespace 
{
    const char* const TextDecoderEncodingFromCDXEncoding(CDXCharSet inCDXEncoding)
    {
        // See https://developer.mozilla.org/en-US/docs/Web/API/TextDecoder/TextDecoder and
        // https://docs.microsoft.com/en-us/windows/win32/intl/code-page-identifiers for supported 
        // encodings and mappings.
        switch (inCDXEncoding)
        {
            case kCDXCharSetUTF8:
                return "utf-8";
            break;

            case kCDXCharSetMacRoman:
                return "macintosh";
            break;

            case kCDXCharSetThai:
                return "windows-874";
            break;

            case kCDXCharSetWin31EasternEuropean:
                return "windows-1250";
            break;

            case kCDXCharSetWin31Cyrillic:
                return "windows-1251";
            break;

            case kCDXCharSetUnknown:
            case kCDXCharSetWin31Latin1:
                return "windows-1252";
            break;

            case kCDXCharSetWin31Greek:
                return "windows-1253";
            break;

            case kCDXCharSetWin31Turkish:
                return "windows-1254";
            break;

            case kCDXCharSetHebrew:
                return "windows-1255";
            break;

            case kCDXCharSetArabic:
                return "windows-1256";
            break;

            case kCDXCharSetBaltic:
                return "windows-1257";
            break;

            case kCDXCharSetVietnamese:
                return "windows-1258";
            break;

            case kCDXCharSetMacCyrillic:
                return "x-mac-cyrillic";
            break;

            case kCDXCharSetChineseSimplified:
                return "gb2312";
            break;

            case kCDXCharSetJapanese:
                return "shift-jis";
            break;    
            
            case kCDXCharSetChineseTraditional:
                return "big5";
            break;           

            case kCDXCharSetKorean:
                return "ks_c_5601-1987";
            break;

            case kCDXCharSetMSDOSRussian:
                return "cp866";
            break;
        }

        // Unsupported encoding
        return nullptr;
    }

    EM_JS(void, TranscodeToUTF8, (const char* inputPtr, size_t inputBufferSize, const char* outputPtr, size_t outputBufferSize, const char* inSourceEncodingUTF8), {
        // Note that the code in this function is written in Javascript! Pointers are all indexes into
        // HEAPU8.
        //
        // This will throw on IE as TextDecoder is not available.

        try 
        {
            // Warning, hackery ahead! We need to get the input string into an ArrayBuffer. To do this
            // we access the Emscripten heap directly and extract a subarray that corresponds to the part
            // we want. See https://github.com/emscripten-core/emscripten/blob/7c3ced64d51be643f7a12bb34e8b875155d49dc7/src/runtime_strings.js#L40
            // for an example of where this is done in Emscripten.
            const sourceBuffer = HEAPU8.subarray(inputPtr, inputPtr+inputBufferSize);

            // Convert from the source encoding to a JS string
            const sourceEncoding = UTF8ToString(inSourceEncodingUTF8);
            const textDecoder = new TextDecoder(sourceEncoding);
            const jsString = textDecoder.decode(sourceBuffer);

            // Now convert from a JS string to UTF-8
            stringToUTF8(jsString, outputPtr, outputBufferSize);
        }
        catch (e)
        {
            // If anything goes wrong then just return, our caller will see an empty string.
            return;
        }
    });
}

#if UTF8_STD_STRING
string TranscodeString(const string& inString, CDXCharSet inSourceEncoding, CDXCharSet inDestinationEncoding)
{
    if ((inSourceEncoding == inDestinationEncoding) || inString.empty())
    {
        // Nothing to do
        return inString;
    }

    const char* const textDecoderDestinationEncoding = TextDecoderEncodingFromCDXEncoding(inDestinationEncoding);
    if (textDecoderDestinationEncoding == nullptr)
    {
        // Unsuppoted destination encoding, just return the original string
        return inString;
    }

    const char* const textDecoderSourceEncoding = TextDecoderEncodingFromCDXEncoding(inSourceEncoding);
    if (textDecoderSourceEncoding == nullptr)
    {
        // Unsupported source encoding, just return the original string
        return inString;
    }

    // Make a buffer to hold the output string. This has a capacity of 4 times the size of the input string
    // plus a null terminator. This will be big enough to hold any encoded output.
    string transcodedString;
    transcodedString.resize((inString.length() * 4) + 1);

    // Call into JS to transcode the string
    TranscodeToUTF8(inString.c_str(), inString.length(), transcodedString.data(), transcodedString.length(), textDecoderSourceEncoding);

    // Truncate to our null-terminated size
    transcodedString.resize(strlen(transcodedString.c_str()));

    // If the string is empty something went wrong, return the original string
    if (transcodedString.empty())
    {
        return inString;
    }

    return transcodedString;
}
#endif
#elif TARGET_OS_LINUX // TARGET_WEB
/**
 *  Convert the CDXCharSet based encoding to an encoding for the iconv backend for boost locale
 * 
 *  All the non-mac based encoding are observed as supported on Linux. For Mac some of the encodings names are guessed based on
 *  their CoreFoundation constant names where an iconv name was not available.
 * 
 *  @param encoding The CDX encoding constant
 *  @return The encoding name to use with boost::locale
 */
const char* const LocaleEncodingFromCDXEncoding(CDXCharSet encoding)
{
    switch (encoding)
    {
        case kCDXCharSetEBCDICOEM:
            return "CP037";
        case kCDXCharSetMSDOSUS:
            return "CP437";
        case kCDXCharSetEBCDIC500V1:
            return "CP500";
        case kCDXCharSetArabicASMO708:
            return "ASMO-708";
        case kCDXCharSetArabicASMO449P:
            return "ASMO_449";
        case kCDXCharSetArabicTransparent:
            return "ASMO-708";
        case kCDXCharSetArabicTransparentASMO:
            return "ASMO-708";
        case kCDXCharSetGreek437G:
            return "CP737";
        case kCDXCharSetBalticOEM:
            return "CP775";
        case kCDXCharSetMSDOSLatin1:
            return "CP850";
        case kCDXCharSetMSDOSLatin2:
            return "CP852";
        case kCDXCharSetIBMCyrillic:
            return "CP855";
        case kCDXCharSetIBMTurkish:
            return "CP857";
        case kCDXCharSetMSDOSPortuguese:
            return "CP860";
        case kCDXCharSetMSDOSIcelandic:
            return "CP861";
        case kCDXCharSetHebrewOEM:
            return "CP862";
        case kCDXCharSetMSDOSCanadianFrench:
            return "CP863";
        case kCDXCharSetArabicOEM:
            return "CP864";
        case kCDXCharSetMSDOSNordic:
            return "CP865";
        case kCDXCharSetMSDOSRussian:
            return "CP866";
        case kCDXCharSetIBMModernGreek:
            return "CP869";
        case kCDXCharSetThai:
            return "CP874";
        case kCDXCharSetEBCDIC:
            return "CP875";
        case kCDXCharSetJapanese:
            return "CP932";
        case kCDXCharSetChineseSimplified:
            return "CP936";
        case kCDXCharSetKorean:
            return "CP949";
        case kCDXCharSetChineseTraditional:
            return "CP950";
        case kCDXCharSetUnicodeISO10646:
            return "UTF-16";
        case kCDXCharSetWin31EasternEuropean:
            return "WINDOWS-1250";
        case kCDXCharSetWin31Cyrillic:
            return "WINDOWS-1251";
        case kCDXCharSetWin31Latin1:
        case kCDXCharSetUnknown:
            return "WINDOWS-1252";
        case kCDXCharSetWin31Greek:
            return "WINDOWS-1253";
        case kCDXCharSetWin31Turkish:
            return "WINDOWS-1254";
        case kCDXCharSetHebrew:
            return "WINDOWS-1255";
        case kCDXCharSetArabic:
            return "WINDOWS-1256";
        case kCDXCharSetBaltic:
            return "WINDOWS-1257";
        case kCDXCharSetVietnamese:
            return "WINDOWS-1258";
        case kCDXCharSetKoreanJohab:
            return "CP1361";
        case kCDXCharSetMacRoman:
            return "MACINTOSH";
        case kCDXCharSetMacJapanese:
            return "MACJAPANESE";
        case kCDXCharSetMacTradChinese:
            return "MACCHINESETRAD";
        case kCDXCharSetMacKorean:
            return "MACKOREAN";
        case kCDXCharSetMacArabic:
            return "MACARABIC";
        case kCDXCharSetMacHebrew:
            return "MACHEBREW";
        case kCDXCharSetMacGreek:
            return "MACGREEK";
        case kCDXCharSetMacCyrillic:
            return "MACCYRILLIC";
        case kCDXCharSetMacDevanagari:
            return "MACDEVANGARI";
        case kCDXCharSetMacGurmukhi:
            return "MACGURMUKHI";
        case kCDXCharSetMacGujarati:
            return "MACGUJARATI";
        case kCDXCharSetMacOriya:
            return "MACORIYA";
        case kCDXCharSetMacBengali:
            return "MACBENGALI";
        case kCDXCharSetMacTamil:
            return "MACTAMIL";
        case kCDXCharSetMacTelugu:
            return "MACTELUGU";
        case kCDXCharSetMacKannada:
            return "MACKANNADA";
        case kCDXCharSetMacMalayalam:
            return "MACMALAYALAM";
        case kCDXCharSetMacSinhalese:
            return "MACSINHALESE";
        case kCDXCharSetMacBurmese:
            return "MACBURMESE";
        case kCDXCharSetMacKhmer:
            return "MACKHMER";
        case kCDXCharSetMacThai:
            return "MACTHAI";
        case kCDXCharSetMacLao:
            return "MACLAO";
        case kCDXCharSetMacGeorgian:
            return "MACGEORGIAN";
        case kCDXCharSetMacArmenian:
            return "MACARMEDIAN";
        case kCDXCharSetMacSimpChinese:
            return "MACCHINESESIMP";
        case kCDXCharSetMacTibetan:
            return "MACTIBETAN";
        case kCDXCharSetMacMongolian:
            return "MACMONGOLIAN";
        case kCDXCharSetMacEthiopic:
            return "MACETHIOPIC";
        case kCDXCharSetMacCentralEuroRoman:
            return "MACCENTRALEUROPE";
        case kCDXCharSetMacVietnamese:
            return "MACVIETNAMESE";
        case kCDXCharSetMacExtArabic:
            return "MACEXTARABIC";        
        case kCDXCharSetMacSymbol:
            return "MACSYMBOL";
        case kCDXCharSetMacDingbats:
            return "MACDINGBATS";
        case kCDXCharSetMacCroatian:
            return "MACCROATIAN";
        case kCDXCharSetMacRomanian:
            return "MACROMANIAN";
        case kCDXCharSetMacCeltic:
            return "MACCELTIC";
        case kCDXCharSetMacGaelic:
            return "MACGAELIC";
        case kCDXCharSetMacKeyboardGlyphs:
        case kCDXCharSetMacIcelandic:
            return "MACICELAND";
        case kCDXCharSetMacTurkish:
            return "MACTURKISH";
        case kCDXCharSetUTF8:
            return "utf-8";

   	    case kCDXCharSetMacReserved:
        case kCDXCharSetMacUninterpreted:
        default:
            break;
    }

    return nullptr;
}

#if UTF8_STD_STRING
string TranscodeString(const string& inString, CDXCharSet inSourceEncoding, CDXCharSet inDestEncoding)
{    
    try
    {
        const char* sourceEncoding = LocaleEncodingFromCDXEncoding(inSourceEncoding);
        const char* destEncoding = LocaleEncodingFromCDXEncoding(inDestEncoding);
        if (sourceEncoding != nullptr && destEncoding != nullptr)
        {
            return boost::locale::conv::between(inString, destEncoding, sourceEncoding);
        }
    }
    catch (...)
    {
        // Fallback to the source string if transcoding fails
    }

    return inString;
}
#endif
#endif // TARGET_OS_LINUX
