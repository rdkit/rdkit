// CommonCS/LibCommon/Src/Base64Converter.cpp
// Copyright (c) 1997-2004, CambridgeSoft Corp., All Rights Reserved

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

#include "cs_univDefs.h"
#include "cs_assert.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <algorithm>
#include "Base64Converter.h"
using namespace std;

#ifdef _AUTHORIZED
static bool bAuth = true;
#else
static bool bAuth = false;
#endif

// sBase64 has 64 characters
const int kBaseSize = 64;
static const unsigned char s_Base64 [kBaseSize + 1] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";


// Checks if the given character string is Net variant of base 64 encoding.
bool  csIsNetVariant (const char *s, long len)
{
	return s && len >= 4 && s[0] == '#' && isdigit(s[1]) && isdigit(s[2]) && isdigit(s[3]);
}

//-------------------------------------------------
// Checks if the given character string is valid base 64 encoding.
// Checks only upto kMaxCheckLength bytes.
const long kMaxCheckLength = 256;
bool  csIsValidBase64 (const char *s, long len)
{
	if (s == NULL  ||  len <= 0)
		return false;
	if (csIsNetVariant (s, len))
	{
		s += 4;
		len -= 4;
	}
	long	checkLen = len;
	if (checkLen > kMaxCheckLength)
		checkLen = kMaxCheckLength;
	for ( ; checkLen > 0;  checkLen--, s++)
		if (!isspace((unsigned char)*s) && !isalnum((unsigned char)*s) && *s != '+' &&
				*s != '/'  &&  *s != '='  &&  *s != '\r'  &&  *s != '\n')
			return false;
	return true;
}
//-------------------------------------------------
static bool BuildBase(const char *s, long len, unsigned char *&pBase)
{
	pBase = NEW unsigned char[kBaseSize];
	memcpy(pBase, s_Base64, kBaseSize);

	if (csIsNetVariant (s, len))
	{
		int i, j, k, l;
		for (i = 1; i <= 3; ++i)
		{
			k = i;
			j = s[i] + k;
			l = 20;
			while (j != i && l-- > 0)
			{
				std::swap (*(pBase + j), *(pBase + k));
				k = j + i;
				j = k + s[i];
				k = k % kBaseSize;
				j = j % kBaseSize;
			}
		}
	}
	return true;
}
//----------------------------------------------------------
// ToBinary() maps a legal ascii character to a 6-bit binary number.
//				Illegal characters are mapped to 255.
//----------------------------------------------------------
static unsigned char ToBinary(const unsigned char base64)
{
	if ((base64 >= 'A') && (base64 <= 'Z'))
		return (unsigned char)(base64 - 'A');
	else if ((base64 >= 'a') && (base64 <= 'z'))
		return (unsigned char)(26 + base64 - 'a');
	else if ((base64 >= '0') && (base64 <= '9' ))
		return (unsigned char)(52 + base64 - '0');
	else if (base64 == '+')
		return (unsigned char)(62);
	else if (base64 == '/')
		return (unsigned char)(63);
	else
		return (unsigned char)(255);
}
//---------------------------------------------------------------
// "pIndex" maps binary 6-bit values to the index (0-63) of the Ascii character assigned to it.
static bool BuildIndex(const char *s, long len, int *&pIndex)
{
	unsigned char *pBase = NULL;
	if (!BuildBase(s, len, pBase))
		return false;
	pIndex = NEW int[kBaseSize];
	unsigned char *pc = pBase;
	for (int i = 0; i < kBaseSize; ++i, ++pc)
		pIndex[ToBinary(*pc)] = i;
	delete [] pBase;
	return true;
}
//--------------------------------------
long  csBase64Encode_Ext (const char	*input,
						  long			inputLen,
						  char			*&output,
						  bool			bNetVariant,
						  unsigned int	lineBreakFrequency /* = 0 */)
{
	long outputLen = inputLen / 3;

	if (inputLen % 3)
		outputLen += 1;

	int iKey = 0;
	if (bNetVariant)
		outputLen += 1;

	outputLen *= 4;
	lineBreakFrequency &= 0xFFFC;	// ensure it is a multiple of 4: clear last two bits
	if (lineBreakFrequency)
	{
		long nlines = (outputLen + lineBreakFrequency - 1) / lineBreakFrequency;
		outputLen += nlines * 2; // add room for cr and lf
	}
	output = NEW char[outputLen + 1];

	const unsigned char *pi = (const unsigned char *) input;
	unsigned char *po = (unsigned char *) output;

	unsigned int nCharsSinceLineBreak = 0;
	if (bNetVariant)	// Write a 4-byte random header
	{
		time_t t;
		srand(static_cast<unsigned int >(time(&t)));
		// Key requires 3 digits.
		while (iKey <= 100) iKey = rand() % 1000;
		*po++ = '#';
		while (iKey > 0)
		{
			*po++ = (unsigned char)(iKey % 10 + '0');
			iKey /= 10;
		}
		nCharsSinceLineBreak = 4;
	}

	unsigned char *pBase = NULL;
	if (!BuildBase(output, bNetVariant ? 4 : 0, pBase))
		return 0;

	long i;
	for (i = inputLen - 3; i >= 0; i -= 3, pi += 3)
	{
		*po++ = pBase[pi[0] / 4];
		*po++ = pBase[((pi[0] & 0x03) << 4) | (pi[1] / 16)];
		*po++ = pBase[((pi[1] & 0x0F) << 2) | (pi[2] / 64)];
		*po++ = pBase[pi[2] & 0x3F];
		nCharsSinceLineBreak += 4;
		if (lineBreakFrequency  &&  nCharsSinceLineBreak >= lineBreakFrequency)
		{
			*po++ = '\015';
			*po++ = '\012';
			nCharsSinceLineBreak = 0;
		}
	}
	if (i == -1 || i == -2)
	{
		*po++ = pBase[pi[0] / 4];
		if (i == -1)
		{
			*po++ = pBase[((pi[0] & 0x03) << 4) | (pi[1] / 16)];
			*po++ = pBase[(pi[1] & 0x0F) << 2];
		}
		else
		{
			*po++ = pBase[(pi[0] & 0x03) << 4];
			*po++ = '=';
		}
		*po++ = '=';
		nCharsSinceLineBreak += 4;
	}
	if ((lineBreakFrequency != 0) && (nCharsSinceLineBreak != 0))
	{
		*po++ = '\015';
		*po++ = '\012';
	}
	*po = '\0';
	ASSERT ((char*)po - output == outputLen);
	delete [] pBase;
	return outputLen;
}

// Decodes the input base64 string.
// If bNetVariantOnly flag is set, decodes only Net-variants.
// Otherwise, decodes any valid string.
// A server handling encoded CDX data (Pro/Std/Net) should not set bNetVariant flag.
// A server handling both plain (Pro/Std) and encoded (Net) Mol or SMILES data should set bNetVariant flag.
long  csBase64Decode(const char *input, long inputLen, char *&output, bool bNetVariantOnly)
{
	// If it is not a valid string (i.e., plain mol file), return immediately.
	if (!csIsValidBase64(input, inputLen))
		return 0;

	// Calculate output length and allocate output string.
	// Note that the output string is not null terminated.
	const unsigned char *pi = (const unsigned char *) input;
	bool isNetVariant = csIsNetVariant (input, inputLen);
	// If it is not a desired type, return immediately.
	if ((bNetVariantOnly && !isNetVariant) || (isNetVariant && !bAuth))
		return 0;
	if (isNetVariant && bAuth)
	{
		inputLen -= 4;
		pi += 4;
	}

	int *pIndex = NULL;
	if (!BuildIndex(input, inputLen, pIndex))
		return 0;
	
	// Since we skip spaces, the input length is no longer a multiple of 4.
	// The output size is set to theoretic maximum size.
	long outputLen = (inputLen + 3) / 4 * 3;
	output = NEW char[outputLen + 1];
	unsigned char *po = (unsigned char *) output;
	unsigned char c[4];
	int j = 0;
	for (long i = 0; i < inputLen; ++i)
	{
		c[j] = ToBinary(*pi++);	// low-order six bits are defined
		// If it is an invalid character, skip it.
		if (c[j] == 255)
			continue;
		c[j] = (unsigned char)(pIndex[c[j]]);
		++j;
		if (j == 4)
		{
			*po++ = (unsigned char)((c[0] << 2) | ((c[1] & 0x30) / 16));		// c0(5-0) plus c1(5-4)
			*po++ = (unsigned char)(((c[1] & 0x0F) << 4) | ((c[2] & 0x3C) / 4));// c1(3-0) plus c2(5-2)
			*po++ = (unsigned char)(((c[2] & 0x03) << 6) | c[3]);				// c2(1-0) plus c3(5-0)
			j = 0;
		}
	}
	if (j != 0)
	{
		if (j > 1)
			*po++ = (unsigned char)((c[0] << 2) | ((c[1] & 0x30) / 16));
		if (j > 2)
			*po++ = (unsigned char)(((c[1] & 0x0F) << 4) | ((c[2] & 0x3C) / 4));
	}
	*po = '\0';

	// Output length does not include the terminating null character.
	outputLen = int (po - (unsigned char *) output);

	// Dispose the index.
	delete [] pIndex;
	return outputLen;
}

// Make the caller authorized to decoded net-encripted b64.
// Obviously, net variant should never call this.
void csSetBase64DecodeAuthorization (bool bAuthorized)
{
	bAuth = bAuthorized;
}
