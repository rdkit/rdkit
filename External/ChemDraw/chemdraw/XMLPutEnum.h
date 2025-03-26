// CommonCS/LibCommon/Hdr/XMLPutEnum.h
// Contains: Class used to convert xml attribute values to binary CDX enum constants
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

#if TARGET_OS_WIN32
// Don't warn about truncation of long symbols
#pragma warning (disable:4786)
// Don't warn about unknown pragmas
#pragma warning (disable:4068)
#endif

#include <algorithm>
#include <functional>
#include <string>

#include <stdexcept>
#ifdef __linux
#include "cs_univDefs.h"
#endif

using std::runtime_error;
using std::ostringstream;
using std::istringstream;
using std::stringstream;
using std::ostream;
using std::ofstream;
using std::streamsize;
using std::exception;
using std::nothrow_t;
using std::bad_alloc;
using std::swap;
using std::make_pair;
using std::fill_n;
using std::ios;
using std::multimap;
using std::find;
using std::streambuf;
using std::vector;
using std::ostream;
using std::min;
using std::max;
using std::string;

// This class is used to implement a lookup table relating an enum to the text names
template <class T>
class XMLPutEnum {
public:
	struct Value { T value; const char *name; };
		// This is used to pass in the values for the tokens

	XMLPutEnum(const Value *v, int sizeOfValues, T deflt) : m_values(v), m_numValues(sizeOfValues / sizeof (Value)), m_default(deflt) {}
		// constructor - note sizeOfValues is the total size (e.g. sizeof valueArray), not the number of values

	const char *lookup(T v);
		// look up the XML token for a specific value of this type

	T lookup(const std::string &s);
		// look up the value for an XML token of this type

	std::string asTokens(T v);
		// look up the tokens for the bits in a bitmask

	T lookupTokens(const std::string &s);
		// look up the combined value for a series of XML tokens, where the values are bits in a bitmask

private:
	struct FindByValue {
		T m_value;
		FindByValue (T val) : m_value(val) {}
		bool operator() (const Value &v) const { return m_value == v.value; }
	};

	struct FindByToken {
		std::string m_value;
		FindByToken (const std::string &val) : m_value(val) {}
		bool operator() (const Value &v) const { return m_value == v.name; }
	};

	const Value *m_values;
	int m_numValues;
	T m_default;
};

template <class T>
const char *XMLPutEnum<T>::lookup(T v)
{
	const Value *result = std::find_if(m_values, m_values + m_numValues, FindByValue(v));
	if (result == m_values + m_numValues)
	{
		result = std::find_if(m_values, m_values + m_numValues, FindByValue(m_default));
		if (result == m_values + m_numValues)
		{
			throw std::runtime_error("Unknown value for attribute");
		}
	}
	return result->name;
}

template <class T>
T XMLPutEnum<T>::lookup(const std::string &s)
{
	const Value *result = std::find_if(m_values, m_values + m_numValues, FindByToken(s));
	if (result == m_values + m_numValues)
		return m_default; // throw std::runtime_error("Unknown value for attribute");
	return result->value;
}

template <class T>
T XMLPutEnum<T>::lookupTokens(const std::string &s)
{
	// Break the attribute value up into space-delimited tokens, look up each one, and combine them into a bitmask
	T result = T(0);
	std::string::size_type startPos = 0;
	while ((startPos = s.find_first_not_of(' ', startPos)) != std::string::npos)
	{
		std::string::size_type	nextSpace = s.find_first_of(' ', startPos);
		result |= lookup(s.substr(startPos, nextSpace-startPos));
		startPos = nextSpace;
	}
	return result;
}

template <class T>
std::string XMLPutEnum<T>::asTokens(T v)
{
	unsigned short intVal = (unsigned short)(v);
	unsigned short i = 1;
	bool first = true;
	std::string result;
	
	// If there are no bits on, and there's a value for that, use that string rather than returning the empty string
	if (v == T(0)
		&& std::find_if(m_values, m_values + m_numValues, FindByValue(T(0))) != m_values + m_numValues)
		return lookup(T(0));
	
	// Go throught the bits one at a time, removing them from the bitmask, until it's empty
	while (intVal != 0)
	{
		if ((intVal & i) != 0)
		{
			if (!first)
				result += " ";
			first = false;
			result += lookup(T(i));
			intVal &= ~i;
		}
		i <<= 1;
	}
	return result;
}

