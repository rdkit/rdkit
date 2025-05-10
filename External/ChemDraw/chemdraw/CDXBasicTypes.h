// CommonCS/LibCommon/Hdr/CDXBasicTypes.h
// Copyright © 1999-2004, CambridgeSoft Corp., All Rights Reserved

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

// disable identifier too long warning
#if TARGET_OS_WIN32
#pragma warning (disable: 4786)
#endif

// IN VC++ 6.0, INT32 and UINT32 are pre-defined in basetsd.h.
// For other compilers, we need to define them.
#if TARGET_OS_WIN32  &&  defined(_MSC_VER)  &&  _MSC_VER >= 1200
#include <basetsd.h>
#else
	typedef long INT32;
	typedef unsigned long UINT32;
#endif // TARGET_OS_WIN32

#if !TARGET_OS_WIN32  ||  !defined(_MSC_VER)  ||  _MSC_VER <= 1200
typedef char			INT8;
#endif

typedef short			INT16;
typedef unsigned char	UINT8;
typedef unsigned short	UINT16;

// I don't think there are any platforms where this is wrong
// There once were Mac compilers where double was 10 or 12 bytes
// and short double is 8.
// This is checked by an assert in CDX.cpp but that may only be
// tripped when reading or writing spectra.
//  - SDR 5/18/98
typedef float FLOAT32;
typedef double FLOAT64;


#ifndef FILELINE
#	define FILELINE
#endif

#ifndef NEW
	#define NEW new
#endif

