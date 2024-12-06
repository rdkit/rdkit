// CommonCS/LibCommon/Hdr/CDXBasicTypes.h
// Copyright © 1999-2004, CambridgeSoft Corp., All Rights Reserved

#pragma once

// disable identifier too long warning
#pragma warning (disable: 4786)

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

