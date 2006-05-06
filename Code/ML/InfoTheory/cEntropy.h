//
//  Copyright (C) 2001,2003 greg Landrum and Rational Discovery LLC
//

// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the CENTROPY_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// CENTROPY_API functions as being imported from a DLL, wheras this DLL sees symbols
// defined with this macro as being exported.

#ifdef WIN32
#ifdef CENTROPY_EXPORTS
#define CENTROPY_API extern "C" __declspec(dllexport)
#else
//#define CENTROPY_API extern "C" __declspec(dllimport)
#define CENTROPY_API extern "C"
#endif
#include <windows.h>
#else  // WIN32
#define CENTROPY_API extern "C"
#endif

CENTROPY_API void initcEntropy(void);
template<class T> extern double InfoEntropy(long int *,long int);
CENTROPY_API double InfoGain(long int *,long int,long int);


#ifdef _DEBUG
#ifndef BOOST_PYTHON_DEBUG
#undef _DEBUG
#endif
#include <Python.h>
#define _DEBUG
#else
#include <Python.h>
#endif
