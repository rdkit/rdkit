
// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the CQUANTIZE_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// CQUANTIZE_API functions as being imported from a DLL, wheras this DLL sees symbols
// defined with this macro as being exported.

#ifdef WIN32
#ifdef CQUANTIZE_EXPORTS
#define CQUANTIZE_API extern "C" __declspec(dllexport)
#else
#define CQUANTIZE_API extern "C" __declspec(dllimport)
#endif
#include <windows.h>
#else  // WIN32
#define CQUANTIZE_API extern "C"
#endif

CQUANTIZE_API void initcQuantize(void);

#ifdef _DEBUG
#ifndef BOOST_PYTHON_DEBUG
#undef _DEBUG
#endif
#include <Python.h>
#define _DEBUG
#else
#include <Python.h>
#endif
