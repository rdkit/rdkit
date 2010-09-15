
// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the CALGORITHMS_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// CALGORITHMS_API functions as being imported from a DLL, wheras this DLL sees symbols
// defined with this macro as being exported.

#ifdef WIN32
#ifdef CALGORITHMS_EXPORTS
#define CALGORITHMS_API extern "C" __declspec(dllexport)
#else
#define CALGORITHMS_API extern "C" __declspec(dllimport)
#endif
#include <windows.h>
#else  // WIN32
#define CALGORITHMS_API extern "C"
#endif

CALGORITHMS_API void initClustering(void);


#ifdef _DEBUG
#undef _DEBUG
#include <Python.h>
#define _DEBUG
#else
#include <Python.h>
#endif

#define PY_ARRAY_UNIQUE_SYMBOL Py_Array_API_Clustering
#ifndef PYTH_FILE_WITH_INIT
  #define NO_IMPORT_ARRAY
#endif



#include <stdlib.h>

