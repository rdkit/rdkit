#pragma once

#ifndef STATIC_MAEPARSER

#ifdef WIN32
#ifdef IN_MAEPARSER
#define EXPORT_MAEPARSER __declspec(dllexport)
#else
#define EXPORT_MAEPARSER __declspec(dllimport)
#endif // IN_MAEPARSER

#else

#define EXPORT_MAEPARSER __attribute__((visibility("default")))
#endif // WIN32

#else

#define EXPORT_MAEPARSER

#endif // STATIC_MAEPARSER
