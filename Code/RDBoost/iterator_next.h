#include <RDGeneral/export.h>
#ifndef RDKIT_ITERATOR_NEXT_INCLUDED
#define RDKIT_ITERATOR_NEXT_INCLUDED

#if PY_MAJOR_VERSION >= 3
const char* const NEXT_METHOD = "__next__";
#else
const char* const NEXT_METHOD = "next";
#endif

#endif
