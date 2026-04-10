//
//  Copyright (C) 2024 Novartis Biomedical Research and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifdef USE_BETTER_ENUMS
#ifdef BETTER_ENUM
#undef BETTER_ENUM
#endif
#ifdef BETTER_ENUM_CLASS
#undef BETTER_ENUM_CLASS
#endif
#include "enum.h"
#define BETTER_ENUM_CLASS BETTER_ENUM
#else
#ifndef BETTER_ENUM
#define BETTER_ENUM(Enum, Underlying, ...) \
  enum Enum : Underlying {                 \
    __VA_ARGS__                            \
  }
#endif

#ifndef BETTER_ENUM_CLASS
#define BETTER_ENUM_CLASS(Enum, Underlying, ...) \
  enum class Enum : Underlying {                 \
    __VA_ARGS__                                  \
  }
#endif

#endif