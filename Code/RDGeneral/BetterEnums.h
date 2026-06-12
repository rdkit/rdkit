//
//  Copyright (C) 2024 Novartis Biomedical Research and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// NOTE: this header intentionally has no include guard / #pragma once.
// What it expands to depends on whether USE_BETTER_ENUMS is defined at the
// point of inclusion. A translation unit may include it once in the default
// (plain-enum) mode -- e.g. transitively via RDKitBase.h, including when that is
// baked into a precompiled header -- and then again after defining
// USE_BETTER_ENUMS to pull in the full better-enums implementation. An include
// guard would pin the file to whichever mode was seen first and break those
// translation units, so it must remain re-includable.

#ifdef USE_BETTER_ENUMS
// Re-entering in better-enums mode: discard any fallback macros left over from
// an earlier default-mode inclusion before enum.h supplies the real ones.
#ifdef BETTER_ENUM
#undef BETTER_ENUM
#endif
#ifdef BETTER_ENUM_CLASS
#undef BETTER_ENUM_CLASS
#endif
#include "enum.h"
#define BETTER_ENUM_CLASS BETTER_ENUM
#else
// Fallback: plain enums. Guarded so a later better-enums inclusion can redefine.
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
