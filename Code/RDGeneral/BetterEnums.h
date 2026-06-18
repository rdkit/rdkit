//
//  Copyright (C) 2024 Novartis Biomedical Research and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// This header intentionally has no include guard / #pragma once: what it expands
// to depends on whether USE_BETTER_ENUMS is defined at the point of inclusion. A
// TU may include it in plain-enum mode (e.g. transitively via RDKitBase.h, even
// when baked into a precompiled header) and then again, after defining
// USE_BETTER_ENUMS, to pull in the full better-enums implementation; a guard
// would pin it to whichever mode was seen first.
//
// Re-includability only un-pins this header, though. The headers that declare
// the enums (e.g. MolOps.h's SanitizeFlags) are #pragma once, so once one is
// baked into the PCH in plain-enum mode its enums stay plain for every consumer.
// A TU that #defines USE_BETTER_ENUMS must therefore skip the PCH; such sources
// are marked SKIP_PRECOMPILE_HEADERS in their CMakeLists.

#ifdef USE_BETTER_ENUMS
// Gate on enum.h's own include guard (it is #pragma once): only the first
// inclusion supplies the real BETTER_ENUM. Without this, a second
// better-enums-mode inclusion in the same TU would #undef the real macro and
// then hit a no-op #include "enum.h", leaving BETTER_ENUM undefined.
#ifndef BETTER_ENUMS_ENUM_H
// Discard any fallback macros from an earlier plain-enum inclusion before
// enum.h supplies the real ones.
#ifdef BETTER_ENUM
#undef BETTER_ENUM
#endif
#ifdef BETTER_ENUM_CLASS
#undef BETTER_ENUM_CLASS
#endif
#include "enum.h"
#define BETTER_ENUM_CLASS BETTER_ENUM
#endif
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
