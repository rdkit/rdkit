#pragma once

#if defined(USE_BETTER_ENUMS) && !defined(SWIG)
#include "BetterEnums.h"
#else
#define BETTER_ENUM(Enum, Underlying, ...) enum Enum : Underlying { __VA_ARGS__ }
#endif
