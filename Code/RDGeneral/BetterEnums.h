#pragma once

#ifdef USE_BETTER_ENUMS
#include "enum.h"
#define BETTER_ENUM_CLASS BETTER_ENUM
#else
#define BETTER_ENUM(Enum, Underlying, ...) enum Enum : Underlying { __VA_ARGS__ }
#define BETTER_ENUM_CLASS(Enum, Underlying, ...) enum class Enum : Underlying { __VA_ARGS__ }
#endif
