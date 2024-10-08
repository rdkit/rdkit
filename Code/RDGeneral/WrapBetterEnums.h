#pragma once

#ifdef USE_BETTER_ENUMS
#include "enum.h"
#else
#define BETTER_ENUM(Enum, Underlying, ...) enum Enum : Underlying { __VA_ARGS__ }
#endif
