
// Copyright 2005-2009 Daniel James.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#if !defined(GBOOST_FUNCTIONAL_HASH_DETAIL_HASH_FLOAT_HEADER)
#define GBOOST_FUNCTIONAL_HASH_DETAIL_HASH_FLOAT_HEADER

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

#if defined(BOOST_MSVC)
#pragma warning(push)
#if BOOST_MSVC >= 1400
#pragma warning(disable:6294) // Ill-defined for-loop: initial condition does
                              // not satisfy test. Loop body not executed 
#endif
#endif

#include <RDGeneral/hash/hash_fwd.hpp>
#include <RDGeneral/hash/detail/float_functions.hpp>
#include <boost/integer/static_log2.hpp>
#include <boost/cstdint.hpp>
#include <boost/limits.hpp>
#include <boost/assert.hpp>

// Select implementation for the current platform.

// Cygwn
#if defined(__CYGWIN__)
#  if defined(__i386__) || defined(_M_IX86)
#    define GBOOST_HASH_USE_x86_BINARY_HASH
#  endif

// STLport
#elif defined(__SGI_STL_PORT) || defined(_STLPORT_VERSION)
// fpclassify aren't good enough on STLport.

// GNU libstdc++ 3
#elif defined(__GLIBCPP__) || defined(__GLIBCXX__)
#  if (defined(__USE_ISOC99) || defined(_GLIBCXX_USE_C99_MATH)) && \
      !(defined(macintosh) || defined(__APPLE__) || defined(__APPLE_CC__))
#    define GBOOST_HASH_USE_FPCLASSIFY
#  endif

// Dinkumware Library, on Visual C++ 
#elif (defined(_YVALS) && !defined(__IBMCPP__)) || defined(_CPPLIB_VER)

// Not using _fpclass because it is only available for double.

#endif

// On OpenBSD, numeric_limits is not reliable for long doubles, but
// the macros defined in <float.h> are.

#if defined(__OpenBSD__)
#include <float.h>
#endif

namespace gboost
{
    namespace hash_detail
    {
        template <class T>
        struct limits : std::numeric_limits<T> {};

#if defined(__OpenBSD__)
        template <>
        struct limits<long double>
             : std::numeric_limits<long double>
        {
            static long double epsilon() {
                return LDBL_EPSILON;
            }

            static long double (max)() {
                return LDBL_MAX;
            }

            static long double (min)() {
                return LDBL_MIN;
            }

            BOOST_STATIC_CONSTANT(int, digits = LDBL_MANT_DIG);
            BOOST_STATIC_CONSTANT(int, max_exponent = LDBL_MAX_EXP);
            BOOST_STATIC_CONSTANT(int, min_exponent = LDBL_MIN_EXP);
        };
#endif // __OpenBSD__

        inline void hash_float_combine(std::hash_result_t& seed, std::hash_result_t value)
        {
            seed ^= value + (seed<<6) + (seed>>2);
        }

// A simple, non-portable hash algorithm for x86.
#if defined(GBOOST_HASH_USE_x86_BINARY_HASH)
        inline std::hash_result_t float_hash_impl(float v)
        {
            boost::uint32_t* ptr = (boost::uint32_t*)&v;
            std::hash_result_t seed = *ptr;
            return seed;
        }

        inline std::hash_result_t float_hash_impl(double v)
        {
            boost::uint32_t* ptr = (boost::uint32_t*)&v;
            std::hash_result_t seed = *ptr++;
            hash_float_combine(seed, *ptr);
            return seed;
        }

        inline std::hash_result_t float_hash_impl(long double v)
        {
            boost::uint32_t* ptr = (boost::uint32_t*)&v;
            std::hash_result_t seed = *ptr++;
            hash_float_combine(seed, *ptr++);
            hash_float_combine(seed, *(boost::uint16_t*)ptr);
            return seed;
        }

#else

        template <class T>
        inline std::hash_result_t float_hash_impl(T v)
        {
            int exp = 0;

            v = gboost::hash_detail::call_frexp(v, &exp);

            // A postive value is easier to hash, so combine the
            // sign with the exponent.
            if(v < 0) {
                v = -v;
                exp += limits<T>::max_exponent -
                    limits<T>::min_exponent;
            }

            // The result of frexp is always between 0.5 and 1, so its
            // top bit will always be 1. Subtract by 0.5 to remove that.
            v -= T(0.5);
            v = gboost::hash_detail::call_ldexp(v,
                    limits<std::hash_result_t>::digits + 1);
            std::hash_result_t seed = static_cast<std::hash_result_t>(v);
            v -= seed;

            // ceiling(digits(T) * log2(radix(T))/ digits(hash_result_t)) - 1;
            std::hash_result_t const length
                = (limits<T>::digits *
                        boost::static_log2<limits<T>::radix>::value - 1)
                / limits<std::hash_result_t>::digits;

            for(std::hash_result_t i = 0; i != length; ++i)
            {
                v = gboost::hash_detail::call_ldexp(v,
                        limits<std::hash_result_t>::digits);
                std::hash_result_t part = static_cast<std::hash_result_t>(v);
                v -= part;
                hash_float_combine(seed, part);
            }

            hash_float_combine(seed, exp);

            return seed;
        }
#endif

        template <class T>
        inline std::hash_result_t float_hash_value(T v)
        {
#if defined(GBOOST_HASH_USE_FPCLASSIFY)
            using namespace std;
            switch (fpclassify(v)) {
            case FP_ZERO:
                return 0;
            case FP_INFINITE:
                return (std::hash_result_t)(v > 0 ? -1 : -2);
            case FP_NAN:
                return (std::hash_result_t)(-3);
            case FP_NORMAL:
            case FP_SUBNORMAL:
                return float_hash_impl(v);
            default:
                BOOST_ASSERT(0);
                return 0;
            }
#else
            return v == 0 ? 0 : float_hash_impl(v);
#endif
        }
    }
}

#if defined(BOOST_MSVC)
#pragma warning(pop)
#endif

#endif
