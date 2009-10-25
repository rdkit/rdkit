
// Copyright 2005-2009 Daniel James.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Modifications by Greg Landrum to get portable hashes across machines
//
//  Based on Peter Dimov's proposal
//  http://www.open-std.org/JTC1/SC22/WG21/docs/papers/2005/n1756.pdf
//  issue 6.18. 

#if !defined(GBOOST_FUNCTIONAL_HASH_HASH_HPP)
#define GBOOST_FUNCTIONAL_HASH_HASH_HPP

#include <RDGeneral/hash/hash_fwd.hpp>
#include <functional>
#include <RDGeneral/hash/detail/hash_float.hpp>
#include <boost/detail/container_fwd.hpp>
#include <string>

#if defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)
#include <boost/type_traits/is_pointer.hpp>
#endif

#if defined(BOOST_NO_FUNCTION_TEMPLATE_ORDERING)
#include <boost/type_traits/is_array.hpp>
#endif

#if BOOST_WORKAROUND(BOOST_MSVC, < 1300)
#include <boost/type_traits/is_const.hpp>
#endif

namespace gboost
{
    std::hash_result_t hash_value(bool);
    std::hash_result_t hash_value(char);
    std::hash_result_t hash_value(unsigned char);
    std::hash_result_t hash_value(signed char);
    std::hash_result_t hash_value(short);
    std::hash_result_t hash_value(unsigned short);
    std::hash_result_t hash_value(int);
    std::hash_result_t hash_value(unsigned int);
    std::hash_result_t hash_value(long);
    std::hash_result_t hash_value(unsigned long);

#if !defined(BOOST_NO_INTRINSIC_WCHAR_T)
    std::hash_result_t hash_value(wchar_t);
#endif
    
#if defined(BOOST_HAS_LONG_LONG)
    std::hash_result_t hash_value(boost::long_long_type);
    std::hash_result_t hash_value(boost::ulong_long_type);
#endif

#if !BOOST_WORKAROUND(__DMC__, <= 0x848)
    template <class T> std::hash_result_t hash_value(T* const&);
#else
    template <class T> std::hash_result_t hash_value(T*);
#endif

#if !defined(BOOST_NO_FUNCTION_TEMPLATE_ORDERING)
    template< class T, unsigned N >
    std::hash_result_t hash_value(const T (&x)[N]);

    template< class T, unsigned N >
    std::hash_result_t hash_value(T (&x)[N]);
#endif

    std::hash_result_t hash_value(float v);
    std::hash_result_t hash_value(double v);
    std::hash_result_t hash_value(long double v);

    template <class A, class B>
    std::hash_result_t hash_value(std::pair<A, B> const&);
    template <class T, class A>
    std::hash_result_t hash_value(std::vector<T, A> const&);
    template <class T, class A>
    std::hash_result_t hash_value(std::list<T, A> const& v);
    template <class T, class A>
    std::hash_result_t hash_value(std::deque<T, A> const& v);
    template <class K, class C, class A>
    std::hash_result_t hash_value(std::set<K, C, A> const& v);
    template <class K, class C, class A>
    std::hash_result_t hash_value(std::multiset<K, C, A> const& v);
    template <class K, class T, class C, class A>
    std::hash_result_t hash_value(std::map<K, T, C, A> const& v);
    template <class K, class T, class C, class A>
    std::hash_result_t hash_value(std::multimap<K, T, C, A> const& v);

    template <class T>
    std::hash_result_t hash_value(std::complex<T> const&);

    // Implementation

    namespace hash_detail
    {
        template <class T>
        inline std::hash_result_t hash_value_signed(T val)
        {
             const int hash_result_t_bits = std::numeric_limits<std::hash_result_t>::digits;
             // ceiling(std::numeric_limits<T>::digits / hash_result_t_bits) - 1
             const int length = (std::numeric_limits<T>::digits - 1)
                 / hash_result_t_bits;

             std::hash_result_t seed = 0;
             T positive = val < 0 ? -1 - val : val;

             // Hopefully, this loop can be unrolled.
             for(unsigned int i = length * hash_result_t_bits; i > 0; i -= hash_result_t_bits)
             {
                 seed ^= (std::hash_result_t) (positive >> i) + (seed<<6) + (seed>>2);
             }
             seed ^= (std::hash_result_t) val + (seed<<6) + (seed>>2);

             return seed;
        }

        template <class T>
        inline std::hash_result_t hash_value_unsigned(T val)
        {
             const int hash_result_t_bits = std::numeric_limits<std::hash_result_t>::digits;
             // ceiling(std::numeric_limits<T>::digits / hash_result_t_bits) - 1
             const int length = (std::numeric_limits<T>::digits - 1)
                 / hash_result_t_bits;

             std::hash_result_t seed = 0;

             // Hopefully, this loop can be unrolled.
             for(unsigned int i = length * hash_result_t_bits; i > 0; i -= hash_result_t_bits)
             {
                 seed ^= (std::hash_result_t) (val >> i) + (seed<<6) + (seed>>2);
             }
             seed ^= (std::hash_result_t) val + (seed<<6) + (seed>>2);

             return seed;
        }
    }

    inline std::hash_result_t hash_value(bool v)
    {
        return static_cast<std::hash_result_t>(v);
    }

    inline std::hash_result_t hash_value(char v)
    {
        return static_cast<std::hash_result_t>(v);
    }

    inline std::hash_result_t hash_value(unsigned char v)
    {
        return static_cast<std::hash_result_t>(v);
    }

    inline std::hash_result_t hash_value(signed char v)
    {
        return static_cast<std::hash_result_t>(v);
    }

    inline std::hash_result_t hash_value(short v)
    {
        return static_cast<std::hash_result_t>(v);
    }

    inline std::hash_result_t hash_value(unsigned short v)
    {
        return static_cast<std::hash_result_t>(v);
    }

    inline std::hash_result_t hash_value(int v)
    {
        return static_cast<std::hash_result_t>(v);
    }

    inline std::hash_result_t hash_value(unsigned int v)
    {
        return static_cast<std::hash_result_t>(v);
    }

    inline std::hash_result_t hash_value(long v)
    {
        return static_cast<std::hash_result_t>(v);
    }

    inline std::hash_result_t hash_value(unsigned long v)
    {
        return static_cast<std::hash_result_t>(v);
    }

#if !defined(BOOST_NO_INTRINSIC_WCHAR_T)
    inline std::hash_result_t hash_value(wchar_t v)
    {
        return static_cast<std::hash_result_t>(v);
    }
#endif

#if defined(BOOST_HAS_LONG_LONG)
    inline std::hash_result_t hash_value(boost::long_long_type v)
    {
        return hash_detail::hash_value_signed(v);
    }

    inline std::hash_result_t hash_value(boost::ulong_long_type v)
    {
        return hash_detail::hash_value_unsigned(v);
    }
#endif

    // Implementation by Alberto Barbati and Dave Harris.
#if !BOOST_WORKAROUND(__DMC__, <= 0x848)
    template <class T> std::hash_result_t hash_value(T* const& v)
#else
    template <class T> std::hash_result_t hash_value(T* v)
#endif
    {
        std::hash_result_t x = static_cast<std::hash_result_t>(
           reinterpret_cast<std::ptrdiff_t>(v));

        return x + (x >> 3);
    }

#if BOOST_WORKAROUND(BOOST_MSVC, < 1300)
    template <class T>
    inline void hash_combine(std::hash_result_t& seed, T& v)
#else
    template <class T>
    inline void hash_combine(std::hash_result_t& seed, T const& v)
#endif
    {
        gboost::hash<T> hasher;
        seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    }

    template <class It>
    inline std::hash_result_t hash_range(It first, It last)
    {
        std::hash_result_t seed = 0;

        for(; first != last; ++first)
        {
            hash_combine(seed, *first);
        }

        return seed;
    }

    template <class It>
    inline void hash_range(std::hash_result_t& seed, It first, It last)
    {
        for(; first != last; ++first)
        {
            hash_combine(seed, *first);
        }
    }

#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x551))
    template <class T>
    inline std::hash_result_t hash_range(T* first, T* last)
    {
        std::hash_result_t seed = 0;

        for(; first != last; ++first)
        {
            gboost::hash<T> hasher;
            seed ^= hasher(*first) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }

        return seed;
    }

    template <class T>
    inline void hash_range(std::hash_result_t& seed, T* first, T* last)
    {
        for(; first != last; ++first)
        {
            gboost::hash<T> hasher;
            seed ^= hasher(*first) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
    }
#endif

#if !defined(BOOST_NO_FUNCTION_TEMPLATE_ORDERING)
    template< class T, unsigned N >
    inline std::hash_result_t hash_value(const T (&x)[N])
    {
        return hash_range(x, x + N);
    }

    template< class T, unsigned N >
    inline std::hash_result_t hash_value(T (&x)[N])
    {
        return hash_range(x, x + N);
    }
#endif

    inline std::hash_result_t hash_value(float v)
    {
        return gboost::hash_detail::float_hash_value(v);
    }

    inline std::hash_result_t hash_value(double v)
    {
        return gboost::hash_detail::float_hash_value(v);
    }

    inline std::hash_result_t hash_value(long double v)
    {
        return gboost::hash_detail::float_hash_value(v);
    }

    template <class A, class B>
    std::hash_result_t hash_value(std::pair<A, B> const& v)
    {
        std::hash_result_t seed = 0;
        hash_combine(seed, v.first);
        hash_combine(seed, v.second);
        return seed;
    }

    template <class T, class A>
    std::hash_result_t hash_value(std::vector<T, A> const& v)
    {
        return hash_range(v.begin(), v.end());
    }

    template <class T, class A>
    std::hash_result_t hash_value(std::list<T, A> const& v)
    {
        return hash_range(v.begin(), v.end());
    }

    template <class T, class A>
    std::hash_result_t hash_value(std::deque<T, A> const& v)
    {
        return hash_range(v.begin(), v.end());
    }

    template <class K, class C, class A>
    std::hash_result_t hash_value(std::set<K, C, A> const& v)
    {
        return hash_range(v.begin(), v.end());
    }

    template <class K, class C, class A>
    std::hash_result_t hash_value(std::multiset<K, C, A> const& v)
    {
        return hash_range(v.begin(), v.end());
    }

    template <class K, class T, class C, class A>
    std::hash_result_t hash_value(std::map<K, T, C, A> const& v)
    {
        return hash_range(v.begin(), v.end());
    }

    template <class K, class T, class C, class A>
    std::hash_result_t hash_value(std::multimap<K, T, C, A> const& v)
    {
        return hash_range(v.begin(), v.end());
    }

    template <class T>
    std::hash_result_t hash_value(std::complex<T> const& v)
    {
        gboost::hash<T> hasher;
        std::hash_result_t seed = hasher(v.imag());
        seed ^= hasher(v.real()) + (seed<<6) + (seed>>2);
        return seed;
    }

    //
    // gboost::hash
    //

#if !BOOST_WORKAROUND(BOOST_MSVC, < 1300)
#define BOOST_HASH_SPECIALIZE(type) \
    template <> struct hash<type> \
         : public std::unary_function<type, std::hash_result_t> \
    { \
        std::hash_result_t operator()(type v) const \
        { \
            return gboost::hash_value(v); \
        } \
    };

#define BOOST_HASH_SPECIALIZE_REF(type) \
    template <> struct hash<type> \
         : public std::unary_function<type, std::hash_result_t> \
    { \
        std::hash_result_t operator()(type const& v) const \
        { \
            return gboost::hash_value(v); \
        } \
    };
#else
#define BOOST_HASH_SPECIALIZE(type) \
    template <> struct hash<type> \
         : public std::unary_function<type, std::hash_result_t> \
    { \
        std::hash_result_t operator()(type v) const \
        { \
            return gboost::hash_value(v); \
        } \
    }; \
    \
    template <> struct hash<const type> \
         : public std::unary_function<const type, std::hash_result_t> \
    { \
        std::hash_result_t operator()(const type v) const \
        { \
            return gboost::hash_value(v); \
        } \
    };

#define BOOST_HASH_SPECIALIZE_REF(type) \
    template <> struct hash<type> \
         : public std::unary_function<type, std::hash_result_t> \
    { \
        std::hash_result_t operator()(type const& v) const \
        { \
            return gboost::hash_value(v); \
        } \
    }; \
    \
    template <> struct hash<const type> \
         : public std::unary_function<const type, std::hash_result_t> \
    { \
        std::hash_result_t operator()(type const& v) const \
        { \
            return gboost::hash_value(v); \
        } \
    };
#endif

    BOOST_HASH_SPECIALIZE(bool)
    BOOST_HASH_SPECIALIZE(char)
    BOOST_HASH_SPECIALIZE(signed char)
    BOOST_HASH_SPECIALIZE(unsigned char)
#if !defined(BOOST_NO_INTRINSIC_WCHAR_T)
    BOOST_HASH_SPECIALIZE(wchar_t)
#endif
    BOOST_HASH_SPECIALIZE(short)
    BOOST_HASH_SPECIALIZE(unsigned short)
    BOOST_HASH_SPECIALIZE(int)
    BOOST_HASH_SPECIALIZE(unsigned int)
    BOOST_HASH_SPECIALIZE(long)
    BOOST_HASH_SPECIALIZE(unsigned long)

    BOOST_HASH_SPECIALIZE(float)
    BOOST_HASH_SPECIALIZE(double)
    BOOST_HASH_SPECIALIZE(long double)

#undef BOOST_HASH_SPECIALIZE
#undef BOOST_HASH_SPECIALIZE_REF

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)
    template <class T>
    struct hash<T*>
        : public std::unary_function<T*, std::hash_result_t>
    {
        std::hash_result_t operator()(T* v) const
        {
#if !BOOST_WORKAROUND(__SUNPRO_CC, <= 0x590)
            return gboost::hash_value(v);
#else
            std::hash_result_t x = static_cast<std::hash_result_t>(
                reinterpret_cast<std::ptrdiff_t>(v));

            return x + (x >> 3);
#endif
        }
    };
#else
    namespace hash_detail
    {
        template <bool IsPointer>
        struct hash_impl;

        template <>
        struct hash_impl<true>
        {
            template <class T>
            struct inner
                : public std::unary_function<T, std::hash_result_t>
            {
                std::hash_result_t operator()(T val) const
                {
#if !BOOST_WORKAROUND(__SUNPRO_CC, <= 590)
                    return gboost::hash_value(val);
#else
                    std::hash_result_t x = static_cast<std::hash_result_t>(
                        reinterpret_cast<std::ptrdiff_t>(val));

                    return x + (x >> 3);
#endif
                }
            };
        };
    }

    template <class T> struct hash
        : public gboost::hash_detail::hash_impl<boost::is_pointer<T>::value>
            ::BOOST_NESTED_TEMPLATE inner<T>
    {
    };
#endif
}

#endif // BOOST_FUNCTIONAL_HASH_HASH_HPP

// Include this outside of the include guards in case the file is included
// twice - once with BOOST_HASH_NO_EXTENSIONS defined, and then with it
// undefined.

#if !defined(GBOOST_HASH_NO_EXTENSIONS) \
    && !defined(GBOOST_FUNCTIONAL_HASH_EXTENSIONS_HPP)
#include <RDGeneral/hash/extensions.hpp>
#endif
