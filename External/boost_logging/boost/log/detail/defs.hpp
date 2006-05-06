// defs.hpp

// Boost Logging Template library
//
// Author: John Torjo
//
// Copyright (C) 2004-2005 John Torjo (john@torjo.com)
//
// Permission to copy, use, sell and distribute this software is granted
// provided this copyright notice appears in all copies.
// Permission to modify the code and to distribute modified code is granted
// provided this copyright notice appears in all copies, and a notice
// that the code was modified is included with the copyright notice.
//
// This software is provided "as is" without express or implied warranty,
// and with no claim as to its suitability for any purpose.
 
// See http://www.boost.org for updates, documentation, and revision history.

#ifndef JT_BOOST_LOG_DEFS_HPP
#define JT_BOOST_LOG_DEFS_HPP


#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif


#include <string>
#include <boost/config.hpp>
#include <boost/shared_ptr.hpp>

#if defined(_WIN32) || defined(__WIN32__) || defined(WIN32) || defined(_WINDOWS) || defined(WINDOWS)
#define BOOST_LOG_WIN32
#endif



#ifdef BOOST_LOG_WIN32
#if (defined(BOOST_ALL_DYN_LINK) || defined(BOOST_LOG_DYN_LINK)) && !defined(BOOST_LOG_NO_LIB)
// export if this is our own source, otherwise import:
#ifdef BOOST_LOG_SOURCE
# define BOOST_LOG_DECL __declspec(dllexport)
#else
# define BOOST_LOG_DECL __declspec(dllimport)
#endif  // BOOST_LOG_SOURCE
#endif  // DYN_LINK
#endif

// if BOOST_LOG_DECL isn't defined yet define it now:
#ifndef BOOST_LOG_DECL
#define BOOST_LOG_DECL
#endif

// for VC
#if defined(_DLL) || defined(_USRDLL)
#ifndef BOOST_LOG_DYN_LINK
#define BOOST_LOG_DYN_LINK
#endif
#endif


#if defined(BOOST_ALL_DYN_LINK) || defined(BOOST_LOG_DYN_LINK)
#  define BOOST_DYN_LINK
#endif


#if !defined(BOOST_LOG_SOURCE) && !defined(BOOST_ALL_NO_LIB) && !defined(BOOST_LOG_NO_LIB)
//
// Set the name of our library, this will get undef'ed by auto_link.hpp
// once it's done with it:
//
#define BOOST_LIB_NAME boost_log
//
// And include the header that does the work:
//
#include <boost/config/auto_link.hpp>
#endif  // auto-linking disabled



#endif

