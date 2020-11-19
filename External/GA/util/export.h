//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#pragma once
#ifndef SWIG
#	ifdef _MSC_VER
#		pragma warning(disable : 4251)
#		pragma warning(disable : 4275)
#	endif

#	include <boost/config.hpp>
#endif

#if defined(BOOST_HAS_DECLSPEC) && defined(RDKIT_DYN_LINK) && !defined(SWIG)
#	define  GA_EXPORT __declspec(dllexport)
#endif
#ifndef GA_EXPORT
#	define GA_EXPORT
#endif


