#ifndef H_EXPATPPLIB
#define H_EXPATPPLIB

/**
	\file expatpplib.h
	ensure expat definitions setup correctly before including expatpp.h
	so client projects don't have lots of hidden settings
*/
#ifdef WIN32
	#define XML_STATIC  // we are using a static lib build, not expat.dll
	#define COMPILED_FROM_DSP
	#define _LIB
#endif

#include "expatpp.h"  




#endif  //  H_EXPATPPLIB
