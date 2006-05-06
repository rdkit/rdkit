//
//  Copyright (C) 2001-2006 Randal Henne, Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_INPUT_FILLER_H
#define _RD_INPUT_FILLER_H
#ifdef WIN32
#include <io.h>
#endif

typedef int (* INPUT_FUNC_TYPE )( char *, int ) ;

void
setInputCharPtr( char * inpStr ); 
void
charPtrCleanup( ); 

void
setInputFilePtr( FILE * fptr );


#endif
