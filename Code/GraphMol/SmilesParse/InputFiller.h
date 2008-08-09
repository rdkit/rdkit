//
//  Copyright (C) 2001-2008 Randal Henne, Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_INPUT_FILLER_H
#define _RD_INPUT_FILLER_H
#include <cstdio>


typedef int (* INPUT_FUNC_TYPE )( char *, int ) ;

void
setInputCharPtr( char * inpStr ); 
void
charPtrCleanup( ); 

void
setInputFilePtr( FILE * fptr );


#endif
