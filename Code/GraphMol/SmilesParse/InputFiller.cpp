// $Id$
//
//  Copyright (C) 2001-2008 Randal Henne, Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "InputFiller.h"
#include <cstring>
#ifdef WIN32
#include <io.h> 	 
#endif
int filePtrFunc( char * , int );

INPUT_FUNC_TYPE gp_myInput = filePtrFunc;

static char *gp_buff=0,* gp_cur=0,* gp_end=0;

static FILE *gp_inputFilePtr = stdin;


int
cStringFunc( char * b, int mx )
{
  int ii = gp_end - gp_cur;
  int nn = ii<mx ? ii: mx;
  if ( nn > 0 ) {
    memcpy( b, gp_cur, nn );
    gp_cur += nn;
  } else {
    delete [] gp_buff;
    gp_buff = gp_cur = gp_end = 0;
  }
  return nn;
}


void
setInputCharPtr( char * inpStr )
{
  int len = strlen( inpStr );
  if(gp_buff) delete [] gp_buff;
  gp_buff = new char[ len + 1 ];
  strcpy( gp_buff, inpStr );
  gp_cur  = gp_buff;
  gp_end  = gp_buff + len;

  gp_myInput = cStringFunc;
}


int
filePtrFunc( char * b, int mx )
{
  return fread( b, 1, mx, gp_inputFilePtr );
}

void
setInputFilePtr( FILE * fptr )
{
  gp_inputFilePtr = fptr;
  gp_myInput      = filePtrFunc; 
}

void charPtrCleanup(){
  delete [] gp_buff;
  gp_buff = 0;
}
