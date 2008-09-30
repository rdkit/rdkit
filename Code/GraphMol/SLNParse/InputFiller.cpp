// $Id: InputFiller.cpp 1066 2008-08-12 06:28:17Z landrgr1 $
//
//
//  Copyright (c) 2008, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met: 
//
//     * Redistributions of source code must retain the above copyright 
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following 
//       disclaimer in the documentation and/or other materials provided 
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
//       nor the names of its contributors may be used to endorse or promote 
//       products derived from this software without specific prior
//       written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  Created by Greg Landrum, September 2006
//
#include <GraphMol/SLNParse/InputFiller.h>
#include <cstring>
#ifdef WIN32
#include <io.h>
#endif


  int filePtrFunc( char * , int );

  INPUT_FUNC_TYPE gp_myInput = filePtrFunc;

  static char 
  * gp_buff=0,
    * gp_cur=0,
    * gp_end=0;

  static FILE 
  * gp_inputFilePtr = stdin;


  int
  cStringFunc( char * b, int mx )
  {
    int
      ii = gp_end - gp_cur,
      nn = ii<mx ? ii: mx;
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

