//
//  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
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
//       products derived from this software without specific prior written permission.
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

/************************************************************************/
/*                                                                      */
/*    File:           forio.h                                           */
/*                                                                      */
/*    Purpose:        This file provide structure definitions for the   */
/*                    C-simulation of line buffered Fortran-I/O.        */
/*                                                                      */
/************************************************************************/

#ifndef _FORIO_
#define _FORIO_      1

#include <stdio.h>

#define MAX_BUFFER      4001

#define FORTRAN_NORMAL       0
#define FORTRAN_EOF        EOF
#define FORTRAN_ERROR   (EOF-1)

typedef
struct fortran_file
   {
      char  buffer[MAX_BUFFER];
      FILE *filep;
      int   status;
      int   in_use;
      int   line_no;
      char *source_string;
      char *string_pointer;
   }Fortran_FILE;

extern
void RemoveTrailingBlanks(char *string);

extern
void GetBuffer(Fortran_FILE *fp);

extern
Fortran_FILE *FortranOpen(char *name, char *mode);

extern
Fortran_FILE *FortranStringOpen(char *string_buffer);

extern
void FortranClose(Fortran_FILE *fp);

extern
int SearchString(Fortran_FILE *fp, const char *string, const char *stop_string);
/*
 * Searches the Fortran_FILE *fp for a line beginning with string[].
 * Returns TRUE if the string was found and FALSE if either stop_string
 * or FORTRAN_EOF was encountered.
 */

extern
FILE * RedirFopen(char * fname, char * mode);
/*
 * works like fopen but if fname = "-" it will return stdin for mode "r" and "R"
 * and stdout for modes wWaA
 */

extern
int RedirFclose(FILE * file);
/*
 * works like fclose but will not close stdin and stdout
 */

#endif

#define SkipTagLine(fp,tag) \
   if (!STRING_BEGINS(fp->buffer,tag)) \
   { \
      fprintf(stderr,"format error looking for '%s' on line\n'%s'\n",\
              tag,fp->buffer); \
      exit (EXIT_FAILURE); \
   } \
   else GetBuffer(fp)
