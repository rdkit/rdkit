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
/*    File:           forio.c                                           */
/*                                                                      */
/*    Purpose:        This file implements the functions used to        */
/*                    implements input/output similar to Fortran in     */
/*                    a line-buffered fashion.                          */
//                                                                      */
//    History:        22-Dec-1992 Added ANSI prototypes.            */
//                    14-2-97     (AG) Added RedirFopen and RedirFclose */
/*                                                                      */
/************************************************************************/

#include "forio.h"

#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "local.h"

#define NFILES 5

/* used to save state of open Fortran_FILE handles */
Fortran_FILE  fortran_files[NFILES] =
{
   {{0}, (FILE *)NULL, FORTRAN_EOF, FALSE},
   {{0}, (FILE *)NULL, FORTRAN_EOF, FALSE},
   {{0}, (FILE *)NULL, FORTRAN_EOF, FALSE},
   {{0}, (FILE *)NULL, FORTRAN_EOF, FALSE},
   {{0}, (FILE *)NULL, FORTRAN_EOF, FALSE}
};

void RemoveTrailingBlanks(char *string)
/*
 * Remove trailing spaces from string to save space when e.g. printing.
 */
{
   register char *cp;

   for (cp=string; *cp!='\0'; cp++)
      if (!isspace(*cp)) string = cp+1;

   *string = '\0';
}

void GetBuffer(Fortran_FILE *fp)
{
   char *cp, *cph;

   if (fp->source_string)
   {
      if (fp->string_pointer[0] == '\0')
      {
	 fp->buffer[0] = '\0';
	 fp->status = FORTRAN_EOF;
      }
      else
      {
         for (cp=fp->string_pointer, cph=fp->buffer;
              (*cp) && (*cp) != '\n' && (*cp) != '\r' &&
                 cp < fp->string_pointer + MAX_BUFFER;
              cp++)
         {
            (*cph) = (*cp); cph++;
         }
         if ((*cp) == '\r') cp++;
         if ((*cp) == '\n') cp++;
         fp->string_pointer = cp;
         (*cph) = '\0';
      }

      RemoveTrailingBlanks(fp->buffer);

      fp->line_no++;
      return;
   }
   if (feof(fp->filep))
   {
      fp->status = FORTRAN_EOF;
      fp->buffer[0] = '\0';
      return;
   }
   cp = fgets(fp->buffer,MAX_BUFFER,fp->filep);
   if (cp == (char *)NULL)
   {
      if (feof(fp->filep)) fp->status = FORTRAN_EOF;
      else                fp->status = FORTRAN_ERROR;
   }

   RemoveTrailingBlanks(fp->buffer);
   fp->line_no++;
}

Fortran_FILE *FortranOpen(char *name, char *mode)
{
   int i;
   Fortran_FILE *fp;

   for (i=0, fp=fortran_files; i<NFILES; i++, fp++)
      if (!fp->in_use) break;
   if (i == NFILES)                                /* no more file slots */
      return ((Fortran_FILE *)NULL);

   fp->filep = RedirFopen(name,mode);
   if (fp->filep == (FILE *)NULL) /* file open failed */
      return((Fortran_FILE *)NULL);

   fp->status    = FORTRAN_NORMAL;
   fp->in_use    = TRUE;
   fp->line_no   = 0;
   fp->buffer[0] = '\0';

   fp->source_string  = (char *)NULL;
   fp->string_pointer = (char *)NULL;

   GetBuffer(fp);

   return (fp);
}

Fortran_FILE *FortranStringOpen(char *string_buffer)
{
   int i;
   Fortran_FILE *fp;

   for (i=0, fp=fortran_files; i<NFILES; i++, fp++)
      if (!fortran_files[i].in_use) break;
   if (i == NFILES)                             /* no more file slots */
      return ((Fortran_FILE *)NULL);

   fp->filep = (FILE *)NULL;
   fp->source_string  = TypeAlloc(strlen(string_buffer)+1, char);
   fp->string_pointer = fp->source_string;

   if (!fp->source_string) /* not enough memory */
      return((Fortran_FILE *)NULL);
   strcpy(fp->source_string, string_buffer);

   fortran_files[i].status    = FORTRAN_NORMAL;
   fortran_files[i].in_use    = TRUE;
   fortran_files[i].line_no   = 0;
   fortran_files[i].buffer[0] = '\0';

   GetBuffer(&fortran_files[i]);

   return(&fortran_files[i]);
}

void FortranClose(Fortran_FILE *fp)
{
   if (fp->source_string)
   {
      MyFree(fp->source_string);
      fp->source_string  = NULL;
      fp->string_pointer = NULL;
   }
   else
   RedirFclose(fp->filep);
   fp->in_use    = FALSE;
   fp->status    = FORTRAN_NORMAL;
   fp->buffer[0] = '\0';
}

int SearchString(Fortran_FILE *fp,
                 const char *string, const char *stop_string)
/*
 * Searches the Fortran_FILE *fp for a line beginning with string[].
 * Returns TRUE if the string was found and FALSE if either stop_string
 * or FORTRAN_EOF was encountered.
 */
{
   while (fp->status != FORTRAN_EOF &&        /* search for header line */
          !STRING_BEGINS(fp->buffer,string))
      if (STRING_BEGINS(fp->buffer,stop_string))
         return(FALSE);
      else
         GetBuffer(fp);
   if (fp->status == FORTRAN_EOF) return(FALSE);
   else                           return(TRUE);
}

FILE *RedirFopen(char * fname, char * mode)
{
   if (strcmp(fname, "-") == 0)
   {
      if (strcmp(mode, "r") == 0 || strcmp(mode, "R") == 0)
      {
         return stdin;
      }
      else if (strcmp(mode, "w") == 0 || strcmp(mode, "a") == 0 ||
               strcmp(mode, "W") == 0 || strcmp(mode, "A") == 0)
      {
         return stdout;
      }
   }
   return fopen(fname, mode);
}

int RedirFclose(FILE * file)
{
   if (file != stdin && file != stdout)
      return fclose(file);

   return 0;
}

