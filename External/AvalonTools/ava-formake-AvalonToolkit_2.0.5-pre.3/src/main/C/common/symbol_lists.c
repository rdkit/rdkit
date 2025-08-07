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
/*    File:           symbol_lists.c                                    */
/*                                                                      */
/*    Author:         B. Rohde                                          */
/*                                                                      */
/*    Purpose:        Implements functions for handling lists of atom   */
/*                    symbol lists.                                     */
/*                                                                      */
/*    History:        12-May-1991     Creation                          */
/*                                                                      */
/************************************************************************/

#include "symbol_lists.h"

#include <string.h>

#include "forio.h"
#include "utilities.h"

void FreeSymbolLists(struct symbol_list_t *symbol_lists)
/*
 * Frees the storage occupied by the symbol lists *symbol_lists.
 */
{
   struct symbol_list_t *slp;

   while (!IsNULL(symbol_lists))
   {
      slp = symbol_lists;
      symbol_lists = slp->next;
      MyFree((char *)slp);
   }
}

struct symbol_list_t *CopySymbolLists(struct symbol_list_t *symbol_lists)
/*
 * Copies the symbol lists *symbol_lists and returns a pointer to the
 * result.
 */
{
   struct symbol_list_t *slp;
   struct symbol_list_t *result;

   result = (struct symbol_list_t *)NULL;
   while (!IsNULL(symbol_lists))
   {
      slp = TypeAlloc(1,struct symbol_list_t);
      *slp = *symbol_lists; slp->next = result;
      result = slp;
      symbol_lists = symbol_lists->next;
   }
   return (result);
}

struct symbol_list_t *ReadSymbolLists(Fortran_FILE *fp, int nlists)
/*
 * Reads nlists symbol list descriptions off the file *fp and returns
 * a pointer to the linked descriptions.
 */
{
   struct symbol_list_t *slp, *result;
   char type_string[4];
   int index, n_elements, type;
   char buffer[MAX_BUFFER+1];
   int i, j;

   result = (struct symbol_list_t *)NULL;
   for (i=0; i<nlists; i++)
   {
      sscanf(fp->buffer,"%d%s%d",&index,type_string,&n_elements);
      for (j=0; j<n_elements; j++)
      {
         if (j==0) strcpy(buffer,"");
         else      strncat(buffer,",",MAX_BUFFER);
         sscanf(&fp->buffer[11+4*j],"%d",&type);
         strncat(buffer,IntToString(periodic_table,type),MAX_BUFFER);
      }
      slp = TypeAlloc(1,struct symbol_list_t);
      slp->next = result; result = slp;
      slp->atom = index;
      strcpy(slp->string, buffer);
      if (0 == strcmp(type_string,"F")) slp->logic = INCLUSIVE;
      else                         slp->logic = EXCLUSIVE;
      GetBuffer(fp);
   }
   return (result);
}

void PrintSymbolLists(FILE *fp, struct symbol_list_t *symbol_lists)
/*
 * Prints the symbol lists pointed to by symbol_lists to the file *fp.
 */
{
   struct symbol_list_t *slp;
   char buffer[MAXSYMBOL+1];
   char *token;
   int n;

   for (slp=symbol_lists; !IsNULL(slp);  slp=slp->next)
   {
      fprintf(fp,"%3d",slp->atom);
      if (slp->logic == INCLUSIVE) fprintf(fp," F ");
      else                         fprintf(fp," T ");
      n = 0;                            /* count symbols */
      strcpy(buffer,slp->string);
      for (token = strtok(buffer," ,");
           !IsNULL(token);
           token = strtok((char *)NULL," ,"))
         n++;
      fprintf(fp," %3d",n);
      strcpy(buffer,slp->string);       /* write symbols */
      for (token = strtok(buffer," ,");
           !IsNULL(token);
           token = strtok((char *)NULL," ,"))
         fprintf(fp," %3d",StringToInt(periodic_table,token));
      fprintf(fp,"\n");
   }
}
