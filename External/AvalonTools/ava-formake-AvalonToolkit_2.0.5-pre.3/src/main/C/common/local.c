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
/*    File:           local.c                                           */
/*                                                                      */
/*    Purpose:        Local utilities                                   */
/*                                                                      */
/************************************************************************/

#define SHOW_MESSAGES 0

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <ctype.h>

#include "local.h"

#define NO_WINDOWS

#ifdef _WINDOWS
#include <windows.h>

#if SHOW_MESSAGES

void ShowMessage(char *msg, char *t)
/*
 * Used to print a message with no parameter.
 */
{
   MessageBox(0, msg, t, MB_OK);
}

void ShowMessageI(char *msg, char *t, int i)
/*
 * Used to print a message with one integer parameter.
 */
{
   char buffer[100];

   wsprintf(buffer, msg, i);

        MessageBox(0, buffer, t, MB_OK);
}

void ShowMessageS(char *msg, char *t, char *s)
/*
 * Used to print a message with one string parameter.
 */
{
        char buffer[100];

        wsprintf(buffer, msg, s);

        MessageBox(0, buffer, t, MB_OK);
}

#else
void ShowMessage(char *msg, char *t)
/*
 * Used to print a message with no parameter.
 */
{
   /* NOP */
}

void ShowMessageI(char *msg, char *t, int i)
/*
 * Used to print a message with one integer parameter.
 */
{
   /* NOP */
}

void ShowMessageS(char *msg, char *t, char *s)
/*
 * Used to print a message with one string parameter.
 */
{
   /* NOP */
}
#endif
#else
void ShowMessage(char *msg, char *t)
/*
 * Used to print a message with no parameter.
 */
{
   fprintf(stderr, "%s: ", t);
   fprintf(stderr, msg);
   fprintf(stderr, "\n");
}

void ShowMessageI(char *msg, char *t, int i)
/*
 * Used to print a message with one integer parameter.
 */
{
   fprintf(stderr, "%s: ", t);
   fprintf(stderr, msg, i);
   fprintf(stderr, "\n");
}

void ShowMessageS(char *msg, char *t, char *s)
/*
 * Used to print a message with one string parameter.
 */
{
   fprintf(stderr, "%s: ", t);
   fprintf(stderr, msg, s);
   fprintf(stderr, "\n");
}
#endif

char *ToLower(char *str)
/*
 * Converts the uppercase characters in str to lowercase and returns
 * str as the result.
 */
{
   char *cp;

   for (cp=str; (*cp) != '\0'; cp++)
      (*cp) = (char)tolower(*cp);
   return (str);
}

char *MyCalloc(unsigned int n, size_t s)
/*
 * Allocates new memory for n objects of size s, and checks if the 
 * attempt was successful.
 */
{
   char *result;
   static long n_alloc = 0;
   static long total_allocation = 0;
   
#ifdef __HEAPCHECK_ENABLED__
   if (heapcheck() != _HEAPOK)
   {
      ShowMessageI("Heap corrupted. heapcheck() returned %d", "MyCalloc", heapcheck());
#ifndef _WINDOWS
      fprintf(stderr,
              "MyCalloc: heap corrupted! heapcheck returned %d\n",
              heapcheck());
      fprintf(stderr, "%ld bytes left in far core\n", farcoreleft());
      fprintf(stderr, "MyCalloc: Stack Trace:\n");
#endif
      exit (EXIT_FAILURE);
   }
#endif
   if (n == 0) return (NULL);
   result = (char *)calloc(n, s);
   if (!IsNULL(result))
   {
      n_alloc++;
      total_allocation += n*s;
      return (result);
   }

#ifndef _WINDOWS
   fprintf(stderr, "Could not allocate %d objects of size %d bytes\n", n, s);
   fprintf(stderr,
          "%ld bytes could be successfully allocated before in %ld calls\n",
           total_allocation, n_alloc);
#else
   MessageBox(0, "MyCalloc()", "Could not allocate more memory!", 0);
#endif

#ifdef __TURBOC__
#ifndef _Windows
#ifndef __WIN32__
        fprintf(stderr, "%ld bytes left in far core\n", farcoreleft());
#endif
#endif
#endif
   abort();
   return (NULL);
}

char *MyRealloc(char *ptr, unsigned int nold, unsigned int nnew, size_t s)
/*
 * Allocates the memory pointed to by ptr to hold nnew objects of size s.
 * Uses nold to clear any newly allocated memory to '\0' and returns the possibly moved pointer
 * to the new memory block.
 */
{
   char *result;
   unsigned int i;
   static long n_alloc = 0;
   static long total_allocation = 0;
   
   if (ptr == NULL) return MyCalloc(nnew, s);

   result = (char *)realloc((void *)ptr, nnew*s);
   if (!IsNULL(result))
   {
       for (i=nold*s; i<nnew*s; i++) result[i] = '\0';
      n_alloc++;
      total_allocation += (nnew-nold)*s;
      return (result);
   }

#ifndef _WINDOWS
   fprintf(stderr, "Could not allocate %d objects of size %d bytes\n", nnew-nold, s);
   fprintf(stderr,
          "%ld bytes could be successfully allocated before in %ld calls\n",
           total_allocation, n_alloc);
#else
   MessageBox(0, "MyCalloc()", "Could not allocate more memory!", 0);
#endif
   
   abort();
   return (NULL);
}

// FORTIFY void CheckFortifyMemory()
// FORTIFY {
// FORTIFY     int nbroken;
// FORTIFY 
// FORTIFY     nbroken = Fortify_CheckAllMemory();
// FORTIFY     if (nbroken > 0)
// FORTIFY     {
// FORTIFY         fprintf(stderr, "fortify memory check failed on %d blocks\n", nbroken);
// FORTIFY         exit (EXIT_FAILURE);
// FORTIFY     }
// FORTIFY }
//
// FORTIFY void logFortifyMessage(const char *msg)
// FORTIFY {
// FORTIFY     fprintf(stderr, "%s", msg);
// FORTIFY }

void MyFree(char *cp)
/*
 * Returns storage pointed to by cp back to free store.
 */
{
#ifdef __HEAPCHECK_ENABLED__
   int ret;

   if ((ret = heapcheck()) != _HEAPOK)
   {
      ShowMessageI("MyFree(1) Heap corrupted! heapcheck() returned %d",
                   "MyFree", ret);
#ifndef _WINDOWS
      fprintf(stderr, "cp = %lX\n", cp);
      fprintf(stderr, "Stack Trace:\n");
#endif
      exit (EXIT_FAILURE);
   }
#endif

   free(cp);
// FORTIFY CheckFortifyMemory();
   return;
}

FILE * CheckedFileOpen(char *name, char *mode)
/*
 * Tries to open file *name with mode *mode. Exits program
 * and writes message if open fails.
 */
{
   FILE *result;

   if (IsNULL(name))
   {
      ShowMessage("<NULL> file name\n", "CheckedFileOpen");
      exit(2);
   }
#ifdef VAX
   result = fopen(name, mode, "rat=cr", "rfm=var");
#else
   result = fopen(name, mode);
#endif

   if (IsNULL(result))
   {
      ShowMessageS("Could not open with mode=%s\n", name, mode);
      exit(2);
   }
   return(result);
}

void BlankToZero(char *string)
/*
 * Changes ' ' characters in *string to '0' to ease
 * FORTRAN like reading of numbers.
 */
{
   for (; *string != '\0'; string++)
      if (*string == ' ') *string = '0';
}

void Squeeze(char *cp)
/*
 * Squezes whitespace out of the string pointed to by cp.
 */
{
   char *cph;

   for (cph=cp; *cp; cp++)
      if (!isspace(*cp))
      {
         *cph = *cp; cph++;
      }
   *cph = *cp;   /* add '\0' */
}

int StringToInt(string_int_table *table, char *key)
/*
 * Returns the integer value assuciated with the string *key
 * in *table. Returns 0 if string is not found in *table.
 */
{
   for ( ; !IsNULL(table->ident); table++)
      if (0 == strcmp(table->ident, key))
       return(table->value);

   return (0);
}

char *IntToString(string_int_table *table, int key)
/*
 * Returns a pointer to the string associated with the integer
 * key in *table. It returns NULL if key was not found in *table.
 */
{
   for ( ; !IsNULL(table->ident); table++)
      if (table->value == key) return(table->ident);

   return ((char *)NULL);
}

string_int_table month_names[] =
   {
      {"???",  0},
      {"JAN",  1},
      {"FEB",  2},
      {"MAR",  3},
      {"APR",  4},
      {"MAI",  5},
      {"JUN",  6},
      {"JUL",  7},
      {"AUG",  8},
      {"SEP",  9},
      {"OCT", 10},
      {"NOV", 11},
      {"DEC", 12},
      {(char *)NULL, 0}
   };

char *FileToCommandName(char *Command, char *File)
/*
 * Extracts the command name *Command from the full path file name
 * *File. It returns the address of Command.
 */
{
   if (!IsNULL(strrchr(File,']')))
      strcpy(Command, strrchr(File,']')+1);
   else if (!IsNULL(strrchr(File,'/')))
      strcpy(Command, strrchr(File,'/')+1);
   else
      strcpy(Command, File);

   if (!IsNULL(strchr(Command,'.'))) strchr(Command,'.')[0] = '\0';

   return(Command);
}

void SearchChar(FILE *fp, int tag)
/*
 * Searches through *fp until character tag is found.
 */
{
   int c;
   
   while (tag != (c=fgetc(fp)) && c != EOF)
      if (c == EOF)
      {
         ShowMessageI("Unexpected end of file while searching '%c'\n",
                      "SearchChar",
                   tag);
         exit (EXIT_FAILURE);
      }
}

FILE *waitopen(char *filename, char *mode,
               int interval, int ntry,
               char *retry_message, char *abort_message)
/*
 * Tries to open the file filename[] with mode mode[]. If not successful and
 * errno == EVMSERR (presumably 'locked by another user'), the functions waits
 * for interval seconds and retries. This cycle is performed ntry times.
 * The file pointer to the opened file is returned or NULL if the file could
 * not be successfully opened.
 * The retry_message[] is shown for each retry and the abort_message[] is shown
 * if all retries failed.
 */
{
   FILE *retval;
   int itry;
   
   itry = 0;
   for (;;)
   {
      errno = 0;
      retval=fopen(filename, mode);
      if (retval                                 /* successfully opened OR */
#ifdef VAX
          || errno != EVMSERR                     /* cannot be locked       */
#endif
       )
         return (retval);
      if (ntry-- > 0)
      {
         itry++;
#ifndef _WINDOWS
#ifndef __WIN32__
         if (retry_message) fprintf(stderr, "\r%d) %s \r", itry, retry_message);
#endif
#endif
        	}
      else break;
   }
#ifndef _WINDOWS
   if (abort_message)
      fprintf(stderr, "file='%s' mode='%s': %s\n", filename, mode, abort_message);
#endif
   return (retval);
}

struct msg_line_t *msg_list = NULL;

char msg_buffer[MAXMSG+1];   /* space saver for message handling */

int MsgsPending()
/*
 * Returns TRUE if a message is pending on the message queue.
 */
{
   return (msg_list != NULL);
}

void FlushMsgs(FILE *fp)
/*
 * Writes the collected messages to the file *fp and frees the storage.
 * If fp == NULL, it just frees the storage.
 */
{
   struct msg_line_t *hp;

   while (msg_list)
   {
      if (fp) fprintf(fp, "%s\n", msg_list->buffer);
      hp = msg_list->next; free(msg_list); msg_list = hp;
   }
}

void PrintMsgs(FILE *fp)
/*
 * Writes the collected messages to the file *fp.
 */
{
   struct msg_line_t *hp, *hp1;

   hp = (struct msg_line_t *)NULL; /* invert order of messages */
   while (msg_list)
   {
      hp1 = msg_list->next; msg_list->next = hp;
      hp = msg_list; msg_list = hp1;
   }
   msg_list = hp;

   for (hp = msg_list; hp; hp = hp->next)
      if (fp) fprintf(fp, "%s\n", hp->buffer);

   hp = (struct msg_line_t *)NULL;   /* restore order of messages */
   while (msg_list)
   {
      hp1 = msg_list->next; msg_list->next = hp;
      hp = msg_list; msg_list = hp1;
   }
   msg_list = hp;
}

void AddMsgToList(char buffer[])
/*
 * Prepends the message stored in buffer[] to *msg_list.
 */
{
   struct msg_line_t *hp;

   hp = TypeAlloc(1, struct msg_line_t);
   strncpy(hp->buffer, buffer, MAXMSG); hp->buffer[MAXMSG] = '\0';
   hp->next = msg_list;
   msg_list = hp;
}

#define PI  3.14159265359

double Angle(double x1, double y1, double x2, double y2)
/*
 * Returns the angle between the two vectors (x1,y1) and (x2,y2) at
 * (0,0).
 */
{
   double l1, l2;
   double cos_alpha, sin_alpha;
   double result;

   l1 = sqrt(x1*x1+y1*y1); l2 = sqrt(x2*x2+y2*y2);
   if (l1 < 0.00001  ||  l2 < 0.00001) return (0.0);

   cos_alpha = (x1*x2 + y1*y2) / (l1*l2);
   if (cos_alpha > 1.0)          /* safeguard against round off erros */
      cos_alpha = 1.0;
   else if (cos_alpha < -1.0)
      cos_alpha = -1.0;
   sin_alpha = (x1*y2 - x2*y1) / (l1*l2);

   result = acos(cos_alpha);
   if (sin_alpha < 0.0) result = 2*PI-result;
   return (result);
}

void FreeMsgList(void)
/*
 * free memory which was allocated for messagelist
 */
{
   struct msg_line_t *hp;

   while (msg_list)
   {
      hp = msg_list->next; free(msg_list); msg_list = hp;
   }
}

char * GetMsgList(void)
/*
 * copy all messages to a char *
 * if no messages are available return NULL
 * non NULL return string values need to be free()
 */
{   struct msg_line_t *hp, *hp1;
    char * buf;
    int    Size;

    if ((buf = (char *)malloc(1)) == NULL)
    {
       fprintf(stderr, "Error allocating memory)\n");
       abort();
    }
    *buf = '\0';
    Size = 0;

    hp = (struct msg_line_t *)NULL; /* invert order of messages */
    while (msg_list)
    {
       hp1 = msg_list->next; msg_list->next = hp;
       hp = msg_list; msg_list = hp1;
    }
    msg_list = hp;

    for (hp = msg_list; hp; hp = hp->next)
    {
       Size += strlen(hp->buffer) + 1;
       if( (buf = (char *)realloc(buf, Size + 1)) == NULL )
       {
          fprintf(stderr, "Error allocating memory(%d)\n", Size + 1);
          abort();
       }
       strcat( buf, hp->buffer );
       buf[Size-1]   = '\n';
       buf[Size] = '\0';
    }

    hp = (struct msg_line_t *)NULL;   /* restore order of messages */
    while (msg_list)
    {
       hp1 = msg_list->next; msg_list->next = hp;
       hp = msg_list; msg_list = hp1;
    }
    msg_list = hp;

    return buf;
}

char *ReplaceOnce(char *str, char *from, char *to)
/*
 * Performs an in-place replacement of from to to[...].
 * The caller needs to make sure that str has enough space for the result.
 * The function returns the result for convenience.
 */
{
   char *cp;
   char buffer[5000];

   cp = strstr(str, from);
   if (!cp) return (str);
   if (strlen(from) >= strlen(to))
   {
      strncpy(cp, to, strlen(to));
      strcpy(cp+strlen(to), cp+strlen(from));
      cp+=strlen(to);
   }
   else
   {
      strncpy(buffer, str, cp-str);
      strcpy(buffer+(cp-str), to);
      strcat(buffer, cp+strlen(from));
      strcpy(str, buffer);
   }
   return (str);
}

char *tmp_dir_name = (char*) NULL;

// FORTIFY long fortifySet()
// FORTIFY {
// FORTIFY     long result;
// FORTIFY     result = 0;
// FORTIFY #ifdef FORTIFY
// FORTIFY    result = Fortify_GetCurrentAllocation();
// FORTIFY #endif
// FORTIFY    return (result);
// FORTIFY }

// FORTIFY  fortifyTest(long bytes_allocated, char* message)
// FORTIFY {
// FORTIFY #ifdef FORTIFY
// FORTIFY    if (bytes_allocated != Fortify_GetCurrentAllocation())
// FORTIFY    {
// FORTIFY        Fortify_ListAllMemory();
// FORTIFY        fprintf(stderr, "fortifyTest(%s): allocation mismatch of size '%ld'\n",
// FORTIFY                message,
// FORTIFY                Fortify_GetCurrentAllocation()-bytes_allocated);
// FORTIFY        exit(1);
// FORTIFY    }
// FORTIFY #endif
// FORTIFY }

