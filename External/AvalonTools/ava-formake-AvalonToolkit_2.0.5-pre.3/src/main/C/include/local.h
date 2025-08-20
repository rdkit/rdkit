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
/*                                                                    */
/*    File:           local.h                                         */
/*                                                                    */
/*    Purpose:        Defines the local constants used to tailor the  */
/*                    code to some machine/compiler.                  */
/*                                                                    */
/************************************************************************/

#ifndef _LOCAL_H_
#define _LOCAL_H_  1

#include <stdlib.h>

// FORTIFY #define FORTIFY
//
// FORTIFY #ifdef FORTIFY
// FORTIFY #ifndef NOFORTIFY
// FORTIFY #include "fortify.h"
// FORTIFY #endif
// FORTIFY /*
// FORTIFY * Returns the currently allocated number of bytes as known to FORTIFY tool.
// FORTIFY */
// FORTIFY extern long fortifySet();
// FORTIFY extern void CheckFortifyMemory();
// FORTIFY extern void logFortifyMessage(const char *msg);
//
// FORTIFY /*
// FORTIFY  * Checks if FORTIFY managed memory has been released to defined level
// FORTIFY  * and exits with message if not.
// FORTIFY  */
// FORTIFY int fortifyTest(long bytes_allocated, char* message);
// FORTIFY #endif

#ifndef FILE
#include <stdio.h>
#endif

#ifndef TRUE
#define TRUE        1
#endif

#ifndef FALSE
#define FALSE       0
#endif

#ifdef __TURBOC__
#define __PROTOTYPES__
extern void sleep(unsigned seconds);
#endif

#ifdef __STDC__
#define __PROTOTYPES__
#endif

#ifdef VAX

#define file_open(name, mode) fopen(name, mode,"rat = cr","rfm = var")
#define __PROTOTYPES__

#else

#define file_open(name, mode) fopen(name, mode)

#endif

#include "mymacros.h"

#define MAXMSG	255

struct msg_line_t
   {
      char buffer[MAXMSG+1];
      struct msg_line_t *next;
   };

typedef struct STRING_INT_TABLE
   {
      char *ident;
      int   value;
   } string_int_table;

char *ToLower(char *str);
/*
 * Converts the uppercase characters in str to lowercase and returns
 * str as the result.
 */

char *MyCalloc(unsigned int n, size_t s);
/*
 * Allocates new memory for n objects of size s, and checks if the 
 * attempt was successful.
 */

char *MyRealloc(char *ptr, unsigned int nold, unsigned int nnew, size_t s);
/*
 * Allocates the memory pointed to by ptr to hold nnew objects of size s.
 * Uses nold to clear any newly allocated memory to '\0' and returns the possibly moved pointer
 * to the new memory block.
 */

void ShowMessage(char *msg, char *title);
/*
 * Used to print a message with no parameter.
 */

void ShowMessageI(char *msg, char *title, int i);
/*
 * Used to print a message with one integer parameter.
 */

void ShowMessageS(char *msg, char *title, char *s);
/*
 * Used to print a message with one string parameter.
 */

FILE *waitopen(char *filename, char *mode,
               int interval, int ntry,
               char *retry_message, char *abort_message);
/*
 * Tries to open the file filename[] with mode mode[]. If not successful and
 * errno == EVMSERR (presumably 'locked by another user'), the functions waits
 * for interval seconds and retries. This cycle is performed ntry times.
 * The file pointer to the opened file is returned or NULL if the file could
 * not be successfully opened.
 * The retry_message[] is shown for each retry and the abort_message[] is shown
 * if all retries failed.
 */

extern FILE * CheckedFileOpen(char *name, char *mode);
/*
 * Tries to open file *name with mode *mode. Exits program
 * and writes message if open fails.
 */

extern int StringToInt(string_int_table *table, char *key);
/*
 * Returns the integer value assuciated with the string *key
 * in *table. Returns 0 if string is not found in *table.
 */

extern char *IntToString(string_int_table *table, int key);
/*
 * Returns a pointer to the string associated with the integer
 * key in *table. It returns NULL if key was not found in *table.
 */

extern void BlankToZero(char *string);
/*
 * Changes ' ' characters in *string to '0' to ease
 * FORTRAN like reading of numbers.
 */

extern void Squeeze(char *cp);
/*
 * Squezes whitespace out of the string pointed to by cp.
 */

extern char *FileToCommandName(char *Command, char *File);
/*
 * Extracts the command name *Command from the full path file name
 * *File. It returns the address of Command.
 */

extern void SearchChar(FILE *fp, int c);
/*
 * Searches through *fp until character c is found.
 */

extern void MyFree(char *cp);
/*
 * Returns storage pointed to by cp back to free store.
 */

extern string_int_table month_names[];

extern FILE *log_file;

extern char msg_buffer[MAXMSG+1];       /* space saver for message handling */

extern int MsgsPending(void);
/*
 * Returns TRUE if a message is pending on the message queue.
 */

extern void FlushMsgs(FILE *fp);
/*
 * Writes the collected messages to the file *fp and frees the storage.
 * If fp == NULL, it just frees the storage.
 */

extern void PrintMsgs(FILE *fp);
/*
 * Writes the collected messages to the file *fp.
 */

extern void AddMsgToList(char buffer[]);
/*
 * Prepends the message stored in buffer[] to *msg_list.
 */

extern double Angle(double x1, double y1, double x2, double y2);
/*
 * Returns the angle between the two vectors (x1,y1) and (x2,y2) at
 * (0,0).
 */

void FreeMsgList(void);
/*
 * free memory which was allocated for messagelist
 */

char *GetMsgList(void);
/*
 * copy all messages to a char *
 * if no messages are available return NULL
 * non NULL return string values need to be free()
 */


/* This is the name to be used by tempnam() to create temporary files.   */
/* However, it is initialized to NULL which means it is ignored. If some */
/* program, however, wants to set it, it can do so.			 */
extern char *tmp_dir_name;

/*
 * Performs an in-place replacement of from to to[...].
 * The caller needs to make sure that str has enough space for the result.
 * The function returns the result for convenience.
 */
extern char *ReplaceOnce(char *str, char *from, char *to);

#endif
