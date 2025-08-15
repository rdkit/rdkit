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
/*
 * File: jni_JNIWinTools.cpp
 *
 * Implements the wrappers for the jni.JNIWinTools class.
 */

#include "jni_JNIWinTools.h"

#include <stdlib.h>
#include <stdio.h>
#include <windows.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "depictutil.h"

FILE *log_file = (FILE *)NULL; /* global file */

/**
 * Puts a MetaFile referenced by hmf to the clipboard with size
 * (xExt,yExt).
 *
 * The HANDLE hmf needs to be system-owned piece of memory allocated by
 * GlobalAlloc(GMEM_SHARE | GMEM_ZEROINIT, size+1) or a straight
 * in-memory metafile handle.
 */
static void SendMetaFileToClipboard(HANDLE hmf, int xExt, int yExt)
{
   GLOBALHANDLE   hGMem;
   LPMETAFILEPICT lpMFP;

   hGMem = GlobalAlloc(GHND, (DWORD) sizeof (METAFILEPICT));
   lpMFP = (LPMETAFILEPICT) GlobalLock(hGMem);

   lpMFP->mm   = MM_ANISOTROPIC;
   lpMFP->xExt = xExt;
   lpMFP->yExt = yExt;
   lpMFP->hMF  = (HMETAFILE__*)hmf;

   GlobalUnlock(hGMem);

   OpenClipboard(NULL);
   SetClipboardData(CF_METAFILEPICT, hGMem);
   CloseClipboard();
}

/* Use #pragma to make sure we get a 2 byte aligned struct (size 22 bytes) */
#pragma pack(push)
#pragma pack(2)
struct placer_t
{
   UINT32    key;
   INT16     hmf; //HANDLE hmf;
   INT16     left;//RECT bbox;
   INT16     top;
   INT16     right;
   INT16     bottom;
   INT16     inch;
   UINT32    reserved;
   INT16     checksum;
};
#pragma pack(pop)

/*
 * Class:     JNIWinTools
 * Method:    putBytesToClipboard
 * Signature: ([BLjava/lang/String;)I
 */
JNIEXPORT jint JNICALL Java_jni_JNIWinTools_putBytesToClipboard
  (JNIEnv *env, jobject thisObj, jbyteArray buffer, jstring format)
{
   jboolean isCopy = JNI_FALSE;
   jbyte *bytes;
   int size;
   int xExt, yExt;

   HANDLE datahnd;
   HANDLE hmf;
   LPSTR databuffer;
   UINT fmt;

   struct placer_t placer;

   const char *fmtstr;
   /* register format if necessary */
   fmtstr = env->GetStringUTFChars(format, &isCopy);
   if (0 == strcmp(fmtstr,"Text"))
      fmt = CF_TEXT;
   else if (0 == strcmp(fmtstr, "Picture"))
      fmt = CF_METAFILEPICT;
   else
      fmt = RegisterClipboardFormat(fmtstr);
   if (isCopy == JNI_TRUE) env->ReleaseStringUTFChars(format, fmtstr);

   /* copy byte array to globally allocated memory block */
   size = env->GetArrayLength(buffer);
   bytes = env->GetByteArrayElements(buffer, &isCopy);
   datahnd = GlobalAlloc(GMEM_SHARE | GMEM_ZEROINIT, size+1);
   databuffer = (char *)GlobalLock(datahnd);
   memcpy(databuffer, bytes, size);
   databuffer[size] = '\0';

   if (fmt == CF_METAFILEPICT)
   {
      memcpy((char *)&placer, databuffer, sizeof(struct placer_t));
      if (placer.inch <= 0) placer.inch = 0;
      xExt = placer.right-placer.left;
      if (xExt < 0) xExt *= -1;
      xExt = (int)(xExt*2540.0/placer.inch);
      yExt = placer.bottom-placer.top;
      if (yExt < 0) yExt *= -1;
      yExt = (int)(yExt*2540.0/placer.inch);
      hmf = SetMetaFileBitsEx(size, (const unsigned char *)(databuffer+sizeof(struct placer_t)));
      SendMetaFileToClipboard(hmf, xExt, yExt);
      GlobalUnlock(datahnd);
      if (isCopy == JNI_TRUE) env->ReleaseByteArrayElements(buffer, bytes, 0);
      return fmt;
   }
   GlobalUnlock(datahnd);
   if (isCopy == JNI_TRUE) env->ReleaseByteArrayElements(buffer, bytes, 0);

   /* Put the data to the clipboard */
   OpenClipboard(NULL);
   EmptyClipboard();
   SetClipboardData(fmt, datahnd);
   CloseClipboard();

   return fmt;
}

static LPSTR data_buffer = NULL;        /* Current data buffer pointer */
static HGLOBAL hbuffer = 0;             /* Current data buffer handle */

static LPSTR GetDataFromClipboard(LPINT sizep, LPSTR format)
/*
 * Return a pointer to the data of the given format or NULL if
 * the format is no available.
 */
{
   HWND hwnd;
   int fmt;
   LPSTR result;
   char buffer[100];
   HGLOBAL hCMem;
   LPSTR lpCMem, cp;

   int data_size, line_size;

   hwnd = GetActiveWindow();
   OpenClipboard(hwnd);
   fmt = 0; result = NULL; (*sizep) = 0;
   for (;;)
   {
      fmt = EnumClipboardFormats(fmt);
      if (fmt == 0) break;
      GetClipboardFormatName(fmt, buffer, 99);
      if (strcmp(buffer, format) == 0)    /* format was found */
      {
         if (data_buffer)
         {
            GlobalUnlock(hbuffer);
            GlobalFree(hbuffer);
            hbuffer = 0; data_buffer = NULL;
         }
         hCMem = GetClipboardData(fmt);
         lpCMem = (char *)GlobalLock(hCMem);    /* make a pointer to handle's memory */
         data_size = (int)GlobalSize(hCMem);
         hbuffer = GlobalAlloc(GMEM_SHARE | GMEM_ZEROINIT, data_size);
         if (hbuffer == 0)
         {
            hbuffer = GlobalAlloc(GMEM_SHARE | GMEM_ZEROINIT, data_size);
         }
         cp = data_buffer = (char *)GlobalLock(hbuffer);
         while (data_size > 0)
         {
            line_size = lpCMem[0];
            lpCMem++; data_size--;
            if (line_size >= 6  &&  0 == strcmp(lpCMem, "M  END"))
            {
               sprintf(cp,"%.*s\n", (int)strlen("M  END"), lpCMem); cp += strlen(cp);
               break;
            }
            sprintf(cp,"%.*s\n",line_size, lpCMem);
            cp += strlen(cp);
            lpCMem += line_size; data_size -= line_size;
         }

         GlobalUnlock(hCMem);
         result = data_buffer; (*sizep) = strlen(data_buffer);
         break;
      }
   }
   CloseClipboard();
   return (result);
}

/**
 * Returns a full path name for a file fname in the directory pointed
 * to by the environment variables TMP or TEMP in that order.
 */
static char *GetTempName(const char *fname)
{
   static char buffer[1000];
   char *varp;
   varp = getenv("TMP");
   if (!varp) varp = getenv("TEMP");
   if (varp)
   {
      strcpy(buffer, varp);
      strcat(buffer, "\\");
   }
   else      strcpy(buffer, "C:\\TEMP\\");
   strcat(buffer, fname);
   return (buffer);
}


static void SaveCTFileFromClipboard(LPSTR fname)
/*
 * Fetches the connection table from the clipboard and places it
 * into file *fname. The file is cleared if there is no MDLCT format
 * available.
 */
{
   HWND hwnd;
   char buffer[100];
   int fmt;
   HGLOBAL hCMem;
   LPSTR lpCMem;
   FILE *fp;

   unsigned int data_size;
   unsigned int line_size;

   fp = fopen(fname,"w");
   if (!fp) return;

   hwnd = GetActiveWindow();
   OpenClipboard(hwnd);
   fmt = 0;
   for (;;)                  /* loop over all available formats */
   {
      fmt = EnumClipboardFormats(fmt);
      if (fmt == 0) break;
      GetClipboardFormatName(fmt, buffer, 99);
      if (strcmp(buffer, "MDLCT") == 0)  /* MDLCT format was found */
      {
         hCMem = GetClipboardData(fmt);
         lpCMem = (char *)GlobalLock(hCMem);    /* make a pointer to handle's memory */
         data_size = (int)GlobalSize(hCMem);

         while (data_size > 1)
         {
            line_size = lpCMem[0];
            lpCMem++; data_size--;
            if (data_size < line_size) break;   /* safeguard */
            strncpy(buffer, lpCMem, line_size); buffer[line_size] = '\0';
            strcat(buffer,"\n");
            fputs(buffer, fp);
            if (data_size > strlen("M  END")+1 &&
                0 == strncmp(lpCMem, "M  END", strlen("M  END"))) break;
            lpCMem += line_size; data_size -= line_size;
         }

         GlobalUnlock(hCMem);
         break;
      }
   }
   CloseClipboard();

   fclose(fp);

   return;
}


/*
 * Class:     JNIWinTools
 * Method:    getSMILESBytesOfClipboardCT
 * Signature: ([B)I
 */
JNIEXPORT
jint JNICALL Java_jni_JNIWinTools_getSMILESBytesOfClipboardCT(JNIEnv *env,
                                                              jobject thisObj,
                                                              jbyteArray buffer)
{
   int size;
   char *smiles;
   jbyte *bytes;
   jboolean isCopy = JNI_FALSE;

   SaveCTFileFromClipboard(GetTempName("tmp.mol"));
   smiles = MOLFileToSMILESString(&size, GetTempName("tmp.mol"));
   if (size <= 0) return (0);


   if (size+1 >= env->GetArrayLength(buffer)) return (-size);

   bytes = env->GetByteArrayElements(buffer, &isCopy);
   strncpy((char *)bytes, smiles, size);
   bytes[size] = '\0';
   if (isCopy == JNI_TRUE) env->ReleaseByteArrayElements(buffer, bytes, 0);
   return (size);
}

/*
 * Class:     JNIWinTools
 * Method:    getSMARTSBytesOfClipboardCT
 * Signature: ([B)I
*/
JNIEXPORT
jint JNICALL Java_jni_JNIWinTools_getSMARTSBytesOfClipboardCT(JNIEnv *env,
                                                              jobject thisObj,
                                                              jbyteArray buffer)
{
   int size;
   jbyte *bytes;
   jboolean isCopy = JNI_FALSE;
   char sbuffer[1000];

   SaveCTFileFromClipboard(GetTempName("tmp.mol"));
   size = 999;
   MOLFileToSMARTSBuffer(&size, sbuffer, GetTempName("tmp.mol"));
   if (size <= 0) return (0);

   if (size+1 >= env->GetArrayLength(buffer)) return (-size);

   bytes = env->GetByteArrayElements(buffer, &isCopy);
   strncpy((char *)bytes, sbuffer, size);
   bytes[size] = '\0';
   if (isCopy == JNI_TRUE) env->ReleaseByteArrayElements(buffer, bytes, 0);
   return (size);
}

/*
 * Class:     JNIWinTools
 * Method:    getMOLFileBytesOfClipboardCT
 * Signature: ([B)I
*/
JNIEXPORT
jint JNICALL Java_jni_JNIWinTools_getMOLFileBytesOfClipboardCT(JNIEnv *env,
                                                             jobject thisObj,
                                                             jbyteArray buffer)
{
   char *ct_string;
   int   size;
   jbyte *bytes;
   jboolean isCopy = JNI_FALSE;

   ct_string = GetDataFromClipboard(&size, (char *)"MDLCT");

   if (!ct_string) return (0);
   if (size <= 0) return (0);
   if (size+1 >= env->GetArrayLength(buffer)) return (-size);

   bytes = env->GetByteArrayElements(buffer, &isCopy);
   strncpy((char *)bytes, ct_string, size);
   bytes[size] = '\0';
   if (isCopy == JNI_TRUE) env->ReleaseByteArrayElements(buffer, bytes, 0);
   return (size);
}

/*
 * Class:     JNIWinTools
 * Method:    emptyClipboard
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_jni_JNIWinTools_emptyClipboard
   (JNIEnv *env, jobject thisObj)
{
   OpenClipboard(NULL);
   EmptyClipboard();
   CloseClipboard();
}

#if 1
static HGLOBAL FileToGlobalMemory(const char * fname, int *fsizep)
/*
 * Reads *fname and returns a handle to a global shared memory block
 * which contains the lines of the file as (line length)+data
 * records. This is the format used by the MDLCT clipboard entry.
 * *fsizep is set to the size of the used part of the allocated memory.
 * The two sizes differ because space has been allocated on the basis
 * of the file size including CR+LF sequences, but only one byte of
 * length information is used in the memory coding of the file.
 */
{
   char buffer[256];      /* line buffer */
   FILE *fp;
   long fsize;             /* size of data block to being allocated */
   int lsize;
   HGLOBAL returnhnd;
   LPSTR molbuffer, cp;
   int data_length;

   fp = fopen(fname,"r");    /* get size of file */
   fseek(fp, 0L, SEEK_END);
   fsize = ftell(fp);
   fseek(fp, 0L, SEEK_SET);

   // allocate one more to catch case where last line does not end with '\n'
   returnhnd = GlobalAlloc(GMEM_SHARE | GMEM_ZEROINIT, fsize+1);
   molbuffer = (char *)GlobalLock(returnhnd);

   data_length = 0;
   cp = molbuffer;
   while (fgets(buffer, 255, fp))
   {
      buffer[255] = '\0';
      lsize = strlen(buffer);
      if (buffer[lsize-1] == '\n')
      {
         lsize--;
         buffer[lsize] = '\0';     /* get rid of '\n' */
      }
      data_length += strlen(buffer)+1;
      (*cp) = lsize; cp++;
      strncpy(cp, buffer, lsize); cp += lsize;
   }

   GlobalUnlock(returnhnd);

   fclose(fp);
   (*fsizep) = data_length;
   (*fsizep) = (int)fsize;

   return (returnhnd);
}
#endif

/*
 * Class:     JNIWinTools
 * Method:    putCTToClipboard
 * Signature: (Ljava/lang/String;)I
 *
 * Reads *fname and posts the contents to the clipboard in
 * MDLCT format.
 */
JNIEXPORT
void JNICALL Java_jni_JNIWinTools_putMOLFileToClipboard(JNIEnv *env,
                                                        jobject thisObj,
					                                  jstring fnameObj)
{
   const char *const_fname;
   jboolean isCopy = JNI_FALSE;

   HWND hwnd;
   int fsize;              // size of data block to be posted
   HANDLE molhnd;
   UINT fmt;

   const_fname = env->GetStringUTFChars(fnameObj, &isCopy);
   molhnd = FileToGlobalMemory(const_fname, &fsize);
   if (isCopy == JNI_TRUE) env->ReleaseStringUTFChars(fnameObj, const_fname);

   fmt = RegisterClipboardFormat("MDLCT");  // post to the clipboard
   hwnd = GetActiveWindow();
   OpenClipboard(hwnd);
   EmptyClipboard();
   SetClipboardData(fmt, molhnd);
   CloseClipboard();

   return;
}

#ifdef __cplusplus
}
#endif
