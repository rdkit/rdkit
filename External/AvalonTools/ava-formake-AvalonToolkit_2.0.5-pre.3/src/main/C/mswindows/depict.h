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
#include <windows.h>

#ifdef __cplusplus
extern "C" {
#endif

/* This set of macros is necessary to make a source file that exports */
/* the DLL entry points in uppercase on both VC++ and Borland C++     */
#define DrawMOLFile 			DRAWMOLFILE
#define SmilesToMWMF			SMILESTOMWMF
#define WMFDepictSmiles			WMFDEPICTSMILES
#define SmilesToWMFBuffer		SMILESTOWMFBUFFER
#define RTFDepictSmiles			RTFDEPICTSMILES
#define WMFDepictCTString		WMFDEPICTCTSTRING
#define GetDataFromClipboard		GETDATAFROMCLIPBOARD
#define RemoveAllButOne			REMOVEALLBUTONE
#define SaveCTFileFromClipboard		SAVECTFILEFROMCLIPBOARD
#define FetchDataFromClipboard		FETCHDATAFROMCLIPBOARD
#define PostMOLFileToClipboard		POSTMOLFILETOCLIPBOARD
#define PostMOLFileEditable		POSTMOLFILEEDITABLE
#define PostDataToClipboard		POSTDATATOCLIPBOARD
#define SMILESStringToMOLFile		SMILESSTRINGTOMOLFILE
#define MOLFileToSMILESString		MOLFILETOSMILESSTRING
#define MOLFileToSMILESExt		MOLFILETOSMILESEXT
#define MOLFileToSMARTSString		MOLFILETOSMARTSSTRING
#define NumberBytesToString		NUMBERBYTESTOSTRING
#define	AttributedSmilesToWMFBuffer	ATTRIBUTEDSMILESTOWMFBUFFER
#define CTStringToSmiles		CTSTRINGTOSMILES
#define CTStringToSmilesExt		CTSTRINGTOSMILESEXT

int FAR PASCAL _export WMFDepictSmiles(LPSTR smiles, LPSTR fname,
                                       int *xextp, int *yextp);
/*
 * Writes a WMF depiction of the structure defined by the SMILES string
 * smiles[] to the file named fname[].
 * The fuction returns TRUE if a depiction was created and FALSE otherwise.
 */

void FAR PASCAL _export SmilesToMWMF(LPSTR smiles,
				     double *mwp,
				     LPSTR mfbuffer, int nbuf);
/*
 * Computes the molecular weight of the molecule defined
 * by smiles and puts it into *mwp. The buffer mfbuffer is filled with
 * a '\0' delimited string of the molecular formula. The caller is
 * responsible for providing enough space for this string.
 */

int FAR PASCAL _export AttributedSmilesToWMFBuffer(LPSTR smiles,
					           LPSTR buffer, int bufsize,
                                                   int *xextp, int *yextp,
						   LPSTR coordinatestring,
						   LPSTR colorstring,
						   LPSTR labelstring);
/*
 * Writes a WMF depiction of the structure defined by the SMILES string
 * smiles[] to the buffer buffer[0..bufsize-1].
 *
 * The fuction returns the number of bytes written to buffer[] or 0
 * in case of failure.
 *
 * The string coordinatestring[] contains a comma separated list of
 * 2D coordinates, i.e. 2 values per atom. They are used to place the
 * the atoms of the SMILES that are not implicit hydrogen atoms.
 * The implicit ones are then positioned with the layout algorithm.
 * If coordinatestring == NULL, all atoms will be layed out.
 * Atoms with missing coordinates are also layed out. Only relative
 * positioning is considered and only connected atoms are layed out
 * as one unit.
 * E.g.: "1.54,1.54,,,3.08,3.08"
 *
 * The string colorstring[] lists the colors to be used for the different
 * atoms of the SMILES. If NULL, then a standard atom symbol dependent
 * scheme is used. Only 16 colors are supported.
 * E.g.: "0,1,1,1,2,0,0,0,3,4,5"
 *
 * The string labelstring[] lists comma separated values for
 * alternate atom labels. If this string is NULL, then standard atom
 * labels will be used. This string can also contain empty values in which
 * case the atom symbol will be used, too.
 * E.g.: "CH3,,X,Y,Z,,"
 *
 * The function prepends the ALDUS placeable header to the buffer.
 *
 * It returns the size of the buffer used if successful, and a number <= 0
 * on failure. 0 means general error, a number < 0 is the size of the necessary
 * buffer to hold all the data. buffer[] is not changed in this case.
 */

int FAR PASCAL _export SmilesToWMFBuffer(LPSTR smiles,
				         LPSTR buffer, int bufsize,
				         int *xextp, int *yextp);
/*
 * Writes a WMF depiction of the structure defined by the SMILES string
 * smiles[] to the buffer buffer[0..bufsize-1].
 *
 * The fuction returns the number of bytes written to buffer[], 0
 * in case of failure, or -n where n is the required size of the buffer
 * if the buffer is too small.
 *
 * The function prepends the ALDUS placeable header to the buffer.
 */

int FAR PASCAL _export RTFDepictSmiles(LPSTR smiles, LPSTR fname,
                                       int *xextp, int *yextp);
/*
 * Writes a WMF depiction of the structure defined by the SMILES string
 * smiles[] in hex-code form to the file named fname[].
 * This form can directly be used for WMF fields in RTF files.
 */

int FAR PASCAL _export WMFDepictCTString(LPSTR ct_string, LPSTR fname,
                                         int *xextp, int *yextp);
/*
 * Writes a WMF depiction of the structure defined by the connection
 * table string ct_string[] to the file named fname[].
 * The fuction returns TRUE if a depiction was created and FALSE otherwise.
 * *xextp and *yextp are set to the dimensions of the resulting picture.
 */

int FAR PASCAL _export CTStringToSmiles(LPSTR ct_string,
                                        LPSTR smiles, int nbuf);
/*
 * Write the isomeric SMILES corresponding to the MOL file in ct_string into
 * the buffer smiles[0..nbuf-1].
 *
 * The function returns the size of the SMILES written, 0 on error, and
 * the negative required size of the buffer if not enough space was provided.
 *
 * Note: The buffer must also provide space for the terminal '\0'.
 */

LPSTR FAR PASCAL _export GetDataFromClipboard(LPINT sizep, LPSTR format);
/*
 * Return a pointer to the data of the given format or NULL if
 * the format is no available.
 */

void FAR PASCAL _export PostMOLFileToClipboard(LPSTR fname);
/*
 * Reads *fname and posts the contents to the clipboard in
 * MDLCT format.
 */

void FAR PASCAL _export PostDataToClipboard(LPSTR data, int size, LPSTR fmtstr);
/*
 * Posts the data contained in data[0..size-1] to the clipboard as format fmtstr.
 * Registers fmt if necessary. The procedure is designed to work with
 * null terminated strings.
 */

void FAR PASCAL _export SaveCTFileFromClipboard(LPSTR fname);
/*
 * Fetches the connection table from the clipboard and places it
 * into file *fname. The file is cleared if there is no MDLCT format
 * available.
 */

int FAR PASCAL _export FetchDataFromClipboard(LPSTR buffer,
                                              int   bufsize,
                                              LPSTR fmtname);
/*
 * Fetches data of the format fmtname from the clipboard and places them
 * into buffer. Returns the actual number of bytes retrieved if successful.
 * If bufsize is not big enough, -1 is returned.
 * If fmtname is not on the clipboard, FetchDataFromClipboard returns 0.
 */

void FAR PASCAL _export SMILESStringToMOLFile(LPSTR smiles, LPSTR fname);
/*
 * Converts the SMILES string *smiles to a MOL-file named *fname.
 *
 * The file is cleared in case of error.
 */

LPSTR FAR PASCAL _export MOLFileToSMILESString(LPINT sizep, LPSTR fname);
/*
 * Converts the MOL file stored in *fname to the corresponding
 * SMILES string and returns a pointer to it. The data becomes
 * private to the DLL. The caller must therefore take care of
 * Copying it to his data area.
 *
 * In case of error, the function returns a NULL pointer.
 */

LPSTR FAR PASCAL _export MOLFileToSMILESExt(LPINT sizep,
				            LPSTR fname,
					    int flags);
/*
 * Converts the MOL file stored in *fname to the corresponding
 * SMILES string and returns a pointer to it. The data becomes
 * private to the DLL. The caller must therefore take care of
 * Copying it to his data area.
 *
 * The parameter flags defines the kind processing required. It
 * can request that stereochemistry be ignored (UNIQUESMILES) and
 * that a comma separated list of coordinates for each atom in the
 * SMILES be appended.
 *
 * In case of error, the function returns a NULL pointer.
 */

LPSTR FAR PASCAL _export MOLFileToSMARTSString(LPINT sizep, LPSTR fname);
/*
 * Converts the MOL file query stored in *fname to the corresponding
 * SMARTS string and returns a pointer to it. The data becomes
 * private to the DLL. The caller must therefore take care of
 * Copying it to his data area.
 *
 * In case of error, the function returns a NULL pointer.
 */

void _export FAR PASCAL MOLFileToSMARTSBuffer(LPINT sizep, LPSTR buffer,
                                              LPSTR fname);
/*
 * Converts the MOL file query stored in *fname to the corresponding
 * SMARTS string and writes the result to buffer[]. Temporary storage is
 * removed. (*sizep) shall contain the buffer size on input and will
 * become the actual size of the copied SMARTS on output.
 *
 * In case of error, the function sets (*sizep) to -1 or -required_buffer_size
 * if the size of the buffer wasn't sufficient.
 */

#ifdef __cplusplus
}
#endif
