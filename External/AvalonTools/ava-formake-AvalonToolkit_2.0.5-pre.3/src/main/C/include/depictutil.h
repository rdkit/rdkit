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
/*   File:           depictutil.h                                       */
/*                                                                      */
/*   Purpose:        Declares utility functions used by depict.c.       */
/*                   Also makes sure the names are all uppercase.       */
/*                                                                      */
/************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

#define SmilesToMWMF SMILESTOMWMF
void SmilesToMWMF(char* smiles, double *mwp, char* mfbuffer, int bufsize);
/*
 * Computes the molecular weight of the molecule defined
 * by smiles and puts it into *mwp. The buffer mfbuffer is filled with
 * a '\0' delimited string of the molecular formula. The caller is
 * responsible for providing enough space for this string.
 */

#define CTStringToSmiles CTSTRINGTOSMILES
int CTStringToSmiles(char* ct_string, char* smiles, int nbuf);
/*
 * Write the isomeric SMILES corresponding to the MOL file in ct_string into
 * the buffer smiles[0..nbuf-1].
 *
 * The function returns the size of the SMILES written, 0 on error, and
 * the negative required size of the buffer if not enough space was provided.
 *
 * Note: The buffer must also provide space for the terminal '\0'.
 */

#define SMILESStringToMOLFile SMILESSTRINGTOMOLFILE
void SMILESStringToMOLFile(char* smiles, char* fname);
/*
 * Converts the SMILES string *smiles to a MOL-file named *fname.
 *
 * The file is cleared in case of error.
 */

#define MOLFileToSMILESString MOLFILETOSMILESSTRING
char* MOLFileToSMILESString(int* sizep, char* fname);
/*
 * Converts the MOL file stored in *fname to the corresponding
 * SMILES string and returns a pointer to it. The data becomes
 * private to the DLL. The caller must therefore take care of
 * Copying it to his data area.
 *
 * In case of error, the function returns a NULL pointer.
 */

#define MOLFileToSMILESExt MOLFILETOSMILESEXT
char* MOLFileToSMILESExt(int* sizep, char* fname, int flags);
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

#define MOLFileToSMARTSString MOLFILETOSMARTSSTRING
char* MOLFileToSMARTSString(int* sizep, char* fname);
/*
 * Converts the MOL file query stored in *fname to the corresponding
 * SMARTS string and returns a pointer to it. The data becomes
 * private to the DLL. The caller must therefore take care of
 * Copying it to his data area.
 *
 * In case of error, the function returns a NULL pointer.
 */

#define MOLFileToSMARTSBuffer MOLFILETOSMARTSBUFFER
int MOLFileToSMARTSBuffer(int* sizep, char* buffer, char* fname);
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
