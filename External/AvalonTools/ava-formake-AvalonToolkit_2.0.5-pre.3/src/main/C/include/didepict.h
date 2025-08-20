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

/**************************************************************************/
/*                                                                        */
/*   File:   	didepict.h                                                */
/*		                                                          */
/*   Purpose:	Prototypes of exported functions of didepict.c which      */
/*              is used to generate device independent graphics commands  */
/*              from an MDL data structure or a SMILES string.            */
/*		These commands are primaryly intended to be intepreted by */
/*              a Java class.                                             */
/*                                                                        */
/*   Histroy:	1996-06-21	Start of development from the Windows-DLL */
/*				version depict.c.                         */
/*                                                                        */
/**************************************************************************/

#define USE_COLORS    0x001

void MoleculeToPlotFile(FILE *fp,
			struct reaccs_molecule_t *mp,
			int xext, int yext,
                        int flags,
                        int colors[],
			char *labels[]);
/*
 * Chemical structure depiction workhorse. Takes the MDL formatted
 * chemical structure *mp and writes device independent drawing
 * commands to *fp.
 * 
 * The flags parameter currently tells if color is to be used.
 *
 * The array colors[] is used to color the different parts of the molecule.
 * colors[i] == 0 means draw in black. 1 <= colors[i] <= 9 means use one of
 * the nine most dark colors to get reasonable readablility. If colors ==
 * NULL, the standard atom label coloring with black bonds is used.
 *
 * xext and yext are used for the desired x- and y- dimensions of the resulting
 * image in twips.
 *
 * The array labels[] contains alternate atom labels to be used instead
 * of the ones defined in the atom_symbol fields.
 * Labels which are "\0" strings are not considered. Use " " to get
 * a blank label.
 */

int DepictSmilesToPlotFile(FILE *fp, char *smiles, int xext, int yext);
/*
 * Writes a device independent metafile of the depiction of the structure
 * defined by the SMILES string smiles[] to the file *fp.
 *
 * The fuction returns TRUE if a depiction was created and FALSE in case of error.
 *
 * xext and yext are desired deminsions of the picture twips or 0 if the caller
 * doesn't care.
 */

int DepictRxnSmilesToPlotFile(FILE *fp, char *smiles);
/*
 * Writes a device independent metafile of the depiction of the reaction
 * defined by the reaction SMILES string smiles[] to the file *fp.
 *
 * The fuction returns TRUE if a depiction was created and FALSE in case of error.
 */

void SmilesToMWMF(char *smiles, double *mwp, char *mfbuffer, int bufsize);
/*
 * Computes the molecular weight of the molecule defined
 * by smiles and puts it into *mwp. The buffer mfbuffer is filled with
 * a '\0' delimited string of the molecular formula. The caller is
 * responsible for providing enough space for this string.
 * bufsize is the usable size of the buffer including the '\0' character.
 */

int AttributedSmilesToPlotFile(FILE *fp,
			       char *smiles,
			       int xext, int yext,
			       char *coordinatestring,
			       char *colorstring,
			       char *labelstring);
/*
 * Writes a device independent depiction of the structure defined by the
 * SMILES string smiles[] to the file *fp.
 *
 * xext and yext are desired deminsions of the picture twips or 0 if the caller
 * doesn't care.
 *
 * The fuction returns TRUE on success and FALSE in case of failure.
 *
 * The string coordinatestring[] contains a comma separated list of
 * 2D coordinates, i.e. 2 values per atom. They are used to place the
 * the atoms of the SMILE that are not implicit hydrogen atoms.
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
 */

