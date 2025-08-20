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
/*                                                                    	*/
/*    File:           reaccsio.h                                      	*/
/*                                                                    	*/
/*    Purpose:        This file contains the prototypes for the       	*/
/*                    functions of reaccsio.c.                        	*/
//                                                                    	*/
//    History:        04-Jan-93       Added some comments.            	*/
//                                                                    	*/
//		      26-Jul-98	      Added declarations of routines	*/
//				      reading/writing molecules from	*/
//				      and to strings (code from AG)	*/
/*                                                                    	*/
/************************************************************************/

#ifndef _REACCSIO_
#define _REACCSIO_        1

#include "forio.h"
#include "reaccs.h"

/*
 * Strip or don't strip trailing zeros on MOL file output.
 * Returns previous setting.
 */
extern
int SetStripZeroes(int do_strip);

extern
int ReadREACCSAtom(Fortran_FILE *fp,struct reaccs_atom_t *ap);

extern
void PrintREACCSAtom(FILE *fp,struct reaccs_atom_t *ap);

extern
int ReadREACCSBond(Fortran_FILE *fp,struct reaccs_bond_t *bp);

extern
void PrintREACCSBond(FILE *fp,struct reaccs_bond_t *bp);

extern
int ReadREACCSMolecule(Fortran_FILE *fp,
                       struct reaccs_molecule_t *mp,
                       const char label[]);
/*
 * Reads the next molecule description in Fortran_FILE *fp starting
 * with a line that begins with label[]. The molecule is put into *mp
 * which must be preallocated by the caller.
 * The molecule must start at the current line if label == "".
 */

extern
void PrintMACCSMolecule(FILE *fp,struct reaccs_molecule_t *mp, char *label);

extern
void PrintREACCSMolecule(FILE *fp,struct reaccs_molecule_t *mp, const char *label);

extern
struct reaccs_reaction_t *ReadREACCSReaction(Fortran_FILE *fp);

extern
void PrintREACCSReaction(FILE *fp,struct reaccs_reaction_t *rp);

extern
struct reaccs_molecule_t *ReadReactionAgents(Fortran_FILE *fp);
/*
 * Reads the reagtion agents (solvents and catalysts) for a
 * reaction entry. It returns a list of those molecules.
 */

struct data_line_t *ReadMACCSDataLines(Fortran_FILE *fp);
/*
 * Reads data lines from FORTRAN file *fp, constructs a data list,
 * and returns this list.
 */

struct data_line_t *ReadREACCSDataLines(Fortran_FILE *fp);
/*
 * Reads data lines from FORTRAN file *fp, constructs a data list,
 * and returns this list.
 */

extern
char *bond_type_names[];


extern
struct reaccs_molecule_t * MolStr2Mol( char * MolStr );
/*
 * Convert the MolFile in MolStr to a molecule struct
 * this should be deallocated with: FreeMolecule()
 */



extern
char * MolToMolStr( struct reaccs_molecule_t * mp );
/*
 * convert a molecule struct to a char * MolFile
 * the returned string should be free()ed
 * In case of problems NULL will be returned and an ErrorMessage appended
 * to the message list
 */

extern
struct symbol_list_t *ParseV30SymbolList(char *symbol,
                                         int iatom,
                                         struct reaccs_molecule_t *mp,
                                         struct symbol_list_t *old_list);
/*
 * Parsed the contents of symbol into a new symbol_list_t structure and
 * prepends it to old_list. atom_symbol is set to 'L' if parsing was OK
 * and to 'Unk' otherwise.
 */
#endif
