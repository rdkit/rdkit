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
/*      File:           utilities.h                                     */
/*                                                                      */
/*      Purpose:        This file contains the prototypes for the       */
/*                      utility functions implemented in utilities.c.   */
/*                                                                      */
/************************************************************************/

#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include "reaccs.h"
#include "symboltable.h"

#define NO_COLOR 0

struct ptable_entry
   {
      char *symbol;             /* atomic symbol */
      int  valence;             /* common valence of element */
      float mass;               /* average mass of atoms */
      float electronegativity;  /* el.neg. according to Pauling */
      float heat_of_formation;  /* h.o.f. of element from atoms */
   };

extern struct ptable_entry ptable[];

#define MAXNEIGHBOURS   20      /* used to be 12, but need to have more space */

typedef struct NEIGHBOURHOOD_T
   {
      short n_ligands;
      unsigned short atoms[MAXNEIGHBOURS];      /* indices into atom_array */
      unsigned short bonds[MAXNEIGHBOURS];      /* indices into bond_array */
   } neighbourhood_t;

extern char *halogenes[];                       /* "X" or "hal" */

extern char *non_metal_hetero_elements[];     /* "Q" */

extern char *metals[];                         /* "M" */

extern char *ONS_table[];                      /* "ONS" or "ons" */

   
extern char *transition_metals[];                /* "trn" */

extern char *lanthanoids[];                  /* "lan" */

extern char *alkali_metals[];                        /* "alk" */

extern char *non_metal_small_solution[];     /* "Qs" */

extern char *amino_acids[];                   /* "Ami" or "ami"*/

extern string_int_table periodic_table[];

                             /* struchk related atom type lists */
extern char *gr2[];
extern char *gr3[];
extern char *gr4[];
extern char *on2[];
extern char *ha2[];

extern char *tra[];
extern char *trb[];
extern char *tm1[];
extern char *tm2[];
extern char *tm3[];
extern char *tm4[];
extern char *tm5[];
extern char *tm6[];
extern char *tm7[];
extern char *tm8[];

extern
struct reaccs_molecule_t *NewMolecule(int n_atoms, int n_bonds);
/*
 * Allocates memory for a molecule with n_atoms atoms and n_bonds bonds.
 */

extern int CountSTextLines(struct stext_line_t *stext_lines);
/*
 * Counts the length of the list stext_lines.
 */

extern void FreeSTextLines(struct stext_line_t *stext_lines);
/*
 * Frees the storage allocated to the STEXT lines stext_lines.
 */

extern int CountPropLines(struct prop_line_t *prop_lines);
/*
 * Counts the length of the list prop_lines.
 */

extern void FreePropLines(struct prop_line_t *prop_lines);
/*
 * Frees the storage allocated to the property lines prop_lines.
 */

void AddNumProperty(struct reaccs_molecule_t *mp,
                    char *tag,
              int atno,
               int prop);
/*
 * Adds the atom property prop for atom atno to the list of property
 * lines in mp using the prefix tag[].
 */

extern int GetNumProperty(struct prop_line_t *prop_lines,
                          char *tag,
                        int atno,
                       int *prop);
/*
 * Fetches the property identified by tag[] for atom atno from the
 * property strings linked at prop_lines. Returns TRUE if a property
 * entry was found for that atom and FALSE otherwise. *prop is set to
 * the value found.
 */

extern
void FreeMolecule(struct reaccs_molecule_t *mp);
/*
 * Deallocates the storage for the REACCS molecule *mp.
 */

void FreeMoleculeChildObjects(struct reaccs_molecule_t *mp);
/*
 * Deallocates the storage for the REACCS molecule *mp leaving the struct itself allocated.
 */

extern
void FreeMoleculeList(struct reaccs_molecule_t *mp);
/*
 * Deallocates the storage occupied by the molecules listed
 * in mp, mp->next, ...
 */

extern
struct reaccs_reaction_t *NewReaction(void);
/*
 * Allocates memory for an empty reaction.
 */

extern
void FreeReaction(struct reaccs_reaction_t *rp);
/*
 * Deallocates the storage for the REACCS reaction *rp.
 */

extern
struct reaccs_molecule_t *CopyMolecule(struct reaccs_molecule_t *mp);
/*
 * Allocates space for a new molecule and copies the contents of *mp
 * to that space. The function returns a pointer to the new molecule.
 */

extern
int AddAtomToMolecule(struct reaccs_molecule_t *mp,
                      double x, double y, double z,
                      char *symbol);
/*
 * Adds a new atom to the structure pointed to by *mp. This atom is
 * not yet connected. It is placed at the given coordinates x, y, z,
 * and is of type *symbol. The function returns the number of the
 * new atom.
 */

extern
void AddBondToMolecule(struct reaccs_molecule_t *mp,
                       int from_atom, int to_atom,
                       int bond_type, int stereo_symbol);
/*
 * Adds a new bond from from_atom to to_atom to the structure pointed
 * to by *mp. The bond get the type bond_type and the stereosymbol
 * stereo_symbol.
 */

extern
struct reaccs_reaction_t *CopyReaction(struct reaccs_reaction_t *rp);
/*
 * Allocates space for a new reaction and copies the contents of *mp
 * to that space. The function returns a pointer to the new reaction.
*/

extern struct reaccs_reaction_t *
ConnectedReactionComponents(struct reaccs_reaction_t *rp);

extern
void RemoveAtomFromMolecule(struct reaccs_molecule_t *mp, int atom);
/*
 * Removes atom from molecule *mp.
 *
 * The atoms are numbered starting from 1!
 * This function has only been tested with isolated atoms although it
 * should also work with multiply connected ones.
 */

extern
struct reaccs_molecule_t *ThreadQueryAtoms(struct reaccs_molecule_t *query);
/*
 * Re-sorts the atoms in query such that (if possible) for all
 * n < query->n_atoms, atom n+1 is connected to one from 1..n.
 *
 * This is needed to improve the recursive phase of substruture matching and
 * shall only be done for atoms that don't have complications like symbol
 * lists and the like. Flag-based stereochemistry is ignored!
 *
 * Note: any neighbourhood arrays may become invalid when this method is called.
 */

extern
void StripMolecule(struct reaccs_molecule_t *mp,
                   int good_atoms[], int good_bonds[]);
/*
 * Deletes from molecule *mp all atoms i and bonds j for which
 * good_atoms[i] == FALSE, and good_bonds == FALSE, respectively.
 * The affected stereosymbols are adjusted.
 * sizeof(good_atoms) = n_atoms+1 !!
 */

extern
void SplitUnconnectedMolecules(struct reaccs_reaction_t *rp);
/*
 * Splits the molecules listed as reactants and products of *rp
 * into their connected components.
 */

extern
struct reaccs_molecule_t *SplitBond(struct reaccs_molecule_t *mp,
                                    int ibond,
                                    char *atsym1, int charge1, int mdiff1,
                                    char *atsym2, int charge2, int mdiff2);
/*
 * Splits the bond with index ibond attaching new atoms with the given
 * characteristics to the residual. The coordinates of the new atoms will be
 * on top of the old ones. It may therefore be needed to convert to SMILES or use
 * SplitMolecule to separate and visualize the result.
 *
 * The change is affected in-place and the modified molecule is returned.
 */

extern
struct reaccs_molecule_t *SplitMolecule(struct reaccs_molecule_t *mp);
/*
 * Checks if mp contains more than one fragment, modifies *mp to be
 * one of them and returns a pointer to a structure containing the others
 */

extern
struct reaccs_molecule_t *SplitMoleculeList(struct reaccs_molecule_t *mp);
/*
 * Scans the molecule list pointed to by mp and returns a list of
 * of all the fragments found therein.
 * The original molecules are destroyed.
 */

extern
struct reaccs_molecule_t *RemoveEmptyMolecules(struct reaccs_molecule_t *mp);
/*
 * Scans the molecule list pointed to by mp, removes empty molecules, and
 * returns the resulting list.
 * The original list of molecules is destroyed.
 */

#define NON_STEREO_HYDROGENS	 1
#define ANCILLARY_STEREO	 2
#define ESSENTIAL_STEREO	 4
#define ISOTOPE_HYDROGENS	 8
#define NO_QUERY_H_COUNT	16

#define ALL_HYDROGENS	(NON_STEREO_HYDROGENS 	|\
			 ANCILLARY_STEREO	|\
			 ESSENTIAL_STEREO	|\
			 ISOTOPE_HYDROGENS)

extern
void MakeSomeHydrogensImplicit(struct reaccs_molecule_t *mp,
			       int which_hydrogens);
/*
 * purpose:     removes explicit hydrogens from molecule *mp;
 *              stereochemistry is preserved as far as possible.
 *		The bit flags in 'which_hydrogens' define which hydrogens
 *		will be candidates for removal.
 */

extern
void MakeHydrogensImplicit(struct reaccs_molecule_t *mp);
/*
 * purpose:     removes all explicit hydrogens from molecule *mp;
 *              stereochemistry is preserved as farr as possible
 */

extern
void RemoveQueryOptions(struct reaccs_molecule_t *mp);
/*
 * Removes all query options found in the molecule *mp.
 * Currently, it only deals with HCounts.
 */

extern
struct reaccs_molecule_t *FuseMoleculeList(struct reaccs_molecule_t *mpl);
/*
 * Merges all molecules on list *mlp into one and returns the result.
 * Global structure components are taken from the first list element.
 */

extern
struct reaccs_molecule_t *LinkMolecules(struct reaccs_molecule_t *mp1,
                                        unsigned int              at1,
                                        int                       bd,
                                        unsigned int              at2,
                                        struct reaccs_molecule_t *mp2);
/*
 * Creates a new molecule which combines *mp1 and *mp2. The two
 * fragments are linked at the former atoms at1 and at2 with a
 * bond of type bd. The return value is a pointer to the constructed
 * molecule. Global structure components are taken from mp1.
 */

extern void ResetColors(struct reaccs_molecule_t *mp);
/*
 *  Purpose:    Sets the color fields of all bond and atom structures in
 *              the molecule *mp to 0.
 */

extern void ResetValues(struct reaccs_molecule_t *mp);
/*
 *  Purpose:    Sets the value fields of all bond and atom structures in
 *              the molecule *mp to 0.0.
 */

extern
int LookupMassDifference(double mass, char *symbol);
/*
 * Converts the mass value for an atomic symbol into a MOL file mass difference
 * value.
 */

extern
int AtomSymbolMatch(char *atsym,char *pattern);
/*
 * Returns TRUE if atsym is in the komma delimited list of atom symbols
 * stored in pattern and FALSE otherwise.
 */

extern
int MapToQueryAtomType(char *symbol,
                       char *source_table[],
                       char *target_symbol);
/*
 * If the string symbol is listed in source_table[], then it is
 * replaced by the string target_symbol.
 * This function is used to map sets of atom types to a query atom
 * type symbol.
 */

void ApplyToAllMolecules(struct reaccs_reaction_t *rp,
                         void (*func)(struct reaccs_molecule_t *));
/*
 *  Purpose:    Executes func(mp) for all molecules, i.e. reactants and
 *              products, in *rp.
 */

extern void ApplyToAllProducts(struct reaccs_reaction_t *rp,
                           void (*func)(struct reaccs_molecule_t *));
/*
 *  Purpose:    Executes func(mp) for all products in *rp.
 */

extern
int CheckNeighbourhood(struct reaccs_molecule_t *mp);
/*
 * Checks if the number of neighbours of the atoms in *mp is within the allowed limits.
 */

extern
int SetupNeighbourhood(struct reaccs_molecule_t *mp,
                       neighbourhood_t *neighbour_array,
                       int nlimit);
/*
 * Computes the array of atom and bond neighbourhoods of the
 * atoms in *mp which have an atom number less than nlimit.
 *
 * Returns FALSE in case of failure.
 */

extern                                          /* contains one entry for */
symbol_entry_t bond_to_string[];                /* each MDL bond type */

extern                                  /* contains one entry per */
symbol_entry_t charge_to_string[];      /* MDL defined charge */

extern                                  /* contains one entry per */
symbol_entry_t radical_to_string[];     /* MDL defined radical */

extern
void SortNeighbourhood(neighbourhood_t *np,
                       struct reaccs_molecule_t *mp);
/*
 * Sorts the neighbourhood *np with respect to the bond orders and
 * atom symbols defined in *mp.
 */

extern
void SortNeighboursByNumbering(struct reaccs_molecule_t *mp,
                               neighbourhood_t *nbp,
                               int numbering[]);
/*
 * Sorts the atom neighbourhood_t objects of nbp[0..mp->n_atoms] using the
 * atom ranking in numbering[].
 */

extern
void SortBondsWithNumbering(struct reaccs_molecule_t *mp,
                            int numbering[]);
/*
 * Sorts the bond_array member of *mp using numbering.
 * Note: the sorting preserves invariants with stereo-bonds.
 */

extern void
ComputeImplicitH(struct reaccs_molecule_t *mp,
                 int H_count[]);
/*
 * Computes the implicit hydrogen counts for the atoms in *mp. The
 * result is returned in H_count[] starting at H_count[1] for atom
 * number 1. H_count[0] is not changed.
 */

extern
int NewParity(unsigned int atn,
              int          old_parity,
              unsigned int n_bonds, struct reaccs_bond_t bonds[],
              int          renumbering[]);
/*
 * Computes the new parity for atom number atn when the structure
 * is renumbered using renumbering[].
 */

extern int RemoveIsolatedAtoms(struct reaccs_molecule_t *mp);
/*
 * Scans *mp for atoms with no bonds attached to them and then deletes
 * those atoms.
 * Renumbers the bonds accordingly.
 * It returns the number of atoms removed.
 */

extern
void RingState(struct reaccs_molecule_t *mp,
               int atom_status[],
            int bond_status[]);
/*
 * Computes how many basis rings each bond shares and how many
 * ring bonds are attached to an atom. The results are stored in
 * atom_status[] and bond_status[] respectively.
 */

extern
void MakeRingsClockwise(struct reaccs_molecule_t *mp);
/*
 * Orients the bonds in *mp such that the small basis rings are
 * traced in a clockwise direction.
 *
 */

extern double MolecularWeight(struct reaccs_molecule_t *mp);
/*
 * Computes and returns the molecular weight of the molecule *mp.
 */

extern int MolecularFormula(struct reaccs_molecule_t *mp, char formula[], int bufsize);
/*
 * Collects the molecular formula of *mp and writes it to formula[].
 * Disconnected fragments are separated by a "." character.
 * It returns the number of fragments found.
 * bufsize is the usable length of the formula buffer including the terminal '\0'.
 */

extern int TotalCharge(struct reaccs_molecule_t *mp);
/*
 * Returns the total charge of the atoms in *mp.
 */

extern int CombineChargeRadical(int charge, int radical);
/*
 * Returns MDL's combined charge/radical code.
 */

extern void SplitChargeRadical(int charge_radical, int *chargep, int *radicalp);
/*
 * Splits MDL's combined charge/radical code into it's components.
 */

extern int IsFieldHeader(char *header, char *field);
/*
 * Returns TRUE if header is a valid SD file data section header
 * for the field *field.
 */

void FlipStereoSymbols(struct reaccs_molecule_t *mp, int color);
/*
 * Flips the stereo symboles of the parts of the molecule *mp that are
 * colored with color.
 */

void FlipMolecule(struct reaccs_molecule_t *mp, int color);
/*
 * Flips the parts of the molecule *mp that are colored with color
 * with respect to the vertical (Y) axis. If color is 'NONE',
 * the whole molecule is flipped.
 */

void FloodFillFragments(struct reaccs_molecule_t *mp);
/*
 * The color of colored atoms is smeared out to the whole fragment
 * containig that atom. If there are more than one distict color,
 * the winning one on each fragment is arbitrary.
 */

int FloodWithColor(struct reaccs_molecule_t *mp,
                   neighbourhood_t *nbp,
	               int aindex, int color);
/*
 * Recursively colors uncolored parts of *mp starting with the atom with
 * index aindex with color. It uses a flood fill algorithm.
 * It returns the number of atoms recolored during the process.
 */

void StripColoredPart(struct reaccs_molecule_t *mp, int color);
/*
 * Removes all atoms and adjacent bonds from *mp that are colored with color.
 */

int
ImplicitHydrogens(char *symbol,
                  int   nsingle,
                  int   naromatic,
                  int   ndouble,
                  int   ntriple,
                  int   radical,
                  int   charge);
/*
 * Computes the number of implicit hydrogens attached to an atom of type
 * symbol with nsingle single bonds, naromatic aromatic bonds, ndouble
 * double bonds, ntriple triple bonds, and radical and charge state
 * radical and charge, resp.
 */
extern void ShowStack(char header[]);

#endif
