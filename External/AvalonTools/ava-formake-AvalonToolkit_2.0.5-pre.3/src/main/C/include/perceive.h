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
/*    File:           perceive.h                                        */
/*                                                                      */
/*    Author:         B. Rohde                                          */
/*                                                                      */
/*    History:        31-May-89       Creation of header file.          */
/*                                                                      */
/*    Purpose:        Declares functions implemented in perceive.c.     */
/*                                                                      */
/************************************************************************/

extern
void PerceiveRingBonds(struct reaccs_molecule_t *mp);
/* Sets the topography of the bonds of *mp to the appropriate ring
 * class and flags the ring size in the dummy field.
 */

extern
void SetRingSizeFlags(struct reaccs_molecule_t *mp, int max_size,
                      neighbourhood_t nbp[]);
/*
 * Sets the perceived ring-sizes as flag bits in the bond rsize_flags field.
 * It used a recursive enumeration algorithm pruned to only ring bonds/atoms.
 * nbp[] is the neighbourhood array to speed up processing.
 */

extern
void PerceiveRingSizes(struct reaccs_molecule_t *mp);
/*
 * Sets the perceived ring-sizes as flag bits in the bond rsize_flags field.
 */

extern
void PerceiveAromaticBonds(struct reaccs_molecule_t *mp);
/*
 * Sets bond types of the bonds of *mp to "AROMATIC" if in
 * six-ring of sp2 atoms.
 */

extern
void PerceiveDYAromaticity(struct reaccs_molecule_t *mp,
			   neighbourhood_t nbp[]);
/*
 * Converts bonds in rings that are perceived as aromatic by Daylight
 * programs into AROMATIC bonds. In addition, hydrogen counts of aromatic
 * non-carbon atoms are remembered in the query_H_count fields.
 */

extern
void PerceiveMarkush(struct reaccs_molecule_t *mp,
		     neighbourhood_t nbp[]);
/*
 * Scans mp for "R" atoms. If there are any, it sets the sub_desc entries
 * of all atoms not attached to an R to s*. This is used to simplify
 * substructure query definition by only specifying free sites.
 */
