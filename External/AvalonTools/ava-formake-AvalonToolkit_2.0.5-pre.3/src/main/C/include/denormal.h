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
/*    File:           denormal.h                                        */
/*                                                                      */
/*    Purpose:        This file declares the functions and data types   */
/*                    for denormalization of connection tables.         */
/*                                                                      */
//    History:        25-Apr-1994     Start of development.             */
//                                                                      */
/************************************************************************/

struct atom_constraint_t
   {
      int  dbl_min;
      int  dbl_max;
      int  is_oxigen;
      int  open;
   };

struct bond_constraint_t
   {
      char can_be_single;
      char can_be_double;
      char scount;        /* total substitution count of atoms of bond */
      char open;
   };

#define MULTIPLE_COMPOSITION   1
#define CONSTRAINT_VIOLATION   2

#define CANDIDATE      1       /* atom/bond is a candidated */
#define TOUCHED        2       /* constraint set might be affected */

extern
int SetupConstraints(struct reaccs_molecule_t *mp,
                     struct atom_constraint_t *acp,
                     struct bond_constraint_t *bcp);
/*
 * This function tries to do its best to setup the constraints for
 * denormalization starting from the minimal atom and bond information
 * stored in *mp. It uses heuristic rules to estimate the hydrogen
 * count range of possibly tautomeric groups.
 * The function returns TRUE when it succeeds and FALSE otherwise.
 */

extern
int Denormalize(struct reaccs_molecule_t *mp,
                struct atom_constraint_t *acp_init,
                struct bond_constraint_t *bcp_init);
/*
 * This function changes the AROMATIC bonds of *mp to single and
 * double bonds. This process is the inverse of the CAS normalization
 * of alternating singe/double bonds and tautomer groups.
 *
 * In the case of tautomeric groups, the hydrogen count fields of the
 * atoms are assumed to contain information on the allowed range in its
 * two bytes.
 *
 * Non-tautomeric groups (i.e. rings) ignore the hydrogen count.
 *
 * The initial constraints are assumed to be set up in *acp and *bcp.
 */

int CCTDenormalize(struct reaccs_molecule_t *mp);
/*
 * Denormalizes the molecule *mp which was the result of decoding
 * a CAS CCT entry. It contains the valence and hydrogen range
 * information in the dummy1 and dummy2 fields of the atom structures.
 *
 * The function returns FALSE if anything went wrong in denormalization,
 */

