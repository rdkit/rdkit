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
/*  File:                     patclean.h                                */
/*                                                                      */
/*  Purpose:                  This module declares the functions        */
/*                            exported by patclean.c used to do a       */
/*                            template driven clean operation.          */
/*                                                                      */
//  History:  29-Sep-92       Created header file.                      */
//                                                                      */
/************************************************************************/

int ForceStereoTemplate(struct reaccs_molecule_t *mp,
                       struct reaccs_molecule_t *ssp);
/*
 * Finds the geometrically closest match of *ssp in *mp and
 * forces the stereosymbols of *ssp onto the matching bonds
 * of *mp.
 */

int TemplateClean(struct reaccs_molecule_t *mp,
                  struct reaccs_molecule_t *ssp);
/*
 * Checks if the substructure *ssp matches *mp, finds the
 * geometrically closest match, and tries to clean the matching
 * portion of *mp to the coordinates of *ssp, while moving the
 * substituents of the match to reasonable positions.
 * The function returns TRUE if the coordinates of *mp have
 * been modified and FALSE if not.
 */

int TemplateRotate(struct reaccs_molecule_t *mp,
                   struct reaccs_molecule_t *ssp);
/*
 * Checks if the substructure *ssp matches *mp, finds the
 * geometrically closest match, and tries to rotate *mp
 * to most closely resemble the coordinates of the pattern.
 */
