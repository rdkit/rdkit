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
/*   File:           geometry.h                                         */
/*                                                                      */
/*  Purpose:         Declares the functions and datatypes needed to     */
/*                        handle geometric objects like point sets.     */
/*                                                                      */
/************************************************************************/

#define USE_INWARDS      1
#define H_GEOMETRY       2
#define NO_RAY_CROSSING  4

extern
void TransformPoints(double r[][2], int n,
                     double p1[2],  double p2[2],
                     double p1p[2], double p2p[2]);
/*
 * Transforms the 2D points in r[...][2] such that the points p1[] and p2[]
 * would be transformed to p1p[] and p2p[], resp. in a best fit manner by
 * just rotating and moving them without scaling.
 */

extern
void PointSetMatchTransformation(double points[][2], unsigned npoints,
                                 double from[][2],
                                 double to[][2],     unsigned nmatch,
                                 int reflection);
/*
 * Computes the transformation which rotates and translates from[0..nmatch-1]
 * such that the points with the same index in from[] and to[] come as
 * close together as possible.  The resulting transformation is applied to
 * points[0..npoints-1]. If reflection is TRUE, the transformation will be
 * improper.
 */

extern
void NextSubstituentPoint(double point[2],
                          double coordinates[][2], unsigned nnodes,
                          unsigned edges[][2], unsigned nedges,
                          unsigned seed,
	                  int flags,
                          int numbers[], int natoms);
/*
 * Finds the next best point to substitute the graph defined by
 * coords[0..nnodes-1][2] and edges[0..nedges-1][2] at node seed.
 * The coordinates found are stored in points[0..1].
 * The parameter flags tells whether candidated positions may 
 * also be put between existing substituents and if a hydrogen needs to be layed out.
 *
 * numbers[0..natoms-1] contains the mapping from molecule to graph numbering.
 */

extern
void DebugGeometry(int flag);
