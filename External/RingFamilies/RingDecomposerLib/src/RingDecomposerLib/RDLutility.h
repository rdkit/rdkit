/*
 * This file is part of the RingDecomposerLib, licensed
 * under BSD New license (see LICENSE in the root directory).
 * Copyright (c) 2016
 * University of Hamburg, ZBH - Center for Bioinformatics
 * Niek Andresen, Florian Flachsenberg, Matthias Rarey
 * 
 * Please cite:
 * 
 * Kolodzik, A.; Urbaczek, S.; Rarey, M.
 * Unique Ring Families: A Chemically Meaningful Description
 * of Molecular Ring Topologies.
 * J. Chem. Inf. Model., 2012, 52 (8), pp 2013-2021
 * 
 * Flachsenberg, F.; Andresen, N.; Rarey, M.
 * RingDecomposerLib: An Open-Source Implementation of
 * Unique Ring Families and Other Cycle Bases.
 * J. Chem. Inf. Model., 2017, 57 (2), pp 122-126
 */

#ifndef RDL_UTILITY_H
#define RDL_UTILITY_H

#include "RingDecomposerLib.h"

unsigned **RDL_alloc2DUIntArray(unsigned n, unsigned m);
char **RDL_alloc2DCharArray(unsigned n, unsigned m);
void RDL_delete2DUIntArray(unsigned **arr, unsigned n);
void RDL_delete2DCharArray(char **arr, unsigned n);

void RDL_print2DUIntArray(unsigned **arr, unsigned n, unsigned m);
void RDL_print2DCharArray(char **arr, unsigned n, unsigned m);

void RDL_writeToStderr(RDL_ERROR_LEVEL level, const char* fmt, ...);
void RDL_writeNothing(RDL_ERROR_LEVEL level, const char* fmt, ...);

#endif
