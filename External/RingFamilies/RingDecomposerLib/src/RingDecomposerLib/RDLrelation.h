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

#ifndef RDL_GAUSSELIM_H
#define RDL_GAUSSELIM_H

#include "RDLforward.h"
#include "RingDecomposerLib.h"

/** Finds potentially URFs-related CFs with the help of Gaussian Elimination,
checks if the relation in fulfilled and stores the result in the URFinfo. In
this step the non-relevant cycle families found before by Vismara's algorithm
are ignored.*/
void RDL_findRelations(RDL_cfURF *, RDL_graph *, RDL_URFinfo *, RDL_sPathInfo*);

#endif

