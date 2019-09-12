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

#ifndef RDL_INFO_H
#define RDL_INFO_H

#include "RDLforward.h"
#include "RingDecomposerLib.h"

struct RDL_URFinfo {
  unsigned nofWeights; /*number of different weights of cycle families occuring*/
  unsigned *nofProtos; /*array storing the number of prototypes for each weight*/
  char ***URFrel; /*array of 2D-matrices. For each weight that a CF can have
                  there is a matrix that stores 1 at position [i,j] if RCFs i
                  and j (of this particular weight) are URF-related and 0
                  otherwise.*/
  unsigned nofURFs; /*number of URFs*/
  RDL_cfam ***URFs; /*array of tuples of RCFs that belong to the same URF*/
  unsigned *nofCFsPerURF; /*array that stores the number of RCFs that belong to each
                     URF*/
};

/** deletes the URFinfo */
void RDL_deleteURFInfo(RDL_URFinfo *);

/** Uses Gaussian elimination to mark potentially URF-related CFs and checks
the previously marked CFs for URF-relation. Fills the URFinfo. */
RDL_URFinfo *RDL_checkURFRelation(RDL_cfURF *, RDL_graph *, RDL_sPathInfo*);

#endif
