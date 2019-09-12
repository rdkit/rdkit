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

#ifndef RDL_CYCLEFAMS_H
#define RDL_CYCLEFAMS_H

#include "RDLforward.h"
#include "RingDecomposerLib.h"

/** a cycle family */
struct RDL_cfam {
  unsigned weight; /*number of edges/vertices in the cf (length of one cycle)*/
  unsigned r, p, q, x; /*vertices defining a cf as in Vismara's definition*/
  char *prototype; /*prototype vector {0,1}^m*/
  int mark; /*mark CF as relevant*/
};

/** contains the Cycle Families */
struct RDL_cfURF {
  RDL_cfam **fams; /*array of CFs sorted by weight*/
  unsigned nofFams; /*number of cycle families*/
  unsigned alloced; /*space for the cycle families*/
};


/** uses Vismara's algorithm to find (Relevant) Cycle Families.
Returns them sorted by weight */
RDL_cfURF *RDL_findCycleFams(RDL_graph *, RDL_sPathInfo *);

/** frees all the memory used by the Cycle Families */
void RDL_deleteCycleFams(RDL_cfURF *);

#endif
