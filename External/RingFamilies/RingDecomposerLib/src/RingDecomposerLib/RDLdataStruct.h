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

#ifndef RDL_DATASTRUCT_H
#define RDL_DATASTRUCT_H

#include "RDLforward.h"
#include "RingDecomposerLib.h"

/** struct that contains all data that was accumulated during the calculation of
the URF-structure. The user does not need to change or read
anything in this structure other than by using the functions provided in
exposed header. */
struct RDL_data {
  RDL_graph *graph;
  unsigned nofURFs; /*number of URFs*/
  unsigned nofRCFs; /*number of RCFS*/

  RDL_BCCGraph* bccGraphs;

  unsigned *nofURFsPerBCC; /*number of URFs*/
  unsigned *nofRCFsPerBCC; /*number of RCFs*/
  RDL_cfURF **CFsPerBCC; /*the cycle families (found by Vismara's algorithm)*/
  RDL_URFinfo **urfInfoPerBCC;/*stores which RCF are URF related and belong to the same URF*/
  RDL_sPathInfo **spiPerBCC; /*shortest paths info*/

  /*
   * this array maps the URFs to BCCs
   * at position 0 the BCC id is stored,
   * at position 1 the internal index of the URF in this BCC
   */
  unsigned (*urf_to_bcc)[2];
  /*
   * this array maps the RCFs to URFs
   * at position 0 the URF id is stored,
   * at position 1 the internal index of the RCF in this URF
   */
  unsigned (*rcf_to_urf)[2];

};

#endif
