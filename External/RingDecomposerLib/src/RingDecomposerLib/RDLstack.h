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
 * J. Chem. Inf. Model., 2012, 52, pp 2013-2021
 * 
 * Flachsenberg, F.; Andresen, N.; Rarey, M.
 * RingDecomposerLib: An Open-Source Implementation of
 * Unique Ring Families and other Cycle Bases.
 * SUBMITTED
 */

#ifndef RDL_STACK_H
#define RDL_STACK_H

typedef struct RDL_stack RDL_stack;

RDL_stack* RDL_stack_new(void);

void* RDL_stack_top(RDL_stack* stack);

int RDL_stack_empty(RDL_stack* stack);

void RDL_stack_pop(RDL_stack* stack);

void RDL_stack_push(RDL_stack* stack, void* element);

void RDL_stack_delete(RDL_stack* stack);

#endif
