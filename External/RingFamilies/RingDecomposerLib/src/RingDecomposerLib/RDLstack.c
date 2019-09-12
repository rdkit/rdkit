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

#include "RDLstack.h"

#include <stdlib.h>

static const unsigned MinSize = 32;

struct RDL_stack {
    void** elements;
    void** top;
    unsigned size;
    unsigned reserved;
};

RDL_stack* RDL_stack_new(void)
{
  RDL_stack* stack;

  stack = malloc(sizeof(*stack));
  stack->size = 0;
  stack->reserved = MinSize;
  stack->elements = malloc(stack->reserved * sizeof(*stack->elements));
  stack->top = stack->elements-1;

  return stack;
}

void* RDL_stack_top(RDL_stack* stack)
{
  if (stack->size > 0) {
    return *stack->top;
  }
  else {
    return NULL;
  }
}

int RDL_stack_empty(RDL_stack* stack)
{
  return (stack->size == 0);
}

void RDL_stack_pop(RDL_stack* stack)
{
  unsigned distance;

  if (stack->size > 0) {
    --stack->size;
    --stack->top;

    if (stack->size * 2 < stack->reserved && stack->reserved > MinSize) {
      distance = stack->top - stack->elements;

      stack->reserved = stack->reserved / 2;

      stack->elements = realloc(stack->elements,
          stack->reserved * sizeof(*stack->elements));

      stack->top = stack->elements + distance;
    }
  }
}

void RDL_stack_push(RDL_stack* stack, void* element)
{
  unsigned distance;

  if (stack->size == stack->reserved) {
    stack->reserved *= 2;

    distance = stack->top - stack->elements;

    stack->elements = realloc(stack->elements,
        stack->reserved * sizeof(*stack->elements));

    stack->top = stack->elements + distance;
  }
  ++stack->top;
  ++stack->size;
  *stack->top = element;
}

void RDL_stack_delete(RDL_stack* stack)
{
  free(stack->elements);
  stack->elements = NULL;
  stack->top = NULL;
  stack->reserved = 0;
  stack->size = 0;
  free(stack);
}
