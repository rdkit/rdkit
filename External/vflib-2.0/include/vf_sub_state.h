/*------------------------------------------------------------
 * vf_sub_state.h
 * Interface of vf_sub_state.cc
 * Definition of a class representing a state of the matching
 * process between two ARGs, for the graph-subgraph isomorphism
 * See: argraph.h state.h
 *
 * Author: P. Foggia
 * $Id: vf_sub_state.h,v 1.1 2001/10/24 01:20:11 glandrum Exp $
 *-----------------------------------------------------------------*/


/*-----------------------------------------------------------------
 * REVISION HISTORY
 *   $Log: vf_sub_state.h,v $
 *   Revision 1.1  2001/10/24 01:20:11  glandrum
 *   added
 *
 *   Revision 1.2  1999/01/05 10:28:44  foggia
 *   Bug fix
 *
 *   Revision 1.1  1998/12/08 13:31:06  foggia
 *   Initial revision
 *
 *----------------------------------------------------------------*/


#ifndef VF_SUB_STATE_H
#define VF_SUB_STATE_H

#include "argraph.h"
#include "state.h"


/*----------------------------------------------------------
 * class VFSubState
 * A representation of the SSR current state
 * for a graph-subgraph isomorphism.
 ---------------------------------------------------------*/
class VFSubState: public State
  { typedef ARGraph_impl Graph;

    private:
      int core_len;
      int t1in_len, t1out_len, t2in_len, t2out_len;
      node_id *core_1;
      node_id *core_2;
      byte *node_flags_1;
      byte *node_flags_2;
      Graph *g1, *g2;
      int n1, n2;
    
    public:
      VFSubState(Graph *g1, Graph *g2);
      VFSubState(const VFSubState &state);
      ~VFSubState(); 
      Graph *GetGraph1() { return g1; }
      Graph *GetGraph2() { return g2; }
      bool NextPair(node_id *pn1, node_id *pn2,
                    node_id prev_n1=NULL_NODE, node_id prev_n2=NULL_NODE);
      bool IsFeasiblePair(node_id n1, node_id n2);
      void AddPair(node_id n1, node_id n2);
      bool IsGoal() { return core_len==n1; };
      bool IsDead() { return n1>n2 || 
                      (t1out_len>t2out_len) ||
                      (t1in_len>t2in_len); 
                    };
      int CoreLen() { return core_len; }
      void GetCoreSet(node_id c1[], node_id c2[]);
      State *Clone();
  };


#endif

