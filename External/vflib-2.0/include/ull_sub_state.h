/*------------------------------------------------------------------
 * ull_sub_state.h
 * Header of ull_sub_state.cc
 * Definition of the class UllSubState
 *
 * Author: P. Foggia
 * $Id: ull_sub_state.h,v 1.1 2001/10/24 01:20:11 glandrum Exp $
 *-----------------------------------------------------------------*/


/*-----------------------------------------------------------------
 * REVISION HISTORY
 *   $Log: ull_sub_state.h,v $
 *   Revision 1.1  2001/10/24 01:20:11  glandrum
 *   added
 *
 *   Revision 1.1  1999/01/06 20:04:04  foggia
 *   Initial revision
 *
 *
 *----------------------------------------------------------------*/


#ifndef ULL_SUB_STATE_H
#define ULL_SUB_STATE_H

#include "argraph.h"
#include "state.h"


/*----------------------------------------------------------
 * class UllSubState
 * A representation of the current search state
 * of the Ullmann's algorithm for graph-subgraph isomorphism
 ---------------------------------------------------------*/
class UllSubState: public State
  { private:
      int core_len;
      node_id *core_1;
      node_id *core_2;
      Graph *g1, *g2;
      int n1, n2;
      byte **M;   // Matrix encoding the compatibility of the nodes

      void refine();
    
    public:
      UllSubState(Graph *g1, Graph *g2);
      UllSubState(const UllSubState &state);
      ~UllSubState(); 
      Graph *GetGraph1() { return g1; }
      Graph *GetGraph2() { return g2; }
      bool NextPair(node_id *pn1, node_id *pn2,
                    node_id prev_n1=NULL_NODE, node_id prev_n2=NULL_NODE);
      bool IsFeasiblePair(node_id n1, node_id n2);
      void AddPair(node_id n1, node_id n2);
      bool IsGoal() { return core_len==n1; };
      bool IsDead() { if (n1>n2) return true;
                      for(int i=core_len; i<n1; i++)
                        { for(int j=0; j<n2; j++)
                            if (M[i][j]!=0) goto next_row;
                          return true;
                      next_row: ;
                        }
                       return false;
                    };
      int CoreLen() { return core_len; }
      void GetCoreSet(node_id c1[], node_id c2[]);
      State *Clone();
  };


#endif
