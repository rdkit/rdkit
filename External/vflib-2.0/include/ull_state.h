/*------------------------------------------------------------------
 * ull_state.h
 * Header of ull_state.cc
 * Definition of the class UllState
 *
 * Author: P. Foggia
 * $Id: ull_state.h,v 1.1 2001/10/24 01:20:11 glandrum Exp $
 *-----------------------------------------------------------------*/


/*-----------------------------------------------------------------
 * REVISION HISTORY
 *   $Log: ull_state.h,v $
 *   Revision 1.1  2001/10/24 01:20:11  glandrum
 *   added
 *
 *   Revision 1.1  1998/09/29 09:50:38  foggia
 *   Initial revision
 *
 *----------------------------------------------------------------*/


#ifndef ULL_STATE_H
#define ULL_STATE_H

#include "argraph.h"
#include "state.h"


/*----------------------------------------------------------
 * class UllState
 * A representation of the current search state
 * of the Ullmann's algorithm
 ---------------------------------------------------------*/
class UllState: public State
  { private:
      int core_len;
      node_id *core_1;
      node_id *core_2;
      Graph *g1, *g2;
      int n1, n2;
      byte **M;   // Matrix encoding the compatibility of the nodes

      void refine();
    
    public:
      UllState(Graph *g1, Graph *g2);
      UllState(const UllState &state);
      ~UllState(); 
      Graph *GetGraph1() { return g1; }
      Graph *GetGraph2() { return g2; }
      bool NextPair(node_id *pn1, node_id *pn2,
                    node_id prev_n1=NULL_NODE, node_id prev_n2=NULL_NODE);
      bool IsFeasiblePair(node_id n1, node_id n2);
      void AddPair(node_id n1, node_id n2);
      bool IsGoal() { return core_len==n1 && core_len==n2; };
      bool IsDead() { if (n1!=n2) return true;
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
