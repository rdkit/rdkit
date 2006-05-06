/*------------------------------------------------------------
 * sd_state.h
 * Interface of sd_state.cc
 * Definition of a class representing a state of the matching
 * process between two ARGs using the Schmidt-Druffel isomorphism
 * algorithm described in:
 * - D. C. Schmidt, L. E. Druffel, "A Fast Backtracking Algorithm
 *   to Test Directed Graphs for Isomorphism Using Distance Matrices",
 *   Journal of ACM, Vol. 23, No. 3, pp. 433-445, 1973.
 *
 * See: argraph.h state.h
 *
 * Author: P. Foggia
 *-----------------------------------------------------------------*/



#ifndef SD_STATE_H
#define SD_STATE_H

#include "argraph.h"
#include "state.h"



/*----------------------------------------------------------
 * class SDState
 * A representation of the SSR current state
 ---------------------------------------------------------*/
class SDState: public State
  { typedef ARGraph_impl Graph;

    private:
	  SDState *parent;
      int core_len, orig_core_len;
      Graph *g1, *g2;
      int n1, n2;
	  node_id *core1, *core2;
	  node_id **dist1, **dist2;
	  node_id *cls1, *cls2;
	  node_id *cnt1, *cnt2;

	  node_id *wrk1, *wrk2;

	  long *share_count;
	  bool dead_end;
    
    public:
      SDState(Graph *g1, Graph *g2);
      SDState(const SDState &state);
      ~SDState(); 
      Graph *GetGraph1() { return g1; }
      Graph *GetGraph2() { return g2; }
      bool NextPair(node_id *pn1, node_id *pn2,
                    node_id prev_n1=NULL_NODE, node_id prev_n2=NULL_NODE);
      bool IsFeasiblePair(node_id n1, node_id n2);
      void AddPair(node_id n1, node_id n2);
      bool IsGoal() { return core_len==n1 && core_len==n2; };
      bool IsDead() { return dead_end; };
      int CoreLen() { return core_len; }
      void GetCoreSet(node_id c1[], node_id c2[]);
      State *Clone();

	  void BackTrack();
  };


#endif

