/*------------------------------------------------------------
 * vf2_mono_state.h
 * Interface of vf2_mono_state.cc
 * Definition of a class representing a state of the matching
 * process between two ARGs.
 * See: argraph.h state.h
 *
 * Author: P. Foggia
 *-----------------------------------------------------------------*/




#ifndef VF2_MONO_STATE_H
#define VF2_MONO_STATE_H

#include "argraph.h"
#include "state.h"



/*----------------------------------------------------------
 * class VF2MonoState
 * A representation of the SSR current state
 * See vf2_sub_state.cc for more details.
 ---------------------------------------------------------*/
class VF2MonoState: public State
  { typedef ARGraph_impl Graph;

    private:
      int core_len, orig_core_len;
      int added_node1;
      int t1both_len, t2both_len, t1in_len, t1out_len, 
          t2in_len, t2out_len; // Core nodes are also counted by these...
      node_id *core_1;
      node_id *core_2;
      node_id *in_1;
      node_id *in_2;
      node_id *out_1;
      node_id *out_2;

      node_id *order;

      Graph *g1, *g2;
      int n1, n2;

	  long *share_count;
    
    public:
      VF2MonoState(Graph *g1, Graph *g2, bool sortNodes=false);
      VF2MonoState(const VF2MonoState &state);
      ~VF2MonoState(); 
      Graph *GetGraph1() { return g1; }
      Graph *GetGraph2() { return g2; }
      bool NextPair(node_id *pn1, node_id *pn2,
                    node_id prev_n1=NULL_NODE, node_id prev_n2=NULL_NODE);
      bool IsFeasiblePair(node_id n1, node_id n2);
      void AddPair(node_id n1, node_id n2);
      bool IsGoal() { 
		  bool goal = core_len==n1;
		  //if(goal) printf("is goal\n");
		  return goal ; 
	  };

      bool IsDead() { 
		  bool dead = n1>n2  ||                      
			  t1both_len>t2both_len ||
              t1out_len>t2out_len ||
              t1in_len>t2in_len;  
		  //if(dead) printf("is dead!\n");
		  return dead;       
	  };

      int CoreLen() { return core_len; }
      void GetCoreSet(node_id c1[], node_id c2[]);
      State *Clone();
	  void printCore() {
		  printf("core length(%d)->", core_len);
		  for(int i=0; i<n1; i++){
			  if(core_1[i] != NULL_NODE)
				printf("(%d,%d),", i, core_1[i]);
		  }
		  printf("\n");
	  }
      virtual void BackTrack();
  };


#endif

