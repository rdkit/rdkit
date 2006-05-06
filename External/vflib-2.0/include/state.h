/*------------------------------------------------------------
 * state.h
 * Definition of an abstract class representing a state of the 
 * matching process between two ARGs.
 * See: argraph.h 
 *
 * Author: P. Foggia
 * $Id: state.h,v 1.1 2001/10/24 01:20:11 glandrum Exp $
 *-----------------------------------------------------------------*/


/*-----------------------------------------------------------------
 * REVISION HISTORY
 *   $Log: state.h,v $
 *   Revision 1.1  2001/10/24 01:20:11  glandrum
 *   added
 *
 *   Revision 1.3  1998/09/29 09:50:16  foggia
 *   Ora definisce una classe astratta State, da cui discendono
 *   VFState e UllState
 *
 *   Revision 1.2  1998/09/26 09:02:32  foggia
 *   minor changes
 *
 *   Revision 1.1  1998/09/19 14:40:35  foggia
 *   Initial revision
 *
 *----------------------------------------------------------------*/


#ifndef STATE_H
#define STATE_H

/*----------------------------
 * Flags used to encode the
 * state
 ---------------------------*/
enum { ST_CORE=0x01,
       ST_TERM_IN=0x02,
       ST_TERM_OUT=0x04
     };

#include "argraph.h"


/*----------------------------------------------------------
 * class State
 * An abstract representation of the SSR current state.
 * NOTE: Respect to pre-2.0 version of the library, class
 *   State assumes explicitly a depth-first search. The
 *   BackTrack method has been added to allow a state 
 *   to clean up things before reverting to its parent. This
 *   can be used, for instance, for sharing resources between
 *   the parent and the child. The BackTrack implementation
 *   can safely assume that at most one AddPair has been
 *   performed on the state.
 ---------------------------------------------------------*/
class State
  { 

    public:
      virtual ~State() {} 
      virtual Graph *GetGraph1()=0;
      virtual Graph *GetGraph2()=0;
      virtual bool NextPair(node_id *pn1, node_id *pn2,
                    node_id prev_n1=NULL_NODE, node_id prev_n2=NULL_NODE)=0;
      virtual bool IsFeasiblePair(node_id n1, node_id n2)=0;
      virtual void AddPair(node_id n1, node_id n2)=0;
      virtual bool IsGoal() =0;
      virtual bool IsDead() =0;
      virtual int CoreLen() =0;
      virtual void GetCoreSet(node_id c1[], node_id c2[]) =0;
      virtual State *Clone() =0;  // Changed clone to Clone for uniformity
     
      virtual void BackTrack() { };
  };


#endif

