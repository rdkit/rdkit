/*------------------------------------------------------------------
 * vf_mono_state.cc
 * Implementation of the class VFMonoState
 *
 * Author: P. Foggia
 * $Id: vf_mono_state.cpp,v 1.1 2001/10/24 01:20:36 glandrum Exp $
 *-----------------------------------------------------------------*/


/*-----------------------------------------------------------------
 * REVISION HISTORY
 *   $Log: vf_mono_state.cpp,v $
 *   Revision 1.1  2001/10/24 01:20:36  glandrum
 *   added
 *
 *   Revision 1.1  1999/01/05 10:29:16  foggia
 *   Initial revision
 *
 *
 *----------------------------------------------------------------*/

/*-----------------------------------------------------------------
 * NOTES: 
 * 1 The roles of G1 and G2 are not symmetric: the whole G1 must
 *   match a subset of G2.
 * 2 The attribute compatibility check (methods CompatibleNode
 *   and CompatibleEdge of ARGraph) is always performed
 *   applying the method to g1, and passing the attribute of
 *   g1 as first argument, and the attribute of g2 as second
 *   argument. This may be important if the compatibility
 *   criterion is not symmetric.
 -----------------------------------------------------------------*/

#include <stddef.h>

#include "vf_mono_state.h"

#include "error.h"


/*----------------------------------------------------------
 * Methods of the class VFMonoState
 ---------------------------------------------------------*/

/*----------------------------------------------------------
 * VFMonoState::VFMonoState(g1, g2)
 * Constructor. Makes an empty state.
 ---------------------------------------------------------*/
VFMonoState::VFMonoState(Graph *ag1, Graph *ag2)
  { g1=ag1;
    g2=ag2;
    n1=g1->NodeCount();
    n2=g2->NodeCount();

    core_len=0;
    t1in_len=t1out_len=0;
    t2in_len=t2out_len=0;

    core_1=new node_id[n1];
    core_2=new node_id[n2];
    node_flags_1=new byte[n1];
    node_flags_2=new byte[n2];
    if (!core_1 || !core_2 || !node_flags_1 || !node_flags_2)
      error("Out of memory");

    int i;
    for(i=0; i<n1; i++)
      { node_flags_1[i]=0;
        core_1[i]=NULL_NODE;
      }
    for(i=0; i<n2; i++)
      { node_flags_2[i]=0;
        core_2[i]=NULL_NODE;
      }
  }


/*----------------------------------------------------------
 * VFMonoState::VFMonoState(state)
 * Copy constructor. 
 ---------------------------------------------------------*/
VFMonoState::VFMonoState(const VFMonoState &state)
  { g1=state.g1;
    g2=state.g2;
    n1=state.n1;
    n2=state.n2;

    core_len=state.core_len;
    t1in_len=state.t1in_len;
    t1out_len=state.t1out_len;
    t2in_len=state.t2in_len;
    t2out_len=state.t2out_len;

    core_1=new node_id[n1];
    core_2=new node_id[n2];
    node_flags_1=new byte[n1];
    node_flags_2=new byte[n2];
    if (!core_1 || !core_2 || !node_flags_1 || !node_flags_2)
      error("Out of memory");

    int i;
    for(i=0; i<n1; i++)
      node_flags_1[i]=state.node_flags_1[i];
    for(i=0; i<n2; i++)
      node_flags_2[i]=state.node_flags_2[i];
    for(i=0; i<n1; i++)
      core_1[i]=state.core_1[i];
    for(i=0; i<n2; i++)
        core_2[i]=state.core_2[i];
  }


/*---------------------------------------------------------------
 * VFMonoState::~VFMonoState()
 * Destructor.
 --------------------------------------------------------------*/
VFMonoState::~VFMonoState() 
  { delete [] core_1;
    delete [] core_2;
    delete [] node_flags_1;
    delete [] node_flags_2;
  }


/*--------------------------------------------------------------------------
 * bool VFMonoState::NextPair(pn1, pn2, prev_n1, prev_n2)
 * Puts in *pn1, *pn2 the next pair of nodes to be tried.
 * prev_n1 and prev_n2 must be the last nodes, or NULL_NODE (default)
 * to start from the first pair.
 * Returns false if no more pairs are available.
 -------------------------------------------------------------------------*/
bool VFMonoState::NextPair(node_id *pn1, node_id *pn2,
              node_id prev_n1, node_id prev_n2)
  { byte cond1=0, cond2=0, cond3;
    if (t1out_len>0 && t2out_len>0)
      cond1=cond2=ST_TERM_OUT;
    else if (t1in_len>0 && t2in_len>0)
      cond1=cond2=ST_TERM_IN;
    else 
      cond1=~0;

    if (prev_n1==NULL_NODE)
      prev_n1=0;
    if (prev_n2==NULL_NODE)
      prev_n2=0;
    else
      prev_n2++;

    while (prev_n1<n1 &&
           (node_flags_1[prev_n1] & cond1)!=cond2)
      { prev_n1++;    
        prev_n2=0;
      }

    cond3=node_flags_1[prev_n1];

    while (prev_n2<n2 &&
               (node_flags_2[prev_n2] & (cond3|ST_CORE))!=cond3)
          prev_n2++;
    if (prev_n1<n1 && prev_n2<n2)
          { *pn1=prev_n1;
            *pn2=prev_n2;
            return true;
          }

    return false;
  }



/*---------------------------------------------------------------
 * bool VFMonoState::IsFeasiblePair(node1, node2)
 * Returns true if (node1, node2) can be added to the state
 * NOTE: 
 *   The attribute compatibility check (methods CompatibleNode
 *   and CompatibleEdge of ARGraph) is always performed
 *   applying the method to g1, and passing the attribute of
 *   g1 as first argument, and the attribute of g2 as second
 *   argument. This may be important if the compatibility
 *   criterion is not symmetric.
 --------------------------------------------------------------*/
bool VFMonoState::IsFeasiblePair(node_id node1, node_id node2)
  { assert(node1<n1);
    assert(node2<n2);
    assert((node_flags_1[node1] & ST_CORE)==0);
    assert((node_flags_2[node2] & ST_CORE)==0);

    if (!g1->CompatibleNode(g1->GetNodeAttr(node1), g2->GetNodeAttr(node2)))
      return false;

    int i, other1, other2, flags;
    void *attr1;
    int termout1=0, termout2=0, termin1=0, termin2=0, new1=0, new2=0;

    // Check the 'out' edges of node1
    for(i=0; i<g1->OutEdgeCount(node1); i++)
      { other1=g1->GetOutEdge(node1, i, &attr1);
        if ((flags=node_flags_1[other1]) & ST_CORE)
          { other2=core_1[other1];
            if (!g2->HasEdge(node2, other2) ||
                !g1->CompatibleEdge(attr1, g2->GetEdgeAttr(node2, other2)))
              return false;
          }
        else 
          { if (flags & ST_TERM_IN)
              termin1++;
            if (flags & ST_TERM_OUT)
              termout1++;
            if (!flags)
              new1++;
          }
      }

    // Check the 'in' edges of node1
    for(i=0; i<g1->InEdgeCount(node1); i++)
      { other1=g1->GetInEdge(node1, i, &attr1);
        if ((flags=node_flags_1[other1]) & ST_CORE)
          { other2=core_1[other1];
            if (!g2->HasEdge(other2, node2) ||
                !g1->CompatibleEdge(attr1, g2->GetEdgeAttr(other2, node2)))
              return false;
          }
        else 
          { if (flags & ST_TERM_IN)
              termin1++;
            if (flags & ST_TERM_OUT)
              termout1++;
            if (!flags)
              new1++;
          }
      }


    // Check the 'out' edges of node2
    for(i=0; i<g2->OutEdgeCount(node2); i++)
      { other2=g2->GetOutEdge(node2, i);
        if ((flags=node_flags_2[other2]) & ST_CORE)
          { /* do nothing */
          }
        else 
          { if (flags & ST_TERM_IN)
              termin2++;
            if (flags & ST_TERM_OUT)
              termout2++;
            if (!flags)
              new2++;
          }
      }

    // Check the 'in' edges of node2
    for(i=0; i<g2->InEdgeCount(node2); i++)
      { other2=g2->GetInEdge(node2, i);
        if ((flags=node_flags_2[other2]) & ST_CORE)
          { /* do nothing */
          }
        else 
          { if (flags & ST_TERM_IN)
              termin2++;
            if (flags & ST_TERM_OUT)
              termout2++;
            if (!flags)
		      new2++;
		  }
	      }

	    return termin1<=termin2 && termout1<=termout2;
  }



/*--------------------------------------------------------------
 * void VFMonoState::AddPair(node1, node2)
 * Adds a pair to the Core set of the state.
 * Precondition: the pair must be feasible
 -------------------------------------------------------------*/
void VFMonoState::AddPair(node_id node1, node_id node2)
  { assert(node1<n1);
    assert(node2<n2);
    assert(core_len<n1);
    assert(core_len<n2);
    int flags;

    flags=node_flags_1[node1];
    if (flags & ST_TERM_IN)
      t1in_len--;
    if (flags & ST_TERM_OUT)
      t1out_len--;

    flags=node_flags_2[node2];
    if (flags & ST_TERM_IN)
      t2in_len--;
    if (flags & ST_TERM_OUT)
      t2out_len--;

    node_flags_1[node1]=ST_CORE;
    node_flags_2[node2]=ST_CORE;

    core_1[node1]=node2;
    core_2[node2]=node1;

    core_len++;

    int i, other;
    for(i=0; i<g1->InEdgeCount(node1); i++)
      { other=g1->GetInEdge(node1, i);
        if (!(node_flags_1[other] & (ST_CORE | ST_TERM_IN)))
          { node_flags_1[other] |= ST_TERM_IN;
            t1in_len++;
          }
        if (!(node_flags_1[other] & (ST_CORE | ST_TERM_OUT)))
          { node_flags_1[other] |= ST_TERM_OUT;
            t1out_len++;
          }
      }

    for(i=0; i<g1->OutEdgeCount(node1); i++)
      { other=g1->GetOutEdge(node1, i);
        if (!(node_flags_1[other] & (ST_CORE | ST_TERM_IN)))
          { node_flags_1[other] |= ST_TERM_IN;
            t1in_len++;
          }
        if (!(node_flags_1[other] & (ST_CORE | ST_TERM_OUT)))
          { node_flags_1[other] |= ST_TERM_OUT;
            t1out_len++;
          }
      }
    
    for(i=0; i<g2->InEdgeCount(node2); i++)
      { other=g2->GetInEdge(node2, i);
        if (!(node_flags_2[other] & (ST_CORE | ST_TERM_IN)))
          { node_flags_2[other] |= ST_TERM_IN;
            t2in_len++;
          }
        if (!(node_flags_2[other] & (ST_CORE | ST_TERM_OUT)))
          { node_flags_2[other] |= ST_TERM_OUT;
            t2out_len++;
          }
      }

    for(i=0; i<g2->OutEdgeCount(node2); i++)
      { other=g2->GetOutEdge(node2, i);
        if (!(node_flags_2[other] & (ST_CORE | ST_TERM_IN)))
          { node_flags_2[other] |= ST_TERM_IN;
            t2in_len++;
          }
        if (!(node_flags_2[other] & (ST_CORE | ST_TERM_OUT)))
          { node_flags_2[other] |= ST_TERM_OUT;
            t2out_len++;
          }
      }

  }



/*--------------------------------------------------------------
 * void VFMonoState::GetCoreSet(c1, c2)
 * Reads the core set of the state into the arrays c1 and c2.
 * The i-th pair of the mapping is (c1[i], c2[i])
 --------------------------------------------------------------*/
void VFMonoState::GetCoreSet(node_id c1[], node_id c2[])
  { int i,j;
    for (i=0,j=0; i<n1; i++)
      if (node_flags_1[i] & ST_CORE)
        { c1[j]=i;
          c2[j]=core_1[i];
          j++;
        }
  }


/*----------------------------------------------------------------
 * Clones a VFMonoState, allocating with new the clone.
 --------------------------------------------------------------*/
State* VFMonoState::Clone()
  { return new VFMonoState(*this);
  }
