/*------------------------------------------------------------------
 * argedit.h
 * Definition of a simple ARG loader which allows graph edit
 * operations
 * See: argraph.h
 *
 * Author: P. Foggia
 * $Id: argedit.h,v 1.1 2001/10/24 01:20:11 glandrum Exp $
 *-----------------------------------------------------------------*/


/*-----------------------------------------------------------------
 * REVISION HISTORY
 *   $Log: argedit.h,v $
 *   Revision 1.1  2001/10/24 01:20:11  glandrum
 *   added
 *
 *   Revision 1.4  1998/12/12 12:16:49  foggia
 *   Added a new constructor
 *
 *   Revision 1.3  1998/12/08 13:30:43  foggia
 *   Minor changs
 *
 *-----------------------------------------------------------------*/

#ifndef ARGEDIT_H
#define ARGEDIT_H

#include "argraph.h"


/*---------------------------------------------------------
 * Class ARGEdit
 * A simple ARGLoader providing graph edit operations.
 * Note: the ARGEdit does not make provisions for the
 *       deallocation of the attributes, which must be
 *       dealt with by the programmer.
 -------------------------------------------------------*/
class ARGEdit: public ARGLoader
  { 
    public:

      ARGEdit();
      ARGEdit(ARGraph_impl &g);
      ARGEdit(ARGLoader &g);
      ~ARGEdit();

      /* Redefined ARGLoader methods */
      virtual int NodeCount();
      virtual void *GetNodeAttr(node_id node);
      virtual int OutEdgeCount(node_id node);
      virtual node_id GetOutEdge(node_id node, int i, void **pattr);

      /* Graph edit operations */
      node_id InsertNode(void *attr);
      void InsertEdge(node_id n1, node_id n2, void *attr);
      void DeleteNode(node_id n);
      void DeleteEdge(node_id n1, node_id n2);

    protected:
      int count;

      struct eNode
        { node_id from;
          node_id to;
          int pos;
          void *attr;
          eNode *next;
        };

      struct nNode 
        { node_id id;
          int count;
          void *attr;
          nNode *next;
          eNode *edges;
        };

      nNode *nodes;
      nNode *lastNode;
      eNode *lastEdge;

      virtual void destroyNodeAttr(void *) {};
      virtual void destroyEdgeAttr(void *) {};
  };

#endif
