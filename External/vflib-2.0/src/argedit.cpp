/*-----------------------------------------------------------
 * argedit.cc
 * Implementation of the class ARGEdit, a simple ARG loader
 * which allows graph edit operations.
 * See: argedit.h, argraph.h
 *
 * Author: P. Foggia
 * $Id: argedit.cpp,v 1.1 2001/10/24 01:20:36 glandrum Exp $
 -----------------------------------------------------------*/

/*-----------------------------------------------------------
 * REVISION HISTORY
 *   $Log: argedit.cpp,v $
 *   Revision 1.1  2001/10/24 01:20:36  glandrum
 *   added
 *
 *   Revision 1.4  1998/12/12 12:17:45  foggia
 *   Added a new constructor
 *
 *   Revision 1.3  1998/12/08 13:31:19  foggia
 *   Minor changes
 *
 ----------------------------------------------------------*/

#include "argedit.h"
#include "error.h"

/*----------------------------------------------------------
 * Constructor
 ---------------------------------------------------------*/
ARGEdit:: ARGEdit() 
  { count=0;
    nodes=NULL;
    lastNode=NULL;
    lastEdge=NULL;
  }

/*----------------------------------------------------------
 * Constructor
 ---------------------------------------------------------*/
ARGEdit:: ARGEdit(ARGraph_impl &g) 
  { count=0;
    nodes=NULL;
    lastNode=NULL;
    lastEdge=NULL;

    node_id n;
    for(n=0; n<g.NodeCount(); n++)
      InsertNode(g.GetNodeAttr(n));
    int i;
    for(n=0; n<g.NodeCount(); n++)
      { for(i=0; i<g.OutEdgeCount(n); i++)
          { void *attr; 
            node_id n2=g.GetOutEdge(n, i, &attr);
            InsertEdge(n, n2, attr);
          }
      }
  }

/*----------------------------------------------------------
 * Constructor
 ---------------------------------------------------------*/
ARGEdit:: ARGEdit(ARGLoader &g) 
  { count=0;
    nodes=NULL;
    lastNode=NULL;
    lastEdge=NULL;

    node_id n;
    for(n=0; n<g.NodeCount(); n++)
      InsertNode(g.GetNodeAttr(n));
    int i;
    for(n=0; n<g.NodeCount(); n++)
      { for(i=0; i<g.OutEdgeCount(n); i++)
          { void *attr; 
            node_id n2=g.GetOutEdge(n, i, &attr);
            InsertEdge(n, n2, attr);
          }
      }
  }


/*----------------------------------------------------------
 * Destructor
 ---------------------------------------------------------*/
ARGEdit::~ARGEdit()
  { nNode *pn, *pn1;
    eNode *pe, *pe1;

    pn=nodes;
    while (pn!=NULL)
      { pe=pn->edges;
        while (pe!=NULL)
          { pe1=pe;
            pe=pe->next;
            destroyEdgeAttr(pe1->attr);
            delete pe1;
          }
        pn1=pn;
        pn=pn->next;
        destroyNodeAttr(pn1->attr);
        delete pn1;
      }
  }


/*------------------------------------------------------
 * Returns the number of nodes
 ----------------------------------------------------*/
int ARGEdit::NodeCount()
  { return count;
  }

/*------------------------------------------------------
 * Returns the attr of a node.
 ----------------------------------------------------*/
void* ARGEdit::GetNodeAttr(node_id id)
  { nNode *n=nodes;

    if (lastNode!=NULL && lastNode->id <= id)
      n=lastNode;
    
    while (n!=NULL && n->id!=id)
      n=n->next;
    if (n==NULL)
      error("Inconsistent data");
    return n->attr;
  }

/*------------------------------------------------------
 * Returns the number of edges coming out of a node.
 ----------------------------------------------------*/
int ARGEdit::OutEdgeCount(node_id id)
  { nNode *n=nodes;

    if (lastNode!=NULL && lastNode->id <= id)
      n=lastNode;
    
    while (n!=NULL && n->id!=id)
      n=n->next;
    if (n==NULL)
      error("Inconsistent data");
    return n->count;
  }

/*------------------------------------------------------
 * Returns an edge
 ----------------------------------------------------*/
node_id ARGEdit::
GetOutEdge(node_id id, int i, void **pattr)
  { nNode *n=nodes;

    if (lastNode!=NULL && lastNode->id <= id)
      n=lastNode;
    
    while (n!=NULL && n->id!=id)
      n=n->next;
    if (n==NULL)
      error("Inconsistent data");

    eNode *e=n->edges;
    int pos=0;

    if (lastEdge!=NULL && lastEdge->from==id && 
        lastEdge->pos>=0 && lastEdge->pos<=i)
      { e=lastEdge;
        pos=e->pos;
      }

    while (e!=NULL && pos<i)
      { e->pos=pos;
        e=e->next;
        pos++;
      }

    if (e==NULL)
      error("Inconsistent data");

    if (pattr!=NULL)
      *pattr = e->attr;
    lastNode = n;
    lastEdge = e;

    return e->to;
  }




/*-------------------------------------------
 * Creates a new node
 ------------------------------------------*/
node_id ARGEdit::InsertNode(void* attr)
  { nNode *n=new nNode;
    if (n==NULL)
      error("Out of memory");
    node_id id=n->id=count++;
    n->attr=attr;
    n->edges=NULL;
    n->count=0;

    nNode *p=nodes, *p0=NULL;
    if (lastNode!=NULL && lastNode->id<id)
      { p0=lastNode;
        p=lastNode->next;
      }
    while(p!=NULL && p->id<id) 
      { p0=p; 
        p=p->next;
      }


    if (p0==NULL)
      { n->next=nodes;
        nodes=n;
      }
    else
      { n->next=p0->next;
        p0->next=n;
      }
    lastNode=n;

     return n->id;
  }



/*-------------------------------------------
 * Creates a new edge
 ------------------------------------------*/
void ARGEdit::InsertEdge(node_id id1, node_id id2, void* attr)
  {
    eNode *e=new eNode;
    if (e==NULL)
      error("Out of memory");
    e->from=id1;
    e->to=id2;
    e->attr=attr;
    e->pos=-1;
    nNode *pn=nodes;
    if (lastNode!=NULL && lastNode->id<=id1)
      { pn=lastNode;
      }
    while(pn!=NULL && pn->id<id1) 
      { pn=pn->next;
      }
    if (pn==NULL || pn->id!=id1)
      error("Bad param 1 in ARGEdit::InsertEdge: %d", (int)id1);

    eNode *p=pn->edges, *p0=NULL;
    if (lastEdge!=NULL && lastEdge->from==id1 && lastEdge->to<id2)
      { p0=lastEdge;
        p=lastEdge->next;
      }
    while (p!=NULL && p->to<id2)
      { p0=p;
        p=p->next;
      }
    if (p!=NULL && p->to==id2)
      error("Bad param 2 in ARGEdit::InsertEdge: %d", (int)id2);

    if (p0==NULL)
      { e->next=pn->edges;
        pn->edges=e;
      }
    else
      { e->next=p0->next;
        p0->next=e;
      }
    pn->count++;
    lastNode=pn;
    lastEdge=e;
}



/*---------------------------------------------------
 * Delete a node
 --------------------------------------------------*/
void ARGEdit::DeleteNode(node_id id)
  { nNode *p=nodes, *p0=NULL;

    if (lastNode!=NULL && lastNode->id<id)
      { p0=lastNode;
        p=lastNode->next;
      }
    while(p!=NULL && p->id<id) 
      { p0=p; 
        p=p->next;
      }

    if (p==NULL || !p->id==id)
      error("Bad param in ARGEdit::DeleteNode");

    while (p->edges != NULL)
      { eNode *ep=p->edges;
        p->edges=ep->next;
        destroyEdgeAttr(ep->attr);
        delete ep;
      }
    lastEdge=NULL;
    destroyNodeAttr(p->attr);
    if (p0==NULL)
      { nodes=p->next;
      }
    else
      { p0->next=p->next;
      }
    delete p;
    count--;

    p=nodes;
    while (p!=NULL)
      { if (p->id > id)
          p->id --;
        eNode *ep0=NULL;
        eNode *ep=p->edges;
        while (ep!=NULL)
          { if (ep->to == id)
              { if (ep0==NULL)
                  p->edges=ep->next;
                else
                  ep0->next=ep->next;
                destroyEdgeAttr(ep->attr);
                delete ep;
		p->count--;
                if (ep0==NULL)
                  ep=p->edges;
                else
                  ep=ep0->next;
                continue;
              }
            if (ep->from > id)
              ep->from --;
            if (ep->to > id)
              ep->to --;
            ep0=ep;
            ep=ep->next;
          }
        p=p->next;
      }
    

    lastNode=p0;
  }

/*------------------------------------
 * Delete an edge
 -----------------------------------*/
void ARGEdit::DeleteEdge(node_id id1, node_id id2)
  { nNode *pn=nodes;

    if (lastNode!=NULL && lastNode->id<=id1)
      pn=lastNode;

    while (pn!=NULL && pn->id < id1)
      pn=pn->next;

    if (pn==NULL || pn->id != id1)
      error("Bad param in ARGEdit::DeleteEdge");

    eNode *pe=pn->edges, *pe0=NULL;

    if (lastEdge!=NULL &&  lastEdge->from==id1 && lastEdge->to < id2)
      { pe0=lastEdge;
        pe=pe0->next;
      }

    while (pe!=NULL && pe->to<id2)
      { pe0=pe;
        pe=pe->next;
      }
    
    if (pe==NULL || pe->to != id2)
      error("Bad param in ARGEdit::DeleteEdge");

    if (pe0==NULL)
      pn->edges=pe->next;
    else
      pe0->next=pe->next;
    destroyEdgeAttr(pe->attr);
    delete pe;
    pn->count--;

    lastNode=pn;
    lastEdge=pe0;
  }
