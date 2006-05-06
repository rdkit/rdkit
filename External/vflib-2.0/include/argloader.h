/*------------------------------------------------------------------
 * argloader.h
 * Interface of argloader.cc
 * Definition of a simple ARG loader based on iostream using text files,
 * and of a binary file unattributed Graph loader.
 * See: argraph.h
 *
 * Author: P. Foggia
 *-----------------------------------------------------------------*/



/*-----------------------------------------------------------------
 * DESCRIPTION OF THE TEXT FILE FORMAT 
 * On the first line there must be the number of nodes;
 * subsequent lines will contain the node attributes, one node per 
 * line, preceded by the node id; node ids must be in the range from
 * 0 to the number of nodes - 1.
 * Then, for each node there is the number of edges coming out of 
 * the node, followed by a line for each edge containing the 
 * ids of the edge ends and the edge attribute.
 * Blank lines, and lines starting with #, are ignored.
 * An example file, where both node and edge attributes are ints, 
 * could be the following:
     # Number of nodes
     3
     # Node attributes
     0 27
     1 42
     2 13

     # Edges coming out of node 0
     2
     0 1  24
     0 2  73

     # Edges coming out of node 1
     1
     1 3  66

     # Edges coming out of node 2
     0
  *-----------------------------------------------------------------*/


/*---------------------------------------------------------------------------
 *   DESCRIPTION OF THE BINARY FILE FORMAT
 * The file is composed by a sequence of 16-bit words; the words are
 * encoded in little-endian format (e.g., LSB first).
 * The first word represents the number of nodes in the graph.
 * Then, for each node, there is a word encoding the number of
 * edges coming out of that node, followed by a sequence of words
 * encoding the endpoints of those edges.
 * An example, represented in hexadecimal, follows:
 *     03 00     Number of nodes (3)
 *     00 00     Number of edges out of node 0 (0)
 *     02 00     Number of edges out of node 1 (2)
 *     00 00     Target of the first edge of node 1 (edge 1 -> 0)
 *     02 00     Target of the second edge of node 1 (edge 1 -> 2)
 *     01 00     Number of edges out of node 2 (1)
 *     00 00     Target of the first (and only) edge of node 2 (edge 2 -> 0)
 ---------------------------------------------------------------------------*/

#ifndef ARGLOADER_H
#define ARGLOADER_H

#include <iostream>
#include <sstream>
//#include <strstream>
#include <ctype.h>


#include "argedit.h"
#include "allocpool.h"

//using namespace std;

template <class Node, class Edge>
class StreamARGLoader: public ARGEdit
  { 
    public:
      typedef Allocator<Node> NodeAllocator;
      typedef Allocator<Edge> EdgeAllocator;

      StreamARGLoader(NodeAllocator *nalloc, EdgeAllocator *ealloc, 
                      std::istream &in);

      static void write(std::ostream &out, ARGraph<Node, Edge> &g);
      static void write(std::ostream &out, ARGLoader &g);

    private:

      enum { MAX_LINE=512 };
      void readLine(std::istream &in, char *line);
      int  readCount(std::istream &in);
      void readNode(NodeAllocator *alloc, std::istream &in);
      void readEdge(EdgeAllocator *alloc, std::istream &in);
  };


class BinaryGraphLoader: public ARGEdit
  { public:
      BinaryGraphLoader(std::istream &in);
      static void write(std::ostream &out, Graph &g);
      static void write(std::ostream &out, ARGLoader &g);

    private:
      static unsigned readWord(std::istream &in);
      static void writeWord(std::ostream &out, unsigned w);
  };


/*------------------------------------------------------------
 * Methods of the class StreamArgLoader
 -----------------------------------------------------------*/

/*----------------------------------------------------------
 * Constructor
 ---------------------------------------------------------*/
template <class Node, class Edge>
StreamARGLoader<Node, Edge>::
StreamARGLoader(Allocator<Node> *nalloc, 
                Allocator<Edge> *ealloc, 
                std::istream &in)
  { 
    int cnt=readCount(in);
    if (cnt<=0)
      { cnt=0;
        return;
      }

    int i;
    for(i=0; i<cnt; i++)
      { readNode(nalloc, in);
      }

    for(i=0; i<cnt; i++)
      { int ecount, j;
        ecount=readCount(in);
        for(j=0; j<ecount; j++)
          readEdge(ealloc, in);
      }
        
  }

/*------------------------------------------------------
 * Reads a line from the input stream
 ----------------------------------------------------*/
template <class Node, class Edge>
void StreamARGLoader<Node, Edge>::
readLine(std::istream &in, char *line)
  { 
    char *p;
    do {
      *line='\0';
      if (!in.good())
        error("End of file or reading error");
      in.getline(line, MAX_LINE);
      for(p=line; isspace(*p); p++)
        ;
    } while (*p=='\0' || *p=='#');
  }

/*------------------------------------------------------
 * Reads an int from a line
 ----------------------------------------------------*/
template <class Node, class Edge>
int StreamARGLoader<Node, Edge>::
readCount(std::istream &in)
  { char line[MAX_LINE+1];
    readLine(in, line);
    
    int i;
    std::istringstream is(line);
    is>>i;

    return i;
  }

/*------------------------------------------------------
 * Reads a node from a line
 ----------------------------------------------------*/
template <class Node, class Edge>
void StreamARGLoader<Node, Edge>::
readNode(Allocator<Node> *alloc, std::istream &in)
  { char line[MAX_LINE+1];
    readLine(in, line);
    std::istringstream is(line);
    
    Node *nattr=alloc->Allocate();
    node_id id;

    is >> id >> *nattr;

    if (id != NodeCount())
      error("File format error\n  Line: %s", line);

    InsertNode(nattr);
  }

/*------------------------------------------------------
 * Reads an edge from a line
 ----------------------------------------------------*/
template <class Node, class Edge>
void StreamARGLoader<Node, Edge>::
readEdge(Allocator<Edge> *alloc, std::istream &in)
  { char line[MAX_LINE+1];
    readLine(in, line);
    std::istringstream is(line);
    
    Edge *eattr=alloc->Allocate();
    node_id id1, id2;

    is >> id1 >> id2 >> *eattr;

    InsertEdge(id1, id2, eattr);
  }


/*-----------------------------------------------------------
 * Writes an ARGraph on a stream in a format 
 * readable by StreamARGLoader.
 * Relies on stream output operators for the
 * Node and Edge types
 ----------------------------------------------------------*/
template <class Node, class Edge>
void StreamARGLoader<Node,Edge>::
write(std::ostream &out, ARGraph<Node, Edge> &g)
  { out << g.NodeCount() << std::endl;

    int i;
    for(i=0; i<g.NodeCount(); i++)
      out << i << ' ' << *g.GetNodeAttr(i) << std::endl;

    int j;
    for(i=0; i<g.NodeCount(); i++)
      { out << g.OutEdgeCount(i) << std::endl;
        for(j=0; j<g.OutEdgeCount(i); j++)
          { int k;
            Edge *attr;
            k=g.GetOutEdge(i, j, &attr);
            out << i << ' ' << k << ' ' << *attr << std::endl;
          }
      }
  }

/*-----------------------------------------------------------
 * Writes an ARGLoader on a stream in a format 
 * readable by StreamARGLoader.
 * Relies on stream output operators for the
 * Node and Edge types
 ----------------------------------------------------------*/
template <class Node, class Edge>
void StreamARGLoader<Node,Edge>::
write(std::ostream &out, ARGLoader &g)
  { out << g.NodeCount() << std::endl;

    int i;
    for(i=0; i<g.NodeCount(); i++)
      out << i << ' ' << *(Node *)g.GetNodeAttr(i) << std::endl;

    int j;
    for(i=0; i<g.NodeCount(); i++)
      { out << g.OutEdgeCount(i) << std::endl;
        for(j=0; j<g.OutEdgeCount(i); j++)
          { int k;
            void *attr;
            k=g.GetOutEdge(i, j, &attr);
            out << i << ' ' << k << ' ' << *(Edge *)attr << std::endl;
          }
      }
  }

#endif
