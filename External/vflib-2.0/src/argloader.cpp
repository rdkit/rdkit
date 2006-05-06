/*------------------------------------------------------------------
 * argloader.cc
 * Implementation of argloader.h
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
 -----------------------------------------------------------------------------*/


#include "argloader.h"



/*------------------------------------------------------------------
 * Create a BinaryGraphLoader reading from a binary istream.
 * NOTE: the input stream must be open with the 
 * ios::binary | ios::in mode.
 -----------------------------------------------------------------*/
BinaryGraphLoader::BinaryGraphLoader(std::istream &in)
  { unsigned n, ne, dest;
    unsigned i,j;

    n = readWord(in);
    for(i=0; i<n; i++)
      InsertNode(NULL);
    for(i=0; i<n; i++)
      { ne = readWord(in);
        for(j=0; j<ne; j++)
	  { dest=readWord(in);
	    InsertEdge(i, dest, NULL);
	  }
       }
   }

 /*--------------------------------------------------------------
  * Save a graph on a file readable by a BinaryGraphLoader.
  * NOTE: the output stream must be open with the 
  * ios::binary | ios::out mode.
  -------------------------------------------------------------*/
void BinaryGraphLoader::write(std::ostream& out, Graph &g)
  { int i,j;
    writeWord(out, g.NodeCount());
    for(i=0; i<g.NodeCount(); i++)
      { writeWord(out, g.OutEdgeCount(i));
        for(j=0; j<g.OutEdgeCount(i); j++)
	  writeWord(out, g.GetOutEdge(i,j,NULL));
      }
  }

 /*-----------------------------------------------------------------
  * Private functions
  ----------------------------------------------------------------*/
unsigned BinaryGraphLoader::readWord(std::istream& in)
{ //unsigned char c1, c2;
  char c1, c2;
    in.get(c1);
    in.get(c2);
    return c1 | (c2 << 8);
  }

void BinaryGraphLoader::writeWord(std::ostream& out, unsigned w)
  { unsigned char c1, c2;

    c1 = w & 0xFF; 
    c2 = w >> 8;
    out << c1 << c2;
  }
