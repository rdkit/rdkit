//
//  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met: 
//
//     * Redistributions of source code must retain the above copyright 
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following 
//       disclaimer in the documentation and/or other materials provided 
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
//       nor the names of its contributors may be used to endorse or promote 
//       products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
package novartis.utilities;

import java.io.PrintStream;
import java.util.Random;
import java.util.Vector;

import novartis.chemistry.molecule.Bond;
import novartis.chemistry.molecule.Layout;

/**
 * local helper class for ring perception.
 */
class TreeCell
{
   static int UNLINKED = (-1);
   static int NO_COLOR = 0;

   int color;
   int link;

   TreeCell()
   {
      color = NO_COLOR;
      link  = UNLINKED;
   }
};

/**
 * This class implements graph algorithms on arrays of edges.
 */
public class Graph
{
   int graph[][] = null;

   /**
    * create a graph from an array of bonds.
    */
   public Graph(Bond bonds[])
   {
      // get largest atom number
      int idummy = 0;
      for (int i=0; i<bonds.length; i++)
      {
         int[] pair = bonds[i].getAtomNumbers();
         if (pair[0] > idummy) idummy = pair[0];
         if (pair[1] > idummy) idummy = pair[1];
      }

      graph = new int[bonds.length][];
      for (int i=0; i<bonds.length; i++)
         if ((bonds[i].getBondFlags() & Layout.RUBBER_BOND) != 0) // rubberband bond
         {
            graph[i] = new int[2];
            idummy++; graph[i][0] = idummy;
            idummy++; graph[i][1] = idummy;
         }
         else
         {
            graph[i] = bonds[i].getAtomNumbers();
         }

   }

   /**
    * create a graph from an array of edges.
    */
   public Graph(int edges[][])
   {
      graph = new int[edges.length][];
      for (int i=0; i<edges.length; i++)
      {
         graph[i] = new int[2];
         graph[i][0] = edges[i][0];
         graph[i][1] = edges[i][1];
      }
   }

   /**
    * Helper function to get bond partner.
    */
   private final static int theOtherAtom(int pair[], int src)
   {
      if (pair[0] == src) return (pair[1]);
      else                return (pair[0]);
   }

   /**
    * Returns an array of basis rings of the graph defined by
    * this.graph[0..this.graph.length-1][0..1]
    */
   public CountedBitSet[] ringList()
   {
      int b;
      int at1, at2;
      int level1, level2;
      int trace, tmp;
      int color, old_color, new_color;
      int tmp_link, new_link;
      CountedBitSet p;
      int natoms;
      TreeCell tree[];

      Vector<CountedBitSet> result = new Vector<CountedBitSet>();

      // get largest node number to construct spanning tree.
      natoms = 0;
      for (int i=0; i<graph.length; i++)
      {
         if (natoms < graph[i][0]) natoms = graph[i][0];
         if (natoms < graph[i][1]) natoms = graph[i][1];
      }

      // allocate and initialize spanning tree
      tree = new TreeCell[natoms+1];
      for (int i=0; i<=natoms; i++)
         tree[i] = new TreeCell();

     /* Add bonds to spanning tree one at a time. If a bond doesn't link */
     /* to an old component, label the tree cells for both atoms with a  */
     /* new (unique) label. If one of the atoms corresponds to a cell    */
     /* that already has a label, label the cell corresponding to the    */
     /* other atom with the same label. If both atoms correspond to      */
     /* different labels, relabel the cells with the higher label to the */
     /* lower one to show, that they belong to the same component. If    */
     /* both cells are label the same, a ring is found. Then the links   */
     /* of the tree cells are followed to trace back to the common parent*/
     /* and the new ring, i.e. the set of bonds, is added to the result. */

      for (b=0,color=TreeCell.NO_COLOR; b<graph.length; b++)
      {
         at1 = graph[b][0]; at2 = graph[b][1];
         if (tree[at1].color == TreeCell.NO_COLOR  &&         /* new component */
             tree[at2].color == TreeCell.NO_COLOR)
         {
            color++; tree[at2].link = b;
            tree[at1].color = tree[at2].color = color;
         }
         else if (tree[at1].color == TreeCell.NO_COLOR)       /* link first atom */
         {
            tree[at1].color = tree[at2].color;
            tree[at1].link  = b;
         }
         else if (tree[at2].color == TreeCell.NO_COLOR)       /* link second atom */
         {
            tree[at2].color = tree[at1].color;
            tree[at2].link  = b;
         }
         else if (tree[at1].color != tree[at2].color)         /* link two compnts. */
         {
            new_color = tree[at1].color; old_color = tree[at2].color;
            for (int i=0; i<=natoms; i++)
               if (tree[i].color == old_color)
                  tree[i].color = new_color;

            tmp_link = tree[at2].link;     /* trace the links of component 2 */
            tree[at2].link = b;            /* and revert linkage             */
            while (tmp_link != TreeCell.UNLINKED)
            {
               at2 = theOtherAtom(graph[tmp_link],at2);
               new_link = tmp_link;
               tmp_link = tree[at2].link;
               tree[at2].link = new_link;
            }
         }
         else                             /* ring found -> add it to list */
         {                                /* trace back to root from both atoms */
            for (trace=at1,level1=0; tree[trace].link != TreeCell.UNLINKED; level1++)
               trace = theOtherAtom(graph[tree[trace].link],trace);
            for (trace=at2,level2=0; tree[trace].link != TreeCell.UNLINKED; level2++)
               trace = theOtherAtom(graph[tree[trace].link],trace);

            if (level1 > level2)   /* make path 1 the shorter one of the two */
            {
               tmp = level1; level1 = level2; level2 = tmp;
               tmp = at1; at1 = at2; at2 = tmp;
            }

            p = new CountedBitSet();
            result.addElement(p);
            p.set(b);

            for (int i=0; i<level2-level1; i++) /* trace back excess of long path */
            {
               p.set(tree[at2].link);
               at2 = theOtherAtom(graph[tree[at2].link],at2);
            }

            while (at1 != at2)     /* simultaneously trace back both paths */
            {
               p.set(tree[at1].link);
               at1 = theOtherAtom(graph[tree[at1].link],at1);
               p.set(tree[at2].link);
               at2 = theOtherAtom(graph[tree[at2].link],at2);
            }
         }       /* else ring found */
      }       /* for all bonds */

      CountedBitSet[] rings = new CountedBitSet[result.size()];
      for (int i=0; i<rings.length; i++)
         rings[i] = (CountedBitSet)result.elementAt(i);

      return(rings);
   }

   /**
    * Sorts list[] into descending order with respect to cardinality.
    * Returns the sorted array object.
    */
   public static CountedBitSet[] SortRings(CountedBitSet[] list)
   {
      // use straight insertion sort
      for (int i=1; i<list.length; i++)
         for (int j=i-1; j>=0; j--)
            if (list[j].getCardinality() < list[j+1].getCardinality()) //swap sets
            {
               CountedBitSet p = list[j];
               list[j] = list[j+1];
               list[j+1] = p;
            }
            else
               break;

      return(list);
   }

   public CountedBitSet[] combineRings(CountedBitSet[] list)
   /*
    * Combines pairs of rings until selfconsistency to get a list of
    * smaller basis rings.
    */
   {
      Random r = new Random(1);   /* make sure that algorithm work reproducibly */

      CountedBitSet p1, p2;
      int size;
      CountedBitSet set, tmp;
      boolean changed;
      int ntoggle;

      if (list == null) return(list);

      ntoggle = 0;         /* safeguard against infinite looping */
      do
      {
         changed = false;
         list = SortRings(list);  /* loop over all pairs of different rings */
         for (int i=0; i<list.length; i++)
         {
            p1 = list[i];
            for (int j=i+1; j<list.length; j++)
            {
               p2 = list[j];
               set = p1.xor(p2);
               size = set.getCardinality();
               if (size > 0 &&
                   (size <= p1.getCardinality() ||  size <= p2.getCardinality()))
               {
                  if (p1.getCardinality() > p2.getCardinality())
                  {
                     if (p1.getCardinality() > size  || ((r.nextInt() / 10) % 2) != 0)
                     {
                        if (p1.getCardinality() == size)
                           ntoggle++;
                        else
                           ntoggle = 0;
                        changed = true;
                        list[i] = set;
                        p1 = set;
                     }
                  }
                  else
                  {
                     if (p2.getCardinality() > size  || ((r.nextInt() / 10) % 2) != 0)
                     {
                        if (p2.getCardinality() == size)
                           ntoggle++;
                        else
                           ntoggle = 0;
                        changed = true;
                        list[j] = set;
                        p2 = set;
                     }
                  }
               }
            }
         }

         if (ntoggle > 4) changed = false;          /* limit to 4 toggles */
      } while(changed);

      return(list);
   }

   public void printRing(PrintStream fp, CountedBitSet ring)
   {
      int nfound = 0;
      fp.print("[" + ring.getCardinality() + "] : ");
      for (int b=0; nfound < ring.getCardinality(); b++)
         if (ring.get(b))
         {
            nfound++;
            fp.print(" " + graph[b][0] + "-" + graph[b][1]);
         }
      fp.println();
   }

   public void printRingList(PrintStream fp, CountedBitSet list[])
   /*
    * Debugging procedure to print a list of rings in a readable way.
    * graph[][] is used to produce a more meaningful output.
    */
   {
      for (int i=0; i<list.length; i++)
         printRing(fp, list[i]);
   }

   /**
    * Test bed for graph functions.
    */
   public static void main(String argv[])
   {
      int edges[][] =
      {
         {21, 22},  // five membered ring
         {22, 23},
         {23, 24},
         {24, 25},
         {25, 21},

         {25, 26},  // sprouted bond

         {1, 2},    // five membered ring
         {2, 3},
         {3, 4},
         {4, 5},
         {5, 1},

         {1, 6},    // fused seven membered ring
         {6, 7},    // bond 1-2 is shared
         {7, 8},
         {8, 9},
         {9, 10},
         {10, 2}
      };

      Graph g = new Graph(edges);
      CountedBitSet[] rings = g.ringList();
      rings = g.combineRings(rings);
      System.err.println();
      g.printRingList(System.err, rings);
   }
}
