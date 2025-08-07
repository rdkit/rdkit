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
package novartis.chemistry.molecule;

import java.util.Random;

import novartis.utilities.CountedBitSet;
import novartis.utilities.HigherMath;

/**
 * Layout of Molecule atoms to get pleasing pictures.
 */
public class Layout
{
   // flags for bonds
   public static final int    DONT_FLIP_BOND = 0x200;
   public static final int    RUBBER_BOND    = 0x100;

   // flags for atoms
   public static final int    ATOM_USED      = 0x1;
   public static final int    KEEP_POSITION  = 0x2;

   public static final double STDBOND		   = 1.514;
   public static final double STRETCH        = 1.4;   // stretch factor for collision resolution

   /**
    * Imposes a new layout onto the molecule m. It assumes that the
    * atoms and bonds of prelayouted fragments have been colored with
    * the same color.
    * The function returns a copy of the molecule with new coordinates.
    * The original molecule is not changed.
    */
   public static Molecule layoutMolecule(Molecule m)
   {
      Molecule result = (Molecule)m.clone();
      randomCoordinates(result);
      result.setupNeighbourhood();
      result.perceiveRingBonds();

   //for (i=0; i<result->n_atoms; i++)
   //   if (nbp[i].n_ligands == 1)/* terminal atoms don't have DB-stereo */
	// if (result->bond_array[nbp[i].bonds[0]].bond_type == DOUBLE)
	//    result->bond_array[nbp[i].bonds[0]].stereo_symbol = NONE;

      boolean is_ring_atom[] = new boolean[result.atoms.length];
      int     ring_size[]    = new int[result.bonds.length];
      int     ring_count[]   = new int[result.bonds.length];
      for (int i=0; i<result.atoms.length; i++) is_ring_atom[i] = false;
      for (int i=0; i<result.bonds.length; i++) ring_size[i]    = 0;
      for (int i=0; i<result.bonds.length; i++) ring_count[i]   = 0;

      for (int i=0; i<result.ringlist.length; i++)
      {
         CountedBitSet plist = result.ringlist[i];
         for (int j=0; j<result.bonds.length; j++)
         {
            Bond b = result.bonds[j];
            if (plist.get(j))
            {
	            is_ring_atom[b.atoms[0]-1] = true;
               is_ring_atom[b.atoms[1]-1] = true;
	            if (ring_size[j] == 0  ||  ring_size[j] > plist.getCardinality())
	               ring_size[j] = plist.getCardinality();
	            ring_count[j] += 1;
            }
         }
      }

      result.removePerceivedInformation();
      hideMetalComplexProblems(result, ring_count);
      result.setupNeighbourhood();
      result.perceiveRingBonds();

      layoutRings(result, ring_count);

      linkRingFragments(result);

      sproutRingSubstituents(result);

      layoutChainAtoms(result);

	   improveSPAtoms(result,0);

   //   linkRemainingFragments(result);

   // improveMoleculeByBondFlip(result, nbp, is_ring_atom, 0);

   // improveCollisions(result);
   /*
   */

   // layoutRubberFragments(result);

      layoutFragments(result, ring_size, ring_count);

      layoutAtomStereo(result);

   // layoutBondStereo(result, nbp, ring_size);

   /* makeLandscape(result); */

   // for (i=0, bp=result->bond_array; i<mp->n_bonds; i++, bp++)
   //    bp->bond_type &= 0xF;	/* upper nibble was used for internal flags */

      return (result);
   }

   /**
    * Rotates and flips product such that atoms mapped to reactant atoms
    * get as close as possible their respective partners.
    */
   public static void adjustReaction(Molecule reactant, Molecule product)
   {
      System.err.println("adjustReaction not yet implemented!");
   }

   /**
    * Sets the X and Y coordinates of the atoms of *mp to random
    * values.
    *
    * This function is not yet completely implemented as it should
    * use the atom colors to position the atoms of each color class relative
    * to each other.
    */
   public static void randomCoordinates(Molecule m)
   {
      Random r = new Random();

      for (int i=0; i<m.atoms.length; i++)
         if ((m.atoms[i].flags & KEEP_POSITION) == 0)
         {
	         m.atoms[i].x =
	            STDBOND*Math.sqrt((double)m.atoms.length)*r.nextDouble();
	         m.atoms[i].y =
	            STDBOND*Math.sqrt((double)m.atoms.length)*r.nextDouble();
	         m.atoms[i].z = 0.0;
         }
   }

   /**
    * Tries to fix the problems with layouting metal complexes by
    * making the bonds to metals RUBBER_BONDs when they are in
    * more than one ring.
    */
   private static void hideMetalComplexProblems(Molecule m, int[] ring_count)
   {
      if (m.bonds.length == 0) return;

      boolean has_metal = false;
      for (int i=0; i<m.atoms.length; i++)
         if (PTable.atomSymbolMatch(m.atoms[i].symbol, "trn"))
	         has_metal = true;
      if (!has_metal) return;

      for (int i=0; i<m.bonds.length; i++)
      {
         Bond b = m.bonds[i];
         if (ring_count[i] > 1)
	      if (PTable.atomSymbolMatch(m.atoms[b.atoms[0]-1].symbol, "trn")             ||
	          PTable.atomSymbolMatch(m.atoms[b.atoms[1]-1].symbol, "trn"))
	         b.flags |= RUBBER_BOND;
	   }
   }

   /**
    * Fills the ring node table rntp[] with the entries derived
    * from the ring's bond set ringp and the info from the molecule
    * m. It returns the number of entries in the table.
    */
   private static int fillRingNodeTable(RingNode[] 	   rnt,
   		                               CountedBitSet  ring,
   		                               Molecule m)
   {
      if (m.bonds.length <= 0) return (0);

      int bonds[][] = new int[m.bonds.length][];
      for (int i=0; i<m.bonds.length; i++)
         bonds[i] = m.bonds[i].getAtomNumbers();

      int i, j, nrnt;
      for (i=0, j=0; i<m.bonds.length; i++)	/* squeeze out irrelevant bonds */
         if (ring.get(i))
         {
	         bonds[j] = bonds[i]; j++;
         }
      nrnt = j;
      for (i=j+1; i<m.bonds.length; i++) bonds[i] = null;

      int h;
      for (i=1; i<nrnt; i++)	/* sort bonds into sequence */
         for (j=i; j<nrnt; j++)
	         if (bonds[j][0] == bonds[i-1][1])
	         {
	             h = bonds[i][0]; bonds[i][0] = bonds[j][0]; bonds[j][0] = h;
	             h = bonds[i][1]; bonds[i][1] = bonds[j][1]; bonds[j][1] = h;
	             break;
	         }
	         else if (bonds[j][1] == bonds[i-1][1])
	         {
	            h = bonds[j][0]; bonds[j][0] = bonds[j][1]; bonds[j][1] = h;
	            h = bonds[i][0]; bonds[i][0] = bonds[j][0]; bonds[j][0] = h;
	            h = bonds[i][1]; bonds[i][1] = bonds[j][1]; bonds[j][1] = h;
	            break;
	         }

      for (i=1; i<nrnt; i++)
         if (m.atoms[bonds[i-1][0]-1].color != m.atoms[bonds[i][0]-1].color) break;
      if (i == nrnt)	/* ring has single color -> already layouted */
         return (0);

      int from, to;
      from = i-1;

      for (i=nrnt-1; i>from; i--)  /* find color sequence of segment */
      {				                 /* that might wrap */
         if (m.atoms[bonds[i][0]-1].color != m.atoms[bonds[from][0]-1].color) break;
      }
      to = i;

      rnt[0] = new RingNode();
      rnt[0].aifirst = bonds[to][1];   rnt[0].ailast  = bonds[from][0];
      rnt[0].color   = m.atoms[bonds[from][0]-1].color;
      nrnt = 1;
      for (i=from+1; i<=to; i++)
         if (m.atoms[bonds[i][0]-1].color != m.atoms[bonds[i-1][0]-1].color)
         {
	         for (j=i; j<to; j++)
	            if (m.atoms[bonds[j][0]-1].color != m.atoms[bonds[j][1]-1].color) break;
            rnt[nrnt] = new RingNode();
	         rnt[nrnt].aifirst = bonds[i][0]; rnt[nrnt].ailast  = bonds[j][0];
	         rnt[nrnt].color   = m.atoms[bonds[i][0]-1].color; nrnt++;
         }

      for (i=0; i<nrnt; i++)
      {
         rnt[i].xfirst = m.atoms[rnt[i].aifirst-1].x;
         rnt[i].yfirst = m.atoms[rnt[i].aifirst-1].y;
         rnt[i].xlast  = m.atoms[rnt[i].ailast-1].x;
         rnt[i].ylast  = m.atoms[rnt[i].ailast-1].y;
      }

      return (nrnt);
   }

   /**
    * Orders the rings in ring_list[] such that six-membered ones
    * come first follwed by the others by increasing size.
    */
   private static void layoutRingSort(CountedBitSet ring_list[])
   {
      for (int i=1; i<ring_list.length; i++)
         for (int j=i-1; j>=0; j--)
         {
            if (ring_list[j].getCardinality() != 6  && ring_list[j+1].getCardinality() == 6)
            {  // prefer six-membered rings
               CountedBitSet bs = ring_list[j]; ring_list[j] = ring_list[j+1]; ring_list[j+1] = bs;
            }
            else if (ring_list[j].getCardinality() != 6  && ring_list[j].getCardinality() > ring_list[j+1].getCardinality())
            {  // prefer small rings over larger ones
               CountedBitSet bs = ring_list[j]; ring_list[j] = ring_list[j+1]; ring_list[j+1] = bs;
            }
            else
               break;
         }
      return;
   }

   /**
    * Places the atoms of *mp that have colors listed in rntp[0..nrnt-1]
    * relative to each other. It uses the perimeter of a circle to for
    * coordinates of the atoms listed as aifirst and ailast in rntp[].
    */
   private static void layoutRingSegment(Molecule m,
		                                   RingNode rnt[],
		                                   int			 nrnt)
   {
      int i, j, n;
      double xc, yc;
      double x1, y1, x2, y2;
      double peri, len;
      double r, angle;
      double p1[]  = new double[2];
      double p2[]  = new double[2];
      double ps[]  = new double[2];
      double p1p[] = new double[2];
      double p2p[] = new double[2];
      double points[][];
      int npts;

      xc = yc = 0.0; n = 0;		/* compute center of circle */
      for (i=1; i<nrnt; i++)
         for (j=0; j<m.atoms.length; j++)
	         if (m.atoms[j].color == rnt[i].color)
	         {
	            xc += m.atoms[j].x;
	            yc += m.atoms[j].y;
	            n++;
	         }
      xc /= n; yc /= n;

      len = 0.0;				/* compute radius of circle */
      for (i=0; i<nrnt; i++)
      {
         len += STDBOND;					/* bond secante */
         rnt[i].d = Math.sqrt(HigherMath.sqr(rnt[i].xfirst-rnt[i].xlast) +
                              HigherMath.sqr(rnt[i].yfirst-rnt[i].ylast));
         len += rnt[i].d;                                  /* fragment secante */
      }
      peri = 0.0;
      for (i=0; i<nrnt; i++)
      {
         peri += STDBOND*(Math.PI*STDBOND/len)/Math.sin(Math.PI*STDBOND/len);
         rnt[i].angle = 2*Math.PI*rnt[i].d/len;
         if (rnt[i].angle > 0.0001)
	         peri += rnt[i].d*(rnt[i].angle/2)/Math.sin(rnt[i].angle/2);
      }
      r = peri/(2*Math.PI);

      points = new double[m.atoms.length][];
      for (i=0; i<m.atoms.length; i++) points[i] = new double[2];

      angle = 0;
      for (i=0; i<nrnt; i++)	/* new coordinates */
      {
         x1 = r*Math.sin(angle)+xc; y1 = r*Math.cos(angle)+yc;
         angle += rnt[i].angle;
         x2 = r*Math.sin(angle)+xc; y2 = r*Math.cos(angle)+yc;
         angle += 2*Math.PI*STDBOND/len;

         p1[0] = rnt[i].xfirst; p1[1] = rnt[i].yfirst;
         p2[0] = rnt[i].xlast;  p2[1] = rnt[i].ylast;
         p1p[0] = x1; p1p[1] = y1; p2p[0] = x2; p2p[1] = y2;
         for (j=0, npts=0; j<m.atoms.length; j++)
         {
	         Atom a = m.atoms[j];
	         if (a.color == rnt[i].color)
	         {
	            points[npts][0] = a.x; points[npts][1] =a.y;
	            npts++;
	         }
	      }

         if (HigherMath.sqr(p1[0]-p2[0])+HigherMath.sqr(p1[1]-p2[1]) > 0.01*HigherMath.sqr(STDBOND))
         {							/* ring fusion */
	         ps[0] = ps[1] = 0.0;	/* for center of gravity of fragment */
	         for (j=0, npts=0; j<m.atoms.length; j++)
	         {
	            Atom a = m.atoms[j];
	            if (a.color == rnt[i].color)
	            {
	               ps[0] += a.x; ps[1] += a.y; npts++;
	            }
	         }
	         ps[0] /= npts; ps[1] /= npts;
	         /* use mirror image of fragment if needed */
	         if ((ps[0]-p1[0])*(p2[1]-p1[1]) - (ps[1]-p1[1])*(p2[0]-p1[0]) > 0)
	         {
	            p1[1] *= (-1); p2[1] *= (-1);
	            for (j=0; j<npts; j++) points[j][1] *= (-1);
	            m.flipStereoSymbols(rnt[i].color);
	         }
         }
         else				/* single atom or spiro */
         /* if (FALSE) */
         {
	         ps[0] = ps[1] = 0.0;	/* compute center of gravity of atoms */
	         npts = 0;		/* attached to link atom */
	         for (j=0; j<m.bonds.length; j++)
	         {
	            Bond b = m.bonds[j];
	            if ((b.atoms[0] == rnt[i].aifirst  &&
	                 m.atoms[b.atoms[1]-1].color == rnt[i].color) ||
	                (b.atoms[1] == rnt[i].aifirst  &&
	                 m.atoms[b.atoms[0]-1].color == rnt[i].color) ||
	                (b.atoms[0] == rnt[i].ailast  &&
	                 m.atoms[b.atoms[1]-1].color == rnt[i].color) ||
	                (b.atoms[1] == rnt[i].ailast  &&
	                 m.atoms[b.atoms[0]-1].color == rnt[i].color))
	            {
	               ps[0] += m.atoms[b.atoms[0]-1].x;
	               ps[1] += m.atoms[b.atoms[0]-1].y;
	               npts++;
	               ps[0] += m.atoms[b.atoms[1]-1].x;
	               ps[1] += m.atoms[b.atoms[1]-1].y;
	               npts++;
	            }
	         }
	         if (npts > 0)
	         {
	            ps[0] /= npts; ps[1] /= npts;
	            p2[0] = ps[0]; p2[1] = ps[1];
	            p2p[0] = p1p[0] +
		             (p1p[0]-xc)*Math.sqrt(HigherMath.sqr(ps[0]-p1[0])+HigherMath.sqr(ps[1]-p1[1]))/
			        	 Math.sqrt(HigherMath.sqr(p1p[0]-xc)+HigherMath.sqr(p1p[1]-yc));
	            p2p[1] = p1p[1] +
		             (p1p[1]-yc)*Math.sqrt(HigherMath.sqr(ps[0]-p1[0])+HigherMath.sqr(ps[1]-p1[1]))/
			        	 Math.sqrt(HigherMath.sqr(p1p[0]-xc)+HigherMath.sqr(p1p[1]-yc));
	         }
         }

         HigherMath.transformPoints(points, m.atoms.length, p1, p2, p1p, p2p);
         for (j=0, npts=0; j<m.atoms.length; j++)
         {
	         Atom a = m.atoms[j];
	         if (a.color == rnt[i].color)
	         {
	            a.x = points[npts][0]; a.y = points[npts][1];
	            npts++;
	         }
	      }
      }
  /*
   System.err.println("peri = " + peri + ", len = " + len + "\n");
   */
   }

   /**
    * Merges the atom colors listed in rntp[0..nrnt-1] into one.
    */
   private static void mergeColors(Molecule m,
		                 RingNode rnt[],
		                 int		  nrnt)
   {
      int i, j;

      for (i=1; i<nrnt; i++)
         for (j=0; j<m.atoms.length; j++)
	    if (m.atoms[j].color == rnt[i].color)
	       m.atoms[j].color = rnt[0].color;
   }

   /**
    * Set coordinates for the atoms within ring systems. It
    * Recolors atoms and bonds which have been placed relative
    * to each other with the same color.
    */
   private static void layoutRings(Molecule m, int[] ring_count)
   {
      RingNode rnt[];
      int nrnt;

      layoutRingSort(m.ringlist);
      for (int i=0; i<m.ringlist.length;)
      {
         CountedBitSet ring = m.ringlist[i];
         rnt = new RingNode[m.atoms.length];
         nrnt = fillRingNodeTable(rnt, ring, m);
         if (nrnt == 0) i++;
         else
         {
	         layoutRingSegment(m, rnt, nrnt);
	         mergeColors(m, rnt, nrnt);
         }
      }
   }

   /**
    * Computes the signed volume of the tetrahedron spanned by the four 3D points
    * tetra[0][0..2] to tetra[3][0..2].
    */
   private static double volume(double[][] tetra)
   {
      double ax, ay, az, bx, by, bz, cx, cy, cz;

      //System.err.println("v1 = " + tetra[0][0] + "," + tetra[0][1] + "," + tetra[0][2]);
      //System.err.println("v2 = " + tetra[1][0] + "," + tetra[1][1] + "," + tetra[1][2]);
      //System.err.println("v3 = " + tetra[2][0] + "," + tetra[2][1] + "," + tetra[2][2]);
      //System.err.println("v4 = " + tetra[3][0] + "," + tetra[3][1] + "," + tetra[3][2]);
      ax = tetra[1][0] - tetra[0][0];
      ay = tetra[1][1] - tetra[0][1];
      az = tetra[1][2] - tetra[0][2];
      bx = tetra[2][0] - tetra[0][0];
      by = tetra[2][1] - tetra[0][1];
      bz = tetra[2][2] - tetra[0][2];
      cx = tetra[3][0] - tetra[0][0];
      cy = tetra[3][1] - tetra[0][1];
      cz = tetra[3][2] - tetra[0][2];

      return (ax*(by*cz-bz*cy) +
              ay*(bz*cx-bx*cz) +
              az*(bx*cy-by*cx));
   }

   // Constants tuning the selection of bonds to bear a stereo symbol.
   static final double  CHAIN_ATOM_Z      = 0.1;    // substituent is a single atom
   static final double  RING_ATOM_Z       = 0.01;   // substituent is a ring atom
   static final double  PREFER_WEDGE      = 1.01;   // filled wedge drawings are preferred
   static final double  STEREO_NEIGHBOUR  = 1.01;   // panelty for stereobond to stereo neighbour atom
   static final int     IS_STEREO         = 1;      // atom color for stereo depictor atoms

   /**
    * Converts the internal atom stereo information into bond symbol depictions.
    *
    * Chain bonds are prefered of ring bonds and bonds to non-stereocenters are
    * preferred over bonds to stereocenters.
    */
   private static void layoutAtomStereo(Molecule m)
   {
      int[] old_colors = m.removeAtomColors();

      if (m.stereodescs == null) return;

      for (int i=0; i<m.stereodescs.length; i++)
      {
         StereoDescription sd = m.stereodescs[i];
         if (sd.depictor instanceof Atom)
          ((Atom)sd.depictor).color |= IS_STEREO;
      }

      for (int i=0; i<m.stereodescs.length; i++)
      {
         StereoDescription sd = m.stereodescs[i];
         if (sd.stereoclass == StereoDescription.TH)
         {
            // collect relevant pieces of the molecule
            Atom a = (Atom)sd.depictor;
            double[][] tetra = new double[4][];
            Bond[] ligands = new Bond[4];
            int it = 0;
            for (int k=0; k<sd.atoms.length; k++)
            {
               if (sd.atoms[k] == StereoDescription.IMPLICIT_H)
               {
                  tetra[it] = new double[3];
                  tetra[it][0] = a.x;
                  tetra[it][1] = a.y;
                  tetra[it][2] = 0;
                  ligands[it] = null;
                  it++;
               }
               else
               {
                  int j;
                  for (j=0; j<a.neighbour_atoms.length; j++)
                  {
                     if (a.neighbour_atoms[j].index == sd.atoms[k]-1) break;
                  }
                  tetra[it] = new double[3];
                  tetra[it][0] = a.neighbour_atoms[j].x;
                  tetra[it][1] = a.neighbour_atoms[j].y;
                  tetra[it][2] = 0;
                  ligands[it] = a.neighbour_bonds[j];
                  it++;
               }
            }

            // find best assignment of wedges
            double best_volume = 0; int best_code = -1;
            for (int code = 1; code < 3*3*3*3; code++)
            {
               int tmp = code;
               int nwedge = 0;
               for (int j=0; j<4; j++)
               {
                  if (sd.atoms[j] == StereoDescription.IMPLICIT_H  ||
                      ligands[j].stereo_symbol != Bond.NORMAL)
                  {
                     tmp = tmp/3;
                     continue;
                  }
                  if (tmp%3 != 0) nwedge++;
                  switch (tmp%3)
                  {
                     case 0: tetra[j][2] = 0; break;
                     case 1: tetra[j][2] = PREFER_WEDGE; break;
                     case 2: tetra[j][2] = -1; break;
                  }
                  if (ligands[j].topography == Bond.RING)
                     tetra[j][2] *= RING_ATOM_Z;
                  else
                     tetra[j][2] *= CHAIN_ATOM_Z;
                  if ((m.atoms[sd.atoms[j]-1].color & IS_STEREO) != 0)
                     tetra[j][2] *= STEREO_NEIGHBOUR;
                  // prefer wedges to small substituents
                  tetra[j][2] /= m.atoms[sd.atoms[j]-1].neighbour_atoms.length;
                  tmp = tmp/3;
               }
               if (volume(tetra) > best_volume*nwedge)
               {
                  best_volume = volume(tetra)/nwedge;
                  best_code = code;
               }
            }
            // assign wedges
            for (int j=0; j<4; j++)
            {
               if (best_code%3 == 0)
               {
                  best_code = best_code/3;
                  continue;
               }
               if (best_code%3 == 1)
                  ligands[j].stereo_symbol = Bond.UP;
               else
                  ligands[j].stereo_symbol = Bond.DOWN;
               best_code = best_code/3;
               if (a.index != ligands[j].atoms[0]-1)
               {
                  int tmp;
                  tmp = ligands[j].atoms[0];
                  ligands[j].atoms[0] = ligands[j].atoms[1];
                  ligands[j].atoms[1] = tmp;
               }
            }
         }
      }
      m.restoreAtomColors(old_colors);
   }

   /**
    * Links directly connected ring fragments.
    */
   private static void linkRingFragments(Molecule m)
   {
      int[][] edges;
      double[][] coords;
      double[] p1, p2, p1p, p2p;
      int col1, col2;
      int n1, n2, h;
      int seed;
      int i, j;

      if (m.bonds == null  ||  m.bonds.length < 1) return;

      int[] numbers = new int[m.atoms.length];

      for (j=0; j<m.bonds.length; j++)
      {
         Bond b = m.bonds[j];
         if (m.atoms[b.atoms[0]-1].color  ==  m.atoms[b.atoms[1]-1].color    ||
	          (b.flags & RUBBER_BOND)  !=  0                                  ||
	          m.atoms[b.atoms[0]-1].topography != Atom.RING                   ||
	          m.atoms[b.atoms[1]-1].topography != Atom.RING)
	      continue;
	      /* b is now pointing to a bond linking two ring fragments */

         col1 = m.atoms[b.atoms[0]-1].color;
         col2 = m.atoms[b.atoms[1]-1].color;

         for (i=0, n1=0; i<m.atoms.length; i++)    /* find reference fragment */
	         if (m.atoms[i].color == col1) n1++;
         for (i=0, n2=0; i<m.atoms.length; i++)
	         if (m.atoms[i].color == col2) n2++;
         if (n1 > n2)
         {
	          h = b.atoms[0]; b.atoms[0] = b.atoms[1]; b.atoms[1] = h;
	          h = col1; col1 = col2; col2 = h;
	          h = n1; n1 = n2; n2 = h;
         }

         edges  = m.getColoredEdges(col2, numbers);
         coords = m.getColoredCoordinates(col2, numbers);
         seed = numbers[b.atoms[1]-1];
         p2p = new double[2];
         p2p[0] = coords[seed][0];
         p2p[1] = coords[seed][1];
         p1p = HigherMath.nextSubstituentPoint(coords, edges, seed, false);

         edges  = m.getColoredEdges(col1, numbers);
         coords = m.getColoredCoordinates(col1, numbers);
         seed = numbers[b.atoms[0]-1];
         p1 = new double[2];
         p1[0] = coords[seed][0];
         p1[1] = coords[seed][1];
         p2 = HigherMath.nextSubstituentPoint(coords, edges, seed, false);

         HigherMath.transformPoints(coords, coords.length, p1, p2, p1p, p2p);

         for (i=0; i<m.atoms.length; i++)
	         if (m.atoms[i].color == col1)
	         {
	            m.atoms[i].x = coords[numbers[i]][0];
	            m.atoms[i].y = coords[numbers[i]][1];
	         }

         improveBondByFlip(m, b, 0.01);
         improveBondByStretch(m, b);

      				/* merge colors */
         for (i=0; i<m.atoms.length; i++)
	         if (m.atoms[i].color == col2) m.atoms[i].color = col1;
      }
   }


   /**
    * Computes a strain measure between the two parts of the molecule m
    * having the colors col1 and col2.
    */
   private static double colorStrain(Molecule m,
		                               int      col1,
		                               int		 col2)
   {
      int i1, i2;
      double result;

      result = 0;
      for (i1=0; i1<m.atoms.length; i1++)
      {
         Atom a1 = m.atoms[i1];
         if (a1.color == col1)
	      for (i2=0; i2<m.atoms.length; i2++)
	      {
            Atom a2 = m.atoms[i2];
	         if (a2.color == col2)
	         result += 1/(0.01 + (a1.x - a2.x)* (a1.x - a2.x)
			                     + (a1.y - a2.y)* (a1.y - a2.y));
			}
      }

      return (result);
   }

   /**
    * Tries to improve atom collisions by flipping one side of bond bp.
    * The function returnes TRUE if i sucessfully improved the coordinates.
    * Changes are made only if strain improves by at least threshold.
    */
   public static boolean improveBondByFlip(Molecule m,
                                           Bond b,
		                                     double threshold)
   {
      double[] p1  = new double[2], p2 = new double[2];
      double[] p   = new double[2], pp = new double[2];
      double[] r12 = new double[2];
      double strain1, strain2, q;
      int i;

      if ((b.flags & DONT_FLIP_BOND) != 0) return (false);

      strain1 = colorStrain(m, m.atoms[b.atoms[0]-1].color,
   			                   m.atoms[b.atoms[1]-1].color);

      p1[0] = m.atoms[b.atoms[0]-1].x;
      p1[1] = m.atoms[b.atoms[0]-1].y;
      p2[0] = m.atoms[b.atoms[1]-1].x;
      p2[1] = m.atoms[b.atoms[1]-1].y;
      r12[0] = p2[0] - p1[0]; r12[1] = p2[1]-p1[1];

      for (i=0; i<m.atoms.length; i++)
      {
         Atom a = m.atoms[i];
         if (a.color == m.atoms[b.atoms[1]-1].color)
         {
	         p[0] = a.x; p[1] = a.y;
	         q = ((p[0]-p1[0])*r12[0] + (p[1]-p1[1])*r12[1])/
	             (r12[0]*r12[0]       + r12[1]*r12[1]);
	         pp[0] = p1[0] - (p[0]-p1[0]) + 2*r12[0]*q;
	         pp[1] = p1[1] - (p[1]-p1[1]) + 2*r12[1]*q;
	         a.x = pp[0]; a.y = pp[1];
         }
      }

      strain2 = colorStrain(m, m.atoms[b.atoms[0]-1].color,
	                            m.atoms[b.atoms[1]-1].color);

      if (strain1 <= strain2+threshold)		/* transform back */
      {
         for (i=0; i<m.atoms.length; i++)
         {
            Atom a = m.atoms[i];
	         if (a.color == m.atoms[b.atoms[1]-1].color)
	         {
	            p[0] = a.x; p[1] = a.y;
	            q = ((p[0]-p1[0])*r12[0] + (p[1]-p1[1])*r12[1])/
		             (r12[0]*r12[0]       + r12[1]*r12[1]);
	            pp[0] = p1[0] - (p[0]-p1[0]) + 2*r12[0]*q;
	            pp[1] = p1[1] - (p[1]-p1[1]) + 2*r12[1]*q;
	            a.x = pp[0]; a.y = pp[1];
	         }
	      }
         return (false);
      }

      return (true);
   }

   /**
    * Tries to improve atom collisions by stretching the bond by 1/3rd.
    */
   public static void improveBondByStretch(Molecule m, Bond b)
   {
      double r12[] = new double[2];
      int i, i1, i2;
      boolean collision;
      int col1, col2;

      col1 = m.atoms[b.atoms[0]-1].color;
      col2 = m.atoms[b.atoms[1]-1].color;

      collision = false;
      for (i1=0; i1<m.atoms.length; i1++)
      {
         Atom a1 = m.atoms[i1];
         if (a1.color == col1)
	         for (i2=0; i2<m.atoms.length; i2++)
	         {
	            Atom a2 = m.atoms[i2];
	            if (a2.color == col2  &&
	                (a1.x - a2.x) * (a1.x - a2.x) +
		             (a1.y - a2.y) * (a1.y - a2.y) <
		             0.05*STDBOND*STDBOND)
	            {
	               collision = true; break;
	            }
	         }
	   }

      if (!collision) return;

      r12[0] = m.atoms[b.atoms[1]-1].x -
               m.atoms[b.atoms[0]-1].x;
      r12[1] = m.atoms[b.atoms[1]-1].y -
               m.atoms[b.atoms[0]-1].y;

      for (i=0; i<m.atoms.length; i++)
      {
         Atom a = m.atoms[i];
         if (a.color == m.atoms[b.atoms[1]-1].color)
         {
	         a.x += 0.2*r12[0]+0.2*r12[1]; a.y += 0.2*r12[1]-0.2*r12[0];
         }
      }
   }

   /**
    * Attaches the first atom of non-ring substituents to already layouted
    * ring fragments in a symmetrical way. This is necessary because
    * step-by-step addition does not produce sufficient results.
    * Only those attachment points are considered which have exclusively
    * non-ring attachments.
    */
   private static void sproutRingSubstituents(Molecule m)
   {
      if (m.bonds == null  ||
          m.bonds.length < 1) return;	/* trivial molecule */

      int[][] edges;
      double[][] coords;
      double[] p1 = new double[2], p2 = new double[2];
      double[] d = new double[2], dp = new double[2], c = new double[2];
      int seed;
      int i, j, k, iatom;
      int[] numbers = new int[m.atoms.length];

      double sina, cosa;
      Atom subst_atoms[] = new Atom[m.atoms.length];  // just to be sure
      int n_subst;
      int n_rbonds, min_rsize;
      boolean problem_substituent;

      /* loop over all ring atoms */
      for (j=0; j<m.atoms.length; j++)
      {
         Atom a = m.atoms[j];
         if (a.topography != Atom.RING) continue;
         n_subst = 0; n_rbonds = 0;
         problem_substituent = false;
         min_rsize = Integer.MAX_VALUE;
         for (i=0; i<a.neighbour_bonds.length; i++)
	         if (0 == (a.neighbour_bonds[i].flags & RUBBER_BOND))
	            if (a.neighbour_bonds[i].smallestRing() != 0)
	            {
	               if (a.neighbour_bonds[i].smallestRing() < min_rsize)
	                  min_rsize = a.neighbour_bonds[i].smallestRing();
	               n_rbonds++;
	            }
	            else
	            {
	               if (a.neighbour_atoms[i].topography == Atom.RING)
		              problem_substituent = true;
	               else
	               {
		               if (m.countColor(a.neighbour_atoms[i].color) != 1)
		                 problem_substituent = true;
		               subst_atoms[n_subst] = a.neighbour_atoms[i];
		               n_subst++;
	               }
	            }
         if (n_subst == 0) continue;	/* nothing to be layed out */
         if (n_rbonds != 2) continue;	/* ring fusion atom */
         if (problem_substituent) continue;/* substituent is part of other */
      					/* ring or already layed out */

         edges  = m.getColoredEdges(a.color, numbers);
         coords = m.getColoredCoordinates(a.color, numbers);
         seed = numbers[j];
         p1[0] = coords[seed][0];
         p1[1] = coords[seed][1];
         if (min_rsize <= 8)
            p2 = HigherMath.nextSubstituentPoint(coords, edges, seed, false);
         else
            p2 = HigherMath.nextSubstituentPoint(coords, edges, seed, true);
         d[0] = p2[0]-p1[0];
         d[1] = p2[1]-p1[1];
         if (n_subst > 1)
         {
   	      for (i=1; i<n_subst; i++)	       /* sort substituents by size */
   	         for (k=i-1; k>=0; k--)
   	            if (subst_atoms[k].neighbour_atoms.length <
   		             subst_atoms[k+1].neighbour_atoms.length)
   	            {
   		            Atom h = subst_atoms[k];
   		            subst_atoms[k] = subst_atoms[k+1];
   		            subst_atoms[k+1] = h;
   	            }
   	            else
   		           break;
   	      cosa = Math.cos(((0.5)*Math.PI/n_subst-Math.PI/2)*80.0/90.0);
   	      sina = Math.sin(((0.5)*Math.PI/n_subst-Math.PI/2)*80.0/90.0);
   	      dp[0] = cosa*d[0] - sina*d[1] + p1[0];
   	      dp[1] = sina*d[0] + cosa*d[1] + p1[1];
   	      c[0] = c[1] = 0.0;
   	      for (i=0; i<coords.length; i++)	/* compute center of ring system */
   	      {
   	         c[0] += coords[i][0]; c[1] += coords[i][1];
   	      }
   	      c[0] /= coords.length; c[1] /= coords.length;
   	      /* if first substituent is more distant from center than others */
   	      if ((dp[0]-c[0])*(dp[0]-c[0]) + (dp[1]-c[1])*(dp[1]-c[1]) <
   	          (p2[0]-c[0])*(p2[0]-c[0]) + (p2[1]-c[1])*(p2[1]-c[1]))
            {	/* swap extreme substituents */
   	         Atom h = subst_atoms[0];
   	         subst_atoms[0] = subst_atoms[n_subst-1];
   	         subst_atoms[n_subst-1] = h;
   	      }
         }
         for (i=0; i<n_subst; i++)
         {			/* scale for 80 deg. with two substituents */
   	      cosa = Math.cos(((i+0.5)*Math.PI/n_subst-Math.PI/2)*80.0/90.0);
   	      sina = Math.sin(((i+0.5)*Math.PI/n_subst-Math.PI/2)*80.0/90.0);
   	      dp[0] = cosa*d[0] - sina*d[1] + p1[0];
   	      dp[1] = sina*d[0] + cosa*d[1] + p1[1];
   	      subst_atoms[i].x = dp[0];
   	      subst_atoms[i].y = dp[1];
   	      subst_atoms[i].color = a.color;
         }
      }
   }

   /**
    * Computes a measure of the branchedness of atom with index ai
    * or atom number (ai+1). The more smaller substituents the better.
    * The quality function strongly prefers atoms connected to already
    * layouted atoms as indicated by used_atoms[].
    */
   private static int branchQuality(Molecule m, Atom a)
   {
      int nbranch[] = new int[a.neighbour_atoms.length];
      int j, k;
      int nsingle, nmore, ntotal;
      boolean can_be_candidate;
      Atom ah;
      // neighbourhood_t *nbph;
      boolean placed_neighbour;
      int already_layouted;

      already_layouted = 0;
      can_be_candidate = false; placed_neighbour = false;
      for (j=0; j<a.neighbour_atoms.length; j++)
      {
         nbranch[j] = 0;
         if ((a.neighbour_bonds[j].flags & RUBBER_BOND) != 0)
	         continue;
         if (a.color != a.neighbour_atoms[j].color)
            can_be_candidate = true;
         else
	         already_layouted++;
         if ((a.neighbour_atoms[j].flags & ATOM_USED) != 0) placed_neighbour = true;
         ah = a.neighbour_atoms[j];
         for (k=0; k<ah.neighbour_bonds.length; k++)
	         if (0 == (ah.neighbour_bonds[k].flags & RUBBER_BOND))
	         {
	            nbranch[j]++;
	         }
      }

      if (!can_be_candidate  ||  already_layouted > 1) return (-1);

      nsingle = nmore = ntotal = 0;
      for (j=0; j<a.neighbour_atoms.length; j++)
      {
         if (nbranch[j] == 1) nsingle++;
         else if (nbranch[j] > 0) nmore++;
         ntotal += nbranch[j];
      }
      if (placed_neighbour)
         return (1000+nsingle*100+nmore*10+ntotal);
      else
         return (nsingle*100+nmore*10+ntotal);
   }

   /**
    * Layouts the environment of those atoms which have only chain
    * substituents.
    */
   private static void layoutChainAtoms(Molecule m)
   {
      if (m.bonds == null  ||  m.bonds.length == 0) return;

      int i, j;
      int quality, ibest;

      int[][] edges;
      double[][] coords;
      double[] p1 = new double[2], p2 = new double[2];
      double[] p1p = new double[2], p2p = new double[2];
      int seed;
      int[] numbers = new int[m.atoms.length];

      for (i=0; i<m.atoms.length; i++) m.atoms[i].flags &= ~ATOM_USED;
      for (;;)	/* get next best chain atom until all have been used */
      {
         quality = (-1); ibest = (-1);
         for (i=0; i<m.atoms.length; i++)
         {
            Atom a = m.atoms[i];
            if (a.topography != Atom.RING)
	         {
	            if (quality < branchQuality(m, a))
	            {
	               ibest = i; quality = branchQuality(m, a);
	            }
	         }
	      }
         if (quality < 0) break;
         m.atoms[ibest].flags |= ATOM_USED;

         Atom abest = m.atoms[ibest];
         int nneigh = 0;
         int col[] = new int[abest.neighbour_atoms.length];
         int oldcolor;
         int size[] = new int[abest.neighbour_atoms.length];
         Atom atoms[] = new Atom[abest.neighbour_atoms.length];
         Bond bonds[] = new Bond[abest.neighbour_atoms.length];
         oldcolor = abest.color;
         abest.color = (-1);
         for (i=0; i<abest.neighbour_bonds.length; i++)
	         if ((abest.neighbour_bonds[i].flags & RUBBER_BOND) == 0)
	         {
	            atoms[nneigh] = abest.neighbour_atoms[i];
	            bonds[nneigh] = abest.neighbour_bonds[i];
	            col[nneigh]  = atoms[nneigh].color;
	            size[nneigh] = 100*m.countColor(col[nneigh]) +
	                           atoms[nneigh].neighbour_atoms.length;
	            nneigh++;
	         }
         for (i=1; i<nneigh; i++)	/* already layed out branch is reference */
	         if (col[i] == oldcolor)
	         {
	            Atom ah = atoms[i]; atoms[i] = atoms[0]; atoms[0] = ah;
	            Bond bh = bonds[i]; bonds[i] = bonds[0]; bonds[0] = bh;
	            int h = col[i]; col[i] = col[0]; col[0] = h;
	            h = size[i]; size[i] = size[0]; size[0] = h;
	         }
         for (;;)			/* optimize for space usage */
         {
	         for (i=1; i<nneigh-1; i++)
	            if (Math.abs(size[i]-size[i-1]) + Math.abs(size[i+1]-size[(i+2)%nneigh]) <
		             Math.abs(size[i]-size[(i+2)%nneigh]) + Math.abs(size[i+1]-size[i-1]))
	            {
	               Atom ah = atoms[i]; atoms[i] = atoms[i+1]; atoms[i+1] = ah;
	               Bond bh = bonds[i]; bonds[i] = bonds[i+1]; bonds[i+1] = bh;
	               int h = col[i]; col[i] = col[i+1]; col[i+1] = h;
	               h = size[i]; size[i] = size[i+1]; size[i+1] = h;
	               break;
	            }
	         if (i >= nneigh-1) break;
         }

         /* layout reference branch */
         edges  = m.getColoredEdges(col[0], numbers);
         coords = m.getColoredCoordinates(col[0], numbers);
         p1p[0] = Math.cos((2*Math.PI*0)/nneigh)*STDBOND;
         p1p[1] = Math.sin((2*Math.PI*0)/nneigh)*STDBOND;
         p2p[0] = 0.0; p2p[1] = 0.0;
         p1[0] = atoms[0].x;
         p1[1] = atoms[0].y;
         if (col[0] == oldcolor)	/* real reference branch */
         {
      	   p2[0] = m.atoms[ibest].x;
      	   p2[1] = m.atoms[ibest].y;
         }
         else			/* no reference branch */
         {
	         // search index of seed atom
	         seed = numbers[0];
	         for (int k=0; k<m.atoms.length; k++)
	            if (m.atoms[k] == atoms[0])
	            {
	               seed = numbers[k];
	               break;
	            }
	         p2 = HigherMath.nextSubstituentPoint(coords, edges, seed, true);
         }
         HigherMath.transformPoints(coords, coords.length, p1, p2, p1p, p2p);
         for (j=0; j<m.atoms.length; j++)
	         if (m.atoms[j].color == col[0])
	         {
	            m.atoms[j].x = coords[numbers[j]][0];
	            m.atoms[j].y = coords[numbers[j]][1];
	         }
         m.atoms[ibest].x = 0;
         m.atoms[ibest].y = 0;

         /* layout other substituents */
         for (i=1; i<nneigh; i++)
         {
            edges  = m.getColoredEdges(col[i], numbers);
            coords = m.getColoredCoordinates(col[i], numbers);
	         p1p[0] = Math.cos((2*Math.PI*i)/Math.max(nneigh,3))*STDBOND;
	         p1p[1] = Math.sin((2*Math.PI*i)/Math.max(nneigh,3))*STDBOND;
	         p2p[0] = 0.0; p2p[1] = 0.0;
	         p1[0] = atoms[i].x;
	         p1[1] = atoms[i].y;
	         // search index of seed atom
	         seed = numbers[0];
	         for (int k=0; k<m.atoms.length; k++)
	            if (m.atoms[k] == atoms[i])
	            {
	               seed = numbers[k];
	               break;
	            }
	         p2 = HigherMath.nextSubstituentPoint(coords, edges, seed, true);
            HigherMath.transformPoints(coords, coords.length, p1, p2, p1p, p2p);
	         for (j=0; j<m.atoms.length; j++)
	            if (m.atoms[j].color == col[i])
	            {
	               m.atoms[j].x = coords[numbers[j]][0];
	               m.atoms[j].y = coords[numbers[j]][1];
	            }
         }
      				/* merge colors */
         for (i=0; i<nneigh; i++)
	         for (j=0; j<m.atoms.length; j++)
	         {
	            Atom a = m.atoms[j];
	            if (a.color == col[i]) a.color = oldcolor;
	         }
         m.atoms[ibest].color = oldcolor;

         // for (i=0; i<nneigh; i++)
	         // if (size[i] > 200)	/* fragment was already layed out */
	            // if (atoms[i].topography != Atom.RING)
	               // makeBondTrans(m, bonds[i]);
      }

      for (i=0; i<m.bonds.length; i++)
      {
         Bond b = m.bonds[i];
         if (b.topography == Bond.CHAIN)
         {
            makeBondTrans(m, b);
         }
      }

      // restore flags
      for (i=0; i<m.atoms.length; i++) m.atoms[i].flags &= ~ATOM_USED;
   }

   /**
    * Makes sp carbon and nitrogen atoms linear.
    * If atom is != 0, only this atom is considered.
    */
   private static void improveSPAtoms(Molecule m, int atom)
   {
      // struct reaccs_atom_t *ap, *aph;
      // struct reaccs_bond_t *bp;
      boolean changed;
      int[][] edges;
      double[][] coords;

      double[] p1 = new double[2], p2 = new double[2], p1p = new double[2], p2p = new double[2];
      int[] numbers;
      int[] colors;

      colors = new int[m.atoms.length];
      for (int i=0; i<m.atoms.length; i++)
         colors[i] = m.atoms[i].color;

      for (int i=0; i<m.atoms.length; i++)
      {
         Atom a = m.atoms[i];
         if (a.neighbour_atoms.length == 2             &&
             (a.symbol.equals("C") ||
	           a.symbol.equals("N"))  &&
	          (a.neighbour_bonds[0].bond_type == Bond.SINGLE  &&
	           a.neighbour_bonds[1].bond_type == Bond.TRIPLE    ||
	           a.neighbour_bonds[1].bond_type == Bond.SINGLE  &&
	           a.neighbour_bonds[0].bond_type == Bond.TRIPLE    ||
	           a.neighbour_bonds[0].bond_type == Bond.DOUBLE  &&
	           a.neighbour_bonds[1].bond_type == Bond.DOUBLE)  &&
	          (atom == 0  ||  atom == i+1))
         {
	         for (int j=0; j<m.atoms.length; j++)
	            m.atoms[j].color = 0;

	         a.color = 1;
	         a.neighbour_atoms[0].color = 2;
	         a.neighbour_atoms[1].color = 3;

	         do
	         {
	            changed = false;
	            for (int j=0; j<m.bonds.length; j++)
	            {
	               Bond b = m.bonds[j];
	               if (m.atoms[b.atoms[0]-1].color > 0  &&
	                   m.atoms[b.atoms[1]-1].color == 0)
	               {
		               changed = true;
		               m.atoms[b.atoms[1]-1].color =
		                  m.atoms[b.atoms[0]-1].color;
	               }
	               else if (m.atoms[b.atoms[1]-1].color > 0  &&
	                        m.atoms[b.atoms[0]-1].color == 0)
	               {
		               changed = true;
		               m.atoms[b.atoms[0]-1].color =
		                  m.atoms[b.atoms[1]-1].color;
	               }
	            }
	         } while (changed);

	         if (a.neighbour_atoms[0].color !=	/* no ring */
	             a.neighbour_atoms[1].color)
            {
	            p1[0] = a.x;
	            p1[1] = a.y;
	            p2[0] = a.neighbour_atoms[0].x;
	            p2[1] = a.neighbour_atoms[0].y;
	            p1p[0] = a.x;
	            p1p[1] = a.y;
	            p2p[0] = 2*p1p[0]-a.neighbour_atoms[1].x;
	            p2p[1] = 2*p1p[1]-a.neighbour_atoms[1].y;

	            numbers = new int[m.atoms.length];
               edges  = m.getColoredEdges(2, numbers);
               coords = m.getColoredCoordinates(2, numbers);
               HigherMath.transformPoints(coords, coords.length, p1, p2, p1p, p2p);

	            for (int j=0; j<m.atoms.length; j++)
	               if (m.atoms[j].color == 2)
	               {
		               m.atoms[j].x = coords[numbers[j]][0];
		               m.atoms[j].y = coords[numbers[j]][1];
	               }
	         }
         }
      }

      for (int i=0; i<m.atoms.length; i++)
         m.atoms[i].color = colors[i] ;
   }

   /**
    * Flips the bond b if necessary to make it look 'trans'.
    */
   public static void makeBondTrans(Molecule m, Bond b)
   {
      int sizea[], na, aia;
      int sizeb[], nb, aib;
      int[] atom_colors = m.removeAtomColors();

      double[] p1 = new double[2], p2 = new double[2];
      double[] p = new double[2], r12 = new double[2];
      double strain1, strain2, q;

   			/* collect ligands and sizes */
      Atom aa = m.atoms[b.atoms[0]-1];
      Atom ab = m.atoms[b.atoms[1]-1];
      aa.color = 1; ab.color = 2;
      na = 0; sizea = new int [aa.neighbour_atoms.length];
      for (int i=0; i<aa.neighbour_atoms.length; i++)
      {
         sizea[na] = aa.neighbour_atoms[i].neighbour_atoms.length;
         na++;
      }
      nb = 0; sizeb = new int [ab.neighbour_atoms.length];
      for (int i=0; i<ab.neighbour_atoms.length; i++)
      {
         sizeb[nb] = ab.neighbour_atoms[i].neighbour_atoms.length;
         nb++;
      }

      strain1 = 0;
      for (int i=0; i<na; i++)
         for (int j=0; j<nb; j++)
	         strain1 += sizea[i]*sizeb[j]/(1 + HigherMath.distSquare(aa.neighbour_atoms[i].x,
	                                                                 aa.neighbour_atoms[i].y,
	                                                                 ab.neighbour_atoms[j].x,
	                                                                 ab.neighbour_atoms[j].y));

      p1[0] = aa.x;	/* compute bond vector */
      p1[1] = aa.y;
      p2[0] = ab.x;
      p2[1] = ab.y;
      r12[0] = p2[0] - p1[0]; r12[1] = p2[1]-p1[1];
						/* transform ligands a */
      for (int i=0; i<na; i++)
      {
         Atom a = aa.neighbour_atoms[i];
         p[0] = a.x; p[1] = a.y;
         q = ((p[0]-p1[0])*r12[0] + (p[1]-p1[1])*r12[1])/
	           (r12[0]*r12[0]       + r12[1]*r12[1]);
         a.x = p1[0] - (p[0]-p1[0]) + 2*r12[0]*q;
         a.y = p1[1] - (p[1]-p1[1]) + 2*r12[1]*q;
      }

      strain2 = 0;
      for (int i=0; i<na; i++)
         for (int j=0; j<nb; j++)
	         strain2 += sizea[i]*sizeb[j]/(1 + HigherMath.distSquare(aa.neighbour_atoms[i].x,
	                                                                 aa.neighbour_atoms[i].y,
	                                                                 ab.neighbour_atoms[j].x,
	                                                                 ab.neighbour_atoms[j].y));
						/* transform back */
      for (int i=0; i<na; i++)
      {
         Atom a = aa.neighbour_atoms[i];
         p[0] = a.x; p[1] = a.y;
         q = ((p[0]-p1[0])*r12[0] + (p[1]-p1[1])*r12[1])/
	           (r12[0]*r12[0]       + r12[1]*r12[1]);
         a.x = p1[0] - (p[0]-p1[0]) + 2*r12[0]*q;
         a.y = p1[1] - (p[1]-p1[1]) + 2*r12[1]*q;
      }

      if (strain1 > strain2)	/* cis -> trans */
      {
         for (int i=0; i<na; i++)
	         m.floodColor(aa, 3);
         for (int i=0; i<m.atoms.length; i++)
         {
	         Atom a = m.atoms[i];
	         if (a.color == 3)
	         {
	            p[0] = a.x; p[1] = a.y;
	            q = ((p[0]-p1[0])*r12[0] + (p[1]-p1[1])*r12[1])/
		              (r12[0]*r12[0]       + r12[1]*r12[1]);
	            a.x = p1[0] - (p[0]-p1[0]) + 2*r12[0]*q;
	            a.y = p1[1] - (p[1]-p1[1]) + 2*r12[1]*q;
	         }
	      }
      }

      m.restoreAtomColors(atom_colors);
   }

/**
 * Links fragments of m not yet layed-out. This can happen because of
 * previous layout rules being too strict to cover all cases.
 */
   private static void linkRemainingFragments(Molecule m)
   {
      if (m.bonds == null  ||  m.bonds.length < 1) return;

      int[] numbers;
      int[][] edges;
      double[][] coords;
      double[] p1 = new double[2],  p2 = new double[2];
      double[] p1p = new double[2], p2p = new double[2];
      int col1, col2;
      int n1, n2, h;
      int seed;
      int i, j;


      for (j=0; j<m.bonds.length; j++)
      {
         Bond b = m.bonds[j];
         if (m.atoms[b.atoms[0]-1].color  ==  m.atoms[b.atoms[1]-1].color  ||
	          (b.flags & RUBBER_BOND) != 0)
	        continue;
	      // b is now pointing to a bond linking to fragments.

         col1 = m.atoms[b.atoms[0]-1].color;
         col2 = m.atoms[b.atoms[1]-1].color;

         for (i=0, n1=0; i<m.atoms.length; i++)    /* find reference fragment */
   	      if (m.atoms[i].color == col1)
   	         n1++;
         for (i=0, n2=0; i<m.atoms.length; i++)
   	      if (m.atoms[i].color == col2)
   	         n2++;
         if (n1 > n2)
         {
	         h = b.atoms[0]; b.atoms[0] = b.atoms[1]; b.atoms[1] = h;
	         h = col1; col1 = col2; col2 = h;
	         h = n1; n1 = n2; n2 = h;
         }

         numbers = new int[m.atoms.length];
         edges  = m.getColoredEdges(col2, numbers);
         coords = m.getColoredCoordinates(col2, numbers);
         seed = numbers[b.atoms[1]-1];
         p2p[0] = coords[seed][0];
         p2p[1] = coords[seed][1];
	      p1p = HigherMath.nextSubstituentPoint(coords, edges, seed, true);

         edges  = m.getColoredEdges(col1, numbers);
         coords = m.getColoredCoordinates(col1, numbers);
         seed = numbers[b.atoms[0]-1];
         p1[0] = coords[seed][0];
         p1[1] = coords[seed][1];
	      p2 = HigherMath.nextSubstituentPoint(coords, edges, seed, true);

         HigherMath.transformPoints(coords, coords.length, p1, p2, p1p, p2p);

         for (i=0; i<m.atoms.length; i++)
	         if (m.atoms[i].color == col1)
	         {
	            m.atoms[i].x = coords[numbers[i]][0];
	            m.atoms[i].y = coords[numbers[i]][1];
	         }

      				/* merge colors */
         for (i=0; i<m.atoms.length; i++)
	         if (m.atoms[i].color == col2) m.atoms[i].color = col1;
      }
   }

   /**
    * Places the fragments of the molecule with respect to
    * each other. A part of the molecule with the same color
    * is considered a fragment even if it consists of more than
    * one disconnected component.
    */
   private static void layoutFragments(Molecule m, int ring_size[], int ring_count[])
   {
      FragmentDesc flist, fp, fph, ftmp;
      int i;
      double xoffset;

   			/* compute fragment descriptions */
      flist = null;
      for (i=0; i<m.atoms.length; i++)
      {
         Atom a = m.atoms[i];
         for (fp=flist; fp != null; fp=fp.next)
	         if (fp.color == a.color) break;
         if (fp != null)	/* fragment found */
         {
	         fp.natoms++;
         }
         else	/* new fragment */
         {
	         fp = new FragmentDesc(1.0e7, 1.0e7, -1.0e7, -1.0e7, a.color, 1);
	         fp.next = flist; flist = fp;
         }
      }

	   /* bubble sort by size of fragment, assuming only natoms and color do matter yet */
      for (fph = flist; fph != null  &&  fph.next != null; fph=fph.next)
         for (fp = flist; fp != null  && fp.next != null; fp=fp.next)
	         if (fp.natoms < fp.next.natoms)
	         {
	            int h;
	            h = fp.next.natoms; fp.next.natoms = fp.natoms; fp.natoms = h;
	            h = fp.next.color;  fp.next.color = fp.color;   fp.color = h;
	         }

   				/* orient fragments */
      // for (fp=flist; fp; fp=fp->next)
      //    if (fp->natoms > 1)
	   //       OrientFragment(mp, fp->color, ring_size, ring_count);

				/* get bounding boxes */
      for (fp=flist; fp != null; fp=fp.next)
         for (i=0; i<m.atoms.length; i++)
         {
	         Atom a = m.atoms[i];
	         if (a.color == fp.color)
	         {
	            if (a.x < fp.xll) fp.xll = a.x;
	            if (a.y < fp.yll) fp.yll = a.y;
	            if (a.x > fp.xur) fp.xur = a.x;
	            if (a.y > fp.yur) fp.yur = a.y;
	         }
	      }

	      /* frame fragment bounding boxes by 0.5 standard bond length */
      for (fp=flist; fp != null; fp=fp.next)
      {
         fp.xll -= STDBOND/2; fp.yll -= STDBOND/2;
         fp.xur += STDBOND/2; fp.yur += STDBOND/2;
         fp.xoffset = -fp.xll; fp.yoffset = -fp.yll;
      }

      // arrange fragments in X direction
      // a future version could treat single atoms separately
      xoffset = 0;
      for (fp=flist; fp != null; fp=fp.next)
      {
         fp.xoffset += xoffset;
         xoffset += fp.xur - fp.xll;
      }

      for (i=0; i<m.atoms.length; i++)
      {
         Atom a = m.atoms[i];
         for (fp=flist; fp != null; fp=fp.next)
	      if (fp.color == a.color) break;

         if (fp != null)
         {
	         a.x += fp.xoffset; a.y += fp.yoffset;
         }
      }

   }
}

/**
 * private helper class to contain information on fragments
 * that are being laid out.
 */
class FragmentDesc
   {
      double xll;
      double yll;
      double xur;
      double yur;
      double xoffset = 0 ;
      double yoffset = 0 ;
      int color;
      int natoms;
      FragmentDesc next = null;        // pointer to next fragment

      FragmentDesc(double xll, double yll, double xur, double yur, int color, int natoms)
      {
         this.xll = xll;
         this.yll = yll;
         this.xur = xur;
         this.yur = yur;
         this.color = color;
         this.natoms = natoms;
      }
   };

/**
 * private helper class to contain information on ring segments
 * that are being laid out.
 */
class RingNode
{
   double xfirst;
   double yfirst;
   double xlast;
   double ylast;
   int aifirst;
   int ailast;
   int color;	   /* color of fragment 			 */
   /* scratch fields */
   double d;		   /* distance of aifirst and ailast */
   double angle;	   /* angle on ring perimeter */
}
