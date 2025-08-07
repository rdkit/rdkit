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

import java.io.IOException;
import java.util.Vector;

import novartis.utilities.CountedBitSet;
import novartis.utilities.FortranInputStream;
import novartis.utilities.Graph;
import novartis.utilities.HigherMath;

/**
 * Class representing a chemical structure model of a molecule.
 *
 * Basically, it consists of an array of atoms, and array of bonds
 * linking those atoms with a defined bond order, an array of
 * stereo properties to describe stereochemistry, and an array
 * of atom set properties which define higher order objects like
 * references to molecule fragments.
 */

public class Molecule
{
   boolean             being_modified = false;
   Atom[]              atoms          = null;
   Bond[]              bonds          = null;
   StereoDescription[] stereodescs    = null;
   SetProperty[]       setprops       = null;

   public int getNAtoms()
   {
       if (atoms == null) return 0;
       return atoms.length;
   }

   public Atom[] getAtoms()
   {
       return atoms;
   }
   public Bond[] getBonds()
   {
       return bonds;
   }

   /**
    * Deep copy of a molecule. Perceived information is NOT copied.
    */
   public Object clone()
   {
      Molecule result = new Molecule();

      if (atoms != null)
      {
         result.atoms = new Atom[atoms.length];
         for (int i=0; i<atoms.length; i++)
            result.atoms[i] = (Atom)this.atoms[i].clone();
      }
      if (bonds != null)
      {
         result.bonds = new Bond[bonds.length];
         for (int i=0; i<bonds.length; i++)
            result.bonds[i] = (Bond)this.bonds[i].clone();
      }
      if (stereodescs != null)
      {
         result.stereodescs = new StereoDescription[stereodescs.length];
         for (int i=0; i<stereodescs.length; i++)
            result.stereodescs[i] = (StereoDescription)this.stereodescs[i].clone(result);
      }
      if (setprops != null)
      {
         result.setprops = new SetProperty[setprops.length];
         for (int i=0; i<setprops.length; i++)
            result.setprops[i] = (SetProperty)this.setprops[i].clone();
      }

      return ((Object)result);
   }

   /**
    * Reset the color fields of all constituents of this molecule.
    */
   public void resetColors()
   {
      for (int i=0; atoms != null  &&  i<atoms.length; i++)
         atoms[i].color = 0;
      for (int i=0; bonds != null  &&  i<bonds.length; i++)
         bonds[i].color = 0;
   }

   /**
    * Recolors the atoms and bonds of this molecule each one
    * with a different color.
    */
   public void recolor()
   {
      int ncolor;

      ncolor = 1;
      for (int i=0; i<atoms.length; i++) atoms[i].color = ncolor++;
      for (int i=0; i<bonds.length; i++) bonds[i].color = ncolor++;
   }

   /**
    * Recursively colors the fragment of this molecule containing the atom a
    * with color. It uses a flood fill algorithm.
    *
    * It returns the number of atoms recolored during the process.
    */
   public int floodColor(Atom a, int color)
   {
      int result;
      int i;

      a.color = color; result = 1;
      for (i=0; i<a.neighbour_atoms.length; i++)
         if (a.neighbour_atoms[i].color == 0  &&
             (a.neighbour_bonds[i].flags & Layout.RUBBER_BOND) == 0)
	    result += floodColor(a.neighbour_atoms[i], color);

      return (result);
   }

   /**
    * Recursively clears the color value of those atoms colored with color.
    * The algorithms recursively goes through all neighbours atom a
    * colored with color. It returns the number of cleared colors.
    */
   public int floodClearColor(Atom a, int color)
   {
      int result;
      int i;

      if (a.color != color) return (0);

      a.color = 0; result = 1;
      for (i=0; i<a.neighbour_atoms.length; i++)
         if (a.neighbour_atoms[i].color == color  &&
             (a.neighbour_bonds[i].flags & Layout.RUBBER_BOND) == 0)
	    result += floodClearColor(a.neighbour_atoms[i], color);

      return (result);
   }

   /**
    * Clears the color fields of this molecule's atoms and returns an array
    * containing the old values.
    */
   public int[] removeAtomColors()
   {
      int[] result = new int[atoms.length];
      for (int i=0; i<atoms.length; i++)
      {
         result[i] = atoms[i].color;
         atoms[i].color = 0;
      }

      return (result);
   } 

   /**
    * Restore the colors of the atoms of this molecule to the state saved
    * in old_colors[].
    */
   public void restoreAtomColors(int old_colors[])
   {
      for (int i=0; i<atoms.length; i++)
         atoms[i].color = old_colors[i];
   }

   /**
    * Removes all perceived information from this molecule.
    */
   public void removePerceivedInformation()
   {
      if (atoms != null)
         for (int i=0; i<atoms.length; i++)
         {
            atoms[i].neighbour_atoms = null;
            atoms[i].neighbour_bonds = null;
         }
      ringlist        = null;
   }

   
   
   /**
    * Collects the neighbourhood information of this molecule. This
    * is a NOP if the information is already present.
    */
   public void setupNeighbourhood()
   {   
      if (atoms == null  || bonds == null) return;

      // count neighbours
      int [] nneighbours = new int[atoms.length];
      int [] nneighbonds = new int[bonds.length];

      for (int i=0; i<bonds.length; i++)
      {   
         nneighbours[bonds[i].atoms[0]-1]++;
         nneighbours[bonds[i].atoms[1]-1]++;
      }
	
      // allocate neighbour arrays
      for (int i=0; i<atoms.length; i++)
      {  
         atoms[i].index = i;
         atoms[i].neighbour_atoms = new Atom[nneighbours[i]];
         atoms[i].neighbour_bonds = new Bond[nneighbours[i]];
         nneighbours[i] = 0;
      }


      // collect neighbours
      for (int i=0; i<bonds.length; i++)
      {
         atoms[bonds[i].atoms[0]-1].neighbour_atoms[nneighbours[bonds[i].atoms[0]-1]] =
            atoms[bonds[i].atoms[1]-1];
         atoms[bonds[i].atoms[0]-1].neighbour_bonds[nneighbours[bonds[i].atoms[0]-1]] =
            bonds[i];
         atoms[bonds[i].atoms[1]-1].neighbour_atoms[nneighbours[bonds[i].atoms[1]-1]] =
            atoms[bonds[i].atoms[0]-1];
         atoms[bonds[i].atoms[1]-1].neighbour_bonds[nneighbours[bonds[i].atoms[1]-1]] =
            bonds[i];

         nneighbours[bonds[i].atoms[0]-1]++;
         nneighbours[bonds[i].atoms[1]-1]++;
      }
       for (int i=0; i<bonds.length; i++)
      {int bond_atom1 = bonds[i].atoms[0]-1;
       int bond_atom2 = bonds[i].atoms[1]-1;
       nneighbonds[i] = atoms[bonds[i].atoms[0]-1].neighbour_atoms.length +
                        atoms[bonds[i].atoms[1]-1].neighbour_atoms.length - 2;
       bonds[i].neighbourbonds = new Bond[nneighbonds[i]]; 
       nneighbonds[i] = 0;
       for (int j=0; j < nneighbours[bonds[i].atoms[0]-1]; j++) {
       if ((bond_atom1 == atoms[bonds[i].atoms[0]-1].neighbour_bonds[j].atoms[0]-1 &&
            bond_atom2 == atoms[bonds[i].atoms[0]-1].neighbour_bonds[j].atoms[1]-1) == false)
            {bonds[i].neighbourbonds[nneighbonds[i]] = atoms[bonds[i].atoms[0]-1].neighbour_bonds[j];
             nneighbonds[i]++;}
      }
        for (int j=0; j < nneighbours[bonds[i].atoms[1]-1]; j++) {
        if ((bond_atom1 == atoms[bonds[i].atoms[1]-1].neighbour_bonds[j].atoms[0]-1 &&
             bond_atom2 == atoms[bonds[i].atoms[1]-1].neighbour_bonds[j].atoms[1]-1) == false)
            {bonds[i].neighbourbonds[nneighbonds[i]] = atoms[bonds[i].atoms[1]-1].neighbour_bonds[j];
             nneighbonds[i]++;}

      }
      }
     
     }

   public void clearTopography()
   {
      for (int i=0; i<atoms.length; i++)
         atoms[i].topography = Atom.CHAIN;

      for (int i=0; i<bonds.length; i++)
      {
         bonds[i].topography = Bond.UNDEFINED;
         bonds[i].ring_sizes = 0;
      }
   }

   // perceived on demand.
   CountedBitSet[] ringlist = null;

   /**
    * Sets the topography of the bonds of this molecule to the appropriate ring
    * class and flags the ring size in the smallest_ring_size field.
    */
   public void perceiveRingBonds()
   {
      Graph g = new Graph(bonds);
      if (ringlist == null)
      {
         ringlist = g.ringList();
         ringlist = g.combineRings(ringlist);
      }

      int i, j;
      CountedBitSet ring;
      /* sort rings such that important ones come first */
      for (i=1; i<ringlist.length; i++)
         for (j=i; j>0; j--)
         {
            if ((ringlist[j].getCardinality() == 6  &&  ringlist[j-1].getCardinality() != 6)  ||
                (ringlist[j].getCardinality() < ringlist[j-1].getCardinality()))
            {
               ring = ringlist[j]; ringlist[j] = ringlist[j-1]; ringlist[j-1] = ring;
            }
            else
               break;
         }

      for (i=0; i<atoms.length; i++)
         atoms[i].topography = Atom.CHAIN;

      for (i=0; i<bonds.length; i++)
      {
         bonds[i].topography = Bond.CHAIN;
         bonds[i].ring_sizes = 0;
      }
      	
      for (int k=0; k<ringlist.length; k++)
      {  
         ring = ringlist[k];
         for (i=0; i<bonds.length; i++)
            if (ring.get(i))
            {  
               bonds[i].topography =  Bond.RING;
               bonds[i].ring_sizes |= 1<<Math.min(15, ring.getCardinality());
       	       atoms[bonds[i].atoms[0]-1].topography = Atom.RING;
               atoms[bonds[i].atoms[1]-1].topography = Atom.RING;
		
            }
      }
   }



   /**
    * Orients the bonds in this molecule such that the small basis rings are
    * traced in a clockwise direction.
    */
   public void makeRingsClockwise()
   {
      Bond b;
      CountedBitSet ring;
      double angle;
      int at1, at2, at3;
      int i, j, k, nbonds;

      perceiveRingBonds();

      if (ringlist == null) return;

      for (k=ringlist.length-1; k>=0; k--)
      {
         try
         {
             ring = ringlist[k];
             int[][] ring_bonds = new int[bonds.length][];
             nbonds = 0;                /* fetch bonds of ring */
             for (i=0; i<bonds.length; i++)
                if (ring.get(i))
                {
                   ring_bonds[nbonds] = new int[2];
                   ring_bonds[nbonds][0] = bonds[i].atoms[0];
                   ring_bonds[nbonds][1] = bonds[i].atoms[1];
                   nbonds++;
                }

             /* line up ring atoms */
             int[] ring_nodes = new int[nbonds];
             ring_nodes[0] = ring_bonds[0][0]; ring_nodes[1] = ring_bonds[0][1];
             for (i=2; i<nbonds; i++)
                for (j=1; j<nbonds; j++)
                {
                   if (ring_bonds[j][0] == ring_nodes[i-1]  &&  ring_bonds[j][1] != ring_nodes[i-2])
                   {
                      ring_nodes[i] = ring_bonds[j][1]; break;
                   }
                   if (ring_bonds[j][1] == ring_nodes[i-1]  &&  ring_bonds[j][0] != ring_nodes[i-2])
                   {
                      ring_nodes[i] = ring_bonds[j][0]; break;
                   }
                }

             angle = 0.0;      /* Sum up angles to test if clockwise */
             for (i=0; i<nbonds; i++)
             {
                at1 = ring_nodes[i]; at2 = ring_nodes[(i+1)%nbonds]; at3 = ring_nodes[(i+2)%nbonds];
                angle += HigherMath.Angle(atoms[at1-1].x - atoms[at2-1].x,
                                          atoms[at1-1].y - atoms[at2-1].y,
                                          atoms[at3-1].x - atoms[at2-1].x,
                                          atoms[at3-1].y - atoms[at2-1].y);
             }
             if (angle > (nbonds-1.5)*Math.PI)
             {                 /* counter clockwise -> swap direction */
                for (i=0; i<nbonds/2; i++)
                {
                   j = ring_nodes[i]; ring_nodes[i] = ring_nodes[nbonds-1-i]; ring_nodes[nbonds-1-i] = j;
                }
             }

             for (i=0; i<nbonds; i++)   // fix bonds
             {
                at1 = ring_nodes[i]; at2 = ring_nodes[(i+1)%nbonds];
                for (j=0; j<bonds.length; j++)
                {
                   b = bonds[j];
                   if (b.atoms[1]      == at1  &&
                       b.atoms[0]      == at2  &&
                       b.stereo_symbol == Bond.NORMAL)  // don't spoil stereo bonds!!!
                   {
                      b.atoms[0] = at1; b.atoms[1] = at2; break;
                   }
                }
             }
         }
         catch (Exception e)
         {
             System.err.println("Warning: makeRingsClockwise failed on ring " + k + " of " + ringlist.length + " rings");
         }
      }
   }


   /**
    * Create an empty molecule;
    */
   public Molecule()
   {
      // This is currently a NOP.
   }

   /**
    * Create molecule from MOL-File to be read from a LineBufferStream.
    * The connection table must start at the current position if sync == "".
    * Otherwise, the input is scanned for sync prior to reading the MOL-File.
    *
    * expect_header tells whether the caller expects a full MOL-File entry
    * or just the CTAB portion.
    */
   public void readMDLCTable(FortranInputStream in) throws IOException
   {
      int natoms, nbonds, natlists,  chiral, nstexts, nprops;
      String version;

      natoms = in.i(3); nbonds = in.i(3); natlists = in.i(3);
      in.x(3); chiral = in.i(3); nstexts = in.i(3);
      in.x(4*3); // ignore CPSS reaction stuff
      nprops = in.i(3);
      version = in.a(6);
      in.getBuffer();

      atoms = new Atom[natoms];
      for (int i=0; i<natoms; i++)
      {
         atoms[i] = parseMDLAtom(in); in.getBuffer();
      }

      bonds = new Bond[nbonds];
      for (int i=0; i<nbonds; i++)
      {
         bonds[i] = parseMDLBond(in); in.getBuffer();
      }

      for (int i=0; i<natlists; i++)
      {
         int iatom = in.i(3);
         boolean notLogic = in.a(2).trim().equals("T");
         int nInList = in.i(5);
         String[] atomList = new String[nInList];
         for (int j=0; j<nInList; j++)
         {
            atomList[j] = PTable.AtomicNumberToSymbol(in.i(4));
         }
         atoms[iatom-1].atomList = atomList;
         atoms[iatom-1].notLogic = notLogic;
         in.getBuffer();   // skip atlists
      }

      for (int i=0; i<nstexts; i++) in.getBuffer();   // skip stexts

      if (!version.equals("")) nprops = Integer.MAX_VALUE;  // there must be an 'M  END'

      readProperties(in, nprops);
   }


    /**
    * Parses the current buffer of in as an MDL atom record and
    * returns a new object containing the corresponding information.
    */
   public Atom parseMDLAtom(FortranInputStream in)
   {  int DEFAULT_HCOUNT = (-1); 
      String symbol;
      int mdiff, charge_radical, mapping;
      double x, y, z;

      x = in.f(10); y = in.f(10); z = in.f(10);
      in.x(1);
      symbol = in.a(3);
      mdiff  = in.i(2);
      charge_radical  = in.i(3);
      in.x(21);
      mapping  = in.i(3);

      return (new Atom(x, y, z,
                       symbol,
                       null, false,
                       DEFAULT_HCOUNT,
                       PTable.mdiffToIsotope(symbol, mdiff),
                       Atom.MDLChargeRadicalToCharge(charge_radical),
                       Atom.MDLChargeRadicalToRadical(charge_radical),
                       mapping));
   }
  
   /**
    * Parses the current buffer of in as an MDL bond record and
    * returns a new object containing the corresponding information.
    */
   public Bond parseMDLBond(FortranInputStream in)
   {
      int at1, at2, type, stereo, topo, mark;

      at1 = in.i(3); at2 = in.i(3);
      type = in.i(3);
      stereo = in.i(3);
      in.x(3);
      topo = in.i(3);
      mark = in.i(3);

      return (new Bond(at1, at2,
                       Bond.interpretMDLBondType(type),
                       Bond.interpretMDLStereoSymbol(stereo),
                       Bond.interpretMDLTopography(topo),
                       Bond.interpretMDLReactionMark(mark)));
   }

   /**
    * Read the property list off the stream in. Up to nprops lines are
    * read or until an "M  END" line is encountered.
    */
   protected void readProperties(FortranInputStream in, int nprops) throws IOException
   {
      boolean charge_found   = false;  // these tags are needed because
      boolean radical_found  = false;  // the property list superceede ALL
      boolean massdiff_found = false;  // charge, radical, and massdiffs in CTable

      int atno, ival;

      while (in.buffer != null  &&  !in.buffer.startsWith("M  END"))
      {
         if (in.buffer.startsWith("G  "))  // ignore CPSS group abbreviations
         {
            nprops -= 2;
            in.getBuffer();
            in.getBuffer();
            continue;
         }
         if (in.buffer.startsWith("A  "))  // ignore CPSS atom texts
         {
            nprops -= 2;
            in.x(3);
            int iatom  = in.i(3);
            in.getBuffer();
            atoms[iatom-1].setStringProperty("A", in.a(80).trim());
            in.getBuffer();
            continue;
         }

	 if (in.buffer.startsWith("V  "))  // ignore CPSS atom values
         {
             // [TODO] Need to use the atom value to drive label and sequence rendering
            nprops -= 1;
            in.getBuffer();
            continue;
         }

         if (in.buffer.startsWith("S  SKP"))  // ignore skip records
         {
            in.x(6);
            int nskip = in.i(3);
            nprops -= nskip;
            for (int i=0; i<nskip; i++) in.getBuffer();
            continue;
         }

         // now we are at the real values
         if (in.buffer.startsWith("M  CHG"))
         {
            if (!charge_found)   // this is the first charge property
            {                    // => clear all atom charges
               charge_found = true;
               for (int i=0; i<atoms.length; i++)
                  atoms[i].charge = 0;
            }
            in.x(6);
            int ncharge = in.i(3);
            for (int i=0; i<ncharge; i++)
            {
               in.x(1); atno = in.i(3);
               in.x(1); ival = in.i(3);
               atoms[atno-1].charge = ival;
            }
            nprops--; in.getBuffer();
         }
         else if (in.buffer.startsWith("M  SUB"))
         {
            in.x(6);
            int nprop = in.i(3);
            for (int i=0; i<nprop; i++)
            {
               in.x(1); atno = in.i(3);
               in.x(1); ival = in.i(3);
               atoms[atno-1].setIntProperty("SUB", ival);
            }
            nprops--; in.getBuffer();
         }
         else if (in.buffer.startsWith("M  ALS"))
         {
            // [TODO] implement object properties with names and make this
            // a String[] valued property of atoms.
            /*
            in.x(6);
            int nprop = in.i(3);
            for (int i=0; i<nprop; i++)
            {
               in.x(1); atno = in.i(3);
               in.x(1); ival = in.i(3);
               atoms[atno-1].setIntProperty("RGP", ival);
            }
            */
            nprops--; in.getBuffer();
         }
         else if (in.buffer.startsWith("M  ISO"))
         {
            if (!massdiff_found)   // this is the first mass_diff property
            {                      // => clear all atom isotopes
               massdiff_found = true;
               for (int i=0; i<atoms.length; i++)
                  atoms[i].isotope = 0;
            }
            in.x(6);
            int nisotope = in.i(3);
            for (int i=0; i<nisotope; i++)
            {
               in.x(1); atno = in.i(3);
               in.x(1); ival = in.i(3);
               atoms[atno-1].isotope = ival;
                  // PTable.mdiffToIsotope(atoms[atno-1].symbol, ival);
            }
            nprops--; in.getBuffer();
         }
         else if (in.buffer.startsWith("M  RGP"))
         {
            in.x(6);
            int nprop = in.i(3);
            for (int i=0; i<nprop; i++)
            {
               in.x(1); atno = in.i(3);
               in.x(1); ival = in.i(3);
               atoms[atno-1].setIntProperty("RGP", ival);
            }
            nprops--; in.getBuffer();
         }
         else if (in.buffer.startsWith("M  RAD"))
         {
            if (!radical_found)   // this is the first radical property
            {                    // => clear all atom radical assignments
               radical_found = true;
               for (int i=0; i<atoms.length; i++)
                  atoms[i].radical = 0;
            }
            in.x(6);
            int nprop = in.i(3);
            for (int i=0; i<nprop; i++)
            {
               in.x(1); atno = in.i(3);
               in.x(1); ival = in.i(3);
               if (ival >= 3) ival = Atom.TRIPLET;     // map to our radical states
               atoms[atno-1].radical = ival;
            }
            nprops--; in.getBuffer();
         }
	 
         else  // skip unknown entries
         {
             // System.err.println("Skipping unknown property '" + in.buffer + "'");
            nprops--; in.getBuffer();
         }
      }
   }

   /**
    * Allows modification of the molecule and disables methods which
    * rely on a consistent state of the data structure.
    */
   public void modifyOn()
   {
      being_modified  = true;
      removePerceivedInformation();
   }

   /**
    * Disallows modification of the molecule and enables methods which
    * rely on a consistent state of the data structure.
    */
   public void modifyOff()
   {
      // The method should fix rapid access links here.
      being_modified = false;
   }

   /**
    * Computes the implicit hydrogen counts for the atoms in *mp. The
    * result is returned as an array of integers. This array has index origin 1,
    * i.e. its element [0] is not used.
    */
   public int[] ComputeImplicitH()
   {
      int[] H_count = new int[atoms.length+1];  // resulting array

      int[] single_bond   = new int[atoms.length+1];      /* <array>[0] is unused */
      int[] aromatic_bond = new int[atoms.length+1];
      int[] double_bond   = new int[atoms.length+1];
      int[] triple_bond   = new int[atoms.length+1];
      int[] radical       = new int[atoms.length+1];
      int[] charge        = new int[atoms.length+1];

      Atom a;
      Bond b;

      for (int i=0; i<atoms.length; i++)
      {
         radical[i+1] = atoms[i].radical;
         charge[i+1]  = atoms[i].charge;
      }

      for (int i=0; i<=atoms.length; i++)
         single_bond[i] = aromatic_bond[i] = double_bond[i] = triple_bond[i] = 0;

      for (int i=0; i<bonds.length; i++)
      {
         b = bonds[i];
         switch (b.bond_type)
         {
            case Bond.SINGLE: single_bond[b.atoms[0]]++;
                              single_bond[b.atoms[1]]++;
                              break;
            case Bond.DOUBLE: double_bond[b.atoms[0]]++;
                              double_bond[b.atoms[1]]++;
                              break;
            case Bond.TRIPLE: triple_bond[b.atoms[0]]++;
                              triple_bond[b.atoms[1]]++;
                              break;
            case Bond.AROMATIC: aromatic_bond[b.atoms[0]]++;
                                aromatic_bond[b.atoms[1]]++;
                                break;
            default :
               single_bond[b.atoms[0]]++;
               single_bond[b.atoms[1]]++;
               break;
         }
      }

      for (int i=0; i<atoms.length; i++)
      {
         H_count[i+1] = PTable.implicitHydrogens(atoms[i].symbol,
                                                 single_bond[i+1],
                                                 aromatic_bond[i+1],
                                                 double_bond[i+1],
                                                 triple_bond[i+1],
                                                 radical[i+1],
                                                 charge[i+1]);
         if (atoms[i].implicit_H_count > H_count[i+1])
            H_count[i+1] = atoms[i].implicit_H_count;
      }

      return (H_count);
   }

   /**
    * Flips the stereo symboles of the parts of this molecule that are
    * colored with color.
    */
   public void flipStereoSymbols(int color)
   {
      for (int i=0; i<bonds.length; i++)
      {
         Bond b = bonds[i];
         if (color == 0  ||
             (atoms[b.atoms[0]-1].color == color &&
              atoms[b.atoms[1]-1].color == color))
            if (b.stereo_symbol == Bond.UP)
               b.stereo_symbol = Bond.DOWN;
            else if (b.stereo_symbol == Bond.DOWN)
               b.stereo_symbol = Bond.UP;
      }
   }

   /**
    * Flips the parts of this molecule that are colored with color
    * with respect to the vertical (Y) axis. If color is 'NONE',
    * the whole molecule is flipped.
    */
   public void flipMolecule(int color)
   {
      int i;
      int natoms;
      double xcenter;

      xcenter = 0.0; natoms = 0;
      for (i=0; i<atoms.length; i++)
         if (color == 0 || atoms[i].color == color)
         {
            xcenter += atoms[i].x;
            natoms++;
         }

      if (natoms == 0) return;

      xcenter /= natoms;

      for (i=0; i<atoms.length; i++)
         if (color == 0 || atoms[i].color == color)
            atoms[i].x = xcenter-atoms[i].x;

      flipStereoSymbols(color);
   }

   /**
    * Fetches the edges from this molecule that are colored with color.
    *
    * numbers[i] is set to the new number of atoms[i] in the returned list
    * edges.
    */
   public int[][] getColoredEdges(int color, int numbers[])
   {
      int[][] edges;

      Vector<int[]> vector = new Vector<int[]>();
      for (int i=0; i<bonds.length; i++)
      {
         Bond b = bonds[i];
         if ((b.flags & Layout.RUBBER_BOND) == 0         &&
             atoms[b.atoms[0]-1].color == color  &&
             atoms[b.atoms[1]-1].color == color)
            vector.addElement(b.getAtomNumbers());
      }
      // collect new numbering
      int nnodes = 0;
      for (int i=0; i<atoms.length; i++)
         if (atoms[i].color == color)
         {
	         numbers[i] = nnodes; nnodes++;
         }

      // recast result and renumber edges
      edges = new int[vector.size()][];
      for (int i=0; i<edges.length; i++)
      {
         edges[i] = (int[]) vector.elementAt(i);
         edges[i][0] = numbers[edges[i][0]-1];
         edges[i][1] = numbers[edges[i][1]-1];
      }

      return (edges);
   }

   /**
    * Fetches the the coordinates of those nodes from this molecule
    * that are colored with color.
    *
    * The array numbers[] is set such that it transforms indices into
    * atoms[] to indices into returned array of coordinates.
    */
   public double[][] getColoredCoordinates(int color, int numbers[])
   {
      Vector<double[]> vector = new Vector<double[]>();
      int nnodes = 0;
      for (int i=0; i<atoms.length; i++)
         if (atoms[i].color == color)
         {
	         double[] point = new double[2];
	         point[0] = atoms[i].x;
	         point[1] = atoms[i].y;
	         vector.addElement(point);
	         numbers[i] = nnodes; nnodes++;
      }
      double [][] result = new double[vector.size()][];
      for (int i=0; i<vector.size(); i++)
         result[i] = (double [])vector.elementAt(i);
      return (result);
   }

   /**
    * Counts the number of atoms in this molecule which are colored in color.
    */
   public int countColor(int color)
   {
      int result = 0;

      for (int i=0; i<atoms.length; i++)
      if (atoms[i].color == color) result++;

      return (result);
   }

   /**
    * Denormalize aromatic bonds to single/double bonds if possible. Gives up
    * on failed denormalized components and returns false. Returns true on success.
    */
   public boolean denormalize()
   {  
      boolean result = true;
      setupNeighbourhood(); perceiveRingBonds();

      AtomConstraint ac[] = perceiveAtomConstraints();
      colorNormalizedFragments();

      for (int i=0; i<atoms.length; i++)
      { 
         Atom a = atoms[i];
         if (ac[i].is_open)
         {   
	    Bond[] tmp_bonds = denormalizeStep(bonds, ac, a.color);
            if (tmp_bonds != null) bonds = tmp_bonds;
            for (int j=0; j<ac.length; j++)
               if (atoms[j].color == a.color) ac[j].is_open = false;
         }
      }
      return (result);
   }
   
   
    /**
    * Refine the bond order assignment bonds based on the given atom constraints ac
    * for the fragment colored by color.
    */
   boolean refineBonds(Bond[] bonds, AtomConstraint[] ac, int color)
   {
      boolean changed;
      int nsingle[]    = new int[atoms.length];
      int ndouble[]    = new int[atoms.length];
      int ntriple[]    = new int[atoms.length];
      int nambiguous[] = new int[atoms.length];

      for (int i=0; i<bonds.length; i++)
      {
         Bond b = bonds[i];
         if (b.bond_type == Bond.SINGLE)
         {
            nsingle[b.atoms[0]-1]++;
            nsingle[b.atoms[1]-1]++;
         }
         else if (b.bond_type == Bond.DOUBLE)
         {
            ndouble[b.atoms[0]-1]++;
            ndouble[b.atoms[1]-1]++;
         }
         else if (b.bond_type == Bond.TRIPLE)
         {
            ntriple[b.atoms[0]-1]++;
            ntriple[b.atoms[1]-1]++;
         }
         else if (b.bond_type == Bond.NOBOND)
         {
            //NOP
         }
         else
         {
            nambiguous[b.atoms[0]-1]++;
            nambiguous[b.atoms[1]-1]++;
         }
      }

      do
      {
         changed = false;
         for (int i=0; i<atoms.length; i++)
         {
            Atom a = atoms[i];
            if (a.color == color)
            {
               // check for constraint violation
               if (nambiguous[i] == 0  &&
                   (ndouble[i] < ac[i].dbl_min  ||
                    ndouble[i] > ac[i].dbl_max))
                  return (false);
               if (nambiguous[i] != 0  &&  ndouble[i] == ac[i].dbl_max)
               {  // all ambiguous bonds must be single
                  changed = true;
                  for (int j=0; j<bonds.length; j++)
                  {
                     Bond b = bonds[j];
                     if ((b.atoms[0]-1 == i  ||  b.atoms[1]-1 == i)  &&
                         (b.bond_type  ==  Bond.AROMATIC  ||
                          b.bond_type  ==  (Bond.SINGLE | Bond.DOUBLE)))
                     {
                        b.bond_type = Bond.SINGLE;
                        nambiguous[b.atoms[0]-1]--;
                        nsingle[b.atoms[0]-1]++;
                        nambiguous[b.atoms[1]-1]--;
                        nsingle[b.atoms[1]-1]++;
                     }
                  }
               }
               else if (nambiguous[i] != 0  &&
                        ndouble[i]+nambiguous[i] == ac[i].dbl_min)
               {  // all ambiguous bonds must be double
                  changed = true;
                  for (int j=0; j<bonds.length; j++)
                  {
                     Bond b = bonds[j];
                     if ((b.atoms[0]-1 == i  ||  b.atoms[1]-1 == i)  &&
                         (b.bond_type  ==  Bond.AROMATIC  ||
                          b.bond_type  ==  (Bond.SINGLE | Bond.DOUBLE)))
                     {
                        b.bond_type = Bond.DOUBLE;
                        nambiguous[b.atoms[1]-1]--;
                        ndouble[b.atoms[1]-1]++;
                        nambiguous[b.atoms[0]-1]--;
                        ndouble[b.atoms[0]-1]++;
                     }
                  }
               }
            }
         }
      } while (changed);

      return (true);
   }
   


  
   /**
    * Compare the two denormalized bond arrays bonds1[] and bonds2[]. It prefers
    * alternating six-membered rings and doubly bonded 'O' atoms.
    */
   int compareDenormalizations(Bond[] bonds1, Bond[] bonds2)
   {
      if (bonds2 == null) return (1);
      else if (bonds1 == null) return (-1);
      // for test purposes just prefer first one
      return (1);
   }

   /**
    * Recursively denormalizes the aromatic and ambiguous bonds in bonds and
    * returns the denormalized bond list.
    */



  Bond[] denormalizeStep(Bond[] bonds, AtomConstraint[] ac, int color)
   {
   if (bonds == null) return (null);
   if (!refineBonds(bonds, ac, color)) // constraint violation
   return(null);

   
    for (int i=0; i<bonds.length; i++)
      {  Bond b = bonds[i];
          if ((b.bond_type == Bond.AROMATIC ||
          b.bond_type == (Bond.SINGLE | Bond.DOUBLE))
          && b.color == color)
          {
          Bond single_case[] = new Bond[bonds.length];
          for (int j=0; j<bonds.length; j++) single_case[j] = (Bond)bonds[j].clone();
          single_case[i].bond_type = Bond.SINGLE;
          single_case = denormalizeStep(single_case, ac, color);
          if (single_case != null) return(single_case);
          Bond double_case[] = new Bond[bonds.length];
          for (int j=0; j<bonds.length; j++) double_case[j] = (Bond)bonds[j].clone();
          double_case[i].bond_type = Bond.DOUBLE;
          double_case = denormalizeStep(double_case, ac, color);
          if (double_case != null) return (double_case);
	  if (single_case == null && double_case == null) return (null);
          }
          }
    return (bonds);
    }

   /**
    * Recolors the atoms and bonds of this molecule such that pieces connected by
    * AROMATIC bonds are colored in mutually distinct colors. Each fragment itself
    * will be colored with the same color. This is necessary to denormalize each

    * independent normalized set of bonds independently and thus reduce combinatorial
    * problems.
    */
   void colorNormalizedFragments()
   {
      recolor();  // all bonds and atoms get a different color.

      for (int i=0; i<bonds.length; i++)
      {
         Bond b = bonds[i];
         if (b.bond_type == Bond.AROMATIC  ||
             b.bond_type == (Bond.SINGLE | Bond.DOUBLE) ||
	     b.bond_type == (Bond.SINGLE | Bond.NOBOND))

         { 
            int col1 = atoms[b.atoms[0]-1].color;
            int col2 = atoms[b.atoms[1]-1].color;
            int colb = b.color;
            int col = Math.max(col1, Math.max(col2, colb));

            if (col1 == col2  &&  col1 == colb) continue;

            for (int j=0; j<atoms.length; j++)
               if (atoms[j].color == col1  ||
                   atoms[j].color == col2  ||
                   atoms[j].color == colb) atoms[j].color = col;
            for (int j=0; j<bonds.length; j++)
               if (bonds[j].color == col1  ||
                   bonds[j].color == col2  ||
                   bonds[j].color == colb) bonds[j].color = col;
         }
      }
   }

   /**
    * Private method perceiving the atom constraints for denormalization.
    */
   AtomConstraint[] perceiveAtomConstraints()
   {
      AtomConstraint[] result = new AtomConstraint[atoms.length];
      int[] nsingle   = new int[atoms.length];
      int[] ndouble   = new int[atoms.length];
      int[] ntriple   = new int[atoms.length];
      int[] naromatic = new int[atoms.length];
      for (int i=0; i<bonds.length; i++)
      {
         Bond b = bonds[i];
         if (b.bond_type == Bond.SINGLE)
         {
            nsingle[b.atoms[0]-1]++; nsingle[b.atoms[1]-1]++;
         }
         else if (b.bond_type == Bond.DOUBLE)
         {
            ndouble[b.atoms[0]-1]++; ndouble[b.atoms[1]-1]++;
         }
         else if (b.bond_type == Bond.TRIPLE)
         {
            ntriple[b.atoms[0]-1]++; ntriple[b.atoms[1]-1]++;
         }
         else if (b.bond_type == Bond.AROMATIC)
         {
            naromatic[b.atoms[0]-1]++; naromatic[b.atoms[1]-1]++;
         }
         else if (b.bond_type == (Bond.DOUBLE | Bond.SINGLE))
         {
            naromatic[b.atoms[0]-1]++; naromatic[b.atoms[1]-1]++;
         }
         else if (b.bond_type == Bond.NOBOND)
         {
            // NOP
         }
         else  // treat as bond with highest order
         {
            if ((b.bond_type & Bond.TRIPLE)  !=  0)
            {
               ntriple[b.atoms[0]-1]++; ntriple[b.atoms[1]-1]++;
            }
            else if ((b.bond_type & Bond.DOUBLE)  !=  0)
            {
               ndouble[b.atoms[0]-1]++; ndouble[b.atoms[1]-1]++;
            }
            else if ((b.bond_type & Bond.AROMATIC)  !=  0)
            {
               naromatic[b.atoms[0]-1]++; naromatic[b.atoms[1]-1]++;
            }
            else
            {
               nsingle[b.atoms[0]-1]++; nsingle[b.atoms[1]-1]++;
            }
         }
      }

      for (int i=0; i<atoms.length; i++)
      {
         Atom a = atoms[i];
         result[i] = new AtomConstraint(a.symbol,
                                        nsingle[i], ndouble[i], ntriple[i], naromatic[i],
                                        a.charge, a.radical,
                                        a.implicit_H_count);
      }
      return (result);
   } // perceiveAtomConstraints()

   /**
    * Perceives reacting bonds in reactant and product.
    */
   public static void perceiveReaction(Molecule reactant, Molecule product)
   {
      for (int i=0; i<reactant.bonds.length; i++)
      {
         Bond br = reactant.bonds[i];
         int atoms[] = br.getAtomNumbers();
         int map1 = reactant.atoms[atoms[0]-1].reaction_mapping;
         int map2 = reactant.atoms[atoms[1]-1].reaction_mapping;
         if (map1 == Atom.NO_MAPPING  &&  map2 != Atom.NO_MAPPING   ||   // bond to removed fragment
             map1 != Atom.NO_MAPPING  &&  map2 == Atom.NO_MAPPING)
         {
            br.reaction_mark = Bond.MAKE_BRAKE;
         }
         else if (map1 != Atom.NO_MAPPING  &&  map2 != Atom.NO_MAPPING)  // mapped bond
         {
            br.reaction_mark = Bond.UNDEFINED;
            Bond bp = null;
            int j;
            for (j=0; j<product.bonds.length; j++)   // search for mapped bond
            {
               bp = product.bonds[j];
               int patoms[] = bp.getAtomNumbers();
               if (map1 == product.atoms[patoms[0]-1].reaction_mapping  &&
                   map2 == product.atoms[patoms[1]-1].reaction_mapping  ||
                   map1 == product.atoms[patoms[1]-1].reaction_mapping  &&
                   map2 == product.atoms[patoms[0]-1].reaction_mapping) break;
            }
            if (j == product.bonds.length)   // no mapped bond found => broken bond
            {
               br.reaction_mark = Bond.MAKE_BRAKE;
            }
            else if (br.bond_type != bp.bond_type)     // this should be done with atom hybridization
            {                                          // to catch aromaticity!!!!
               br.reaction_mark = Bond.CHANGED;
            }
         }
      }

      for (int i=0; i<product.bonds.length; i++)
      {
         Bond bp = product.bonds[i];
         int atoms[] = bp.getAtomNumbers();
         int map1 = product.atoms[atoms[0]-1].reaction_mapping;
         int map2 = product.atoms[atoms[1]-1].reaction_mapping;
         if (map1 == Atom.NO_MAPPING  &&  map2 != Atom.NO_MAPPING   ||   // bond to removed fragment
             map1 != Atom.NO_MAPPING  &&  map2 == Atom.NO_MAPPING)
         {
            bp.reaction_mark = Bond.MAKE_BRAKE;
         }
         else if (map1 != Atom.NO_MAPPING  &&  map2 != Atom.NO_MAPPING)  // mapped bond
         {
            bp.reaction_mark = Bond.UNDEFINED;
            Bond br = null;
            int j;
            for (j=0; j<reactant.bonds.length; j++)   // search for mapped bond
            {
               br = reactant.bonds[j];
               int ratoms[] = br.getAtomNumbers();
               if (map1 == reactant.atoms[ratoms[0]-1].reaction_mapping  &&
                   map2 == reactant.atoms[ratoms[1]-1].reaction_mapping  ||
                   map1 == reactant.atoms[ratoms[1]-1].reaction_mapping  &&
                   map2 == reactant.atoms[ratoms[0]-1].reaction_mapping) break;
            }
            if (j == reactant.bonds.length)   // no mapped bond found => broken bond
            {
               bp.reaction_mark = Bond.MAKE_BRAKE;
            }
            else if (br.bond_type != bp.bond_type)     // this should be done with atom hybridization
            {                                          // to catch aromaticity!!!!
               bp.reaction_mark = Bond.CHANGED;
            }
         }
      }
   }

}  // Molecule
