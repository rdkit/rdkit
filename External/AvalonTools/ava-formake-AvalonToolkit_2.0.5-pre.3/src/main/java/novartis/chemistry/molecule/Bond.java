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

/**
 * Represents a bond of a connection table.
 */
public class Bond
{
   // instance variables
   int atoms[];   // Numbers of atoms connected by bond. Index origin is 1; 

   public int[] getAtomNumbers()
   {
      int result[] = new int[2];
      result[0] = atoms[0];
      result[1] = atoms[1];
      return (result);
   }

   int bond_type; // bond order type of bond.
   // Bond order constants. Note: They are NOT 1, 2, 3, but 2, 4, 8!!
   // This is to allow for bond type sets in the same field.
   public final static int NOBOND   = (1 << 0);
   public final static int SINGLE   = (1 << 1);
   public final static int DOUBLE   = (1 << 2);
   public final static int TRIPLE   = (1 << 3);
   public final static int AROMATIC = (1 << 4);

   /**
    * Accessor function for bond type
    */
   public int getBondType()
   {
      return (this.bond_type);
   }

   /**
    * Accessor function for bond type
    */
   public int getBondFlags()
   {
      return (this.flags);
   }

   /**
    * Converts the MDL CTable bond type type into our bond_type convention.
    */
   public static int interpretMDLBondType(int type)
   {
      if      (type == 1) return (SINGLE);
      else if (type == 2) return (DOUBLE);
      else if (type == 3) return (TRIPLE);
      else if (type == 4) return (AROMATIC);
      else if (type == 5) return (SINGLE | DOUBLE);
      else if (type == 6) return (SINGLE | AROMATIC);
      else if (type == 7) return (DOUBLE | AROMATIC);
      else if (type == 8) return (NOBOND | SINGLE | DOUBLE | TRIPLE | AROMATIC);
      else                return (NOBOND);
   }

   // common (usually perceived) but non-essential bond attributes
   int stereo_symbol;
   // Bond stereo symbol types. Note: the values are not identical to
   // the MDL bond_stereo types.
   public final static int NORMAL     = 0; // use symbol defined for given bond type
   public final static int UP         = 1; // up-pointing wedge
   public final static int DOWN       = 2; // down pointing hash
   public final static int EITHER     = 3; // wiggle bond
   public final static int DBL_EITHER = 4; // unknown DB configuration, '?' over bond
   public final static int CIS        = 5; // ligands with lower numbers should be on same side
   public final static int TRANS      = 6; // ligands with lower numbers should be on opposite sides

   /**
    * Convert the MDL CTable stereo symbol to our bond symbol convention
    */
   public static int interpretMDLStereoSymbol(int stereo)
   {
      if      (stereo == 1) return (UP);
      else if (stereo == 6) return (DOWN);
      else if (stereo == 4) return (EITHER);
      else if (stereo == 3) return (DBL_EITHER);
      else                  return (NORMAL);
   }

   int topography;
   // Bond ring indicator values. Note: the values are not identical to
   // the MDL topography types.
   public final static int UNDEFINED  = 0; // no topography defined
   public final static int RING       = 1; // bond is defined as a ring bond
   public final static int CHAIN      = 2; // bond is defined as a chain bond

   /**
    * Convert the MDL CTable stereo symbol to our bond symbol convention
    */
   public static int interpretMDLTopography(int topo)
   {
      if      (topo == 1) return (RING);
      else if (topo == 2) return (CHAIN);
      else                return (UNDEFINED);
   }

   int reaction_mark;
   // Bond reaction indicators. Note: the values are not identical to
   // the MDL reaction mark types.
 //public final static int UNDEFINED   = 0; // no reaction indocator defined
   public final static int CHANGED     = 1; // bond order changed, but bond wasn't broken
   public final static int MAKE_BRAKE  = 2; // bond was made or broken
   public final static int UNCHANGED   = 4; // bond order unchanged
   public final static int ANY_CHANGE  = 3; // changed, made, or broken
   public final static int UNKNOWN     = 7; // unknown change

   /**
    * Convert the MDL CTable stereo symbol to our bond symbol convention
    */
   public static int interpretMDLReactionMark(int mark)
   {
      if      (mark ==  1) return (ANY_CHANGE);
      else if (mark ==  2) return (UNCHANGED);
      else if (mark ==  4) return (MAKE_BRAKE);
      else if (mark ==  6) return (UNKNOWN);
      else if (mark ==  8) return (CHANGED);
      else if (mark == -1) return (UNCHANGED);
      else                 return (UNDEFINED);
   }

   // perceived information. not always present
   int ring_sizes; // bit set with one bit per ring size up to size 15.
                   // 0 for chain bonds

   /**
    * Returns the size of the smallest ring (up to 15) of which this bond
    * is a member or 0 if it is a chain bond.
    */
   public int smallestRing()
   {
      if (ring_sizes == 0) return (0);
      for (int i=3; i<=15; i++)
         if ((ring_sizes & (1<<i)) != 0) return (i);
      return (15);
   }

   // convenience fields used for graph algorithms
   float value  = 0;
   int   color  = 0;
   int   flags  = 0;
   int   number = 0;

   public Bond(int at1, int at2, int type, int stereo, int topo, int mark)
   {
      atoms = new int[2]; atoms[0] = at1; atoms[1] = at2;
      bond_type = type; stereo_symbol = stereo;
      topography = topo;
      reaction_mark = mark;
   }

   public Object clone()
   {
      Bond result = new Bond(atoms[0], atoms[1],
                             bond_type, stereo_symbol, topography, reaction_mark);
      result.value  = this.value;
      result.color  = this.color;
      result.flags  = this.flags;
      result.number = this.number;

      return (result);
   }


   Bond[] neighbourbonds = null;    // filled in Molecule.setupNeighbourBonds;  

}
