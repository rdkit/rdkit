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

import novartis.combinatorics.Permutation;

/**
 * Stereodescription class
 */
public class StereoDescription
{
   Molecule parent;     // parent molecule to which the constituents of this StereoDescription refer.
                        // this link is needed to make the correct references for a cloned object.
   Object depictor;     // Reference to the object used for depicting the stereo info
                        // This can be an atom, a bond, or a special depiction object.
   int    atoms[];      // Array of atoms which this stereodescription defines (restricts)

   int    stereoclass;  // stereo description class. One of the above constants

   public StereoDescription(Molecule parent, Object depictor, int atoms[], int stereoclass)
   {
      this.parent      = parent;
      this.depictor    = depictor;
      this.atoms       = new int[atoms.length];
      System.arraycopy(atoms, 0, this.atoms, 0, atoms.length);
      this.stereoclass = stereoclass;
   }

   public StereoDescription(Molecule parent, StereoDescription source)
   {
      this.parent      = parent;
      this.depictor    = null;
      this.atoms       = new int[source.atoms.length];
      System.arraycopy(source.atoms, 0, this.atoms, 0, source.atoms.length);
      this.stereoclass = source.stereoclass;
      if (source.parent.atoms != null)
         for (int i=0; i<source.parent.atoms.length; i++)
            if (source.depictor == source.parent.atoms[i])
            {
               this.depictor = parent.atoms[i];
               return;
            }
      if (source.parent.bonds != null)
         for (int i=0; i<source.parent.bonds.length; i++)
            if (source.depictor == source.parent.bonds[i])
            {
               this.depictor = parent.bonds[i];
               return;
            }
   }

   // atom number of implicit hydrogen
   public final static int IMPLICIT_H = Integer.MAX_VALUE;

   // stereo class tags
   public final static int NONE = 0;   // no definition

   public final static int TH = 1;     // tetrahedral center
   public final Permutation THPermutations[] =  // alternating group of 4 elements
   {
      new Permutation("1 2 3 4"),
      new Permutation("1 3 4 2"),
      new Permutation("1 4 2 3"),
      new Permutation("2 3 1 4"),
      new Permutation("2 1 4 3"),
      new Permutation("2 4 3 1"),
      new Permutation("3 1 2 4"),
      new Permutation("3 2 4 1"),
      new Permutation("3 4 1 2"),
      new Permutation("4 2 1 3"),
      new Permutation("4 1 3 2"),
      new Permutation("4 3 2 1"),
   };

   public final static int DB = 2;     // cis/trans double bond
   public final Permutation DBPermutations[] =
   {
      new Permutation("1 2 3 4"),
      new Permutation("2 1 4 3"),
      new Permutation("3 4 1 2"),
      new Permutation("4 3 2 1"),
   };

   public Object clone(Molecule new_parent)
   {
      StereoDescription result = new StereoDescription(new_parent, this);
      return (result);
   }
}
