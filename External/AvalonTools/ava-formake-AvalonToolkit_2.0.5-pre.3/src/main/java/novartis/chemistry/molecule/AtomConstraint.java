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
 * Constraint bond order set for an atom. To be used during denormalization.
 */
public class AtomConstraint
{
   int dbl_min;      // Minimum number of double bonds attached to this atom.
   int dbl_max;      // Maximum number of double bonds attached to this atom.
   boolean is_open;  // Flag used during processing indicating that the bonding
                     // pattern of this atom is to be defined.


   public AtomConstraint(AtomConstraint source)
   {
      this.dbl_min = source.dbl_min;
      this.dbl_max = source.dbl_max;
      this.is_open = source.is_open;
   }
   

   /**
    * Create an atom constraint based of the number of bonds of each type connected to
    * an atom and its atom symbol. This method uses a set of quick rules for common
    * cases and a lookup in PTable otherwise.
    */
   public AtomConstraint(String symbol,
                         int nsingle, int ndouble, int ntriple, int naromatic,
                         int charge, int radical,
                         int hcount)
   {
      this.is_open = (naromatic > 0);
      if (naromatic == 0)
     { // normal non-aromatic atom
      this.dbl_min = ndouble; this.dbl_max = ndouble; this.is_open = false;
      }

      else if (ndouble >  0  &&  naromatic > 0  &&
               charge  == 0  &&  radical == 0   &&
               symbol.equals("C"))
     { // pseudo-aromatic carbon used in Daylight conventions
                this.dbl_min = ndouble; this.dbl_max = ndouble; this.is_open = true;}

      else if (nsingle == 0  &&  ndouble == 0  &&  ntriple == 0  && naromatic == 2  &&
               charge  == 0  &&  radical == 0                                       &&
               symbol.equals("C")  &&
               (hcount == Atom.DEFAULT_HCOUNT  ||  hcount == 1))
      {  // unsubstituted aromatic carbon
         this.dbl_min = 1; this.dbl_max = 1; this.is_open = true;
      }

      else if (nsingle == 0  &&  ndouble == 0  &&  ntriple == 0  && naromatic == 2  &&
               charge  == 0  &&  radical == 0  &&
               symbol.equals("N")  &&
               (hcount == Atom.DEFAULT_HCOUNT  ||  hcount == 0))
      {  // unsubstituted aromatic nitrogen, like in pyridine
         this.dbl_min = 1; this.dbl_max = 1; this.is_open = true;
      }

      else if (ndouble == 0  &&  ntriple == 0  && naromatic == 2  &&
               charge  == 0  &&  radical == 0                     &&
               symbol.equals("N")  &&
               ((nsingle == 0  &&  hcount == 1)  ||  nsingle == 1))
      {  // nitrogen like in pyrole
         this.dbl_min = 0; this.dbl_max = 0; this.is_open = true;
      }

      else if (nsingle == 0  &&  ndouble == 0  &&  ntriple == 0  && naromatic == 2  &&
               charge  == 0  &&  radical == 0                                       &&
               symbol.equals("O")  &&
               (hcount == Atom.DEFAULT_HCOUNT  ||  hcount == 0))
      {  // oxygen like in furan
         this.dbl_min = 0; this.dbl_max = 0; this.is_open = true;
      }
      else if (nsingle == 0  &&  ndouble == 0  &&  ntriple == 0  && naromatic == 1  &&
               charge  == 0  &&  radical == 0                                       &&
               symbol.equals("O")  &&
               (hcount == Atom.DEFAULT_HCOUNT  ||  hcount == 0))
      {  // oxygen of phosphate groups
         this.dbl_min = 0; this.dbl_max = 1; this.is_open = true;
       }

      else if (nsingle == 0  &&  ndouble == 0  &&  ntriple == 0  && naromatic == 1  &&
             charge  == 0  &&  radical == 0                                       &&
            symbol.equals("N")  &&
           (hcount == Atom.DEFAULT_HCOUNT  ||  hcount == 0))
                         {  // terminal nitrogen guanidine
                          this.dbl_min = 0; this.dbl_max = 1; this.is_open = true;
                         }


      else if (nsingle == 0  &&  ndouble == 0  &&  ntriple == 0  && naromatic == 3  &&
               charge  == 0  &&  radical == 0  &&
               symbol.equals("C")                                                   &&
               (hcount == Atom.DEFAULT_HCOUNT  ||  hcount == 0))
      {  // substituted aromatic carbon
         this.dbl_min = 1; this.dbl_max = 1; this.is_open = true;
      }
      else if (nsingle == 0  &&  ndouble == 0  &&  ntriple == 0  && naromatic == 4  &&
               charge  == 0  &&  radical == 0                                       &&
               symbol.equals("S")                                                   &&
               (hcount == Atom.DEFAULT_HCOUNT  ||  hcount == 0))
      {  // sulfuric acid derivatives
         this.dbl_min = 2; this.dbl_max = 2; this.is_open = true;
      }
      else if (hcount == Atom.DEFAULT_HCOUNT)
      {
         this.dbl_min = ndouble + (naromatic / 2);
         this.dbl_max = ndouble + (naromatic / 2);
         this.is_open = true;
      }
      else
      {
         int tmp_hcount = PTable.implicitHydrogens(symbol,
                                                   nsingle, naromatic, ndouble, ntriple,
                                                   radical, charge);
         if (hcount >= 0  &&  tmp_hcount > hcount)
         {
            this.dbl_min = ndouble + (naromatic / 2) - (tmp_hcount-hcount);
            this.dbl_max = ndouble + (naromatic / 2) - (tmp_hcount-hcount);
         }
         else
         {
            this.dbl_min = ndouble + (naromatic / 2);
            this.dbl_max = ndouble + (naromatic / 2);
         }
      }
   }


   public Object clone()
   {
      AtomConstraint result = new AtomConstraint(this);
      return (result);
   }
}
