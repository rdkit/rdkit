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

import java.util.BitSet;

/**
 * Adds a cardinality member function to java.util.BitSet.
 */
public class CountedBitSet
{
   int cardinality = 0;  // this member is updated whenever it might change.
   BitSet set;

   /**
    * Create a bit set with a quick cardinality.
    */
   public CountedBitSet()
   {
      set = new BitSet();
      cardinality = 0;
   }

   private int countBits()
   {
      int count = 0;

      for (int i=0; i<=set.size(); i++)
         if (set.get(i)) count++;

      return (count);
   }

   public CountedBitSet and(CountedBitSet set)
   {
      CountedBitSet result;

      if (this.set.size() > set.set.size())
      {
         result = (CountedBitSet)this.clone();
         result.set.and(set.set);
      }
      else
      {
         result = (CountedBitSet)set.clone();
         result.set.and(this.set);
      }
      result.cardinality = result.countBits();
      return (result);
   }

   public CountedBitSet or(CountedBitSet set)
   {
      CountedBitSet result;

      if (this.set.size() > set.set.size())
      {
         result = (CountedBitSet)this.clone();
         result.set.or(set.set);
      }
      else
      {
         result = (CountedBitSet)set.clone();
         result.set.or(this.set);
      }
      result.cardinality = result.countBits();
      return (result);
   }

   public CountedBitSet xor(CountedBitSet set)
   {
      CountedBitSet result;

      if (this.set.size() > set.set.size())
      {
         result = (CountedBitSet)this.clone();
         result.set.xor(set.set);
      }
      else
      {
         result = (CountedBitSet)set.clone();
         result.set.xor(this.set);
      }
      result.cardinality = result.countBits();
      return (result);
   }

   public void clear(int bit)
   {
      this.set.clear(bit);
      cardinality = countBits();
   }

   public void set(int bit)
   {
      this.set.set(bit);
      cardinality = countBits();
   }

   public boolean get(int bit)
   {
      return (this.set.get(bit));
   }

   public Object clone()
   {
      CountedBitSet bs = new CountedBitSet();

      bs.set = (BitSet)set.clone();
      bs.cardinality = cardinality;

      /*
      for (int i=0; i<set.size(); i++)
         if (set.get(i))
         {
            //bs.set.set(i);
            bs.cardinality++;
         }
      */

      return (bs);
   }

   public int getCardinality()
   {
      return (cardinality);
   }
}
