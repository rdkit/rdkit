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
package novartis.combinatorics;

import java.io.IOException;
import java.io.StreamTokenizer;
import java.io.StringReader;

/**
 * This class implements the standard group theoretic operations on
 * parmutations. It also implements creating permuations from their
 * string representations, either from permuted numbers or in cycle
 * notation.
 *
 * Permutations operate on the numbers 1..n.
 *
 * Note: The implementation allocates n+1 integers, ignoring perm[0].
 *       This is done to make the code more readable.
 *
 * Permutation objects are always compatible even if they operate on
 * different size basis sets {1..n1} and {1..n2}. The resulting permutation
 * is always promoted to the bigger basis set by assuming stabilization for
 * the missing members.
 */
public class Permutation extends Object
{
   final static int MAXN = 16000;   // maximum size of permutation base
   final static int ident_perm[] = {0,  // element zero is ignored!!
                                    1};

   public static final Permutation identity = new Permutation(ident_perm, 1);

   protected int[] perm;          // vector of permutation images
   protected int   n;             // permutation operates on {1..n}

   /**
    * Create Permutation from integer array.
    */
   Permutation(int[] perm, int n)
   {
      this.perm = new int[n+1];

      this.n = n;

      for(int i=1; i<=n; i++)
         this.perm[i] = perm[i];
   }

   /**
    * Computes the inverse of this permutation.
    */
   Permutation inverse()
   {
      int perm[] = new int[n+1];

      for (int i=1; i<=n; i++)
         perm[this.perm[i]] = i;

      return (new Permutation(perm, this.n));
   }

   /**
    * Creates a copy of this permutation.
    */
   public Object clone()
   {
      return (new Permutation(perm, n));
   }

   /**
    * Create Permutation from string. The constructor assumes lists of mapped integers
    * as the representation.
    */
   public Permutation(String s)
   {
      Permutation tperm;

      int[] tmp, result;
      int tsize, rsize, tok;
      StreamTokenizer st = new StreamTokenizer(new StringReader(s));
      st.parseNumbers();

      result = new int[16];
      for (int i=1; i<16; i++) result[i] = i;
      rsize = 0;

      try
      {
         tok = st.nextToken();

         if (tok == StreamTokenizer.TT_NUMBER)
         {
            int i = 1;
            do
            {
               if (tok == '_') // treat as blank => NOP
               {
                  tok = st.nextToken();
                  continue;
               }

               if (st.nval > MAXN)
                  throw new IllegalArgumentException("Permutation: Illegal token '" + st.sval + "' found");
               else if (i >= result.length  ||  st.nval >= result.length)
               {
                  tmp = result;
                  result = new int[Math.max(i+16, (int)st.nval+16)];
                  for (int j=1; j<tmp.length; j++)             result[j] = tmp[j];
                  for (int j=tmp.length; j<result.length; j++) result[j] = j;
               }

               result[i] = (int)st.nval;
               if (rsize < i) rsize = i;
               if (rsize < st.nval) rsize = (int)st.nval;

               i++;
               tok = st.nextToken();
            } while (tok == StreamTokenizer.TT_NUMBER  ||  tok == '_');
         }
         else
            throw new IllegalArgumentException("Permutation: Illegal token '" + st.sval + "' found");
      }
      catch (IOException e)
      {
         throw new IllegalArgumentException("Permutation: IOException for '" + s + "'");
      }

      if (tok != StreamTokenizer.TT_EOF && tok != StreamTokenizer.TT_EOL)
         throw new IllegalArgumentException("Permutation: Illegal token '" + st.sval + "' found");

      perm = new int[rsize+1];
      n = rsize;
      for (int i=1; i<=rsize; i++) this.perm[i] = result[i];

      if (!permOK())
         throw new IllegalArgumentException("Permutation: Syntax check failed for '" + s + "'");
   }

   /**
    * Checks if this permutation is identical to p.
    *
    * Basis sets are adapted if needed.
    */
   public boolean equals(Object q)
   {
        if (!(q instanceof Permutation)) return false;
        Permutation p = (Permutation)q;
        int i;
        int minn = Math.min(n, p.n);

        for (i=1; i<=minn; i++)
            if (perm[i] != p.perm[i]) return (false);

        // check if extension is identity
        if (n > p.n)
        {
            for (i=minn+1; i<=n; i++) if (perm[i] != i) return (false);
        }
        else
        {
            for (i=minn+1; i<=p.n; i++) if (p.perm[i] != i) return (false);
        }

        return (true);
   }

    /**
     * Provide an equals-compatible hash code.
     */
    public int hashCode()
    {
        int result = this.getClass().getName().hashCode();
        result <<= 2; result ^= toString().hashCode();
        return result;
    }

   /**
    * Returns a string that contains the images of the integer {1..n} under
    * Permutation this.
    *
    * Suppresses identities at end of string. 
    */
   public String toString()
   {
      StringBuffer tmp;

      tmp = new StringBuffer(Integer.toString(perm[1]));
      int lastMismatch = 1;
      for (int i=2; i<=n; i++)
          if (i != perm[i]) lastMismatch = i;
      for (int i=2; i<=lastMismatch; i++)
         tmp.append(" "+perm[i]);

      return (new String(tmp));
   }

   /**
    * Returns the cycle representation of Permutation this. Identities are suppressed.
    */
   String toCycleString()
   {
      StringBuffer tmp;
      boolean visited[] = new boolean[n+1];
      int i, j;

      tmp = null;
      i=1;
      do
      {
         if (perm[i] == i  ||  visited[i])
            i++;
         else
         {
            if (tmp == null)
               tmp = new StringBuffer("(" + i);
            else
               tmp.append(" (" + i);
            visited[i] = true;

            j = perm[i];
            while (!visited[j])
            {
               visited[j] = true;
               tmp.append(" " + j);
               j = perm[j];
            };
            tmp.append(")");
            i++;
         }
      } while (i<n);

      if (tmp == null)  // identity
         tmp = new StringBuffer("()");

      // System.err.println("toCycleString(" + toString() + ") = " + tmp);
      return (new String(tmp));
   }

   /**
    * This method implements the right multiply of this with right_perm
    * as used in C. M. Hoffmann, Lecture Notes in Computer Science, Vol. 136.
    */
   public Permutation times(Permutation right_perm)
   {
      int maxn, tmp;
      maxn = Math.max(this.n, right_perm.n);
      int result_array[] = new int[maxn+1];

      for (int i=1; i<=maxn; i++)
      {
         if (this.n < i)   // extend right_perm
            tmp = i;
         else
            tmp = this.perm[i];

         if (right_perm.n < tmp)       // extend this
            result_array[i] = tmp;
         else
            result_array[i] = right_perm.perm[tmp];
      }

      return (new Permutation(result_array, maxn));
   }

   /**
    * This method implements the right multiply of this with right_perm
    * as shown in "dtv-Atlas zur Mathematik". There seems to be no standard
    * way of doing this.
    */
   public Permutation dtv_times(Permutation right_perm)
   {
      int maxn, tmp;
      maxn = Math.max(this.n, right_perm.n);
      int result_array[] = new int[maxn+1];

      for (int i=1; i<=maxn; i++)
      {
         if (right_perm.n < i)   // extend right_perm
            tmp = i;
         else
            tmp = right_perm.perm[i];

         if (this.n < tmp)       // extend this
            result_array[i] = tmp;
         else
            result_array[i] = this.perm[tmp];
      }

      return (new Permutation(result_array, maxn));
   }

   /**
    * Tests this permutation for internal consistency. Returns true if checks were passed.
    */
   public boolean permOK()
   {
      int counts[] = new int[n+1];
      int i;

      for (i=1; i<=n; i++)
         if (perm[i] > n) return (false);
         else             counts[perm[i]]++;

      for (i=1; i<=n; i++)
         if (counts[i] != 1) return (false);

      return (true);
   }


   /**
    * Test driver program.
    */
   public static void main(String argv[])
   {
      Permutation p1, p2;

      if (argv.length == 2)
      {
         System.out.println(argv[0]);
         System.out.println(argv[1]);
         System.out.println("");
         p1 = new Permutation(argv[0]);
         p2 = new Permutation(argv[1]);
         System.out.println(p1);
         System.out.println(p2);
         System.out.println(p1.times(p2));
         System.out.println("");
         System.out.println(p1.toCycleString());
         System.out.println(p2.toCycleString());
         System.out.println(p1.times(p2).toCycleString());
      }
      else
         System.out.println("sample usage: java Permutation '1 2 3' '2 3 1'");
   }
}
