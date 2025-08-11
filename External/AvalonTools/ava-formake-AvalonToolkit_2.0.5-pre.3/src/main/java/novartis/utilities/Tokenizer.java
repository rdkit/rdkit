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

/**
 * This class tokenizes its input into tokens from some set of
 * predefined token classes. Currently, it works on String objects,
 * but it could be extended to run also from Streams.
 *
 * This is the basis of the SMILES parser in the Novartis.Chemistry package.
 */
public class Tokenizer
{
   String source = null; // input string
   int pos = 0;          // start of next token

   public String sval = "";
   public int    ival = 0;

   public int ttype  = BOF;
   public final static int BOF            =   0;  // state after constructor was called
   public final static int SPECIAL_IDENT  =   1;  // identifier starting with upper case
   public final static int IDENT          =   2;  // one from a list of identifiers, e.g. 'c', 'n', ..
   public final static int RING_NUMBER    =   4;  // single digit or %<digit><digit>
   public final static int UNSIGNED       =   8;  // integer number
   public final static int ATTACHMENT     =  16;  // '&<digit>' string
   public final static int BRACE_STRING   =  32;  // string enclosed in '{}', e.g. {c1ccccc1&1}
   public final static int BRACKET_STRING =  64;  // string enclosed in '[]', e.g. [NH2] or [Ala;Gly]
   public final static int CHARACTER      = 128;  // any character, kind of a catch all
   public final static int UCIDENT        = 256;  // identifier made of upper case characters
   public final static int EOF            = 512;  // string was exhausted
   public final static int ERROR          =  -1;  // an error or unexpected character was detected.

   /**
    * Peek at character i+pos of source.
    * Returns (-1) if beyond the end of the source string.
    */
   private final char charAt(int i)
   {
      if (pos+i >= source.length()) return ('\0');
      else                          return (source.charAt(pos+i));
   }

   /**
    * Create and initialize the Tokenizer from a String in.
    */
   public Tokenizer(String source)
   {
      this.source = source;
      this.pos    = 0;
      this.ttype  = BOF;
      this.sval   = "";
      this.ival   = 0;
   }

   /**
    * Retrieve the next token from source. expecting is a bit mask telling
    * which token types are expected. specials contains special symbols
    * if one of those is expected.
    */
   public int nextToken(int expecting, String[] specials)
   {
      sval = ""; ival = 0;

      if (pos >= source.length())
      {
         ttype = EOF;
         return (EOF);
      }

      if (0 != (SPECIAL_IDENT  & expecting))
      {
         if (specials == null) return (ERROR);

         for (int i=0; i<specials.length; i++)
            if (charAt(0) == specials[i].charAt(0)  &&
                source.substring(pos).startsWith(specials[i]))
            {
               sval = specials[i];
               pos += specials[i].length();
               return (ttype = SPECIAL_IDENT);
            }
      }
      if (0 != (UCIDENT          & expecting)  &&
          Character.isUpperCase(charAt(0)))
      {
         StringBuffer result = new StringBuffer();
         result.append(charAt(0));
         pos++;
         while (Character.isUpperCase(charAt(0)))
         {
            result.append(charAt(0));
            pos++;
         }
         sval = result.toString();
         return (ttype = UCIDENT);
      }
      else if (0 != (IDENT          & expecting)  &&
               Character.isUpperCase(charAt(0)))
      {
         StringBuffer result = new StringBuffer();
         result.append(charAt(0));
         pos++;
         while (Character.isLowerCase(charAt(0)))
         {
            result.append(charAt(0));
            pos++;
         }
         sval = result.toString();
         return (ttype = IDENT);
      }
      if (0 != (RING_NUMBER & expecting)  &&
          (charAt(0) == '%'  &&  Character.isDigit(charAt(1))   ||
           Character.isDigit(charAt(0))))
      {
         if (charAt(0) == '%')
         {
            pos++;
            if (Character.isDigit(charAt(1)))
            {
               ival = 10*(charAt(0)-'0') + (charAt(1)-'0');
               pos += 2;
            }
            else
            {
               ival = (charAt(0)-'0');
               pos += 1;
            }
         }
         else
         {
            ival = (charAt(0)-'0');
            pos += 1;
         }
         sval = Integer.toString(ival);
         return (ttype = RING_NUMBER);
      }
      if (0 != (UNSIGNED & expecting)  &&
          Character.isDigit(charAt(0)))
      {
         ival = charAt(0) - '0';
         pos++;
         while (Character.isDigit(charAt(0)))
         {
            ival *= 10; ival += charAt(0) - '0'; pos++;
         }
         sval = Integer.toString(ival);
         return (ttype = UNSIGNED);
      }
      if (0 != (ATTACHMENT & expecting)  &&
          charAt(0) == '&'               &&
          Character.isDigit(charAt(1)))
      {
         ival = charAt(1)-'0';
         pos += 2;
         sval = Integer.toString(ival);
         return (ttype = ATTACHMENT);
      }
      if (0 != (BRACE_STRING & expecting)  &&
          charAt(0) == '{')
      {
         int ibrace = 1;
         StringBuffer result = new StringBuffer();
         pos++;

         while (pos < source.length())
         {
            if (charAt(0) == '}') ibrace--;
            if (charAt(0) == '{') ibrace++;
            if (ibrace == 0) break;
            result.append(charAt(0));
            pos++;
         }
         if (charAt(0) != '}')
         {
            pos++;
            return (ERROR);
         }
         else
         {
            pos++;
            sval = result.toString();
            return (ttype = BRACE_STRING);
         }
      }
      if (0 != (BRACKET_STRING & expecting)  &&
          charAt(0) == '[')
      {
         int ibracket = 1;
         StringBuffer result = new StringBuffer();
         pos++;

         while (pos < source.length())
         {
            if (charAt(0) == ']') ibracket--;
            if (charAt(0) == '[') ibracket++;
            if (ibracket == 0) break;
            result.append(charAt(0));
            pos++;
         }
         if (charAt(0) != ']')
         {
            pos++;
            return (ERROR);
         }
         else
         {
            pos++;
            sval = result.toString();
            return (ttype = BRACKET_STRING);
         }
      }
      if (0 != (CHARACTER & expecting))
      {
         sval = source.substring(pos, pos+1);
         pos++;
         return (ttype = CHARACTER);
      }
      return (ttype = ERROR);
   }

   public static void main(String argv[])
   {
      Tokenizer t = new Tokenizer("c1ccccc1[NH2]");
      //Tokenizer t = new Tokenizer("C");

      while (true)
      {
         String specials[] = {"c", "n", "o", "s", "p"};
         t.nextToken(IDENT      | SPECIAL_IDENT | UNSIGNED |
                     ATTACHMENT | RING_NUMBER   | BRACKET_STRING | CHARACTER,
                     specials);

         System.err.println("ttype = " + t.ttype +
                         ",  sval = '" + t.sval +
                         "', ival = "  + t.ival);
         if (t.ttype == EOF  ||  t.ttype == ERROR) break;
      }
   }
}
