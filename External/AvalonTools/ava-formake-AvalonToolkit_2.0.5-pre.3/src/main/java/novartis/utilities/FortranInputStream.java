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

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * Line buffered input class. Also contains utilities to ease column
 * based format conversion as done in Fortran.
 */
public class FortranInputStream
{
   public    String buffer = null;
   int              pos    = 0;
   BufferedReader   in     = null;

   /**
    * Create from a DataInput stream in. The class only uses the
    * readLine method of DataInputStream to access the data of in.
    */
   public FortranInputStream(BufferedReader in) throws IOException
   {
      this.in = in;
      buffer = in.readLine();
   }

   /**
    * Fetch next input record.
    */
   public void getBuffer() throws IOException
   {
      buffer = in.readLine();
      pos = 0;
   }

   public String retrieveBuffer() throws IOException
   {
   buffer = in.readLine();
   return(buffer); 
   }

   /**
    * Reset character position in line buffer.
    */
   public void resetPos()
   {
      pos = 0;
   }

   /**
    * Fortran X-format.
    */
   public void x(int width)
   {
      pos += width;
      if (pos > buffer.length()) pos = buffer.length();
   }

   /**
    * Fortran A-Format.
    *
    * Returns the String value represented by the next
    * width characters of this.buffer starting at this.pos.
    *
    * Trailing blanks are stripped.
    */
   public String a(int width)
   {
      StringBuffer result = new StringBuffer();
      int last_nonblank;
      int i;

      last_nonblank = pos-1;
      for (i=pos; i<buffer.length()  &&  i-pos < width; i++)
      {
         result.append(buffer.charAt(i));
         if (buffer.charAt(i) != ' ') last_nonblank = i;
      }

      result.setLength(last_nonblank - pos + 1);
      pos = i;
      return (result.toString());
   }

   /**
    * Fortran F-Format.
    *
    * Returns the double value represented by the next
    * width characters of this.buffer starting at this.pos.
    */
   public double f(int width)
   {
      StringBuffer str = new StringBuffer(this.a(width));
      if (str.length() == 0) return (0);

      return (Double.valueOf(str.toString()).doubleValue());
   }

   /**
    * Fortran I-Format.
    *
    * Returns the int value represented by the next
    * width characters of this.buffer starting at this.pos.
    */
   public int i(int width)
   {
      StringBuffer str = new StringBuffer(this.a(width));
      if (str.length() == 0) return (0);

      // Blank == 0 in Fortran
      if (str.toString().trim().equals("")) str.setCharAt(str.length()-1, '0');

      return (Integer.parseInt(str.toString().trim()));
   }

   public static void main(String[] argv) throws IOException
   {
      FortranInputStream in =
         new FortranInputStream(
            new BufferedReader(
               new FileReader("tmp.mol")));
      System.err.println("'"+in.buffer+"'");
      in.getBuffer();
      System.err.println("'"+in.buffer+"'");
      in.getBuffer();
      System.err.println("'"+in.buffer+"'");
      in.getBuffer();
      System.err.println(in.i(3));
      System.err.println(in.i(3));
      in.getBuffer();
      System.err.println(in.f(10));
      System.err.println(in.f(10));
      System.err.println(in.f(10));
      in.x(1);
      System.err.println(in.a(3));
   }
}
