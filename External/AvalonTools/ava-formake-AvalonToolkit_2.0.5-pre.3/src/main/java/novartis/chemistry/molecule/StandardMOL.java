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

import novartis.utilities.FortranInputStream;
import java.io.*;

/**
 * Class representing a standard MDL MOL-file structure.
 *
 * It extends Novartis.Chemistry.Molecule by adding the header information.
 */
public class StandardMOL extends Molecule
{
   String molname;

   String user_initials;   // only 2 characters are used
   String program;
   String date;            // placeholder for a date string of the form MMDDYYhhmm
   int    dimensionality;  // 2 for 2D and 3 for 3D
   int    int_scale;       // defined but not used
   double double_scale;    // defined but not used
   double energy;          // normally not used
   int    regno;           // numeric registry number

   String comment;

   public StandardMOL()
   {
      super();

      molname = "";
      user_initials = "";
      program = "";
      date = "";
      dimensionality = 2;
      int_scale = 1;
      double_scale = 1.0;
      energy = 0.0;
      regno = 0;
   }

   public void readMOLFile(FortranInputStream in) throws IOException
   {
      this.modifyOn();

      if (in.buffer.startsWith(">") && in.buffer.indexOf('<') >=0)
      {
         return;
      }

      molname = in.buffer; in.getBuffer();

      user_initials = in.a(2);
      program = in.a(8);
      date    = in.a(10);
      String dimtmp = in.a(2);
      if (dimtmp.equals("3D")) dimensionality = 3;
      else                     dimensionality = 2;
      try
      {
          int_scale    = in.i(2);
          double_scale = in.f(10);
          energy       = in.f(12);
          regno        = in.i(6);
      }
      catch (Throwable t)
      {
          System.err.println(t); 
          // ignore parsing errors
      }
      if (int_scale    == 0)   int_scale = 1;
      if (double_scale == 0.0) double_scale = 1.0;
      in.getBuffer();

      comment = in.buffer;
      in.getBuffer();

      this.readMDLCTable(in);

      this.modifyOff();
   }

   public static void main(String[] argv) throws IOException
   {
      FortranInputStream in =
         new FortranInputStream(
            new BufferedReader(
               new FileReader(argv[0])));
      StandardMOL m = new StandardMOL();
      m.readMOLFile(in);
      String[] d = new MoleculeDepicter().computeDepiction((Molecule)m,
                                                          0, 0,
                                                          MoleculeDepicter.USE_COLORS,
                                                          (int [])null,
                                                          (String[])null);
      for (int i=0; i<d.length; i++)
         System.err.println(d[i]);
   }
}
