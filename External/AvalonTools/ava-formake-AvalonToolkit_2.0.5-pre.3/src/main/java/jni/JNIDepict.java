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
package jni;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;


import avalon.tools.LibraryToolbox;

/**
 * The class JNIDepict provides wrappers for functions to access the
 * DEPICT C-code.
 *
 * It is intended to be used from Java applications and trusted applets.
 * @History:    H. Plakties 22-Dec-2000 added checking for writability of
 *              tmpMolFile.
 */
public class JNIDepict
{

   static String tmpMolFilename = System.getProperty("java.io.tmpdir") + System.getProperty("file.separator") + "tmp.mol";
   static JNIDepict depictor = null;

   private JNIDepict()
   {
   }

   static
   {
     loadLibraries();
     depictor = new JNIDepict();
   }

   public static JNIDepict getDepictor()
   {
      return depictor;
   }

   static void loadLibraries()
   {
	  String libName = "JNIDepict";
      // try to load from directory next to source of this class
      if (LibraryToolbox.loadLibrary(JNIDepict.class, libName)) return;
	  // use system defined conventions to load library (possibly distributed via WebStart in a jar file)
      try {
    	  System.loadLibrary(libName);
		  return;
	  } catch (SecurityException e) {
    	  System.err.println(libName + ": Security exception, can't load native library from class path, "+e.getMessage());
      } catch (Throwable e) {
    	  System.err.println(libName + ": Can't load native library via loadLibrary(), " + e.getMessage());
          (new Exception("TRACE")).printStackTrace();
	  }
	  
      // Make sure that we get the links to the needed DLLs.
      String JNIDepictDLL;

      JNIDepictDLL = System.getProperty("user.dir") + "/JNIDepict.dll";
      System.err.println("trying to load '" + JNIDepictDLL + "'");
      if ((new File(JNIDepictDLL)).exists())
      {
          try
          {
             System.load(JNIDepictDLL);
             return;
          }
          catch (Exception e)
          {
             System.err.println("Can't load native library '" +
                                JNIDepictDLL + "'!");
             e.printStackTrace();
          }
          catch (Error err)
          {
            System.err.println("Can't load native library '" +
                               JNIDepictDLL + "'!");
            err.printStackTrace();
          }
      }
      }

   /**
    * This method wraps the corresponding native method call to calculate
    * molecular weights.
    */
   public synchronized double MolecularWeightFromSmiles(String smiles)
   {
      if (smiles == null  ||  smiles.trim().equals("")) return 0.0;
      System.err.println("Entering wrapper mwFromSmiles");
      System.err.println("SMILES = " + smiles);
      double mw = mwFromSmiles(smiles);
      System.err.println("mw = " + mw);
      return (mw);
   }

   private void assureTmpMoleFileIsWritable(String tmpMolFilename)
   {
      try {
        FileWriter f = new FileWriter(tmpMolFilename);
        f.close();
      } catch (IOException e) {
        System.err.println("Error opening "+tmpMolFilename+" for writing! "+e);
        LibraryToolbox.showError(
                 "Error opening "+tmpMolFilename+" for writing! "+
                 "Avalon will exit now.",
                 "File Error", false);
        System.exit(1);
      }

   }

   public synchronized String smiToMOL(String smiles) throws IOException
   {
      int size;
      if (smiles == null  ||  smiles.trim().equals("")) return null;
      assureTmpMoleFileIsWritable(tmpMolFilename);
      // Call native convertion to write actual MOL file.
      smilesToMOLFile(smiles, tmpMolFilename);
      // Read the file into a StringBuffer.
      BufferedReader r = new BufferedReader(new FileReader(tmpMolFilename));
      StringBuffer mol = new StringBuffer();
      String line = null;
      while (null != (line = r.readLine()))
      {
         mol.append(line); mol.append('\n');
      }
      r.close();

      if (mol.length() < 10) return (null);	// Cannot be a real molecule.
      else                   return (mol.toString());
   }

   /**
    * Wrapper for template driven smilesToMOLFileWithTemplate.
    */
   public synchronized String smiToMOLWithTemplate(String smiles,
                                                   String tplCT)
      throws IOException
   {
      int size;
      if (smiles == null  ||  smiles.trim().equals("")) return null;
      assureTmpMoleFileIsWritable(tmpMolFilename);
      // Call native convertion to write actual MOL file.
      smilesToMOLFileWithTemplate(smiles, tmpMolFilename, tplCT);
      // Read the file into a StringBuffer.
      BufferedReader r = new BufferedReader(new FileReader(tmpMolFilename));
      StringBuffer mol = new StringBuffer();
      String line = null;
      while (null != (line = r.readLine()))
      {
         mol.append(line); mol.append('\n');
      }
      r.close();

      if (mol.length() < 10) return (null);	// Cannot be a real molecule.
      else                   return (mol.toString());
   }

   /**
    * Method to be used in substructure matching.
    */
   public synchronized boolean smilesMatchesQueryCT(String smiles,
                                                    String queryCT)
   {
      if (smiles == null)
          return false;
      if (queryCT == null  ||     // null query match everything
          queryCT.length() < 3)
          return true;
      return smilesMatchesQueryCTNative(smiles, queryCT);
   }

    /**
     * Tests if the given smiles matches queryCT.
     *
     * The string queryCT is coded as a MOLFile. It may be preprocessed
     * be better suited for aromaticity matching.
     */
    public synchronized native boolean smilesMatchesQueryCTNative
        (String smiles, String queryCT);

   /**
    * This function writes the MOL file corresponding to smiles to the
    * file fname.
    *
    * It returns 0 if the operation was successful and non-0 otherwise.
    */
   public synchronized native int smilesToMOLFile(String smiles, String fname);

   /**
    * This function writes the MOL file corresponding to smiles to the
    * file fname using the connection table in tplCT to clean it.
    *
    * It returns 0 if the operation was successful and non-0 otherwise.
    */
   public synchronized native int smilesToMOLFileWithTemplate(String smiles,
                                                              String fname,
                                                              String tplCT);

   private static String fixCT(String molfile)
   {
      // Do some checking and normalization
      if (molfile == null) return null;
      // terminate after "M  END" line, if any
      int mEndIndex = molfile.indexOf("M  END");
      if (mEndIndex > 0)
          molfile = molfile.substring(0, mEndIndex) + "M  END\n";
      if (molfile.trim().equals("")) return null;
      // make sure there are complete lines
      if (molfile.charAt(molfile.length()-1) != '\n')
          molfile += "\n";
      return molfile;
   }

   /**
    * Convert a MOL file formatted connection table <code>ct</code> into
    * the corresponding SMILES.
    */
   public synchronized String getSMILESFromCT(String ct)
   {
       return getSMILESFromCT(ct, tmpMolFilename);
   }

   /**
    * Convert a MOL file formatted connection table <code>ct</code> into
    * the corresponding SMILES.
    */
   public synchronized String getSMILESFromCT(String ct, String tmpFileName)
   {
      assureTmpMoleFileIsWritable(tmpFileName);

      try
      {
         // First, we write the clensed ct to a temporary file.
         ct = fixCT(ct);
         PrintWriter w = new PrintWriter(new FileWriter(tmpFileName));
         w.print(ct);
         w.close();

         // Now, we do the real thing
         byte buffer[] = new byte[1000];
         String smiles;
         int size;

         if (0 < (size = MOLFileToSmilesBytes(buffer, tmpFileName)))
            smiles = new String(buffer, 0, size);
         else
            smiles = null;
         return (smiles);
      }
      catch (Exception e)
      {
         e.printStackTrace();
      }
      return null;
   }

   /**
    * Convert a MOL file formatted connection table <code>ct</code> into
    * the corresponding SMARTS.
    */
   public synchronized String getSMARTSFromCT(String ct)
   {
      assureTmpMoleFileIsWritable(tmpMolFilename);

      try
      {
         // First, we write the clensed ct to a temporary file.
         ct = fixCT(ct);
         PrintWriter w = new PrintWriter(new FileWriter(tmpMolFilename));
         w.print(ct);
         w.close();

         // Now, we do the real thing
         byte buffer[] = new byte[1000];
         String smarts;
         int size;

         if (0 < (size = MOLFileToSMARTSBytes(buffer, tmpMolFilename)))
            smarts = new String(buffer, 0, size);
         else
            smarts = null;
         return (smarts);
      }
      catch (Exception e)
      {
         e.printStackTrace();
      }
      return null;
   }

   /**
    * This method returns the molecular weight corresponding to
    * the chemical structure that has the SMILES given as the parameter.
    * It uses the DEPICT32.DLL to actually do the calculations.
    */
   public synchronized native double mwFromSmiles(String smiles);

   /**
     * Converts a MOL file <code>fname</code> to the corresponding
     * SMILES.
     *
     * The method returns the used bytes of <code>buffer</code> or a negative
     * number (-n) that would indicate that n bytes would actually be needed
     * to store the result.
     */
    public synchronized native int MOLFileToSmilesBytes(byte buffer[],
                                                        String fname);

    /**
     * Converts a MOL file <code>fname</code> to the corresponding
     * SMILES.
     *
     * The method returns the used bytes of <code>buffer</code> or a negative
     * number (-n) that would indicate that n bytes would actually be needed
     * to store the result.
     */
    public synchronized native int MOLFileToSMARTSBytes(byte buffer[],
                                                        String fname);

    static public void main(String argv[])
        throws IOException
    {
        loadLibraries();

        JNIDepict jni = new JNIDepict();
        StringBuffer molfile = new StringBuffer();
        if (argv.length == 1) 
        {
            BufferedReader reader =
                new BufferedReader(new InputStreamReader(System.in));
            String line;
            while (null != (line = reader.readLine()))
            {
                molfile.append(line); molfile.append("\n");
            }
        }
        else
        {
            System.err.println("usage: java JNIDepict 'smiles' <molfile.mol");
            return;
        }
        BufferedReader reader = new BufferedReader(new FileReader(argv[0]));
        String line;
        while (null != (line = reader.readLine()))
        {
            System.err.println("smilesMatchesQueryCT() yields " +
                    jni.smilesMatchesQueryCT(line, molfile.toString()));
        }
        reader.close();
    }
}
