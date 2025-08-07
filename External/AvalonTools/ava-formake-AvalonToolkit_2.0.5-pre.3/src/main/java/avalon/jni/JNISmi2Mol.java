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
package avalon.jni;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;


import avalon.tools.LibraryToolbox;

/**
 * The class JNISmi2Mol provides wrappers for the non-Windows functions
 * in DEPICTUTIL.C, which implements chemical structure depiction and
 * conversion fom SMILES to MOL fileand back.
 *
 * It is intended to be used from Java applications and trusted applets.
 * @History:    B. Rohde 30-NOV-2001
 */
public class JNISmi2Mol
{
   static JNISmi2Mol smi2mol = null;

    static String tmpMolFileName = null;

    static  // make sure we get a unique file that is cleared on exit
    {
        try
        {
            File tmpFile = File.createTempFile("JNI", ".mol");
            tmpFile.deleteOnExit();
            tmpMolFileName = tmpFile.getCanonicalPath();
        }
        catch (IOException ioe)
        {
            ioe.printStackTrace();
            System.exit(1);
        }
    }

   private JNISmi2Mol()
   {
   }

   private static void init()
   {
     if (smi2mol != null) return;
     loadLibraries();
     smi2mol = new JNISmi2Mol();
   }

   public static JNISmi2Mol getSmi2Mol()
   {
      init();
      return smi2mol;
   }

   static boolean loadFromDirectory(String dir)
   {
      String JNISmi2MolLib = dir + "/avalon_jni_JNISmi2Mol.dll";
      if (System.getProperty("os.name").toLowerCase().indexOf("windows") < 0)
          JNISmi2MolLib = dir + "/libJNISmi2Mol.so";

      System.err.println("trying to load '" + JNISmi2MolLib + "'");
      if ((new File(JNISmi2MolLib)).exists())
      {
          try
          {
             System.load(JNISmi2MolLib);
             return true;
          }
          catch (Exception e)
          {
             System.err.println("Can't load native library '" + JNISmi2MolLib + "'!");
             e.printStackTrace();
          }
          catch (Error err)
          {
            System.err.println("Can't load native library '" + JNISmi2MolLib + "'!");
            err.printStackTrace();
          }
      }
      return false;
   }

   static void loadLibraries()
   {
     // Make sure that we get the links to the needed binary code.
     String libName = "JNISmi2Mol";
     if (System.getProperty("os.name").toLowerCase().indexOf("windows") >= 0) libName ="avalon_jni_JNISmi2Mol";
     if (LibraryToolbox.loadLibrary(JNISmi2Mol.class, libName)) return;

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

      if (loadFromDirectory(System.getProperty("user.dir"))) return;

      String JNISmi2MolLib = null;
      try
      {
          System.err.println("looking up bundle");
          java.util.ResourceBundle bundle =
              java.util.ResourceBundle.getBundle("libraries");

          System.err.println("trying to load key '" + "jnismi2mol" +
                             "' from bundle " + bundle);
          JNISmi2MolLib = bundle.getString("jnismi2mol");
          System.err.println("Loading library " + JNISmi2MolLib);

         // load from name provided in resource bundle
         // This must work eventually and is the method to be used on Unix/Linux 
         System.load(JNISmi2MolLib);
      }
      catch (Exception e)
      {
         e.printStackTrace();
         LibraryToolbox.showError(
                 "Error loading "+JNISmi2MolLib+"! \n\n"+
                 "Avalon will exit now.\n",
                 "Configuration Error", false);
        System.exit(1);
      }
      catch (Error err)
      {
        err.printStackTrace();
        System.err.println("Can't load native library!");
        LibraryToolbox.showError(
             "Error loading "+JNISmi2MolLib+"! \n\n"+
             "Avalon will exit now.\n",
             "Configuration Error", false);
        System.exit(1);
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

   private void assureTmpMoleFileIsWritable(String tmpMolFileName)
   {
      try {
        FileWriter f = new FileWriter(tmpMolFileName);
        f.close();
      } catch (IOException e) {
        System.err.println("Error opening "+tmpMolFileName+" for writing! "+e);
        LibraryToolbox.showError(
                 "Error opening "+tmpMolFileName+" for writing! "+
                 "Avalon will exit now.",
                 "File Error", false);
        System.exit(1);
      }
      //System.err.println(tmpMolFileName+" writable!");
   }

   /**
    * Transforms a structure represented as smiles into MOL format.
    *
    * @param smiles - String with smiles according to Daylight convention 
    * @return String containing the connection table for the structure,
    *         represented by input smiles or null if input smiles is invalid.
    * @throws IOException - the resulting connection table is temporarily
    *          written to a MOL file. If writing is not successfull,
    *          an IOException is thrown.
    * Convenience method forwarding to smiToMOLWithFlags().
    */
   public synchronized String smiToMOL(String smiles) throws IOException
   {
       return smiToMOLWithFlags(smiles, NO_FLAGS);
   }

   public static int NO_FLAGS = 0;
   public final static int APPLY_SHORTCUTS        = 0x000000FF;
   public final static int APPLY_AMINO_ACIDS      = 0x00000001;
   public final static int APPLY_PROTECTING_GROUPS= 0x00000002;
   public final static int APPLY_RNA_SHORTCUTS    = 0x00000004;
   public final static int EXTENDED_SHORTCUTS     = 0x00000100;
   public final static int NON_STANDARD_SHORTCUTS = 0x00000200;
   public final static int CATCH_ALL_SHORTCUTS    = 0x00008000;

   public final static int EXPECT_SMARTS          = 0x00010000;

   /**
    * Transforms a structure represented as smiles into MOL format.
    *
    * @param smiles - String with smiles according to Daylight convention 
    * @param flags - int with flags indicating some additional processing.
    *                Currently, only APPLY_SHORTCUTS is defined.
    * @return String containing the connection table for the structure,
    *         represented by input smiles or null if input smiles is invalid.
    * @throws IOException - the resulting connection table is temporarily
    *          written to a MOL file. If writing is not successfull,
    *          an IOException is thrown.
    */
   public synchronized String smiToMOLWithFlags(String smiles, int flags) throws IOException

   {
      int size;
      if (smiles == null  ||  smiles.trim().equals("")) return null;
      assureTmpMoleFileIsWritable(tmpMolFileName);
      // Call native convertion to write actual MOL file.
      smilesToMOLFileWithFlags(smiles, tmpMolFileName, flags);
      // Read the file into a StringBuffer.
      BufferedReader r = new BufferedReader(new FileReader(tmpMolFileName));
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
    *
    * The MOL file generated from the SMILES has the atoms highlighted that match
    * to an atom in the template MOL file tplCT.
    *
    * The MOL file is not highlighted of the template does not match and the method returns
    * null on error. 
    */
   public synchronized String smiToMOLWithTemplate(String smiles, String tplCT)
      throws IOException
   {
      int size;
      if (smiles == null  ||  smiles.trim().equals("")) return null;
      assureTmpMoleFileIsWritable(tmpMolFileName);
      // Call native conversion to write actual MOL file.
      smilesToMOLFileWithTemplate(smiles, tmpMolFileName, tplCT);
      // Read the file into a StringBuffer.
      BufferedReader r = new BufferedReader(new FileReader(tmpMolFileName));
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
    * This function writes the MOL file corresponding to smiles to the
    * file fname.
    *
    * It returns 0 if the operation was successful and non-0 otherwise.
    */
   public synchronized native int smilesToMOLFile(String smiles, String fname);

   /**
    * This function writes the MOL file corresponding to smiles to the
    * file fname using the options provided as the flags bitmap.
    *
    * It returns 0 if the operation was successful and non-0 otherwise.
    */
   public synchronized native int smilesToMOLFileWithFlags(String smiles, String fname, int flags);

   /**
    * This function writes the MOL file corresponding to smiles to the
    * file fname using the connection table in tplCT to clean it.
    *
    * It returns 0 if the operation was successful and non-0 otherwise.
    */
   public synchronized native int smilesToMOLFileWithTemplate(String smiles,
                                                              String fname,
                                                              String tplCT);

   /**
    * This method tries to detect and fix some common problems in input CTs.
    *
    * This includes superflous lines and wrong line endings.
    */
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
    *
    * Uses default temporary MOL file. 
    * 
    * If connection table is invalid, returns null
    */
   public synchronized String getSMILESFromCT(String ct)
   {
       return getSMILESFromCT(ct,tmpMolFileName);
   }

   /**
    * Convert a MOL file formatted connection table <code>ct</code> into
    * the corresponding SMILES.
    */
   public synchronized String getSMILESFromCT(String ct, String tmpFilename)
   {
      assureTmpMoleFileIsWritable(tmpFilename);

      try
      {
         // First, we write the clensed ct to a temporary file.
         ct = fixCT(ct);
         PrintWriter w = new PrintWriter(new FileWriter(tmpFilename));
         w.print(ct);
         w.close();

         // Now, we do the real thing
         byte buffer[] = new byte[1000];
         String smiles;
         int size;

         if (0 < (size = MOLFileToSmilesBytes(buffer, tmpFilename)))
            smiles = new String(buffer, 0, size);
         else
         {
            // try again with enlarged buffer
            int requiredSize = -size;
            if (requiredSize > 999  &&  requiredSize < 64000)
            {
                buffer = new byte[10+requiredSize];
                size = MOLFileToSmilesBytes(buffer, tmpFilename);
            }
            if (0 < size)
                smiles = new String(buffer, 0, size);
            else
                smiles = null;
         }
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
      assureTmpMoleFileIsWritable(tmpMolFileName);

      try
      {
         // First, we write the clensed ct to a temporary file.
         ct = fixCT(ct);
         PrintWriter w = new PrintWriter(new FileWriter(tmpMolFileName));
         w.print(ct);
         w.close();

         // Now, we do the real thing
         byte buffer[] = new byte[2000];        // SMARTS are not expected to be really large
         String smarts;
         int size;

         if (0 < (size = MOLFileToSMARTSBytes(buffer, tmpMolFileName)))
            smarts = new String(buffer, 0, size);
         else
            smarts = null;
System.err.println("SMARTS = " + smarts);
System.err.println("size = " + size);
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
       String smiles = "CCO";
       System.err.println("MW("+smiles+") = " +
           JNISmi2Mol.getSmi2Mol().MolecularWeightFromSmiles("CCO"));
   }
}
