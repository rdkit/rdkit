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
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;

import avalon.tools.LibraryToolbox;

/**
 * The class JNIMatch provides wrappers for functions to access the
 * DEPICT C-code.
 *
 * It is intended to be used from Java applications and trusted applets.
 * @History:    H. Plakties 22-Dec-2000 added checking for writability of
 *              tmpMolFile.
 */
public class JNIMatch
{
   static JNIMatch matcher = null;

   public final static int USE_RING_PATTERN      = 0x000001;
   public final static int USE_RING_PATH         = 0x000002;
   public final static int USE_ATOM_SYMBOL_PATH  = 0x000004;
   public final static int USE_ATOM_CLASS_PATH   = 0x000008;
   public final static int USE_ATOM_COUNT        = 0x000010;
   public final static int USE_AUGMENTED_ATOM    = 0x000020;
   public final static int USE_HCOUNT_PATH       = 0x000040;
   public final static int USE_HCOUNT_CLASS_PATH = 0x000080;
   public final static int USE_HCOUNT_PAIR       = 0x000100;
   public final static int USE_BOND_PATH         = 0x000200;
   public final static int USE_AUGMENTED_BOND    = 0x000400;
   public final static int USE_RING_SIZE_COUNTS  = 0x000800;
   public final static int USE_DEGREE_PATH       = 0x001000;
   public final static int USE_CLASS_SPIDERS     = 0x002000;
   public final static int USE_FEATURE_PAIRS     = 0x004000;

   public final static int USE_ALL_BITS          = 0x007FFF;

   public final static int USE_SCAFFOLD_IDS      = 0x100000;
   public final static int USE_SCAFFOLD_COLORS   = 0x200000;
   public final static int USE_SCAFFOLD_LINKS    = 0x400000;

   public final static int USE_NON_SSS_BITS      = 0xF00000;

   static String[] patternNames =
   {
       "RING_PATTERN",
       "RING_PATH",
       "ATOM_SYMBOL_PATH",
       "ATOM_CLASS_PATH",
       "ATOM_COUNT",
       "AUGMENTED_ATOM",
       "HCOUNT_PATH",
       "HCOUNT_CLASS_PATH",
       "HCOUNT_PAIR",
       "BOND_PATH",
       "AUGMENTED_BOND",
       "RING_SIZE_COUNTS",
       "DEGREE_PATH",
       "CLASS_SPIDERS",
       "FEATURE_PAIRS",
       "",
       "","","","",
       "SCAFFOLD_IDS",
       "SCAFFOLD_COLORS",
       "SCAFFOLD_LINKS",
   };

   static String[] shortPatternNames =
   {
       "RPT",
       "RPH",
       "ASP",
       "ACP",
       "AC",
       "AA",
       "HPH",
       "HCP",
       "HP",
       "BP",
       "AB",
       "RSC",
       "DP",
       "CS",
       "FP",
       "",
       "","","","",
       "SI",
       "SC",
       "SL",
   };

   public static String maskToString(int mask)
   {
	   String result = "";
	   for (int i=0; i<patternNames.length; i++)
	   {
		   if (0 == ((1<<i)&mask)) continue;
		   if (result.length() > 0) result += "|";
		   result += patternNames[i];
	   }
	   return result;
   }

   public static String maskToShortString(int mask)
   {
	   String result = "";
	   for (int i=0; i<patternNames.length; i++)
	   {
		   if (0 == ((1<<i)&mask)) continue;
		   if (result.length() > 0) result += "|";
		   result += shortPatternNames[i];
	   }
	   return result;
   }

   public static int stringToMask(String s)
   {
	   int result = 0;
	   for (int i=0; i<patternNames.length; i++)
		   if (patternNames[i].length() > 0  &&
			   0 <= s.indexOf(patternNames[i])) result |= (1<<i);
	   return result;
   }

   // Constants defined in canonzer.h
   /* use double-bond stereochemistry */
   private final static int DB_STEREO     = 0x01;
   /* use tetrahedral center stereochemistry */
   private final static int CENTER_STEREO = 0x02;
   /* list every component only once */
   private final static int COMPONENT_SET = 0x08;
   /* only keep non-trivial pieces */
   private final static int MAIN_PART     = 0x10;
   /* only keep non-trivial pieces */
   private final static int FORCE_OCTET   = 0x20;

   // Surrogate constants to be used in Java code
   public final static int CANONIZE_AS_IS = FORCE_OCTET|
                                            DB_STEREO|CENTER_STEREO;
   public final static int TO_COMPONENT_SET = FORCE_OCTET|
                                              COMPONENT_SET|
                                              DB_STEREO|CENTER_STEREO;
   public final static int TO_MAIN_PARTS = FORCE_OCTET|
                                           COMPONENT_SET|MAIN_PART|
                                           DB_STEREO|CENTER_STEREO;
   public final static int TO_NON_STEREO_MAIN = FORCE_OCTET|
                                                COMPONENT_SET|MAIN_PART;

    /**
     * This function computes the canonical SMILES from inSmiles
     * processing it according to flags and putting the resulting bytes
     * into outBuffer.
     *
     * The method is implemented as a JNI wrapper to the code in canonizer.c.
     */
    public native
    synchronized int canonicalSmilesNative(String inSmiles,
                                           byte[] outBuffer,
                                           int flags);

    /**
     * Non-native wrapper to canonicalSmilesNative to ensure synchronization.
     *
     * It returns the canonical SMILES or null if something was wrong.
     */
    public synchronized String canonicalSmiles(String smiles, int flags)
    {
        // heapCheck();
        if (smiles == null) return null;
        byte buffer[] = new byte[smiles.length()*2];
        int canLength = canonicalSmilesNative(smiles, buffer, flags);
        if (canLength <= 0) return null;
        return new String(buffer, 0, canLength);
    }

   static int[] bits_per_byte = null;

   static
   {
      // Initialize number of bits for the bit patterns of each byte
      bits_per_byte = new int[256];
      for (int i=0; i<256; i++)
      {
         int n=0;
         for (int j=0; j<8; j++)
             if ((i&(1<<j)) != 0) n++;  // count bits
         bits_per_byte[i] = n;
     }
     loadLibraries();
     matcher = new JNIMatch();
   }

   private JNIMatch()
   {
   }

   /**
    * Static method that returns singleton instance of this class.
    *
    * This is necessary to synchronize access to certain resources.
    *
    * @return - singleton instance of JNIMatch
    */
   public static JNIMatch getMatcher()
   {
      return matcher;
   }

   static boolean loadFromDirectory(String dir)
   {
      String libfile = dir + "/JNIMatch.dll";
      System.err.println("trying to load '" + libfile + "'");
      try
      {
          if ((new File(libfile)).exists())
          {
             System.load(libfile);
             System.err.println("Loaded from file '" + libfile + "'");
             return true;
          }
          else
              return false;
      }
      catch (Exception e)
      {
          System.err.println("Can't load native library '" +
                             libfile + "'!");
          e.printStackTrace();
      }
      catch (Error err)
      {
          System.err.println("Can't load native library '" +
                             libfile + "'!");
          err.printStackTrace();
      }
      return false;
   }

   static boolean loadFromClasspath(String library)
   {
       String classpath = System.getProperty("java.class.path");
       // System.err.println("CLASSPATH = " + classpath);
       String separator = System.getProperty("path.separator");
       String fileSeparator = System.getProperty("file.separator");
       String[] tokens = avalon.tools.StringTools.stringToArray(classpath,separator.charAt(0));
       String oldLibFile = "#NONE#";
       for (int i=0; i<tokens.length; i++)
       {
           if (tokens[i].lastIndexOf(fileSeparator) < 0) continue;
           if (tokens[i].lastIndexOf(fileSeparator) > 0  &&
               (tokens[i].toLowerCase().endsWith(".jar")  ||  tokens[i].toLowerCase().endsWith(".zip")))
               tokens[i] = tokens[i].substring(0,tokens[i].lastIndexOf(fileSeparator));
    	   String libfile = tokens[i] + "/" + library;
           if (oldLibFile.equals(libfile)) continue;
           oldLibFile = libfile;
    	   // System.err.println("loadFromClasspath(): Trying to load " + libfile);
    	   try
    	   {
               // System.err.println("Checking file '" + libfile + "'");
    		   if (!(new File(libfile)).isFile()) continue;
               System.load(libfile);
               System.err.println("Loaded from file '" + libfile + "'");
    		   return true;
    	   }
    	   catch (Exception e)
    	   {
    	          System.err.println("Can't load native library '" +
    	                             libfile + "'!");
    	   }
    	   catch (Error err)
    	   {
    	       System.err.println("Can't load native library '" +
    	                          libfile + "'!");
    	   }
       }
       return false;
   }

   static synchronized void loadLibraries()
   {
      // Make sure that we get the links to the needed DLLs.
      String JNIMatchDLL = null;

      if (LibraryToolbox.loadLibrary(JNIMatch.class, "JNIMatch")) return;
      if (loadFromClasspath("JNIMatch.dll")) return;
      if (loadFromDirectory(System.getProperty("user.dir"))) return;

      String JNIMatchLib = null;
      try
      {
          java.util.ResourceBundle bundle =
              java.util.ResourceBundle.getBundle("libraries");
          JNIMatchLib = bundle.getString("jnimatch");
          System.err.println("loading from bundle: trying to load " + JNIMatchLib);
         // Note that the absolute path is required here.
         System.load(JNIMatchLib);
         System.err.println("Loaded from file '" + JNIMatchLib + "'");
         return;
      }
      catch (Throwable t)
      {
            System.err.println("Could not find " + JNIMatchLib);
      }

System.err.println("java.library.path: " + System.getProperty("java.library.path"));
      JNIMatchDLL = "JNIMatch";
System.err.println("JNIMatch(1): " + JNIMatchDLL);

      try
      {
         // Note that no absolute path is required here.
         System.loadLibrary(JNIMatchDLL);
      }
      catch (Exception e)
      {
         e.printStackTrace();
         LibraryToolbox.showError(
                 "Error loading "+JNIMatchDLL+"! \n\n"+
                 "Avalon will exit now.\n",
                 "Configuration Error", false);
         System.exit(1);
      }
      catch (Error err)
      {
        err.printStackTrace();
        System.err.println("Can't load native library!");
         LibraryToolbox.showError(
                 "Error loading "+JNIMatchDLL+"! \n\n"+
                 "Avalon will exit now.\n",
                 "Configuration Error", false);
        System.exit(1);
      }
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

   static int DEFAULT_SIZE = 64*1;
   /**
    * A JAVA wrapper for native method computeFingerprintFromCT, that computes
    * 512-bit fingerprint for a ConnectionTable 
    * @param ct - String representing a connection table according to MDL conventions
    *             if connection table is invalid, null is returned
    * @param asQuery - boolean; if true, the input connection table is treated
    *                  as query fragment, i.e., open valencies in connection table are 
    *                  left open; if false - the input connection table is treated as 
    *                  a whole molecule, with open valencies filled with hydrogens 
    * @param whichBits - a hexadecimal integer that masks certain bits in resulting 
    *                    fingerprint. The list of valid masking constants is given in
    *                    the beginning of this class. For normal purposes, all bits must be computed,
    *                    i.e., USE_ALL_BITS, or 0xFFFF, which results in the same fingerprints.
    * @return - returns byte array of DEFAULT_SIZE bytes, representing a 8*DEFAULT_SIZE-bit fingerprint,
    *           or an array of all zeros in case of error.
    */
   public synchronized byte[] getFingerprintFromCT(String ct,
                                                   boolean asQuery,
                                                   int whichBits)
   {
      byte[] fingerprint = new byte[DEFAULT_SIZE];
      matcher.computeFingerprintFromCT(ct, fingerprint, asQuery, whichBits);
      return fingerprint;
   }

   public synchronized byte[] getFingerprintFromCT(String ct,
                                                   boolean asQuery,
                                                   int whichBits,
                                                   int size)
   {
      byte[] fingerprint = new byte[size];
      matcher.computeFingerprintFromCT(ct, fingerprint, asQuery, whichBits);
      return fingerprint;
   }

   /**
	* A JAVA wrapper for native method computeFingerprintFromSmiles, that computes
	* 512-bit fingerprint for a ConnectionTable 
	* @param smiles - String representing smiles according to Daylight conventions
	*                 if smiles string is invalid, null is returned
	* @param asQuery - boolean; if true, the input connection table is treated
	*                  as query fragment, i.e., open valencies in connection table are 
	*                  left open; if false - the input connection table is treated as 
	*                  a whole molecule, with open valencies filled with hydrogens 
	* @param whichBits - a hexadecimal integer that masks certain bits in resulting 
	*                    fingerprint. The list of valid masking constants is given in
	*                    the beginning of this class. For normal purposes, all bits must be computed,
	*                    i.e., USE_ALL_BITS, or 0xFFFF, which results in the same fingerprints.
	* @return - returns byte array of DEFAULT_SIZE bytes, representing a 8*DEFAULT_SIZE-bit fingerprint,
	*           or null if smiles is invalid
	*/
   public synchronized byte[] getFingerprintFromSmiles(String smiles,
                                                       boolean asQuery,
                                                       int whichBits)
   {
      byte[] fingerprint = new byte[DEFAULT_SIZE];
      matcher.computeFingerprintFromSmiles(smiles, fingerprint, asQuery, whichBits);
      return fingerprint;
   }

   public synchronized byte[] getFingerprintFromSmiles(String smiles,
                                                       boolean asQuery,
                                                       int whichBits,
                                                       int size)
   {
      byte[] fingerprint = new byte[size];
      matcher.computeFingerprintFromSmiles(smiles, fingerprint, asQuery, whichBits);
      return fingerprint;
   }

    /**
     * Parses the given query SMILES smiles and returns a handle
     * (pointer) to the corresponding C data structure.
     *
     * The string smiles is may contain explicit hydrogens and R atoms
     * to block free sites.
     * It will be prepared for substructure matching.
     * 
     * Returns 0 if the connection table queryCT is invalid
     */
    public native
    synchronized long querySmilesToHandleNative(String smiles);

    /**
     * Java wrapper for querySmilesToHandle().
     */
    public synchronized long querySmilesToHandle(String smiles)
    {
        return querySmilesToHandleNative(smiles);
    }

    /**
     * Parses the given connection table queryCT and returns a handle
     * (pointer) to the corresponding C data structure.
     *
     * The string queryCT is coded as a MOLFile. It will be prepared for
     * substructure matching.
     * 
     * Returns 0 if the connection table queryCT is invalid
     */
    public native
    synchronized long queryCTToHandleNative(String queryCT);

    /**
     * Java wrapper for queryCTToHandleNative().
     */
    public synchronized long queryCTToHandle(String queryCT)
    {
        return queryCTToHandleNative(queryCT);
    }

    /**
     * Returns the connection table C data structure back to the
     * non-garbage collected heap.
     */
    public native
    synchronized void disposeCTHandleNative(long handle);

    /**
     * Java wrapper for disposeCTHandleNative().
     */
    public synchronized void disposeCTHandle(long handle)
    {
        disposeCTHandleNative(handle);
    }

    /**
     * Tests if the given smiles matches the query referred to by queryHandle.
     * 
     * String smiles represents smiles according to Daylight convention;
     * if smiles is invalid, false is returned
     */
    public native
    synchronized boolean smilesMatchesHandleNative(String smiles, long queryHandle);

    /**
     * Java wrapper for smilesMatchesHandleNative().
     */
    public synchronized boolean smilesMatchesHandle(String smiles,
                                                    long queryHandle)
    {
        return smilesMatchesHandleNative(smiles, queryHandle);
    }

    /**
     * Tests if the given smiles matches queryCT.
     *
     * The string queryCT is coded as a MOLFile. It may be preprocessed
     * be better suited for aromaticity matching.
     */
    public native
    synchronized boolean smilesMatchesQueryCTNative
        (String smiles, String queryCT);

    /**
     * Tests if the heap is still OK.
     */
    public synchronized native boolean heapCheck();

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
    * This function computes the fingerprint of the structure or
    * query represented by molCT and sets the corresponding bits
    * in fingerprint[].
    *
    * It returns the number of fingerprint components enumerated during
    * that process. The fingerprint is normally DEFAULT_SIZE bytes wide.
    */
   public native
   synchronized int computeFingerprintFromCT(String molCT,
                                             byte[] fingerprint,
                                             boolean asQuery,
                                             int whichBits);

   /**
    * This function sets the colors[] to non-zero if an atom is part of
    * a match of queryCT. Returns true if a match was found.
    *
    * Note: colors[] needs to be already allocated.
    */
   public native
   synchronized boolean computeMatchColorsForCT(String molCT,
                                      String queryCT,
                                      int[] colors);

   /**
    * This function sets the colors indicating the involvement of
    * the atoms in setting bits of fingerprint[].
    *
    * It updates colors[] and returns the number of atoms that are
    * required to set all the bits of fingerprint[].
    * Note: colors[] needs to be already allocated.
    */
   public native
   synchronized int computeFingerprintColorsForCT(String molCT,
                                            byte[] fingerprint,
                                            boolean asQuery,
                                            int whichBits,
                                            int[] colors);

   public synchronized int setFingerprintColorsForCT(String molCT,
                                                     byte[] fingerprint,
                                                     boolean asQuery,
                                                     int whichBits,
                                                     int[] colors)
   {
       return computeFingerprintColorsForCT(molCT, fingerprint, asQuery, whichBits, colors);
   }

    /**
     * One off testing procedure to be enables by editing the source.
     */
    private static void test(String smiles, String queryCT)
    {
        System.err.println("===== Testing ============");
        JNIMatch matcher = JNIMatch.getMatcher();
        boolean match = matcher.smilesMatchesQueryCT(smiles,queryCT);
        System.err.println("test() returns " + match);
        System.exit(1);
    }

   /**
    * This function computes the fingerprint of the structure
    * represented by smiles and sets the corresponding bits
    * in fingerprint[].
    *
    * It returns the number of fingerprint components enumerated during
    * that process. The fingerprint is normally DEFAULT_SIZE bytes wide.
    */
   public native
   synchronized int computeFingerprintFromSmiles(String smiles,
                                           byte[] fingerprint,
                                           boolean asQuery,
                                           int whichBits);

    /**
     * Computes an array of color indices with 0, i.e. black, for non-matching atoms
     * and matchColorIndex for matching atoms. 
     */
    public static int[] getMatchColoring(String molfile, String tpl, int matchColorIndex)
    {
        if (molfile == null) return null;
        // allocate array to be eventually passed to a JNI call
        int[] result = new int[molfile.split("\\r?\\n").length];
        int[] tmp = new int[result.length];
        matcher.computeMatchColorsForCT(molfile, tpl, tmp);
        for (int i=0; i<result.length; i++)
            if (tmp[i] != 0) result[i] = matchColorIndex;
            else             result[i] = 0;

        return result;
    }

    static public void main(String argv[])
        throws IOException
    {
        loadLibraries();

        JNIMatch matcher = JNIMatch.getMatcher();
        StringBuffer molfile = new StringBuffer();
        long queryHandle = 0;
        boolean doFingerprint=true;
        if (argv.length > 2) doFingerprint = false;
        // Test code for querySmilesToHandle
        if (argv.length > 3)
        {
            queryHandle = matcher.querySmilesToHandle(argv[3]);
        }
        if (argv.length >= 1) 
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
            System.err.println("usage: java jni.JNIMatch " +
                               "smiles-file label [no-detail] <query.mol");
            return;
        }
        if (argv.length > 1) System.err.println("Query = " +argv[1]);
        byte[] qfp = null;
        String line;
        int ntest= 0;
        int nfpmatch = 0;
        int nmatch = 0;
        int nconflict = 0;
        double totalBits = 0;
        double totalPassFail = 0;
        double totalMatch = 0;
        double nPassFail = 0;
        int byteSelectivity[] = null;
        byte qfpByFeatureTmp[][] = new byte[16][];
        byte qfpByFeature[][] = new byte[16][];
        byte[] emptyFP = null;
        if (molfile.length() > 10)
        {
            qfp = matcher.getFingerprintFromCT(molfile.toString(), true, USE_ALL_BITS);
            if (queryHandle == 0)
                queryHandle = matcher.queryCTToHandle(molfile.toString());
            byteSelectivity = new int[qfp.length];
            emptyFP = new byte[qfp.length];

            for (int i=0; i<patternNames.length  &&  i<16; i++)
            {
                int mask = 1<<i;
                qfpByFeature[i] =
                    matcher.getFingerprintFromCT(molfile.toString(), true, mask);
                mask = USE_ALL_BITS - (1<<i);
                qfpByFeatureTmp[i] =
                    matcher.getFingerprintFromCT(molfile.toString(), true, mask);
            }
        }
        BufferedReader reader = new BufferedReader(new FileReader(argv[0]));
        int nConflictByFeature[] = new int[16];
        double nFPMatchByFeature[] = new double[16];
        double nFPPassFailByFeature[] = new double[16];
        double bitSumByFeature[] = new double[16];
        double bitMatchSumByFeature[] = new double[16];
        double bitPassFailSumByFeature[] = new double[16];
        int modulus = 0;
        System.err.println("doFingerprint = " + doFingerprint);
        while (null != (line = reader.readLine()))
        {
            // test(line, molfile.toString());
            // if ((modulus++) % 20 != 0) continue;
            if (queryHandle == 0) // no query => test canonicalization
            {
                String newSmiles =
                    matcher.canonicalSmiles(line, JNIMatch.CANONIZE_AS_IS);
                if (!line.equalsIgnoreCase(newSmiles))
                    System.out.println(line + "\t" + newSmiles);
                ntest++;
                if (ntest%100 == 0) System.err.print('.');
                if (ntest%5000 == 0) System.err.println(" "+ntest);
                continue;
            }
            byte[] fp;
            if (doFingerprint)
                // fp = matcher.getFingerprintFromSmiles(line, false, 0xFFFF);
                fp = matcher.getFingerprintFromSmiles(line, false, JNIMatch.USE_ALL_BITS);
            else
                fp = emptyFP;
            int count = countBits(fp);
            totalBits += count;
            if (false && !matcher.heapCheck())
            {
                System.err.println(
                    "getFingerprintFromSmiles corrupted heap with smiles '"+
                    line + "'");
                break;
            }
            boolean fpmatch = fingerTest(fp,qfp);
            for (int i=0; i<qfp.length; i++)
               if ((fp[i]&qfp[i]) == qfp[i])
                  byteSelectivity[i]++;
            // boolean match = matcher.smilesMatchesQueryCT(line,molfile.toString());
            long startMillis = System.currentTimeMillis();
            boolean match = matcher.smilesMatchesHandle(line,queryHandle);
            if (System.currentTimeMillis()-startMillis > 1000)
                System.out.println((System.currentTimeMillis()-startMillis) +
                        "\t" + line);
            if (doFingerprint)
                for (int i=0; i<patternNames.length  &&  i<16; i++)
                {
                    // int mask = USE_ALL_BITS - (1<<i);
                    int mask = (1<<i);
                    byte[] fptmp = matcher.getFingerprintFromSmiles(line, false, mask);
                    // byte[] fptmp = fp;
                    int countTmp = countBits(fptmp);
                    bitSumByFeature[i] += countTmp;
                    if (fingerTest(fptmp, qfpByFeature[i]))
                        nFPMatchByFeature[i]++;
                    else
                    {
                        if (!match)
                        {
                            nFPPassFailByFeature[i]++;
                            bitPassFailSumByFeature[i] += countTmp;
                        }
                        else
                            nConflictByFeature[i]++;
                    }
                    if (match) bitMatchSumByFeature[i] += countTmp;
                }
            else
            {
                if (match) System.out.println(argv[1] + "\t" + line);
            }
            if (fpmatch && ! match)
            {
                totalPassFail += count;
                nPassFail++;
            }
            if (false && !matcher.heapCheck())
            {
                System.err.println(
                    "smilesMatchesQueryCT corrupted heap with smiles '"+
                    line + "'");
                break;
            }
            ntest++;
            if (ntest%100 == 0) System.err.print('.');
            if (ntest%5000 == 0) System.err.println(" "+ntest);
            if (fpmatch) nfpmatch++;
            if (match)
            {
                nmatch++;
                totalMatch += count;
            }
            if (!fpmatch && match)
            {
                if (doFingerprint)
                {
                    nconflict++;
                    PrintWriter pw = new PrintWriter(new FileOutputStream("conflict.smi", true), true);
                    pw.println(argv[0]+"\t"+argv[1].replaceAll("_.*","")+"\t"+line);
                    pw.close();
                }
            }
            /*
            if (fpmatch)
               System.out.println(
                  matcher.smilesMatchesQueryCT(line,molfile.toString()) + ": " + line);
            */
            /*
            else
            {
               System.out.println(fingerTest(fp, qfp) + line);
            System.err.println("fingerTest(fp, qfp) = " + fingerTest(fp,qfp));
            System.err.println(
                "smilesMatchesQueryCT(line, molfile.toString()) = " +
                matcher.smilesMatchesQueryCT(line,molfile.toString()));
            }
            */
        }
        if (queryHandle == 0) return;
        matcher.disposeCTHandle(queryHandle);
        reader.close();
        int nBestByteMatch = ntest;
        for (int i=0; i<qfp.length; i++)
           if (byteSelectivity[i]<nBestByteMatch)
              nBestByteMatch = byteSelectivity[i];
        System.err.println(" "+ntest);
        System.err.println("countQuery = " + countBits(qfp));
        System.err.println("ntest = " + ntest);
        System.err.println("nBestByteMatch = " + nBestByteMatch);
        System.err.println("nfpmatch = " + nfpmatch);
        System.err.println("nmatch = " + nmatch);
        System.err.println("nconflict = " + nconflict);
        System.err.println("avgBits = " + totalBits/ntest);
        System.err.println("avgFailBits = " + totalPassFail/nPassFail);
        System.err.println("avgMatchBits = " + totalMatch/nmatch);
        if (doFingerprint)
        {
            if (argv.length > 1)
            {
                System.out.print("Query");
                System.out.print("\t");
            }
            System.out.print("Pattern");
            System.out.print("\t");
            System.out.print("Bits in Query");
            System.out.print("\t");
            System.out.print("ntest");
            System.out.print("\t");
            System.out.print("nmatch");
            System.out.print("\t");
            System.out.print("nFPMatchAll");
            System.out.print("\t");
            System.out.print("nBestByteMatch");
            System.out.print("\t");
            System.out.print("nFPMatch");
            System.out.print("\t");
            System.out.print("nFPConflict");
            System.out.print("\t");
            System.out.print("avgBits");
            System.out.print("\t");
            System.out.print("yield");
            System.out.println();
            for (int i=0; i<patternNames.length  &&  i<16; i++)
            {
                if (patternNames[i].length() < 1) continue;
                int cntTmp = countBits(qfpByFeature[i]);
                // if (cntTmp == 0) continue;
                if (argv.length > 1)
                {
                    System.out.print(argv[1]);
                    System.out.print("\t");
                }
                System.out.print(patternNames[i]);
                System.out.print("\t");
                System.out.print(cntTmp);
                System.out.print("\t");
                System.out.print(ntest);
                System.out.print("\t");
                System.out.print(nmatch);
                System.out.print("\t");
                System.out.print(nfpmatch);
                System.out.print("\t");
                System.out.print(nBestByteMatch);
                System.out.print("\t");
                System.out.print(nFPMatchByFeature[i]);
                System.out.print("\t");
                System.out.print(nConflictByFeature[i]);
                System.out.print("\t");
                System.out.print((bitSumByFeature[i])/ntest);
                System.out.print("\t");
                if (cntTmp > 0  &&  nfpmatch > 0)
                {
                    double yield = (double)nFPMatchByFeature[i]/nfpmatch;
                    yield = countBits(qfp)*Math.log(yield)/cntTmp;
                    System.out.print(yield);
                }
                System.out.println();
            }
            if (argv.length > 1)
            {
                System.out.print(argv[1]);
                System.out.print("\t");
            }
            System.out.print("ALL_PATTERNS");
            System.out.print("\t");
            System.out.print(countBits(qfp));
            System.out.print("\t");
            System.out.print(ntest);
            System.out.print("\t");
            System.out.print(nmatch);
            System.out.print("\t");
            System.out.print(nfpmatch);
            System.out.print("\t");
            System.out.print(nBestByteMatch);
            System.out.print("\t");
            System.out.print(nfpmatch);
            System.out.print("\t");
            System.out.print(nconflict);
            System.out.print("\t");
            System.out.print(totalBits/ntest);
            System.out.print("\t1");
            System.out.println();
        }

	/*
        System.err.println("\nbyteSelectivity:");
        for (int i=0; i<qfp.length; i++)
        {
           System.err.print((1000*byteSelectivity[i])/ntest);
           if ((i+1) % 10 == 0) System.err.println();
           else                 System.err.print('\t');
        }
	*/
        System.err.println();
    }

    static public boolean fingerTest(byte[] fp, byte[] qfp)
    {
        if (fp == null  ||  qfp == null) return false;
        if (fp.length != qfp.length)     return false;
        for (int i=0; i<fp.length; i++)
            if ((fp[i] & qfp[i]) != qfp[i]) return false;
        return true;
    }

    /**
     * Returns the number of bits in the byte array fp.
     */
    public final static int countBits(byte[] fp)
    {
        int n=0;
        int length = fp.length;
 
        for (int i=0; i<length; i++)
        {  if(fp[i] != 0)     // optimization: for oracle without jit
              n += bits_per_byte[0x00FF&fp[i]];
        }

        return n;
    }

    /**
     * Returns the Tanimoto similarity of the two fingerprints fp1 and fp2.
     */
    public final static double similarity(byte[] fp1, byte[] fp2)
    {
        int nfp1 = countBits(fp1);
        int nfp2 = countBits(fp2);
        int nAnd = 0;
        for (int i=0; i<fp1.length && i<fp2.length; i++)
            nAnd += bits_per_byte[0xFF&(fp1[i]&fp2[i])];
        return nAnd/((double)(nfp1+nfp2-nAnd));
    }
}
