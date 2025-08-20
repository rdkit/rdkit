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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;

import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.SystemFlavorMap;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;

import avalon.tools.LibraryToolbox;

/**
 * The class JNIWinTools provides wrappers for Windows services.
 * It is intended to be used from Java applications and trusted applets.
 * @History:    B. Rohde 30-Aug-2001 created based on JNIDepict.java
 */
public class JNIWinTools implements ClipboardTools
{
    static ClipboardTools toolHost = null;
       
    static String tmpMolFilename = System.getProperty("java.io.tmpdir") + System.getProperty("file.separator") + "tmp.mol";

    private JNIWinTools()
    {
    }
 
    static
    {
      loadLibraries();
      toolHost = new JNIWinTools();
    }

    public static ClipboardTools getToolHost()
    {
        return toolHost;
    }

    static void loadLibraries()
    {
  	  // try to load from classpath, including jar files
		String libName = "JNIWinTools";
        if (LibraryToolbox.loadLibrary(JNIWinTools.class, libName)) return;
		try {
			System.loadLibrary(libName);
			return;
		} catch (SecurityException e) {
			System.err.println(libName
                      + ": Security exception, can't load native library from jar, "
                      + e.getMessage());
		} catch (UnsatisfiedLinkError e) {
			System.err.println(libName
                      + ": Not found, can't load native library from class path, "
                      + e.getMessage());
		}

		// Make sure that we get the links to the needed DLLs.

      String JNIWinTools = System.getProperty("user.dir") + "/JNIWinTools.dll";
      System.err.println("trying to load '" + JNIWinTools + "'");
      if ((new File(JNIWinTools)).exists())
      {
          try
          {
             System.load(JNIWinTools);
             return;
          }
          catch (Exception e)
          {
             System.err.println("Can't load native library '" +
                                JNIWinTools + "'!");
             e.printStackTrace();
          }
          catch (Error err)
          {
            System.err.println("Can't load native library '" +
                               JNIWinTools + "'!");
            err.printStackTrace();
          }
      }

    }

    /* (non-Javadoc)
	 * @see jni.ClipboardTools#putDataToClipboard(byte[], java.lang.String)
	 */
    public synchronized boolean putDataToClipboard(byte[] data, String format)
    {
        return false;
    }

    /* (non-Javadoc)
	 * @see jni.ClipboardTools#putWMFToClipboard(java.io.InputStream)
	 */
    public synchronized boolean putWMFToClipboard(InputStream file)
        throws IOException
    {
        // skip placeable header
        byte[] header = new byte[22];
        if (file.read(header) < header.length) return false;

        byte[] bytes = new byte[file.available()];
        if (file.read(bytes) < bytes.length) return false;
        int mtSize =
            ((0xFF&bytes[6])<<0) + ((0xFF&bytes[7])<<8)+
            ((0xFF&bytes[8])<<16)+ ((0xFF&bytes[9])<<24);
        byte[] newBytes = new byte[header.length+18+mtSize*2];
System.err.println("newBytes.length = " + newBytes.length);
        // Put placer in front of bytes
        for (int i=0; i<header.length; i++)
            newBytes[i] = header[i];
        // Copy real metafile bytes
        for (int i=0; i<mtSize*2; i++)
            newBytes[i+header.length] = bytes[i];
        System.err.println("putBytesToClipboard returned " +
            this.putBytesToClipboard(newBytes, "Picture"));
        return true;
    }

    // make the MDLCT format known
    static DataFlavor ctFlavor = null;
    static
    {
        ctFlavor = new DataFlavor("chemical/x-mdlct-molfile", "MDLCT");
        SystemFlavorMap sfm = (SystemFlavorMap) SystemFlavorMap.getDefaultFlavorMap();
        sfm.setFlavorsForNative("MDLCT", new DataFlavor[]{ctFlavor});
        sfm.setNativesForFlavor(ctFlavor, new String[]{"MDLCT"});
    }

    /* (non-Javadoc)
     * @see jni.ClipboardTools#emptyClipboard()
     */
    public native void emptyClipboard();

    /* (non-Javadoc)
     * @see jni.ClipboardTools#putBytesToClipboard(byte[], java.lang.String)
     */
    public native int putBytesToClipboard(byte buffer[], String format);

    /* (non-Javadoc)
     * @see jni.ClipboardTools#getMOLFileBytesOfClipboardCT(byte[])
     */
    public synchronized native int getMOLFileBytesOfClipboardCT(byte buffer[]);

    /* (non-Javadoc)
     * @see jni.ClipboardTools#getMOLFileOfClipboardCT()
     */
    public synchronized String getMOLFileOfClipboardCT()
    {
        StringBuffer result = new StringBuffer();
        try
        {
            // get the clipboard contents
            Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
            if (clipboard.isDataFlavorAvailable(ctFlavor))
            {
                Transferable transferable = clipboard.getContents(null);
                InputStream stream = (InputStream) transferable.getTransferData(ctFlavor);
                int len;
                while (0 <= (len = stream.read()))
                {
                    for (int j=0; j<len; j++)
                    {
                        int cByte = stream.read();
                        if (cByte > 0)
                        {
                            result.append((char)cByte);
                        }
                        else       break;
                    }
                    result.append("\n");
                }
                return result.toString();
            }
        }
        catch (UnsupportedFlavorException ufe)
        {
            ufe.printStackTrace();
            System.err.println("fallback to JNI method");

            byte buffer[] = new byte[100000];
            String molfile;
            int size;

            size = getMOLFileBytesOfClipboardCT(buffer);
            if (0 < size)
                return (new String(buffer, 0, size));
            else
            {
                buffer = new byte[(-1)*size+1];
                size = getMOLFileBytesOfClipboardCT(buffer);
            }

            if (0 < size) return (new String(buffer, 0, size));

            return null;
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
        return null;

    }

    /* (non-Javadoc)
     * @see jni.ClipboardTools#getSMARTSBytesOfClipboardCT(byte[])
     */
    public synchronized native int getSMARTSBytesOfClipboardCT(byte buffer[]);

    /* (non-Javadoc)
     * @see jni.ClipboardTools#getSMARTSOfClipboardCT()
     */
    public synchronized String getSMARTSOfClipboardCT()
    {
	byte buffer[] = new byte[1000];
	String smarts;
	int size;
	
	if (0 < (size = getSMARTSBytesOfClipboardCT(buffer)))
	{
            smarts = new String(buffer, 0, size);
	}
	else
            smarts = null;
	System.err.println("size = " + size);
	return (smarts);
    }

	/* (non-Javadoc)
	 * @see jni.ClipboardTools#getSMILESBytesOfClipboardCT(byte[])
	 */
	   public synchronized native int getSMILESBytesOfClipboardCT(byte buffer[]);

	/* (non-Javadoc)
	 * @see jni.ClipboardTools#getSMILESOfClipboardCT()
	 */
	   public synchronized String getSMILESOfClipboardCT()
	   {
	      byte buffer[] = new byte[1000];
	      String smiles;
	      int size;
	
	      if (0 < (size = getSMILESBytesOfClipboardCT(buffer)))
	      {
	        smiles = new String(buffer, 0, size);
	      }
	      else
	        smiles = null;
	      return (smiles);
	   }

	/* (non-Javadoc)
	 * @see jni.ClipboardTools#putCTToClipboard(java.lang.String)
	 */
	public synchronized void putCTToClipboard(String molfile)
	   {
	      assureTmpMoleFileIsWritable(tmpMolFilename);
	
	      try
	      {
	         molfile = fixCT(molfile);
	         PrintWriter w = new PrintWriter(new FileWriter(tmpMolFilename));
	         w.print(molfile);
	         w.close();
	         putMOLFileToClipboard(tmpMolFilename);
	      }
	      catch (Exception e)
	      {
	         e.printStackTrace();
	      }
	
	   }

	/* (non-Javadoc)
	 * @see jni.ClipboardTools#putMOLFileToClipboard(java.lang.String)
	 */
	   public synchronized native void putMOLFileToClipboard(String fname);

	/* (non-Javadoc)
	 * @see jni.ClipboardTools#showClipboardCT()
	 */
	   public synchronized native void showClipboardCT();

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
	      //System.err.println(tmpMolFilename+" writable!");
	
	   }

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

	static public void main(String argv[])
        throws IOException
    {
        loadLibraries();

        byte[] bytes = null;

        ClipboardTools jni = getToolHost();
        if (argv.length >= 1)   // there is a file => use it
        {
            FileInputStream file = new FileInputStream(argv[0]);
            System.err.println("file " + argv[0] + " has " + 
                               file.available() + " bytes");
            bytes = new byte[22];
            file.read(bytes);    // skip placeable header
            int left   = ((0xFF&bytes[6])<<0)  + ((0xFF&bytes[7])<<8);
            int top    = ((0xFF&bytes[8])<<0)  + ((0xFF&bytes[9])<<8);
            int right  = ((0xFF&bytes[10])<<0) + ((0xFF&bytes[11])<<8);
            int bottom = ((0xFF&bytes[12])<<0) + ((0xFF&bytes[13])<<8);
            int inch = ((0xFF&bytes[14])<<0) + ((0xFF&bytes[15])<<8);
System.err.println("left = " + left);
System.err.println("top = " + top);
System.err.println("right = " + right);
System.err.println("bottom = " + bottom);
System.err.println("inch = " + inch);
            int xExt = Math.abs((right-left)*2540/inch);
            int yExt = Math.abs((bottom-top)*2540/inch);
System.err.println("xExt = " + xExt);
System.err.println("yExt = " + yExt);
            bytes = new byte[file.available()];
            if (0 >= file.read(bytes))
            {
                file.close();
                return;
            }
            System.err.println("mtType = " +
                Integer.toHexString(((0xFFFF)&(bytes[1]*256+bytes[0]))));
            System.err.println("mtHeaderSize = " +
                Integer.toHexString(((0xFFFF)&(bytes[3]*256+bytes[2]))));
            System.err.println("mtVersion = " +
                Integer.toHexString(((0xFFFF)&(bytes[5]*256+bytes[4]))));
            int mtSize =
                ((0xFF&bytes[6])<<0) + ((0xFF&bytes[7])<<8)+
                ((0xFF&bytes[8])<<16)+ ((0xFF&bytes[9])<<24);
            System.err.println("mtSize " + Integer.toHexString(mtSize));
            System.err.println("mtNoObjects = " +
                Integer.toHexString(((0xFFFF)&(bytes[11]*256+bytes[10]))));
            System.err.println("mtMaxRecords " +
                Integer.toHexString(0xFF&bytes[12]+
                                    ((0xFF&bytes[13])<<8)+
                                    ((0xFF&bytes[14])<<16)+
                                    ((0xFF&bytes[15])<<24)));
            System.err.println("mtNoParameters = " +
                Integer.toHexString(((0xFFFF)&(bytes[17]*256+bytes[16]))));
            // bytes = new byte[mtSize*2];
            byte[] newBytes = new byte[18+mtSize*2];
System.err.println("totalLength = " + (18+mtSize*2));
System.err.println("bytes.length = " + bytes.length);
System.err.println("newBytes.length = " + newBytes.length);
            for (int i=0; i<mtSize*2; i++)
                newBytes[i] = bytes[i];
            System.err.println("file " + argv[0] + " has " + bytes.length + " bytes");
            System.err.println("putBytesToClipboard returned " +
                jni.putBytesToClipboard(newBytes, "Picture"));
            file.close();
            return;
        }
        else
        {
            bytes = "Hello brave new world!".getBytes();
            System.err.println("putBytesToClipboard returned " +
                jni.putBytesToClipboard(bytes, "Text"));
        }
    }
}
