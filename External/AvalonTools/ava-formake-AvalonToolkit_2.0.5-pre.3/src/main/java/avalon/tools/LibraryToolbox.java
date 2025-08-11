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
package avalon.tools;

import java.awt.Toolkit;
import java.io.File;
import java.util.StringTokenizer;
import javax.swing.JOptionPane;

public class LibraryToolbox {
	/**
	 * Scans the current CLASSPATH for avalon.jar and tries to load libName from
	 * the same directory as avalon.jar is located.
	 */
	public static boolean loadLibraryFromClasspath(String libName) {
		String classpath = System.getProperty("java.class.path");
		String separator = System.getProperty("path.separator");
		libName = System.mapLibraryName(libName);

		StringTokenizer tokens = new StringTokenizer(classpath, separator);
		String libPath = null;
		while (tokens.hasMoreTokens() && libPath == null) {
			String path = tokens.nextToken();
			if (path.toLowerCase().endsWith("avalon.jar"))
				libPath = path.substring(0, path.length()
						- "avalon.jar".length())
						+ libName;
		}
		try {
			System.load(libPath);
			return true;
		} catch (Throwable t) {
			System.err.println("Could not load library via CLASSPATH.");
			System.err.println("   CLASSPATH is: "
					+ System.getProperty("java.class.path"));
			System.err.println("   Library Name is: " + libName);
			System.err.println("   Library should be in: " + libPath);
		}
		return false;
	}

    /**
     * Loads a library with the file root fileRootName from the directory that contains the JAR file of this class.
     */
    static public boolean loadLibrary(Class refClass, String fileRootName)
    {
        try
        { 
            // trying to load file from dirictory containing current JAR file
            System.err.println("Java library path: "+System.getProperty("java.library.path"));
            // System.err.println("refClass = " + refClass);
            System.err.println("using " + refClass + ".getProtectionDomain().getCodeSource().getLocation().getPath()");
            // System.err.println("Protection Domain: " + refClass.getProtectionDomain());
            // System.err.println("Code Source: " + refClass.getProtectionDomain().getCodeSource());
            // System.err.println("Location: " + refClass.getProtectionDomain().getCodeSource().getLocation());
            System.err.println("ClassLoader path: " + refClass.getProtectionDomain().getCodeSource().getLocation().getPath());
            String libDir = (new File(refClass.getProtectionDomain().getCodeSource().getLocation().getPath())).getParent()+File.separatorChar;
            // try different versions for different architectures
            if (File.separatorChar == '\\')
            {
                // only try DLL on Windows-like OSes
                try
                {
                    System.load(libDir+fileRootName+".DLL");
                    System.err.println("Loaded '" + libDir+fileRootName+".DLL'");
                    return true;
                }
                catch (Throwable t)
                {
                    System.err.println("Could not load '" + libDir+fileRootName+".DLL'");
                    t.printStackTrace();
                    return false;
                }
            }
            else
            {
                try
                {
                    System.load(libDir+fileRootName+".so");
                    System.err.println("Loaded '" + libDir+fileRootName+".so'");
                    return true;
                }
                catch (Throwable t) {}
                try
                {
                    System.load(libDir+"lib"+fileRootName+".so");
                    System.load(libDir+"lib"+fileRootName+".so");
                    return true;
                }
                catch (Throwable t) {}
                try
                {
                    System.load(libDir+fileRootName+".jnilib");
                    System.load(libDir+fileRootName+".jnilib");
                    return true;
                }
                catch (Throwable t) {}
                System.err.println("Could not load library '" + fileRootName + "' from directory '" + libDir + "'");
                return false;
            }
        }
        catch (Exception e)
        {
            e.printStackTrace();
            return false;
        }
    }

    /**
     * Encapsulates showing an error message.
     *
     * The idea is to avoid trying to open a dialog when the program is
     * executed headlessly. It's replicated from the Toolbox class to be used from
     * the WebServer's boot classpath.
     */
    public static void showError(String message, String title, boolean doBeep)
    {
        if (doBeep) Toolkit.getDefaultToolkit().beep();
        if (!java.awt.GraphicsEnvironment.isHeadless())
            JOptionPane.showMessageDialog(null, message, title, JOptionPane.ERROR_MESSAGE);
        else
            System.err.println(title +": " + message);
    }
}
