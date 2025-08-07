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

import java.io.IOException;
import java.io.InputStream;

public class UnixTools implements ClipboardTools {
	
    static ClipboardTools toolHost = null;
       
    static String tmpMolFilename = System.getProperty("java.io.tmpdir") + System.getProperty("file.separator") + "tmp.mol";

    private UnixTools()
    {
    }
 
    public static ClipboardTools getToolHost()
    {
    	if (toolHost == null) {
      	  toolHost = new UnixTools(); 
        }
    	return toolHost;
    }

	public void emptyClipboard() {
		// TODO Auto-generated method stub

	}

	public int getMOLFileBytesOfClipboardCT(byte[] buffer) {
		// TODO Auto-generated method stub
		return 0;
	}

	public String getMOLFileOfClipboardCT() {
		// TODO Auto-generated method stub
		return null;
	}

	public int getSMARTSBytesOfClipboardCT(byte[] buffer) {
		// TODO Auto-generated method stub
		return 0;
	}

	public String getSMARTSOfClipboardCT() {
		// TODO Auto-generated method stub
		return null;
	}

	public int getSMILESBytesOfClipboardCT(byte[] buffer) {
		// TODO Auto-generated method stub
		return 0;
	}

	public String getSMILESOfClipboardCT() {
		// TODO Auto-generated method stub
		return null;
	}

	public int putBytesToClipboard(byte[] buffer, String format) {
		// TODO Auto-generated method stub
		return 0;
	}

	public void putCTToClipboard(String molfile) {
		// TODO Auto-generated method stub

	}

	public boolean putDataToClipboard(byte[] data, String format) {
		// TODO Auto-generated method stub
		return false;
	}

	public void putMOLFileToClipboard(String fname) {
		// TODO Auto-generated method stub

	}

	public boolean putWMFToClipboard(InputStream file) throws IOException {
		// TODO Auto-generated method stub
		return false;
	}

	public void showClipboardCT() {
		// TODO Auto-generated method stub

	}

}
