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

public interface ClipboardTools {

	/**
	 * This method wraps the corresponding native method call to place
	 * binary <code>data</code> to the system clipboard in a
	 * particular <code>format</code>.
	 */
	public abstract boolean putDataToClipboard(byte[] data, String format);

	/**
	 * This method reads the InputStream in for a WMF (can be placeable
	 * but may be not) and puts the "Picture" format on the clipboard.
	 */
	public abstract boolean putWMFToClipboard(InputStream file)
			throws IOException;

	/**
	 * Removes all formats from the Windows sytem clipboard.
	 */
	public abstract void emptyClipboard();

	/**
	 * Native call to put data to the system clipboard using format.
	 */
	public abstract int putBytesToClipboard(byte buffer[], String format);

	/**
	 * This function fills the buffer with the MOL-File on the clipboard.
	 * It returns the actual number of bytes written to the buffer or (-1)
	 * times the number required for a subsequent call with a larger buffer.
	 */
	public abstract int getMOLFileBytesOfClipboardCT(byte buffer[]);

	/**
	 * This function wraps the native function that returns the MOL-File string
	 * corresponding to the structure on the System Clipboard.
	 */
	public abstract String getMOLFileOfClipboardCT();

	/**
	 * This function fills the SMARTS string corresponding to the
	 * structure on the System Clipboard into a byte array given as parameter.
	 * It returns the actual number of bytes written or (-1) time the
	 * number required for a subsequent call with a larger buffer.
	 */
	public abstract int getSMARTSBytesOfClipboardCT(byte buffer[]);

	/**
	 * This function wraps the native function that returns the SMARTS string
	 * corresponding to the structure on the System Clipboard.
	 */
	public abstract String getSMARTSOfClipboardCT();

	/**
	 * This function fills the SMILES string corresponding to the
	 * structure on the System Clipboard into a byte array given as parameter.
	 * It returns the actual number of bytes written or (-1) time the
	 * number required for a subsequent call with a larger buffer.
	 */
	public abstract int getSMILESBytesOfClipboardCT(byte buffer[]);

	/**
	 * This function directly returns the SMILES string corresponding
	 * to the structure on the System Clipboard.
	 */
	public abstract String getSMILESOfClipboardCT();

	public abstract void putCTToClipboard(String molfile);

	/**
	 * This native function sends the contents of file to the Clipboard
	 * tagged as a MOL file.
	 */
	public abstract void putMOLFileToClipboard(String fname);

	/**
	 * This method prints the connection table that is on the System
	 * Clipboard or an error message if there is none. The function
	 * is used for learning purposes only.
	 */
	public abstract void showClipboardCT();

}
