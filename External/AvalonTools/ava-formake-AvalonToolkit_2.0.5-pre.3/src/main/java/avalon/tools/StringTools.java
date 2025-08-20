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

import java.util.Vector;

/**
 * Implements a few general pupose String utilities that may be useful
 * in other places.
 *
 * Note, some of these tools may be available elsewhere, but it was easier
 * to reimplement them than to trace and use another source.
 */
public class StringTools
{
    /**
     * Utility function to parse a string into an array of its separator
     * delimited components.
     *
     * treats two consecutive separators as two tokens 
     * (this is different from StringTokenizer )
     */
    public static String[] stringToArray(String line, char separator)
    {
 
        Vector<String> tokens = new Vector<String>();
        while (!line.equals(""))
        {
            int index = line.indexOf(separator);
            if (index < 0)
	        {
	            tokens.add(line);
                line = "";
            }
	        else
	        {
	            String token = line.substring(0, index);
	            line = line.substring(line.indexOf(separator)+1);
	            tokens.add(token);
	        }
        }
        String result[] = new String[tokens.size()];
        for (int i=0; i<result.length; i++)
	        result[i] = (String)tokens.elementAt(i);

        return result;
    }

    /**
     * Replaces all non-overlapping occurrences of <code>from</code> in
     * <code>source</code> with <code>to</code>
     *
     * It returns a modified string.
     */
    public static String replace(String source, String from, String to)
    {
        if (source == null  ||
            from   == null  ||
            source.indexOf(from,0) < 0) return source;

        if (to == null) to = "";        // Safeguard
        StringBuffer buffer = new StringBuffer();
        int istart=0;
        while (true)
        {
            int imatch = source.indexOf(from, istart);
            if (imatch < 0)      // no match in rest of string => just copy it
            {
                buffer.append(source.substring(istart));
                break;
            }
            else                 // match => replace and continue scan
            {
                 buffer.append(source.substring(istart, imatch));
                 buffer.append(to);
                 istart = imatch+from.length();
            }
        }
        return buffer.toString();
    }

}
