/*
 *
 *  Copyright (c) 2016, Greg Landrum
 *  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *     * Neither the name of Novartis Institutes for BioMedical Research Inc.
 *       nor the names of its contributors may be used to endorse or promote
 *       products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package org.RDKit;
import static org.junit.Assert.*;

import org.junit.Test;

public class BitOpsTests extends GraphMolTest {
    @Test public void testFPS() {
      ExplicitBitVect bv;
	    bv = new ExplicitBitVect(32);
      String fps;
      fps = RDKFuncs.BitVectToFPSText(bv);
      assertEquals(fps,"00000000");
      bv.setBit(0);
      bv.setBit(1);
      bv.setBit(17);
      bv.setBit(23);
      bv.setBit(31);

      fps = RDKFuncs.BitVectToFPSText(bv);
      assertEquals(fps,"03008280");
    }
    @Test public void testFPS2() {
      ExplicitBitVect bv;
	    bv = new ExplicitBitVect(32);
      String fps="03008280";
      RDKFuncs.UpdateBitVectFromFPSText(bv,fps);
      assertEquals(bv.getNumOnBits(),5);
      assertTrue(bv.getBit(0));
      assertTrue(bv.getBit(1));
      assertTrue(bv.getBit(17));
      assertTrue(bv.getBit(23));
      assertTrue(bv.getBit(31));
    }

    public static void main(String args[]) {
        org.junit.runner.JUnitCore.main("org.RDKit.BitOpsTests");
    }

}
