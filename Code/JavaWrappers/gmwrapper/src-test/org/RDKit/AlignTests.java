/* 
 * $Id$
 *
 *  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
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

import java.io.File;

import org.junit.Test;

public class AlignTests extends GraphMolTest {


	@Test
	public void testO3ABasic () {
            String fname = new File(getRdBase(),
                                    "Code/GraphMol/MolAlign/test_data/ref_e2.sdf").getPath();
            SDMolSupplier sdsup = new SDMolSupplier(fname);
            ROMol m1 = sdsup.next();
            ROMol m2 = sdsup.next();
            Double_Pair res = m1.O3AAlignMol(m2);
            assertEquals(res.getFirst(),0.049,.001);
            assertEquals(res.getSecond(),119.98,.01);
            
	}
    
	@Test
	public void testgetAlignmentTransform(){
		//Test alignment with Transform base on GraphMol/MolAlign/testMolAlign.cpp#test1MolAlign()
		String fname0 = new File(getRdBase(),
					 "Code/GraphMol/MolAlign/test_data/1oir.mol").getPath();
		String fname1 = new File(getRdBase(),
					 "Code/GraphMol/MolAlign/test_data/1oir_conf.mol").getPath();
		ROMol m0 = RWMol.MolFromMolFile(fname0);
		ROMol m1 = RWMol.MolFromMolFile(fname1);
		Transform3D trans = new Transform3D();
		double res = m0.getAlignmentTransform(m1, trans);
		assertEquals(res, 0.6578, 0.001);
		m0.delete();
		m1.delete();
	}
	
	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.AlignTests");
	}

}
