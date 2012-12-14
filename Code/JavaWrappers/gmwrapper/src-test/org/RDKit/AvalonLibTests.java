/* 
 * $Id$
 *
 *  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
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

public class AvalonLibTests extends GraphMolTest {
    @Test public void testAvalonTools1() {
        ROMol m1;
        m1 = RWMol.MolFromSmiles("n1ccccc1");
        String newsmi=RDKFuncs.getCanonSmiles(m1);
        assertEquals(newsmi,"c1ccncc1",newsmi);

        newsmi=RDKFuncs.getCanonSmiles("n1ccccc1",true);
        assertEquals(newsmi,"c1ccncc1",newsmi);
    }
    
    @Test public void testStereochem() {
	ROMol m1;
	m1 = RWMol.MolFromSmiles("C[C@H](F)C(F)Cl");
	RDKFuncs.assignStereochemistry(m1,true,true,true);
	int nKnown=0;
	int nPossible=0;
	for(int i=0;i<m1.getNumAtoms();i++){
	    if(m1.getAtomWithIdx(i).hasProp("_CIPCode")){
		nKnown++;
		nPossible++;
	    } else if(m1.getAtomWithIdx(i).hasProp("_ChiralityPossible")){
		nPossible++;
	    }
	}
	assertEquals(nKnown,1);
	assertEquals(nPossible,2);
    }

    @Test public void testAvalonTools2() {
	StringInt_Pair res;
	res=RDKFuncs.checkMolString("c1ccccn1",true);
	assertEquals(res.getSecond(),0,res.getSecond());
	assertFalse(res.getFirst(),""==res.getFirst());

	res=RDKFuncs.checkMolString("c1c(R)cccc1C1(CC-C(C)C1)",true);
	assertEquals(1,res.getSecond());
	assertEquals(res.getFirst(),"",res.getFirst());
    }

    @Test public void testAvalonTools3() {
	ROMol m1,m2;
	m1 = RWMol.MolFromSmiles("c1ccccc1");
	m2 = RWMol.MolFromSmiles("c1ccccc1C");
	ExplicitBitVect bv1,bv2;
	bv1 = new ExplicitBitVect(1024);
	bv2 = new ExplicitBitVect(1024);
	RDKFuncs.getAvalonFP(m1,bv1,1024,false,false,RDKFuncs.getAvalonSSSBits());
	RDKFuncs.getAvalonFP(m2,bv2,1024,false,false,RDKFuncs.getAvalonSSSBits());
	assertEquals(RDKFuncs.TanimotoSimilarity(bv1,bv1),1.0,0.001);
	assertEquals(RDKFuncs.TanimotoSimilarity(bv2,bv2),1.0,0.001);
	assertTrue(RDKFuncs.TanimotoSimilarity(bv1,bv2)<1.0);
	assertTrue(bv1.getNumOnBits()>0);
	assertTrue(bv2.getNumOnBits()>0);
	assertTrue(bv1.getNumOnBits()<bv2.getNumOnBits());
	assertEquals(RDKFuncs.NumOnBitsInCommon(bv1,bv2),bv1.getNumOnBits());

    }
    @Test public void testAvalonTools4() {
	ROMol m1,m2;
	m1 = RWMol.MolFromSmiles("c1ccccc1");
	m2 = RWMol.MolFromSmiles("c1ccccc1C");
	ExplicitBitVect bv1,bv2;
	bv1 = new ExplicitBitVect(1024);
	bv2 = new ExplicitBitVect(1024);
	RDKFuncs.getAvalonFP(m1,bv1,1024,false,false,RDKFuncs.getAvalonSimilarityBits());
	RDKFuncs.getAvalonFP(m2,bv2,1024,false,false,RDKFuncs.getAvalonSimilarityBits()
);
	assertEquals(RDKFuncs.TanimotoSimilarity(bv1,bv1),1.0,0.001);
	assertEquals(RDKFuncs.TanimotoSimilarity(bv2,bv2),1.0,0.001);
	assertTrue(RDKFuncs.TanimotoSimilarity(bv1,bv2)<1.0);
	assertTrue(bv1.getNumOnBits()>0);
	assertTrue(bv2.getNumOnBits()>0);
	assertTrue(bv1.getNumOnBits()<bv2.getNumOnBits());
    }

    public static void main(String args[]) {
        org.junit.runner.JUnitCore.main("org.RDKit.AvalonLibTests");
    }

}
