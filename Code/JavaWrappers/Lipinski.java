/* 
 * $Id: Lipinski.java 120 2011-01-18 06:24:55Z bill.smith $
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

// Class attempts to provide functionality in
// RDBASE/rdkit/Chem/Lipinski.py

public class Lipinski {

	static ROMol HDonorSmarts = RWMol.MolFromSmarts("[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]");
	// changes log for HAcceptorSmarts:
	//  v2, 1-Nov-2008, GL : fix amide-N exclusion; remove Fs from definition
	static ROMol HAcceptorSmarts = RWMol.MolFromSmarts("[$([O,S;H1;v2]-[!$(*=[O,N,P,S])])," +
	"$([O,S;H0;v2]),$([O,S;-])," +
	"$([N;v3;!$(N-*=!@[O,N,P,S])])," +
	"$([nH0,o,s;+0])" +
	"]");
	static ROMol HeteroatomSmarts = RWMol.MolFromSmarts("[!#6;!#1]");
	//  NOTE: the Rotatable bond smarts here doesn"t treat deuteriums (which are left in the graph
	//  and therefore contribute to the degree of a carbon) the same as hydrogens (which are removed
	//  from the graph). So the bond in [2H]C([2H])([2H])C([2H])([2H])[2H] *is* considered
	//  rotatable.
	static ROMol RotatableBondSmarts = RWMol.MolFromSmarts("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]");
	static ROMol StrictRotatableBondSmarts = RWMol.MolFromSmarts("[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-!@[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])]");
	static ROMol NHOHSmarts = RWMol.MolFromSmarts("[#8H1,#7H1,#7H2,#7H3]");
	static ROMol NOCountSmarts = RWMol.MolFromSmarts("[#7,#8]");

	private static int getMatchCount(ROMol test, ROMol smarts) {
		
		return (int) test.getSubstructMatches(smarts, true).size();
	}
	public static int getHDonorCount(ROMol mol) {
		return getMatchCount(mol, HDonorSmarts);
	}
	public static int getHAcceptorCount(ROMol mol) {
		return getMatchCount(mol, HAcceptorSmarts);
	}
	public static int getHeteroatomCount(ROMol mol) {
		return getMatchCount(mol, HeteroatomSmarts);
	}
	public static int getRotatableBondCount(ROMol mol) {
		return getMatchCount(mol, RotatableBondSmarts);
	}
	public static int getStrictRotatableBondCount(ROMol mol) {
		return getMatchCount(mol, StrictRotatableBondSmarts);
	}
	public static int getNHOHCount(ROMol mol) {
		return getMatchCount(mol, NHOHSmarts);
	}
	public static int getNOCount(ROMol mol) {
		return getMatchCount(mol, NOCountSmarts);
	}

}
