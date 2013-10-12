/* 
 * $Id: ForceFieldsTests.java 131 2011-01-20 22:01:29Z ebakke $
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

import java.io.File;

import org.junit.Test;

public class ForceFieldsTests extends GraphMolTest {

	File testDataDir = new File(getRdBase(), 
	"Code/GraphMol/ForceFieldHelpers/UFF/test_data");

	// From GraphMol/ForceFieldHelpers/Wrap/testHelpers.py
	@Test
	public void testForceFieldOptimizationBasic () {
		File molFile = new File(testDataDir, "benzene.mol");
		ROMol m = RWMol.MolFromMolFile(molFile.getPath());
		assertEquals(0, ForceField.UFFOptimizeMolecule(m));
	}
	@Test
	public void testForceFieldOptimizationBasicWithArgs() {
		File molFile = new File(testDataDir, "benzene.mol");

		ROMol m = RWMol.MolFromMolFile(molFile.getPath());
		assertFalse(0 == ForceField.UFFOptimizeMolecule(m, 1));

		m = RWMol.MolFromMolFile(molFile.getPath());
		assertEquals(0, ForceField.UFFOptimizeMolecule(m, 200, 2.0)); 
		// 200 is default maxIts
	
		m = RWMol.MolFromMolFile(molFile.getPath());
		assertEquals(0, ForceField.UFFOptimizeMolecule(m, 200, 10.0, -1));
		// 10.0 is default vdwThresh
	}
	@Test(expected=ConformerException.class)
	public void testForceFieldOptimizationRaisingError() {
		File molFile = new File(testDataDir, "benzene.mol");

		ROMol m = RWMol.MolFromMolFile(molFile.getPath());
		assertFalse(0 == ForceField.UFFOptimizeMolecule(m, 200, 10.0, 1));
	}	
	@Test
	public void testForceFieldCalculationsBasic() {
		File molFile = new File(testDataDir, "benzene.mol");

		ROMol m = RWMol.MolFromMolFile(molFile.getPath());
		ForceField ff = ForceField.UFFGetMoleculeForceField(m);
		assertNotNull(ff);
		double e1 = ff.calcEnergy();
		int r = ff.minimize();
		assertEquals(0, r);
		double e2 = ff.calcEnergy();
		assertTrue(e2 < e1);
	}
	@Test
	public void testForceFieldCalculationsWithArgs() {
		File molFile = new File(testDataDir, "benzene.mol");

		ROMol m = RWMol.MolFromMolFile(molFile.getPath());
		ForceField ff = ForceField.UFFGetMoleculeForceField(m);
		int r = ff.minimize(200, 1e-8); //  200 is default maxIts
		assertEquals(0, r);
		r = ff.minimize(200, 1e-4, 1e-3); // 1e-4 is default forceTol
		assertEquals(0, r);
	}
	@Test
	public void testForceFieldCalculationsOnMolFromMolBlock() {
		String molB = "" + "\n" + 
		"" + "\n" + 
		"" + "\n" + 
		"  4  4  0  0  0  0  0  0  0  0999 V2000" + "\n" + 
		"   -0.8500    0.4512   -0.6671 C   0  0  0  0  0  0  0  0  0  0  0  0" + "\n" + 
		"   -0.3307   -0.9436   -0.3641 C   0  0  0  0  0  0  0  0  0  0  0  0" + "\n" + 
		"    0.6796   -0.4074    0.5894 C   0  0  0  0  0  0  0  0  0  0  0  0" + "\n" + 
		"    0.5011    0.8998   -0.1231 C   0  0  0  0  0  0  0  0  0  0  0  0" + "\n" + 
		"  1  2  1  0" + "\n" + 
		"  2  3  1  0" + "\n" + 
		"  3  4  1  0" + "\n" + 
		"  1  4  1  0" + "\n" + 
		"M  END" + "\n";
		ROMol m = RWMol.MolFromMolBlock(molB);

	    ForceField ff = ForceField.UFFGetMoleculeForceField(m);
	    assertNotNull(ff);
	    double e1 = ff.calcEnergy();
	    int r = ff.minimize();
	    assertEquals(r, 0);
	    double e2 = ff.calcEnergy();
	    assertTrue(e2 < e1);

	}
	@Test
	public void testForceFieldMoleculeParamChecking() {
		ROMol m = RWMol.MolFromSmiles("[Cu](C)(C)(C)(C)C");
		assertFalse(ForceField.UFFHasAllMoleculeParams(m));
		m = RWMol.MolFromSmiles("C(C)(C)(C)C");
		assertTrue(ForceField.UFFHasAllMoleculeParams(m));
	}

	
	// from testHelpers.cpp
	@Test
	public void testUFFTyper1(){

		ROMol mol;
		mol = RWMol.MolFromSmiles("[SiH3]CC(=O)NC");
		assertNotNull(mol);	
		assertEquals("Si3", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(0)));
		assertEquals("C_3", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(1)));
		assertEquals("C_R", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(2)));
		assertEquals("O_R", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(3)));
		assertEquals("N_R", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(4)));

		mol = RWMol.MolFromSmiles("CC(=O)C");

		assertNotNull(mol);

		assertEquals("C_3", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(0)));
		assertEquals("C_2", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(1)));
		assertEquals("O_2", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(2)));
		assertEquals("C_3", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(3)));



		mol = RWMol.MolFromSmiles("C(=O)S");
		assertNotNull(mol);

		assertEquals("C_2", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(0)));
		assertEquals("O_2", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(1)));
		assertEquals("S_3+2", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(2)));


		mol = RWMol.MolFromSmiles("SCS(=O)S(=O)(=O)O");
		assertNotNull(mol);

		assertEquals("S_3+2", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(0)));
		assertEquals("C_3", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(1)));
		assertEquals("S_3+4", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(2)));
		assertEquals("S_3+6", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(4)));


		mol = RWMol.MolFromSmiles("PCP(O)CP(=O)(=O)");
		assertNotNull(mol);

		assertEquals("P_3+3", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(0)));
		assertEquals("C_3", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(1)));
		assertEquals("P_3+3", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(2)));
		assertEquals("P_3+5", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(5)));



		mol = RWMol.MolFromSmiles("C(F)(Cl)(Br)I");
		assertNotNull(mol);

		assertEquals("C_3", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(0)));
		assertEquals("F_", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(1)));
		assertEquals("Cl", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(2)));
		assertEquals("Br", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(3)));
		assertEquals("I_", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(4)));


		mol = RWMol.MolFromSmiles("[Li].[Na].[K].[Rb].[Cs]");
		assertNotNull(mol);

		assertEquals("Li", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(0)));
		assertEquals("Na", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(1)));
		assertEquals("K_", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(2)));
		assertEquals("Rb", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(3)));
		assertEquals("Cs", ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(4)));
	}

	// test1 from testForceField.cpp
	@Test
	public void testFFBasics() {

		ForceField ff = new ForceField();
		assertEquals(3, ff.dimension());


		Point3D_Vect ps = ff.positions3D();
		ps.add(new Point3D(0,0,0));
		ps.add(new Point3D(1,0,0));
		ps.add(new Point3D(2,0,0));
		ps.add(new Point3D(0,1,0));

		assertEquals(4, ff.positions3D().size());

		ff.initialize();

		assertEquals(1.0, ff.distance(0,1),defaultDoubleTol);
		assertEquals(1.0, ff.distance(1,0),defaultDoubleTol);
		assertEquals(0.0, ff.distance(0,0),defaultDoubleTol);
		assertEquals(2.0, ff.distance(0,2),defaultDoubleTol);
		assertEquals(2.0, ff.distance(2,0),defaultDoubleTol);
		assertEquals(1.0, ff.distance(0,3),defaultDoubleTol);
		assertEquals(1.0, ff.distance(3,0),defaultDoubleTol);
		assertEquals(0.0, ff.distance(3,3),defaultDoubleTol);
		assertEquals(1.0, ff.distance(1,2),defaultDoubleTol);
		assertEquals(1.0, ff.distance(2,1),defaultDoubleTol);

	}
	// testUFF1 from testForceField.cpp
	@Test
	public void testBasicsOfUFFBondStretchTerms () {

		AtomicParams p1 = new AtomicParams();
		AtomicParams p2 = new AtomicParams();

		double restLen;
		double forceConstant;

		// sp3 carbon:
		p1.setR1(.757);
		p1.setZ1(1.912);
		p1.setGMP_Xi(5.343);

		// sp3 - sp3: checks basics
		restLen = RDKFuncs.calcBondRestLength(1.0, p1, p1);
		assertEquals(1.514, restLen, defaultDoubleTol);

		forceConstant = RDKFuncs.calcBondForceConstant(restLen, p1, p1);
		assertEquals(699.5918, forceConstant, defaultDoubleTol);

		// sp2 carbon:
		p2.setR1(.732);
		p2.setZ1(1.912);
		p2.setGMP_Xi(5.343);
		// sp2 - sp2: checks rBO
		restLen=RDKFuncs.calcBondRestLength(2.0, p2, p2);
		assertEquals(1.32883,restLen,defaultDoubleTol);

		forceConstant=RDKFuncs.calcBondForceConstant(restLen, p2, p2);
		assertEquals(1034.69, forceConstant, 1e-2);
	}
	// testUFF2 from testForceField.cpp
	@Test
	public void testUFFBondStretchTerms() {
		ForceField ff = new ForceField();
		Point3D_Vect ps = ff.positions3D();
		ps.add(new Point3D(0,0,0));
		ps.add(new Point3D(1.514,0,0));
		AtomicParams param1 = new AtomicParams();
		// sp3 carbon:
		param1.setR1(.757);
		param1.setZ1(1.912);
		param1.setGMP_Xi(5.343);

		// C_3 - C_3, r0=1.514, k01=699.5918
		ForceFieldContrib bs = new BondStretchContrib(ff,0,1,1, param1, param1);
		ff.contribs().add(bs);
		ff.initialize();

		Double_Array p = new Double_Array(6);
		Double_Array g = new Double_Array(6);
		for (int i = 0; i < 6; i++) {
			p.setitem(i, 0.0);
			g.setitem(i, 0.0);
		}

		double E;
		// edge case: zero bond length:
		E = bs.getEnergy(p.cast());
		assertTrue(E > 0.0);
		bs.getGrad(p.cast(), g.cast());
		for (int i=0;i<6;i++){
			assertTrue(Math.abs(g.getitem(i)) > 0.0);
		}

		p.setitem(0, 0);
		p.setitem(3, 1.514);
		for(int i=0;i<6;i++){
			g.setitem(i, 0.0);
		}
		ff.initialize();
		E = bs.getEnergy(p.cast());
		assertEquals(0.0,E,defaultDoubleTol);
	}
	// testUFF3 from testForceField.cpp
	@Test
	public void testBasicsOfUFFAngleTerms () {
		AtomicParams p1 = new AtomicParams();
		AtomicParams p2 = new AtomicParams();
		AtomicParams p3 = new AtomicParams();

		double restLen;
		double forceConstant;
		// sp3 carbon:
		p3.setR1(.757);
		p3.setZ1(1.912);
		p3.setGMP_Xi(5.343);
		p3.setTheta0(109.47 * Math.PI /180.0);

		// sp3 - sp3: checks basics
		restLen=RDKFuncs.calcBondRestLength(1.0, p3, p3);
		assertEquals(1.514,restLen,defaultDoubleTol);

		// amide bond bend:
		// C_R - N_R - C_3
		// C_R:
		p1.setR1(.729);
		p1.setZ1(1.912);
		p1.setGMP_Xi(5.343);
		// N_R:
		p2.setR1(.699);
		p2.setZ1(2.544);
		p2.setGMP_Xi(6.899);
		p2.setTheta0(120.0* Math.PI / 180.);
		restLen = RDKFuncs.calcBondRestLength(RDKFuncs.getAmideBondOrder(), p1, p2);
		assertEquals(1.357,restLen,1e-3);
		restLen=RDKFuncs.calcBondRestLength(1.0, p2, p3);
		assertEquals(1.450, restLen, 1e-3);

		forceConstant=RDKFuncs.calcAngleForceConstant(p2.getTheta0(), RDKFuncs.getAmideBondOrder(),1,
				p1, p2, p3);
		assertEquals(211.0, forceConstant, 1e-1); //  paper has 105.5



	}
	// testUFF4 from testForceFields.cpp
	@Test
	public void testForUFFAngleBendTerms() {
		ForceField ff = new ForceField();
		Point3D_Vect ps = ff.positions3D();
		Point3D p1 = new Point3D(1.514,0,0);
		Point3D p2 = new Point3D(0,0,0); 
		Point3D p3 = new Point3D(0.1,1.5,0);
		ps.add(p1);
		ps.add(p2);
		ps.add(p3);
		AtomicParams param1 = new AtomicParams();
		// sp3 carbon:
		param1.setR1(.757);
		param1.setZ1(1.912);
		param1.setGMP_Xi(5.343);
		// cheat to get the angle to 90 so that testing is easier:
		param1.setTheta0(90.0* Math.PI / 180.);

		// C_3 - C_3, r0=1.514, k01=699.5918
		ff.contribs().add(new BondStretchContrib(ff,0,1,1, param1, param1));
		ff.contribs().add(new BondStretchContrib(ff,1,2,1, param1, param1));
		ff.contribs().add(new AngleBendContrib(ff,0,1,2,1,1, param1, param1, param1));


		Point3D v1,v2;
		double theta;
		// ------- ------- ------- ------- ------- ------- -------
		// try a bit of minimization
		ff.initialize();
		ff.minimize(10,1e-8,1e-8);

		v1 = ff.positions3D().get(0).minus(ff.positions3D().get(1));
		v2 = ff.positions3D().get(1).minus(ff.positions3D().get(2));
		theta = v1.angleTo(v2);

		assertEquals(1.514,v1.length(),1e-3);
		assertEquals(1.514,v2.length(),1e-3);
		assertEquals(theta,90* Math.PI / 180.,1e-4);

		// ------- ------- ------- ------- ------- ------- -------
		// more complicated atomic coords:

		p1.setX(1.3);
		p1.setY(0.1);
		p1.setZ(0.1);
		p2.setX(-0.1);
		p2.setY(0.05);
		p2.setZ(-0.05);
		p3.setX(0.1);
		p3.setY(1.5);
		p3.setZ(0.05);
		ff.initialize();
		ff.minimize(10,1e-8,1e-8);

		v1 = ff.positions3D().get(0).minus(ff.positions3D().get(1));
		v2 = ff.positions3D().get(1).minus(ff.positions3D().get(2));
		theta = v1.angleTo(v2);

		assertEquals(1.514,v1.length(),1e-3);
		assertEquals(1.514,v2.length(),1e-3);
		assertEquals(theta,90* Math.PI / 180., 1e-4);

		// ------- ------- ------- ------- ------- ------- -------
		// try for the tetrahedral angle instead of 90:
		param1.setTheta0(109.47* Math.PI / 180.);
		ff.contribs().set((int) (ff.contribs().size() - 1), 
				new AngleBendContrib(ff,0,1,2,1,1, param1, param1, param1));

		p1.setX(1.3);
		p1.setY(0.1);
		p1.setZ(0.1);
		p2.setX(-0.1);
		p2.setY(0.05);
		p2.setZ(-0.05);
		p3.setX(0.1);
		p3.setY(1.5);
		p3.setZ(0.05);

		ff.initialize();
		ff.minimize(100,1e-8,1e-8);
		v1 = ff.positions3D().get(0).minus(ff.positions3D().get(1));
		v2 = ff.positions3D().get(2).minus(ff.positions3D().get(1));
		theta = v1.angleTo(v2);
		assertEquals(1.514,v1.length(),1e-3);
		assertEquals(1.514,v2.length(),1e-3);
		assertEquals(theta,param1.getTheta0(),1e-4);


		// ------- ------- ------- ------- ------- ------- -------
		//
		// Do a series of "special cases" (i.e. test the functional forms
		// for linear, trigonal planar, square planar and octahedral)
		//
		// ------- ------- ------- ------- ------- ------- -------

		// ------- ------- ------- ------- ------- ------- -------
		// test a linear molecule:
		param1.setTheta0(Math.PI);
		ff.contribs().set((int) (ff.contribs().size() - 1), 
				new AngleBendContrib(ff,0,1,2,1,1, param1, param1, param1, 2));
		p1.setX(1.3);
		p1.setY(0.1);
		p1.setZ(0.0);
		p2.setX(0.0);
		p2.setY(0.0);
		p2.setZ(0.0);
		p3.setX(-1.3);
		p3.setY(0.1);
		p3.setZ(0.00);
		ff.initialize();
		ff.minimize(100,1e-8,1e-8);

		v1 = ff.positions3D().get(0).minus(ff.positions3D().get(1));
		v2 = ff.positions3D().get(2).minus(ff.positions3D().get(1));
		theta = v1.angleTo(v2);

		assertEquals(1.514,v1.length(),1e-3);
		assertEquals(1.514,v2.length(),1e-3);
		assertEquals(theta,param1.getTheta0(),1e-4);


		// ------- ------- ------- ------- ------- ------- -------
		// test n=3:
		param1.setTheta0(120.* Math.PI / 180.0);
		ff.contribs().set((int) (ff.contribs().size() - 1), 
				new AngleBendContrib(ff,0,1,2,1,1, param1, param1, param1, 3));

		p1.setX(1.3);
		p1.setY(0.1);
		p1.setZ(0.0);
		p2.setX(0.0);
		p2.setY(0.0);
		p2.setZ(0.0);
		p3.setX(-.3);
		p3.setY(-1.3);
		p3.setZ(0.00);

		ff.initialize();
		ff.minimize(100,1e-8,1e-8);

		v1 = ff.positions3D().get(0).minus(ff.positions3D().get(1));
		v2 = ff.positions3D().get(2).minus(ff.positions3D().get(1));
		theta = v1.angleTo(v2);

		assertEquals(1.514,v1.length(),1e-3);
		assertEquals(1.514,v2.length(),1e-3);
		assertEquals(theta,param1.getTheta0(),1e-4);


	}

	// part of testUFF6 of testForceField.cpp
	@Test
	public void testVDWContribs () {
		ForceField ff = new ForceField();
		Point3D_Vect ps = ff.positions3D();
		Point3D p1 = new Point3D(0,0,0);
		Point3D p2 = new Point3D(0,0,0); 
		ps.add(p1);
		ps.add(p2);
		AtomicParams param1 = new AtomicParams();
		// sp3 carbon:
		param1.setR1(.757);
		param1.setZ1(1.912);
		param1.setGMP_Xi(5.343);
		param1.setX1(3.851);
		param1.setD1(0.105);

		ff.initialize();
		ForceFieldContrib contrib = new vdWContrib (ff,0,1, param1, param1);
		ff.contribs().add(contrib);

		// try a bit of minimization
		ff.initialize();

		// edge case: our energy at zero length should be zero:
		double E;
		E=ff.calcEnergy();
		assertEquals(0.0,E,defaultDoubleTol);
	}
	// from void testUFF7 of testForceFields.cpp
	@Test
	public void testUFFTorsionalTerms () {
		ForceField ff = new ForceField();
		Point3D p1 = new Point3D();
		Point3D p2 = new Point3D();
		Point3D p3 = new Point3D();
		Point3D p4 = new Point3D();;
		ff.positions3D().add(p1);
		ff.positions3D().add(p2);
		ff.positions3D().add(p3);
		ff.positions3D().add(p4);

		AtomicParams param1 = new AtomicParams();
		AtomicParams param2 = new AtomicParams();
		// sp3 carbon:
		param1.setR1(.757);
		param1.setZ1(1.912);
		param1.setGMP_Xi(5.343);
		param1.setX1(3.851);
		param1.setD1(0.105);
		param1.setV1(2.119);
		param1.setU1(2.0);

		// H_1:
		param2.setR1(0.354);
		param2.setZ1(0.712);
		param2.setGMP_Xi(4.528);

		double cosPhi;

		ForceFieldContrib contrib;
		// ------- ------- ------- ------- ------- ------- -------
		// Basic SP3 - SP3
		// ------- ------- ------- ------- ------- ------- -------
		contrib = new TorsionAngleContrib(ff,0,1,2,3,1,
				6,6,
				Atom.HybridizationType.SP3, Atom.HybridizationType.SP3,
				param1, param1);
		ff.contribs().add(contrib);

		p1.setX(0);
		p1.setY(1.5);
		p1.setZ(0);

		p2.setX(0.0);
		p2.setY(0.0);
		p2.setZ(0.0);

		p3.setX(1.5);
		p3.setY(0.0);
		p3.setZ(0.0);

		p4.setX(1.5);
		p4.setY(0.0);
		p4.setZ(1.5);

		ff.initialize();
		ff.minimize(10,1e-8,1e-8);
		cosPhi = RDKFuncs.calculateCosTorsion(ff.positions3D().get(0),
				ff.positions3D().get(1),
				ff.positions3D().get(2),
				ff.positions3D().get(3));
		assertEquals(0.5,cosPhi,1e-4);

		// ------- ------- ------- ------- ------- ------- -------
		// Basic SP2 - SP2
		// ------- ------- ------- ------- ------- ------- -------
		contrib = new TorsionAngleContrib(ff,0,1,2,3,1,
				6,6,
				Atom.HybridizationType.SP2,Atom.HybridizationType.SP2,
				param1, param1);
		ff.contribs().set((int) (ff.contribs().size() - 1), contrib);
		p1.setX(0);
		p1.setY(1.5);
		p1.setZ(0.1);

		p2.setX(0.0);
		p2.setY(0.0);
		p2.setZ(0.0);

		p3.setX(1.5);
		p3.setY(0.0);
		p3.setZ(0.0);

		p4.setX(1.5);
		p4.setY(0.2);
		p4.setZ(1.5);

		ff.initialize();
		ff.minimize(10,1e-8,1e-8);
		cosPhi = RDKFuncs.calculateCosTorsion(ff.positions3D().get(0),
				ff.positions3D().get(1),
				ff.positions3D().get(2),
				ff.positions3D().get(3));
		assertEquals(1.0,cosPhi,1e-4);

		// ------- ------- ------- ------- ------- ------- -------
		// Basic SP2 - SP3
		// ------- ------- ------- ------- ------- ------- -------
		contrib = new TorsionAngleContrib(ff,0,1,2,3,1,
				6,6,
				Atom.HybridizationType.SP2,Atom.HybridizationType.SP3,
				param1, param1);
		ff.contribs().set((int) (ff.contribs().size() - 1), contrib);
		p1.setX(0);
		p1.setY(1.5);
		p1.setZ(0.1);

		p2.setX(0.0);
		p2.setY(0.0);
		p2.setZ(0.0);

		p3.setX(1.5);
		p3.setY(0.0);
		p3.setZ(0.0);

		p4.setX(1.5);
		p4.setY(0.2);
		p4.setZ(1.5);

		ff.initialize();
		ff.minimize(100,1e-8,1e-8);
		cosPhi = RDKFuncs.calculateCosTorsion(ff.positions3D().get(0),
				ff.positions3D().get(1),
				ff.positions3D().get(2),
				ff.positions3D().get(3));
		assertEquals(0.5,cosPhi,1e-4);

		// ------- ------- ------- ------- ------- ------- -------
		// special case for group 6 - group 6 bonds:
		// ------- ------- ------- ------- ------- ------- -------
		contrib = new TorsionAngleContrib(ff,0,1,2,3,1,
				8,8,
				Atom.HybridizationType.SP3,Atom.HybridizationType.SP3,
				param1, param1);
		ff.contribs().set((int) (ff.contribs().size() - 1), contrib);
		p1.setX(0);
		p1.setY(1.5);
		p1.setZ(0.1);

		p2.setX(0.0);
		p2.setY(0.0);
		p2.setZ(0.0);

		p3.setX(1.5);
		p3.setY(0.0);
		p3.setZ(0.0);

		p4.setX(1.5);
		p4.setY(0.2);
		p4.setZ(1.5);

		ff.initialize();
		ff.minimize(100,1e-8,1e-8);
		cosPhi = RDKFuncs.calculateCosTorsion(ff.positions3D().get(0),
				ff.positions3D().get(1),
				ff.positions3D().get(2),
				ff.positions3D().get(3));
		assertEquals(0.0,cosPhi,1e-4);

		// ------- ------- ------- ------- ------- ------- -------
		// special case for SP3 group 6 - SP2 other group
		// ------- ------- ------- ------- ------- ------- -------
		contrib = new TorsionAngleContrib(ff,0,1,2,3,1,
				8,6,
				Atom.HybridizationType.SP3,Atom.HybridizationType.SP2,
				param1, param1);
		ff.contribs().set((int) (ff.contribs().size() - 1), contrib);
		p1.setX(0);
		p1.setY(1.5);
		p1.setZ(0.1);

		p2.setX(0.0);
		p2.setY(0.0);
		p2.setZ(0.0);

		p3.setX(1.5);
		p3.setY(0.0);
		p3.setZ(0.0);

		p4.setX(1.5);
		p4.setY(0.2);
		p4.setZ(1.5);

		ff.initialize();
		ff.minimize(100,1e-8,1e-8);
		cosPhi = RDKFuncs.calculateCosTorsion(ff.positions3D().get(0),
				ff.positions3D().get(1),
				ff.positions3D().get(2),
				ff.positions3D().get(3));
		assertEquals(0.0,cosPhi,1e-4);

		// ------- ------- ------- ------- ------- ------- -------
		// special case for (SP2 -) SP2 - SP3
		// ------- ------- ------- ------- ------- ------- -------
		contrib = new TorsionAngleContrib(ff,0,1,2,3,1,
				6,6,
				Atom.HybridizationType.SP2,Atom.HybridizationType.SP3,
				param1, param1,true);
		ff.contribs().set((int) (ff.contribs().size() - 1), contrib);
		p1.setX(0);
		p1.setY(1.5);
		p1.setZ(0.1);

		p2.setX(0.0);
		p2.setY(0.0);
		p2.setZ(0.0);

		p3.setX(1.5);
		p3.setY(0.0);
		p3.setZ(0.0);

		p4.setX(1.5);
		p4.setY(0.2);
		p4.setZ(1.5);

		ff.initialize();
		ff.minimize(100,1e-8,1e-8);
		cosPhi = RDKFuncs.calculateCosTorsion(ff.positions3D().get(0),
				ff.positions3D().get(1),
				ff.positions3D().get(2),
				ff.positions3D().get(3));
		assertEquals(0.5,cosPhi,1e-4);
	}


	@Test
	public void testUFFParameterObjects () {
		ParamCollection params = ParamCollection.getParams();
		assertNotNull(params);

		AtomicParams param;
		param = params.get("C_3");
		assertNotNull(param);
		assertEquals(0.757, param.getR1(), defaultDoubleTol);
		assertEquals(109.47* Math.PI / 180., param.getTheta0(), defaultDoubleTol);
		assertEquals(3.851, param.getX1(), defaultDoubleTol);
		assertEquals(0.105, param.getD1(), defaultDoubleTol);
		assertEquals(12.73, param.getZeta(), defaultDoubleTol);
		assertEquals(1.912, param.getZ1(), defaultDoubleTol);
		assertEquals(2.119, param.getV1(), defaultDoubleTol);
		assertEquals(5.343, param.getGMP_Xi(), defaultDoubleTol);
		assertEquals(5.063, param.getGMP_Hardness(), defaultDoubleTol);
		assertEquals(0.759, param.getGMP_Radius(), defaultDoubleTol);

		param = params.get("N_3");
		assertNotNull(param);

		param = params.get("C_5");
		assertNull(param);


	}
	// from testUFF*
	@Test
	public void testSimpleUFFMoleculeOptimization () {

		ForceField ff = new ForceField();
		Point3D p1 = new Point3D();
		Point3D p2 = new Point3D();
		Point3D p3 = new Point3D();
		Point3D p4 = new Point3D();
		Point3D p5 = new Point3D();
		Point3D p6 = new Point3D();
		ff.positions3D().add(p1);
		ff.positions3D().add(p2);
		ff.positions3D().add(p3);
		ff.positions3D().add(p4);
		ff.positions3D().add(p5);
		ff.positions3D().add(p6);

		ParamCollection params = ParamCollection.getParams();
		AtomicParams param1,param2;

		// C_2 (sp2 carbon):
		param1 = params.get("C_2");
		assertNotNull(param1);
		// H_:
		param2 = params.get("H_");
		assertNotNull(param2);
		ForceFieldContrib contrib;

		// build ethylene:
		// BONDS:
		contrib = new BondStretchContrib(ff,0,1,2,param1,param1);
		ff.contribs().add(contrib);
		contrib = new BondStretchContrib(ff,0,2,1,param1,param2);
		ff.contribs().add(contrib);
		contrib = new BondStretchContrib(ff,0,3,1,param1,param2);
		ff.contribs().add(contrib);
		contrib = new BondStretchContrib(ff,1,4,1,param1,param2);
		ff.contribs().add(contrib);
		contrib = new BondStretchContrib(ff,1,5,1,param1,param2);
		ff.contribs().add(contrib);

		// ANGLES:
		contrib = new AngleBendContrib(ff,1,0,2,2,1,param1,param1,param2,3);
		ff.contribs().add(contrib);
		contrib = new AngleBendContrib(ff,1,0,3,2,1,param1,param1,param2,3);
		ff.contribs().add(contrib);
		contrib = new AngleBendContrib(ff,2,0,3,1,1,param2,param1,param2,3);
		ff.contribs().add(contrib);
		contrib = new AngleBendContrib(ff,0,1,4,2,1,param1,param1,param2,3);
		ff.contribs().add(contrib);
		contrib = new AngleBendContrib(ff,0,1,5,2,1,param1,param1,param2,3);
		ff.contribs().add(contrib);
		contrib = new AngleBendContrib(ff,4,1,5,1,1,param2,param1,param2,3);
		ff.contribs().add(contrib);

		// DIHEDRALS:
		contrib = new TorsionAngleContrib(ff,2,0,1,4,2,
				6,6,
				Atom.HybridizationType.SP3,Atom.HybridizationType.SP3,
				param1,param1);
		ff.contribs().add(contrib);
		contrib = new TorsionAngleContrib(ff,2,0,1,5,2,
				6,6,
				Atom.HybridizationType.SP3,Atom.HybridizationType.SP3,
				param1,param1);
		ff.contribs().add(contrib);
		contrib = new TorsionAngleContrib(ff,3,0,1,4,2,
				6,6,
				Atom.HybridizationType.SP3,Atom.HybridizationType.SP3,
				param1,param1);
		ff.contribs().add(contrib);
		contrib = new TorsionAngleContrib(ff,3,0,1,5,2,
				6,6,
				Atom.HybridizationType.SP3,Atom.HybridizationType.SP3,
				param1,param1);
		ff.contribs().add(contrib);


		p1.setX(-0.58);
		p1.setY(-0.33);
		p1.setZ(0.1);

		p2.setX(0.58);
		p2.setY(0.33);
		p2.setZ(0.1);

		p3.setX(-0.61);
		p3.setY(-1.43);
		p3.setZ(0.0);

		p4.setX(-1.54);
		p4.setY(0.20);
		p4.setZ(0.0);

		p5.setX(0.61);
		p5.setY(1.43);
		p5.setZ(0.0);

		p6.setX(1.54);
		p6.setY(-0.20);
		p6.setZ(0.0);

		Point3D v1,v2;
		double theta;
		// ------- ------- ------- ------- ------- ------- -------
		// try a bit of minimization
		ff.initialize();
		ff.minimize(100,1e-8,1e-8);

		double CCDblBondLen=RDKFuncs.calcBondRestLength(2,param1,param1);
		double CHBondLen=RDKFuncs.calcBondRestLength(1,param1,param2);

		v1=ff.positions3D().get(0).minus(ff.positions3D().get(1));
		v2=ff.positions3D().get(0).minus(ff.positions3D().get(2));
		theta = v1.angleTo(v2);
		assertEquals(CCDblBondLen,v1.length(),1e-3);
		assertEquals(CHBondLen, v2.length(),1e-3);
		assertEquals(theta,param1.getTheta0(),1e-4);
		v2=ff.positions3D().get(0).minus(ff.positions3D().get(3));
		theta = v1.angleTo(v2);
		assertEquals(v2.length(),CHBondLen,1e-3);
		assertEquals(theta,param1.getTheta0(),1e-4);

		v1=ff.positions3D().get(0).minus(ff.positions3D().get(2));
		theta = v1.angleTo(v2);
		assertEquals(theta,param1.getTheta0(),1e-4);

	}
	// from testHelpers.cpp
	@Test
	public void testUFFAtomTyper2 () {

		ROMol mol,mol2;
		String key;

		mol = RWMol.MolFromSmiles("[SiH3]CC(=O)NC");
		assertNotNull(mol);

		Atomic_Params_Vect types;
		Flagged_Atomic_Params_Vect flaggedParams = ForceField.UFFGetAtomTypes(mol);
		assertTrue(flaggedParams.getSecond());
		types = flaggedParams.getFirst();
		assertEquals(mol.getNumAtoms(), types.size());
		for (int i = 0; i < types.size(); i++ )
			assertNotNull(types.get(i));
		mol2 = mol.addHs(false);

		flaggedParams = ForceField.UFFGetAtomTypes(mol2);
		assertTrue(flaggedParams.getSecond());
		types = flaggedParams.getFirst();
		assertEquals(mol2.getNumAtoms(), types.size());
		for (int i = 0; i < types.size(); i++ )
			assertNotNull(types.get(i));

		// connected with sf.net bug 2094445 
		mol = RWMol.MolFromSmiles("[SiH2]=C");
		assertNotNull(mol);
		key = ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(0));
		assertEquals("Si3", key);
		key = ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(1));
		assertEquals("C_2", key);

		mol = RWMol.MolFromSmiles("[AlH]=C");
		assertNotNull(mol);
		key = ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(0));
		assertEquals("Al3", key);
		key = ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(1));
		assertEquals("C_2", key);

		mol = RWMol.MolFromSmiles("[Mg]=C");
		assertNotNull(mol);
		key = ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(0));
		assertEquals("Mg3+2", key);
		key = ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(1));
		assertEquals("C_2", key);

		mol = RWMol.MolFromSmiles("[SiH3][Si]([SiH3])=C");
		assertNotNull(mol);
		key = ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(0));
		assertEquals("Si3", key);
		key = ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(1));
		assertEquals("Si3", key);
		key = ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(2));
		assertEquals("Si3", key);
		key = ForceField.UFFGetAtomLabel(mol.getAtomWithIdx(3));
		assertEquals("C_2", key);


	}
	@Test
	public void testUFFBuilderSpecialCases(){
		int needMore;
		Point3D v1,v2;

		String basePath = new File(getRdBase(), 
		"/Code/GraphMol/ForceFieldHelpers/UFF/test_data").getPath();
		// ----------
		//  Trigonal bipyramid
		// ----------
		RWMol mol = RWMol.MolFromMolFile(new File(basePath, "tbp.mol").getPath(),false);
		assertNotNull(mol);
		mol.sanitizeMol();

		Conformer conf = mol.getConformer();
		ForceField field = ForceField.UFFGetMoleculeForceField((ROMol)mol);
		assertNotNull(field);
		field.initialize();
		needMore = field.minimize(200,1e-8,1e-4);
		assertEquals(0, needMore);
		v1 = conf.getAtomPos(0).directionVector(conf.getAtomPos(1));
		v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(2));
		assertEquals(v1.dotProduct(v2), -1.0, 1e-3);
		v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(3));
		assertEquals(v1.dotProduct(v2), 0.0, 1e-3);
		v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(4));
		assertEquals(v1.dotProduct(v2), 0.0, 1e-3);
		v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(5));
		assertEquals(v1.dotProduct(v2), 0.0, 1e-3);

		v1 = conf.getAtomPos(0).directionVector(conf.getAtomPos(2));
		v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(3));
		assertEquals(v1.dotProduct(v2), 0.0, 1e-3);
		v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(4));
		assertEquals(v1.dotProduct(v2), 0.0, 1e-3);
		v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(5));
		assertEquals(v1.dotProduct(v2), 0.0, 1e-3);

		v1 = conf.getAtomPos(0).directionVector(conf.getAtomPos(3));
		v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(4));
		assertEquals(v1.dotProduct(v2), -0.5, 1e-3);
		v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(5));
		assertEquals(v1.dotProduct(v2), -0.5, 1e-3);

		v1 = conf.getAtomPos(0).directionVector(conf.getAtomPos(4));
		v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(5));
		assertEquals(v1.dotProduct(v2), -0.5, 1e-3);
	}

	// testIssue239 from testForceFields.cpp
	@Test
	public void testCalcEnergy() {
		@SuppressWarnings("unused")
		int needMore;
		ForceField field;
		double e1,e2;

		String pathName = new File(getRdBase(), 
		"/Code/GraphMol/ForceFieldHelpers/UFF/test_data/Issue239.mol").getPath();
		RWMol mol = RWMol.MolFromMolFile(pathName, true);
		assertNotNull(mol);

		field = ForceField.UFFGetMoleculeForceField((ROMol)mol);
		assertNotNull(field);
		field.initialize();
		needMore = field.minimize(200,1e-6,1e-3);
		e1 = field.calcEnergy();
		needMore = field.minimize(200,1e-6,1e-3);
		e2 = field.calcEnergy();
		assertEquals(e2,e1,0.1);

	}
	@Test
	public void testCalcEnergyWithInputs() {
		@SuppressWarnings("unused")
		int needMore;
		ForceField field;
		@SuppressWarnings("unused")
		double e1,e2;

		String pathName = new File(getRdBase(), 
		"/Code/GraphMol/ForceFieldHelpers/UFF/test_data/Issue239.mol").getPath();
		RWMol mol = RWMol.MolFromMolFile(pathName, true);
		assertNotNull(mol);

		field = ForceField.UFFGetMoleculeForceField((ROMol)mol);
		assertNotNull(field);
		field.initialize();
		needMore = field.minimize(200,1e-6,1e-3);
		e1 = field.calcEnergy();
		needMore = field.minimize(200,1e-6,1e-3);
		e2 = field.calcEnergy();

		// A very basic test of calcEnergy(double[])
		Double_Array p = new Double_Array((int) field.numPoints() * 3);
		for (int i = 0; i < (int) field.numPoints() * 3; i++) {
			p.setitem(i, 0.1 * (i + 3));
		}
		// This is supposed to change the distance matrix
		double dOld = field.distance(0, 1);
		e1 = field.calcEnergy(p.cast());
		assertFalse(Math.abs(dOld - field.distance(0, 1)) < defaultPointTol);
	}
	// from testHelpers.cpp
	@Test
	public void testMissingParams() {
		Atomic_Params_Vect types;
		boolean foundAll;

		RWMol mol = RWMol.MolFromSmiles("[Cu](C)(C)(C)(C)C");
		assertNotNull(mol);

		ROMol mol2 = mol.addHs(false);
		assertTrue(DistanceGeom.EmbedMolecule(mol2) >= 0);

		Flagged_Atomic_Params_Vect flaggedParams = ForceField.UFFGetAtomTypes(mol2);
		foundAll = flaggedParams.getSecond();
		assertTrue(!foundAll);
		types = flaggedParams.getFirst();
		assertEquals(types.size(), mol2.getNumAtoms());
		assertNull(types.get(0));

		// make sure we can optimize anyway:
		ForceField field = ForceField.UFFGetMoleculeForceField(mol2);
		assertNotNull(field);
		field.initialize();
		double e1=field.calcEnergy();
		field.minimize();
		double e2 = field.calcEnergy();
		assertTrue(e2<e1);
	}


	@Test
	public void testMMMFFBasics1() {
            @SuppressWarnings("unused")
		int needMore;

            String pathName = new File(getRdBase(), 
                                       "/Code/GraphMol/ForceFieldHelpers/MMFF/test_data/benzene.mol").getPath();
            RWMol mol = RWMol.MolFromMolFile(pathName, false);
            assertNotNull(mol);
            mol.sanitizeMol();

            ForceField field = ForceField.MMFFGetMoleculeForceField((ROMol)mol);
            assertNotNull(field);
            field.initialize();
            needMore = field.minimize(200);
            double e1 = field.calcEnergy();
            needMore = field.minimize(200);
            double e2 = field.calcEnergy();
            assertTrue(e2<e1);
	}

	@Test
	public void testMMMFFBasics2() {
            String pathName = new File(getRdBase(), 
                                       "/Code/GraphMol/ForceFieldHelpers/MMFF/test_data/benzene.mol").getPath();
            RWMol mol = RWMol.MolFromMolFile(pathName, false);
            assertNotNull(mol);
            mol.sanitizeMol();

            int needMore = ForceField.MMFFOptimizeMolecule((ROMol)mol);
            assertTrue(needMore==0);
	}
    

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.ForceFieldsTests");
	}

}
