/*
*  Copyright (C) 2016 Sereina Riniker, Paolo Tosco
*
*    @@ All Rights Reserved @@
*   This file is part of the RDKit.
*   The contents are covered by the terms of the BSD license
*   which is included in the file license.txt, found at the root
*   of the RDKit source tree.
*/

package org.RDKit;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;

import org.junit.Test;

public class TrajectoryTests extends GraphMolTest {
    @Test
    public void testBasicInstantiation_Snapshot() {
        Shared_Double_Array pos = new Shared_Double_Array();
        Snapshot s = new Snapshot(pos);
        assertNotNull(s);
    }

    @Test
    public void testSnapshot() {
        Shared_Double_Array pos = new Shared_Double_Array();
        Snapshot s = new Snapshot(pos);
        boolean e = false;
        try {
            s.getPoint2D(12);
        }
        catch (GenericRDKitException ex) {
            e = true;
        }
        assertEquals(true, e);
    }

    @Test
    public void testTrajectory2D() {
        final int dim = 2;
        final int np = 10;
        final int ns = 5;
        Trajectory traj = new Trajectory(dim, np);
        assertEquals(traj.dimension(), dim);
        assertEquals(traj.numPoints(), np);
        final int posLen = np * dim;
        Double_Vect dv = new Double_Vect(posLen);
        for (int i = 0; i < posLen; ++i)
            dv.set(i, (double)i);
        for (int i = 0; i < ns; ++i)
            traj.addSnapshot(new Snapshot(dv, (double)i));
        assertEquals(traj.size(), ns);
        boolean e = false;
        try {
            traj.getSnapshot(ns);
        }
        catch (GenericRDKitException ex) {
            e = true;
        }
        assertEquals(true, e);
        e = false;
        try {
            traj.getSnapshot(0).getPoint2D(np);
        }
        catch (GenericRDKitException ex) {
            e = true;
        }
        assertEquals(true, e);
        for (int i = 0; i < np; ++i) {
            assertEquals(traj.getSnapshot(0).getPoint2D(i).getX(), (double)(i * dim), 0.001);
            assertEquals(traj.getSnapshot(0).getPoint2D(i).getY(), (double)(i * dim + 1), 0.001);
            e = false;
            try {
                assertEquals(traj.getSnapshot(0).getPoint3D(i).getZ(), 0.0, 0.001);
            }
            catch (GenericRDKitException ex) {
                e = true;
            }
            assertEquals(false, e);
        }
        for (int i = 0; i < ns; ++i)
            assertEquals((double)(i), traj.getSnapshot(i).getEnergy(), 0.001);
        traj.removeSnapshot(0);
        assertEquals(ns - 1, traj.size());
        for (int i = 0; i < (ns - 1); ++i)
            assertEquals((double)(i + 1), traj.getSnapshot(i).getEnergy(), 0.001);
        traj.insertSnapshot(0, new Snapshot(dv, 999.0));
        assertEquals(ns, traj.size());
        Snapshot copySnapshot = new Snapshot(traj.getSnapshot(0));
        traj.addSnapshot(copySnapshot);
        assertEquals(ns + 1, traj.size());
        assertEquals(999.0, traj.getSnapshot(0).getEnergy(), 0.001);
        assertEquals(1.0, traj.getSnapshot(1).getEnergy(), 0.001);
        assertEquals(999.0, traj.getSnapshot(traj.size() - 1).getEnergy(), 0.001);
        Trajectory traj2 = new Trajectory(traj);
        assertEquals(traj2.size(), traj.size());
    }

    @Test
    public void testTrajectory3D() {
        final int dim = 3;
        final int np = 10;
        final int ns = 5;
        Trajectory traj = new Trajectory(dim, np);
        assertEquals(dim, traj.dimension());
        assertEquals(np, traj.numPoints());
        final int posLen = np * dim;
        Double_Vect dv = new Double_Vect(posLen);
        for (int i = 0; i < posLen; ++i)
            dv.set(i, (double)i);
        for (int i = 0; i < ns; ++i)
            traj.addSnapshot(new Snapshot(dv, (double)i));
        assertEquals(ns, traj.size());
        boolean e = false;
        try {
            traj.getSnapshot(ns);
        }
        catch (GenericRDKitException ex) {
            e = true;
        }
        assertEquals(true, e);
        e = false;
        try {
            traj.getSnapshot(0).getPoint2D(np);
        }
        catch (GenericRDKitException ex) {
            e = true;
        }
        assertEquals(true, e);
        for (int i = 0; i < np; ++i) {
            assertEquals((double)(i * dim), traj.getSnapshot(0).getPoint3D(i).getX(), 0.001);
            assertEquals((double)(i * dim + 1), traj.getSnapshot(0).getPoint3D(i).getY(), 0.001);
            assertEquals((double)(i * dim + 2), traj.getSnapshot(0).getPoint3D(i).getZ(), 0.001);
            if (i == 0) {
                e = false;
                try {
                    traj.getSnapshot(0).getPoint2D(i);
                }
                catch (GenericRDKitException ex) {
                    e = true;
                }
                assertEquals(true, e);
            }
        }
        for (int i = 0; i < ns; ++i)
            assertEquals((double)(i), traj.getSnapshot(i).getEnergy(), 0.001);
        traj.removeSnapshot(0);
        assertEquals(ns - 1, traj.size());
        for (int i = 0; i < (ns - 1); ++i)
            assertEquals((double)(i + 1), traj.getSnapshot(i).getEnergy(), 0.001);
        traj.insertSnapshot(0, new Snapshot(dv, 999.0));
        assertEquals(ns, traj.size());
        Snapshot copySnapshot = new Snapshot(traj.getSnapshot(0));
        traj.addSnapshot(copySnapshot);
        assertEquals(ns + 1, traj.size());
        assertEquals(999.0, traj.getSnapshot(0).getEnergy(), 0.001);
        assertEquals(1.0, traj.getSnapshot(1).getEnergy(), 0.001);
        assertEquals(999.0, traj.getSnapshot(traj.size() - 1).getEnergy(), 0.001);
        Trajectory traj2 = new Trajectory(traj);
        assertEquals(traj2.size(), traj.size());
    }

    @Test
    public void testReadAmber() {
		String rdpath = System.getenv("RDBASE");
		if (rdpath == null)
			org.junit.Assert.fail("No definition for RDBASE");
		File base = new File(rdpath);
		File testFile = new File(base, "Code" + File.separator + "GraphMol"
				+ File.separator + "test_data" + File.separator + "water_coords_bad.trx");
		String fName = testFile.getAbsolutePath();
        Trajectory traj = new Trajectory(2, 0);
        boolean e = false;
        try {
            RDKFuncs.readAmberTrajectory(fName, traj);
        }
        catch (GenericRDKitException ex) {
            e = true;
        }
        assertEquals(true, e);
        traj = new Trajectory(3, 3);
        e = false;
        try {
            RDKFuncs.readAmberTrajectory(fName, traj);
        }
        catch (GenericRDKitException ex) {
            e = true;
        }
        assertEquals(true, e);
		testFile = new File(base, "Code" + File.separator + "GraphMol"
				+ File.separator + "test_data" + File.separator + "water_coords_bad2.trx");
		fName = testFile.getAbsolutePath();
        e = false;
        traj = new Trajectory(3, 3);
        try {
            RDKFuncs.readAmberTrajectory(fName, traj);
        }
        catch (GenericRDKitException ex) {
            e = true;
        }
        assertEquals(true, e);
		testFile = new File(base, "Code" + File.separator + "GraphMol"
				+ File.separator + "test_data" + File.separator + "water_coords.trx");
		fName = testFile.getAbsolutePath();
        traj = new Trajectory(3, 3);
        RDKFuncs.readAmberTrajectory(fName, traj);
        assertEquals(traj.size(), 1);
		testFile = new File(base, "Code" + File.separator + "GraphMol"
				+ File.separator + "test_data" + File.separator + "water_coords2.trx");
		fName = testFile.getAbsolutePath();
        traj = new Trajectory(3, 3);
        RDKFuncs.readAmberTrajectory(fName, traj);
        assertEquals(traj.size(), 2);
    }

    @Test
    public void testReadAmberJava() {
        /*
        reimplemented the Amber trajectory reader in Java
        let's check we get the same data as the C++ reader
        (test for building a trajectory out of Snapshots from Java)
        */
		String rdpath = System.getenv("RDBASE");
		if (rdpath == null)
			org.junit.Assert.fail("No definition for RDBASE");
		File base = new File(rdpath);
		File testFile = new File(base, "Code" + File.separator + "GraphMol"
				+ File.separator + "test_data" + File.separator + "water_coords2.trx");
		String fName = testFile.getAbsolutePath();
        Trajectory traj = new Trajectory(3, 3);
        long nCoords = traj.numPoints() * 3;
        int nSnapshots = 0;
        BufferedReader r = null;
        try {
            r = new BufferedReader(new FileReader(testFile));
        }
        catch (Exception e) {
            org.junit.Assert.fail("Could not open " + fName);
        }
        int lineNum = 0;
        Double_Vect dv = new Double_Vect();
        int i = 0;
        String line = null;
        try {
            line = r.readLine();
        }
        catch (Exception e) {
            org.junit.Assert.fail("Could not read from " + fName);
        }
        while (line != null) {
            ++lineNum;
            if (lineNum > 1) {
                String[] tok = line.split(" +");
                int j = 0;
                for (; (i < nCoords) && (j < tok.length); ++j) {
                    if (tok[j].length() == 0) continue;
                    dv.add(Double.parseDouble(tok[j]));
                    ++i;
                }
                line = "";
                if (i == nCoords) {
                    ++nSnapshots;
                    traj.addSnapshot(new Snapshot(dv));
                    dv.clear();
                    i = 0;
                    for (; j < tok.length; ++j) {
                        if (tok[j].length() > 0)
                            line += tok[j] + " ";
                    }
                }
            }
            else {
                line = "";
            }
            try {
                String lineNew = r.readLine();
                if (lineNew != null)
                    line += lineNew;
                else if (line.length() == 0)
                    line = null;
            }
            catch (Exception e) {
                org.junit.Assert.fail("Could not read from " + fName);
            }
        }
        try {
            r.close();
        }
        catch (Exception e) {
            org.junit.Assert.fail("Could not close " + fName);
        }
        assertEquals(0, i);
        assertEquals(2, nSnapshots);
        Trajectory traj2 = new Trajectory(3, 3);
        RDKFuncs.readAmberTrajectory(fName, traj2);
        assertEquals(traj2.size(), traj.size());
        assertEquals(traj2.numPoints(), traj.numPoints());
        for (int snapshotNum = 0; snapshotNum < traj.size(); ++snapshotNum) {
            for (int pointNum = 0; pointNum < traj.numPoints(); ++pointNum) {
                assertEquals(traj2.getSnapshot(snapshotNum).getPoint3D(pointNum).getX(),
                    traj.getSnapshot(snapshotNum).getPoint3D(pointNum).getX(), 0.001);
                assertEquals(traj2.getSnapshot(snapshotNum).getPoint3D(pointNum).getY(),
                    traj.getSnapshot(snapshotNum).getPoint3D(pointNum).getY(), 0.001);
                assertEquals(traj2.getSnapshot(snapshotNum).getPoint3D(pointNum).getZ(),
                    traj.getSnapshot(snapshotNum).getPoint3D(pointNum).getZ(), 0.001);
            }
        }
    }

    @Test
    public void testReadGromos() {
		String rdpath = System.getenv("RDBASE");
		if (rdpath == null)
			org.junit.Assert.fail("No definition for RDBASE");
		File base = new File(rdpath);
		File testFile = new File(base, "Code" + File.separator + "GraphMol"
				+ File.separator + "test_data" + File.separator + "water_coords_bad.trc");
		String fName = testFile.getAbsolutePath();
        Trajectory traj = new Trajectory(2, 0);
        boolean e = false;
        try {
            RDKFuncs.readGromosTrajectory(fName, traj);
        }
        catch (GenericRDKitException ex) {
            e = true;
        }
        assertEquals(true, e);
        traj = new Trajectory(3, 3);
        e = false;
        try {
            RDKFuncs.readGromosTrajectory(fName, traj);
        }
        catch (GenericRDKitException ex) {
            e = true;
        }
        assertEquals(true, e);
		testFile = new File(base, "Code" + File.separator + "GraphMol"
				+ File.separator + "test_data" + File.separator + "water_coords_bad2.trc");
		fName = testFile.getAbsolutePath();
        e = false;
        traj = new Trajectory(3, 3);
        try {
            RDKFuncs.readGromosTrajectory(fName, traj);
        }
        catch (GenericRDKitException ex) {
            e = true;
        }
        assertEquals(true, e);
		testFile = new File(base, "Code" + File.separator + "GraphMol"
				+ File.separator + "test_data" + File.separator + "water_coords.trc");
		fName = testFile.getAbsolutePath();
        traj = new Trajectory(3, 3);
        RDKFuncs.readGromosTrajectory(fName, traj);
        assertEquals(traj.size(), 1);
		testFile = new File(base, "Code" + File.separator + "GraphMol"
				+ File.separator + "test_data" + File.separator + "water_coords2.trc");
		fName = testFile.getAbsolutePath();
        traj = new Trajectory(3, 3);
        RDKFuncs.readGromosTrajectory(fName, traj);
        assertEquals(traj.size(), 2);
    }

    @Test
    public void testAddConformersFromTrajectory() {
        String molBlock =
            "\n" +
            "     RDKit          3D\n" +
            "\n" +
            " 71 74  0  0  0  0  0  0  0  0999 V2000\n" +
            "    8.2543    3.1901   -0.3005 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    7.4558    1.9712    0.0938 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    7.3934    1.0441   -0.9483 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    6.6660   -0.0533   -0.4641 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    5.1928    0.2346   -0.4609 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    4.3713   -0.9410   -0.5770 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    3.1852   -1.0034   -1.2291 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    2.2914    0.1276   -1.6316 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    0.9308   -0.4468   -1.9908 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    0.1417   -0.7821   -0.7545 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -0.1848    0.3695    0.0456 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -1.5661    0.7686   -0.0745 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -2.4768   -0.0640    0.8206 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -3.8874    0.1143    0.3941 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -4.6333   -0.9984    0.0264 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -6.0127   -0.9516   -0.0400 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -6.7062    0.1599    0.3963 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -8.0408    0.4828   -0.1977 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -7.7914    1.1180   -1.5591 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -8.7622    1.4403    0.7265 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -8.8409   -0.7397   -0.4395 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -8.9121   -1.6637    0.4258 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -9.7414   -0.7636   -1.5059 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -5.9736    1.2357    0.8565 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -4.5843    1.2252    0.8530 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    0.6263    1.4884   -0.3942 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    2.0541    1.0258   -0.4230 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    2.9225   -2.3317   -1.2963 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    3.6061   -2.9745   -0.3180 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    3.3554   -4.1536    0.3735 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    3.7653   -4.2712    1.6948 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    4.8254   -3.4613    2.0796 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    5.1978   -2.3436    1.3419 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    4.5694   -2.0799    0.1305 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    9.3138    3.1372    0.0031 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    7.8117    4.0754    0.1798 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    8.2358    3.3535   -1.4074 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    6.4027    2.2146    0.3634 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    7.9270    1.5444    1.0040 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    7.0677   -0.2415    0.5615 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    6.9530   -0.9105   -1.1025 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    4.9578    0.7259    0.5137 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    4.9985    0.9430   -1.3033 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    2.7171    0.7264   -2.4494 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    0.3994    0.2339   -2.6810 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    1.1342   -1.4171   -2.5076 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -0.7632   -1.3370   -1.0391 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    0.7845   -1.4394   -0.1311 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    0.0125    0.1989    1.0673 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -1.6672    1.8215    0.2925 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -1.8705    0.7271   -1.1337 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -2.3045    0.3159    1.8590 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -2.1980   -1.1367    0.7635 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -4.1513   -1.9468   -0.2114 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -6.6138   -1.7460   -0.4718 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -7.0727    0.4399   -2.0858 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -7.3144    2.1076   -1.4482 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -8.7609    1.1720   -2.1135 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -8.3137    2.4504    0.5729 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -8.6170    1.0817    1.7580 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -9.8244    1.4444    0.4200 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -6.4629    2.0541    1.3719 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   -4.0445    2.0563    1.3058 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    0.3329    1.8224   -1.3991 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    0.4920    2.3164    0.3160 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    2.2025    0.3766    0.4766 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    2.7945    1.8369   -0.3969 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    2.4404   -4.6964    0.1303 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    3.3157   -5.0055    2.3587 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    5.4272   -3.7654    2.9380 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    5.5668   -1.5069    1.9380 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "  1  2  1  0\n" +
            "  2  3  1  0\n" +
            "  3  4  1  0\n" +
            "  4  5  1  0\n" +
            "  5  6  1  0\n" +
            "  6  7  1  0\n" +
            "  7  8  1  0\n" +
            "  8  9  1  0\n" +
            "  9 10  1  0\n" +
            " 10 11  1  0\n" +
            " 11 12  1  0\n" +
            " 12 13  1  0\n" +
            " 13 14  1  0\n" +
            " 14 15  2  0\n" +
            " 15 16  1  0\n" +
            " 16 17  2  0\n" +
            " 17 18  1  0\n" +
            " 18 19  1  0\n" +
            " 18 20  1  0\n" +
            " 18 21  1  0\n" +
            " 21 22  2  0\n" +
            " 21 23  1  0\n" +
            " 17 24  1  0\n" +
            " 24 25  2  0\n" +
            " 11 26  1  0\n" +
            " 26 27  1  0\n" +
            "  7 28  2  0\n" +
            " 28 29  1  0\n" +
            " 29 30  2  0\n" +
            " 30 31  1  0\n" +
            " 31 32  2  0\n" +
            " 32 33  1  0\n" +
            " 33 34  2  0\n" +
            " 34  6  1  0\n" +
            " 27  8  1  0\n" +
            " 34 29  1  0\n" +
            " 25 14  1  0\n" +
            "  1 35  1  0\n" +
            "  1 36  1  0\n" +
            "  1 37  1  0\n" +
            "  2 38  1  0\n" +
            "  2 39  1  0\n" +
            "  4 40  1  0\n" +
            "  4 41  1  0\n" +
            "  5 42  1  0\n" +
            "  5 43  1  0\n" +
            "  8 44  1  0\n" +
            "  9 45  1  0\n" +
            "  9 46  1  0\n" +
            " 10 47  1  0\n" +
            " 10 48  1  0\n" +
            " 11 49  1  0\n" +
            " 12 50  1  0\n" +
            " 12 51  1  0\n" +
            " 13 52  1  0\n" +
            " 13 53  1  0\n" +
            " 15 54  1  0\n" +
            " 16 55  1  0\n" +
            " 19 56  1  0\n" +
            " 19 57  1  0\n" +
            " 19 58  1  0\n" +
            " 20 59  1  0\n" +
            " 20 60  1  0\n" +
            " 20 61  1  0\n" +
            " 24 62  1  0\n" +
            " 25 63  1  0\n" +
            " 26 64  1  0\n" +
            " 26 65  1  0\n" +
            " 27 66  1  0\n" +
            " 27 67  1  0\n" +
            " 30 68  1  0\n" +
            " 31 69  1  0\n" +
            " 32 70  1  0\n" +
            " 33 71  1  0\n" +
            "M  CHG  2  11   1  23  -1\n" +
            "M  END\n";
        ROMol mol = RWMol.MolFromMolBlock(molBlock, true, false);
        final int everySteps = 10;
        final int maxIts = 1000;
        final double gradTol = 0.01;
		String rdpath = System.getenv("RDBASE");
		if (rdpath == null)
			org.junit.Assert.fail("No definition for RDBASE");
		File base = new File(rdpath);
		File testFile = new File(base, "Code" + File.separator + "GraphMol"
				+ File.separator + "test_data" + File.separator + "bilastine_trajectory_java.sdf");
		String fName = testFile.getAbsolutePath();
        SDWriter w = new SDWriter(fName);
        ForceField field = ForceField.MMFFGetMoleculeForceField(mol);
        field.initialize();
        SWIGTYPE_p_std__vectorT_RDKit__Snapshot_t sv = Snapshot.SnapshotVect();
        int res = field.minimize(everySteps, sv, maxIts, gradTol);
        assertEquals(0, res);
        Trajectory traj = new Trajectory(3, mol.getNumAtoms(), sv);
        mol.removeConformer(0);
        traj.addConformersToMol(mol);
        for (int nConf = 0; nConf < mol.getNumConformers(); ++nConf) {
            String ss = String.format("%.4f", traj.getSnapshot(nConf).getEnergy());
            mol.setProp("ENERGY", ss, false);
            w.write(mol, nConf);
        }
        w.close();
        traj.clear();
        long n1 = mol.getNumConformers();
        traj.addConformersToMol(mol);
        long n2 = mol.getNumConformers();
        assertEquals(n2, n1);
        // getSnapshot should raise exception after Clear()
        boolean e = false;
        try {
            traj.getSnapshot(0);
        }
        catch (GenericRDKitException ex) {
            e = true;
        }
        assertEquals(true, e);
    }

    @Test
    public void testAddConformersFromAmberTrajectory() {
        ROMol mol = RWMol.MolFromSmiles("CCC");
		String rdpath = System.getenv("RDBASE");
		if (rdpath == null)
			org.junit.Assert.fail("No definition for RDBASE");
		File base = new File(rdpath);
		File testFile = new File(base, "Code" + File.separator + "GraphMol"
				+ File.separator + "test_data" + File.separator + "water_coords.trx");
		String fName = testFile.getAbsolutePath();
        {
            Trajectory traj = new Trajectory(3, mol.getNumAtoms());
            RDKFuncs.readAmberTrajectory(fName, traj);
            assertEquals(1, traj.size());
            for (int i = 0; i < 2; ++i) {
                traj.addConformersToMol(mol);
                assertEquals(i + 1, mol.getNumConformers());
                assertEquals(3, mol.getConformer(i).getNumAtoms());
                assertEquals(0.1941767, mol.getConformer(i).getAtomPos(0).getX(), 0.001);
                assertEquals(-0.4088006, mol.getConformer(i).getAtomPos(2).getZ(), 0.001);
            }
            mol.clearConformers();
            boolean e = false;
            try {
                traj.addConformersToMol(mol, 1);
            }
            catch (GenericRDKitException ex) {
                e = true;
            }
            assertEquals(true, e);
            assertEquals(0, mol.getNumConformers());
		}
        testFile = new File(base, "Code" + File.separator + "GraphMol"
				+ File.separator + "test_data" + File.separator + "water_coords2.trx");
        fName = testFile.getAbsolutePath();
        {
            Trajectory traj = new Trajectory(3, mol.getNumAtoms());
            RDKFuncs.readAmberTrajectory(fName, traj);
            assertEquals(2, traj.size());
            traj.addConformersToMol(mol);
            assertEquals(2, mol.getNumConformers());
            mol.clearConformers();
            traj.addConformersToMol(mol, 0, 0);
            assertEquals(1, mol.getNumConformers());
            traj.addConformersToMol(mol, 1);
            assertEquals(2, mol.getNumConformers());
        }
    }

    @Test
    public void testAddConformersFromGromosTrajectory() {
        ROMol mol = RWMol.MolFromSmiles("CCC");
		String rdpath = System.getenv("RDBASE");
		if (rdpath == null)
			org.junit.Assert.fail("No definition for RDBASE");
		File base = new File(rdpath);
		File testFile = new File(base, "Code" + File.separator + "GraphMol"
				+ File.separator + "test_data" + File.separator + "water_coords.trc");
		String fName = testFile.getAbsolutePath();
        {
            Trajectory traj = new Trajectory(3, mol.getNumAtoms());
            RDKFuncs.readGromosTrajectory(fName, traj);
            assertEquals(1, traj.size());
            for (int i = 0; i < 2; ++i) {
                traj.addConformersToMol(mol);
                assertEquals(i + 1, mol.getNumConformers());
                assertEquals(3, mol.getConformer(i).getNumAtoms());
                assertEquals(1.941767, mol.getConformer(i).getAtomPos(0).getX(), 0.001);
                assertEquals(-4.088006, mol.getConformer(i).getAtomPos(2).getZ(), 0.001);
            }
            mol.clearConformers();
            boolean e = false;
            try {
                traj.addConformersToMol(mol, 1);
            }
            catch (GenericRDKitException ex) {
                e = true;
            }
            assertEquals(true, e);
            assertEquals(0, mol.getNumConformers());
		}
        testFile = new File(base, "Code" + File.separator + "GraphMol"
				+ File.separator + "test_data" + File.separator + "water_coords2.trc");
        fName = testFile.getAbsolutePath();
        {
            Trajectory traj = new Trajectory(3, mol.getNumAtoms());
            RDKFuncs.readGromosTrajectory(fName, traj);
            assertEquals(2, traj.size());
            traj.addConformersToMol(mol);
            assertEquals(2, mol.getNumConformers());
            mol.clearConformers();
            traj.addConformersToMol(mol, 0, 0);
            assertEquals(1, mol.getNumConformers());
            traj.addConformersToMol(mol, 1);
            assertEquals(2, mol.getNumConformers());
        }
    }

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.TrajectoryTests");
	}
}
