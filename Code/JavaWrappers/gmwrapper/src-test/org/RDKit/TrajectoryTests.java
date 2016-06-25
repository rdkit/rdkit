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
        final int ns = 1;
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
    }

    @Test
    public void testTrajectory3D() {
        final int dim = 3;
        final int np = 10;
        final int ns = 1;
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

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.TrajectoryTests");
	}
}
