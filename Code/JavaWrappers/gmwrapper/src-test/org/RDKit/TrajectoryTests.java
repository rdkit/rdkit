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

import org.RDKit.RDKFuncsJNI;

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
        System.out.println("ok1");
        Double_Array da = new Double_Array(posLen);
        System.out.println("ok2, da = " + da);
        System.out.println("ok3");
        for (int i = 0; i < posLen; ++i) {
            da.setitem(i, (double)i);
            System.out.println("i = " + i + ", v = " + da.getitem(i));
        }
        System.out.println("ok4");
        for (int i = 0; i < ns; ++i) {
            System.out.println("ok5, i = " + i);
            traj.addSnapshot(new Snapshot(da.cast(), (double)i));
            System.out.println("ok6");
        }
        /*
        System.out.println("ok5");
        assertEquals(traj.size(), ns);
        System.out.println("ok6");
        boolean e = false;
        try {
            traj.getSnapshot(ns);
        }
        catch (GenericRDKitException ex) {
            e = true;
        }
        System.out.println("ok7");
        assertEquals(true, e);
        System.out.println("ok8");
        e = false;
        try {
            traj.getSnapshot(0).getPoint2D(np);
        }
        catch (GenericRDKitException ex) {
            e = true;
        }
        System.out.println("ok9");
        assertEquals(true, e);
        System.out.println("ok10");
        */
    }

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.TrajectoryTests");
	}
}
