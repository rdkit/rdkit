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
        SWIGTYPE_p_double pos = RDKFuncs.new_SWIGArrayUtility(0);
        Snapshot s = new Snapshot(pos);
        assertNotNull(s);
    }

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.TrajectoryTests");
	}
}
