package org.RDKit;

import static org.junit.Assert.*;

import java.io.*;
import java.util.ArrayList;

import org.junit.*;

public class GzStreamTests extends GraphMolTest {

	@Test
	public void test11GZstream() {
		// NCI_aids_few.sdf.gz
		File base = getRdBase();
		File gzpath  = new File(base, "Code" + File.separator + "GraphMol" + File.separator +
								"FileParsers" + File.separator + "test_data");
		File fileN = new File(gzpath, "NCI_aids_few.sdf.gz");
		assertTrue(fileN.exists());
		gzstream stream = new gzstream(fileN.getPath());
		ForwardSDMolSupplier suppl = new ForwardSDMolSupplier(stream);
		assertFalse(suppl.atEnd());
		ArrayList<ROMol> ms = new ArrayList<ROMol>();
		ROMol m;
		do {
		m = suppl.next();
		if (m != null)
			ms.add(m);
		} while (!suppl.atEnd());
		assertEquals(16, ms.size());
	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.GzStreamTests");
	}

}
