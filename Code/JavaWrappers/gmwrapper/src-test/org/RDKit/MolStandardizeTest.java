package org.RDKit;

import static org.junit.Assert.*;

import java.io.File;

import org.junit.*;

/**
 * JUnit tests for MolStandardize wrappers.
 */
 public class MolStandardizeTest extends GraphMolTest {
	@Test
	public void testStandardize1() {
		assertEquals("fail", RDKFuncs.standardizeSmiles("[Na]OC(=O)c1ccccc1"),"O=C([O-])c1ccccc1.[Na+]");
	}
	@Test
	public void testPipelineBadInput() {
		Pipeline pipeline = new Pipeline();

		PipelineResult result = pipeline.run(
			"\n" +
            "             sldfj;ldskfj sldkjfsd;lkf\n" +
            "M  V30 BEGIN CTAB"
		);

	    assertEquals(result.getStage(), PipelineStage.PARSING_INPUT);
        assertFalse(result.getStatus() == PipelineStatus.NO_EVENT);
        assertTrue((result.getStatus() & PipelineStatus.INPUT_ERROR) != PipelineStatus.NO_EVENT);

		result.delete();
		pipeline.delete();
	}
	@Test
	public void testPipelineUnsupportedFeatures() {
		Pipeline pipeline = new Pipeline();

		PipelineResult result = pipeline.run(
			"\n" +
            "  Mrv2311 01162413552D          \n" +
            "\n" +
            "  0  0  0     0  0            999 V3000\n" +
            "M  V30 BEGIN CTAB\n" +
            "M  V30 COUNTS 2 1 0 0 0\n" +
            "M  V30 BEGIN ATOM\n" +
            "M  V30 1 R# -17.3747 6.9367 0 0 RGROUPS=(1 0)\n" +
            "M  V30 2 C -18.7083 6.1667 0 0\n" +
            "M  V30 END ATOM\n" +
            "M  V30 BEGIN BOND\n" +
            "M  V30 1 1 2 1\n" +
            "M  V30 END BOND\n" +
            "M  V30 END CTAB\n" +
            "M  END"
		);

	    assertEquals(result.getStage(), PipelineStage.COMPLETED);
        assertFalse(result.getStatus() == PipelineStatus.NO_EVENT);
        assertTrue((result.getStatus() & PipelineStatus.VALIDATION_ERROR) != PipelineStatus.NO_EVENT);
        assertTrue((result.getStatus() & PipelineStatus.FEATURES_VALIDATION_ERROR) != PipelineStatus.NO_EVENT);

		result.delete();
		pipeline.delete();
	}
	public static void main(String args[]) {
						org.junit.runner.JUnitCore.main("org.RDKit.MolStandardizeTest");
	}
}
