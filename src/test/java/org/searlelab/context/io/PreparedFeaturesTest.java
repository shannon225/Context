package org.searlelab.context.io;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collections;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import org.searlelab.context.datastructures.FeatureTable;
import org.searlelab.context.datastructures.PreparedFeatures;

public class PreparedFeaturesTest {

	@Rule
	public TemporaryFolder folder = new TemporaryFolder();

	private static final String HEADER = "id\tLabel\tScanNr\tfeatA\tconstFeat\tsequence\tProteins";

	private File background() throws IOException {
		return write("run_background.features.txt", HEADER,
				"bg1\t1\t10\t1.5\t7\tPEPTIDEK\tsp|P1|A",
				"bg2\t-1\t11\t0.5\t7\tKEDITPEP\tDECOY_sp|P1|A",
				"bg3\t1\t12\t2.5\t7\tOTHERK\tsp|P2|B",
				"bg4\t-1\t13\t0.2\t7\tKREHTO\tDECOY_sp|P2|B");
	}

	private File reference() throws IOException {
		return write("run_reference.features.txt", HEADER,
				"ref1\t1\t20\t3.5\t7\tTARGETK\tsp|P3|C",
				"ref2\t-1\t21\t0.1\t7\tKTEGRAT\tDECOY_sp|P3|C");
	}

	@Test
	public void prunesTheSameColumnsFromBothHalves() throws IOException {
		PreparedFeatures prepared = PreparedFeatures.prepare(background(), reference(), folder.getRoot(), "run");

		assertEquals(Collections.singletonList("constFeat"), prepared.getDroppedFeatures());

		assertEquals(Arrays.asList("featA"), FeatureTable.read(prepared.getBackground()).getFeatureNames());
		assertEquals(Arrays.asList("featA"), prepared.getReferenceTable().getFeatureNames());
	}

	@Test
	public void constantColumnsAreJudgedOnTheBackgroundOnly() throws IOException {
		File reference = write("run_reference.features.txt", HEADER,
				"ref1\t1\t20\t9.9\t7\tTARGETK\tsp|P3|C",
				"ref2\t-1\t21\t9.9\t7\tKTEGRAT\tDECOY_sp|P3|C");

		PreparedFeatures prepared = PreparedFeatures.prepare(background(), reference, folder.getRoot(), "run");

		assertEquals(Arrays.asList("featA"), prepared.getReferenceTable().getFeatureNames());
	}

	@Test
	public void countsBothHalves() throws IOException {
		PreparedFeatures prepared = PreparedFeatures.prepare(background(), reference(), folder.getRoot(), "run");

		assertEquals(2, prepared.getBackgroundTargets());
		assertEquals(2, prepared.getBackgroundDecoys());
		assertEquals(1, prepared.getReferenceTargets());
		assertEquals(1, prepared.getReferenceDecoys());
	}

	@Test
	public void writesThePrunedTablesWhereItSaysItDoes() throws IOException {
		PreparedFeatures prepared = PreparedFeatures.prepare(background(), reference(), folder.getRoot(), "run");

		assertEquals(new File(folder.getRoot(), "run.background.pin"), prepared.getBackground());
		assertEquals(new File(folder.getRoot(), "run.reference.pin"), prepared.getReference());
		assertTrue(prepared.getBackground().isFile());
		assertTrue(prepared.getReference().isFile());
	}

	@Test
	public void rejectsAReferencePanelWithNoDecoys() throws IOException {
		File reference = write("run_reference.features.txt", HEADER,
				"ref1\t1\t20\t3.5\t7\tTARGETK\tsp|P3|C");

		try {
			PreparedFeatures.prepare(background(), reference, folder.getRoot(), "run");
			fail("expected a targets-only reference panel to be rejected");
		} catch (IOException expected) {
			assertTrue(expected.getMessage().contains("no decoys"));
		}
	}

	@Test
	public void rejectsABackgroundWithOnlyOneClass() throws IOException {
		File background = write("run_background.features.txt", HEADER,
				"bg1\t1\t10\t1.5\t7\tPEPTIDEK\tsp|P1|A",
				"bg3\t1\t12\t2.5\t7\tOTHERK\tsp|P2|B");

		try {
			PreparedFeatures.prepare(background, reference(), folder.getRoot(), "run");
			fail("expected a background without decoys to be rejected");
		} catch (IOException expected) {
			assertTrue(expected.getMessage().contains("training needs both"));
		}
	}

	@Test
	public void rejectsHalvesFromDifferentRuns() throws IOException {
		File reference = write("run_reference.features.txt",
				"id\tLabel\tScanNr\tfeatZ\tsequence\tProteins",
				"ref1\t1\t20\t3.5\tTARGETK\tsp|P3|C",
				"ref2\t-1\t21\t0.1\tKTEGRAT\tDECOY_sp|P3|C");

		try {
			PreparedFeatures.prepare(background(), reference, folder.getRoot(), "run");
			fail("expected mismatched headers to be rejected");
		} catch (IOException expected) {
			assertTrue(expected.getMessage().contains("same columns"));
		}
	}

	private File write(String name, String... lines) throws IOException {
		File file = new File(folder.getRoot(), name);
		try (PrintWriter out = new PrintWriter(file, "UTF-8")) {
			for (String line : lines) {
				out.println(line);
			}
		}
		return file;
	}
}
