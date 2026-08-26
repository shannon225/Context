package org.searlelab.context.io;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import org.searlelab.context.datastructures.EncyclopediaFeatures;

public class EncyclopediaFeaturesTest {
	private static final double EPSILON = 1e-9;

	@Rule
	public TemporaryFolder folder = new TemporaryFolder();

	private File standardTable() throws IOException {
		return write("features.pin",
				"id\tLabel\tScanNr\tfeatA\tconstFeat\tsequence\tProteins",
				"psm1\t1\t10\t1.5\t7\tPEPTIDEK\tsp|P1|A",
				"psm2\t-1\t11\t0.5\t7\tDECOYK\tDECOY_sp|P1|A",
				"psm3\t1\t12\t2.5\t7\tOTHERK\tsp|P2|B");
	}

	@Test
	public void identifiesNonFeatureColumnsByPosition() throws IOException {
		EncyclopediaFeatures table = EncyclopediaFeatures.read(standardTable());

		assertEquals(3, table.size());
		assertEquals(Arrays.asList("featA", "constFeat"), table.getFeatureNames());
		assertEquals(2, table.getScanIndex());
		assertEquals(5, table.getPeptideIndex());
		assertEquals(6, table.getProteinIndex());
	}

	@Test
	public void readsPercolatorLabels() throws IOException {
		EncyclopediaFeatures table = EncyclopediaFeatures.read(standardTable());

		assertFalse(table.isDecoy(table.getRows().get(0)));
		assertTrue(table.isDecoy(table.getRows().get(1)));
	}

	@Test
	public void mergesExtraTrailingFieldsIntoTheProteinColumn() throws IOException {
		File file = write("features.pin",
				"id\tLabel\tScanNr\tfeatA\tsequence\tProteins",
				"psm1\t1\t10\t1.5\tPEPTIDEK\tsp|P1|A\tsp|P2|B\tsp|P3|C");

		EncyclopediaFeatures table = EncyclopediaFeatures.read(file);

		assertEquals("sp|P1|A,sp|P2|B,sp|P3|C", table.getRows().get(0)[table.getProteinIndex()]);
	}

	@Test
	public void rejectsRowsThatAreTooShort() throws IOException {
		File file = write("features.pin",
				"id\tLabel\tScanNr\tfeatA\tsequence\tProteins",
				"psm1\t1\t10");

		try {
			EncyclopediaFeatures.read(file);
			fail("expected a failure on a truncated row");
		} catch (IOException expected) {
			assertTrue(expected.getMessage().contains("header declares"));
		}
	}

	@Test
	public void nonNumericFeatureValuesBecomeZero() throws IOException {
		File file = write("features.pin",
				"id\tLabel\tScanNr\tfeatA\tfeatB\tsequence\tProteins",
				"psm1\t1\t10\tNaN\t\tPEPTIDEK\tsp|P1|A");

		EncyclopediaFeatures table = EncyclopediaFeatures.read(file);
		double[] values = table.getFeatureValues(table.getRows().get(0), table.getFeatureIndices());

		assertEquals(0.0, values[0], EPSILON);
		assertEquals(0.0, values[1], EPSILON);
	}

	@Test
	public void findsConstantFeatureColumns() throws IOException {
		EncyclopediaFeatures table = EncyclopediaFeatures.read(standardTable());

		assertEquals(Collections.singletonList("constFeat"), table.getNearConstantFeatureNames());
	}

	@Test
	public void droppedColumnsDisappearFromTheWrittenFile() throws IOException {
		EncyclopediaFeatures table = EncyclopediaFeatures.read(standardTable());
		File pruned = new File(folder.getRoot(), "pruned.pin");

		table.writeWithout(pruned, new HashSet<>(Collections.singletonList("constFeat")));

		EncyclopediaFeatures reread = EncyclopediaFeatures.read(pruned);
		assertEquals(Collections.singletonList("featA"), reread.getFeatureNames());
		assertEquals(3, reread.size());
		assertEquals(1.5, reread.getFeatureValues(reread.getRows().get(0), reread.getFeatureIndices())[0], EPSILON);
	}

	@Test
	public void writesUnixLineEndings() throws IOException {
		EncyclopediaFeatures table = EncyclopediaFeatures.read(standardTable());
		File pruned = new File(folder.getRoot(), "pruned.pin");

		table.writeWithout(pruned, Collections.<String>emptySet());

		String written = new String(Files.readAllBytes(pruned.toPath()), StandardCharsets.UTF_8);
		assertFalse(written.contains("\r"));
	}

	@Test
	public void headersMustMatchBetweenBackgroundAndReference() throws IOException {
		File background = standardTable();
		File reference = write("other.pin", "id\tLabel\tScanNr\tfeatZ\tsequence\tProteins");

		try {
			EncyclopediaFeatures.requireMatchingHeaders(background, reference);
			fail("expected mismatched headers to be rejected");
		} catch (IOException expected) {
			assertTrue(expected.getMessage().contains("same columns"));
		}

		EncyclopediaFeatures.requireMatchingHeaders(background, background);
	}

	private File write(String name, String... lines) throws IOException {
		File file = folder.newFile(name);
		try (PrintWriter out = new PrintWriter(file, "UTF-8")) {
			for (String line : lines) {
				out.println(line);
			}
		}
		return file;
	}
}
