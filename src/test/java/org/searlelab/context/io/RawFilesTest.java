package org.searlelab.context.io;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;
import org.searlelab.context.io.RawFiles;

public class RawFilesTest {

	@Test
	public void stripsEveryFormatMSRawJavaCanOpen() {
		assertEquals("/data/run01", RawFiles.stripExtension("/data/run01.raw"));
		assertEquals("/data/run01", RawFiles.stripExtension("/data/run01.d"));
		assertEquals("/data/run01", RawFiles.stripExtension("/data/run01.dia"));
		assertEquals("/data/run01", RawFiles.stripExtension("/data/run01.mzML"));
	}

	@Test
	public void isCaseInsensitive() {
		assertEquals("/data/run01", RawFiles.stripExtension("/data/run01.RAW"));
		assertEquals("/data/run01", RawFiles.stripExtension("/data/run01.MzMl"));
		assertEquals("/data/run01", RawFiles.stripExtension("/data/run01.DIA"));
	}

	@Test
	public void prefersTheLongestMatchingExtension() {
		assertEquals("/data/sample.mz", RawFiles.stripExtension("/data/sample.mz.d"));
		assertEquals("/data/run01", RawFiles.stripExtension("/data/run01.mzml"));
	}

	@Test
	public void leavesUnrecognisedNamesAlone() {
		assertEquals("/data/run01.wiff", RawFiles.stripExtension("/data/run01.wiff"));
		assertEquals("/data/run01", RawFiles.stripExtension("/data/run01"));
		assertEquals("/my.data/run01.txt", RawFiles.stripExtension("/my.data/run01.txt"));
	}

	@Test
	public void handlesNull() {
		assertEquals(null, RawFiles.stripExtension(null));
		assertFalse(RawFiles.isRecognised(null));
	}

	@Test
	public void recognisesSupportedFormats() {
		assertTrue(RawFiles.isRecognised("run01.raw"));
		assertTrue(RawFiles.isRecognised("run01.d"));
		assertTrue(RawFiles.isRecognised("run01.dia"));
		assertTrue(RawFiles.isRecognised("run01.mzML"));
		assertFalse(RawFiles.isRecognised("run01.wiff"));
		assertFalse(RawFiles.isRecognised("notes.txt"));
	}

	@Test
	public void namesTheFormatsInHelpText() {
		String supported = RawFiles.supportedExtensions();
		assertTrue(supported.contains(".raw"));
		assertTrue(supported.contains(".d"));
		assertTrue(supported.contains(".dia"));
		assertTrue(supported.contains(".mzML"));
	}
}
