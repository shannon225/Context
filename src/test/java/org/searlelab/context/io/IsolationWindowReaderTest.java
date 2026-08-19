package org.searlelab.context.io;

import static org.junit.Assert.*;

import java.net.URL;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Objects;

import org.junit.Test;
import org.searlelab.context.datastructures.IsolationWindow;

public class IsolationWindowReaderTest {

	@Test
	public void testThatFileWorks() throws Throwable {

		// Does this correctly read a mass list? Read the mass list, and determine if
		// the first two precursor values and retention times are the same.
		// For the file used in the test case, this is true because its a target-decoy
		// pair/
		URL massListFile = Objects.requireNonNull(getClass().getResource("/org/searlelab/context/io/IL2A_GPFDIA_0combined_masked0_assay.txt"));

		String massListPath = Paths.get(massListFile.toURI()).toString();
		System.out.println("Mass list test file is " + massListPath);
		ArrayList<IsolationWindow> massList = IsolationWindowReader.parseMassList(massListPath);

		IsolationWindow firstWindow = massList.get(0);
		IsolationWindow secondWindow = massList.get(1);

		System.out.println("First precursor m/z: " + firstWindow.getTargetMz());
		System.out.println("Second precursor m/z: " + secondWindow.getTargetMz());

		assertEquals(firstWindow.getTargetMz(), secondWindow.getTargetMz(), 0.000001);
	}


}
