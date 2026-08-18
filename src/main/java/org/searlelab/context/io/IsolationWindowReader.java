package org.searlelab.context.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import org.searlelab.context.datastructures.IsolationWindow;

public class IsolationWindowReader {

	// formatted as a mass list that is output when generating targeted assays with
	// encyclopedia
	public static ArrayList<IsolationWindow> parseMassList(String massListFile) {

		// Variables to fill in with the assay.csv entries
		ArrayList<IsolationWindow> isolationWindows = new ArrayList<>();
		File massList = new File(massListFile);
		boolean hasPrintedDebugInfo = false;
		boolean hasPrintedAddingPrecursor = false;

		try (BufferedReader br = new BufferedReader(new FileReader(massList))) {
			String delim = getDelimiter(massListFile);

			String header = br.readLine();

			String line;
			while ((line = br.readLine()) != null) {

				String[] columns = line.split(delim, -1);

				if (!hasPrintedDebugInfo) {
					hasPrintedDebugInfo = true;
				}

				System.out.println(line); // Console will print what the data looks like as its read in

				String peptide = columns[0];
				double targetMz = Double.parseDouble(columns[3]);
				byte charge = Byte.parseByte(columns[4]);
				float rtCenter = Float.parseFloat(columns[5]);
				float rtWindow = Float.parseFloat(columns[6]);

				float rtMin = (rtCenter - (rtWindow / 2)) * 60;
				float rtMax = (rtCenter + (rtWindow / 2)) * 60;

				boolean isDecoy = Boolean.parseBoolean(columns[7]);

				// Assemble each window
				IsolationWindow window = new IsolationWindow(peptide, targetMz, charge, rtMin, rtMax, isDecoy);
				isolationWindows.add(window);

				if (!hasPrintedAddingPrecursor) {

					hasPrintedAddingPrecursor = true;
				}

			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		return isolationWindows;

	}

	// Detect the delim - Thermo Mass Lists are usually in .csv, but this program
	// accepts .csv or .txt
	private static String getDelimiter(String filePath) {
		String lowerPath = filePath.toLowerCase();

		if (lowerPath.endsWith(".csv")) {
			return ",";
		}

		if (lowerPath.endsWith(".txt")) {
			return "\t";
		}

		if (lowerPath.endsWith(".tsv")) {
			return "\t";
		}

		throw new IllegalArgumentException(
				"Error in reading Isolation Windows. Mass list file must be a .csv or .tsv file: " + filePath);
	}
}
