package org.searlelab.context.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import org.searlelab.context.mprophet.IsolationWindow;

public class IsolationWindowReader {

	public static ArrayList<IsolationWindow> parseMassList(String massListFile) {

		// Variables to fill in with the assay.csv entries
		ArrayList<IsolationWindow> isolationWindows = new ArrayList<>();
		File massList = new File(massListFile);

		try (BufferedReader br = new BufferedReader(new FileReader(massList))) {

			String header = br.readLine();
			if (header == null) {
				throw new IOException("Mass list is empty: " + massListFile);
			}

			String delim = getDelimiter(massListFile, header);

			String line;
			int lineNumber = 1;
			while ((line = br.readLine()) != null) {
				lineNumber++;
				if (line.trim().isEmpty()) continue;

				String[] columns = line.split(delim, -1);

				// a row without an RT window cannot describe an isolation window
				if (columns.length < 7) {
					System.err.println("Skipping mass-list line " + lineNumber + " (expected at least 7 fields, got "
							+ columns.length + "): " + line);
					continue;
				}

				String peptide = columns[0];
				double targetMz = Double.parseDouble(columns[3]);
				byte charge = Byte.parseByte(columns[4]);
				float rtCenter = Float.parseFloat(columns[5]);
				float rtWindow = Float.parseFloat(columns[6]);

				float rtMin = (rtCenter - (rtWindow / 2)) * 60;
				float rtMax = (rtCenter + (rtWindow / 2)) * 60;

				// Skyline exports stop at the RT window; only our own assays carry the flag
				boolean isDecoy = columns.length > 7 && Boolean.parseBoolean(columns[7]);

				// Assemble each window
				isolationWindows.add(new IsolationWindow(peptide, targetMz, charge, rtMin, rtMax, isDecoy));
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		return isolationWindows;

	}

	private static String getDelimiter(String filePath, String header) {
		int tabs = countChar(header, '\t');
		int commas = countChar(header, ',');
		if (tabs > 0 || commas > 0) {
			return tabs >= commas ? "\t" : ",";
		}

		String lowerPath = filePath.toLowerCase();
		if (lowerPath.endsWith(".csv")) return ",";
		if (lowerPath.endsWith(".txt") || lowerPath.endsWith(".tsv")) return "\t";

		throw new IllegalArgumentException(
				"Error in reading Isolation Windows. Mass list file must be a .csv or .tsv file: " + filePath);
	}

	private static int countChar(String s, char c) {
		int n = 0;
		for (int i = 0; i < s.length(); i++) {
			if (s.charAt(i) == c) n++;
		}
		return n;
	}
}
