package org.searlelab.context.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import org.searlelab.context.mprophet.IsolationWindow;

public class IsolationWindowReader {


	// formatted as a mass list that is output when generating targeted assays with
	// encyclopedia
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

		boolean hasPrintedDebugInfo = false;
		boolean hasPrintedAddingPrecursor = false;

		try (BufferedReader br = new BufferedReader(new FileReader(massList))) {
			String delim = getDelimiter(massListFile);

			@SuppressWarnings("unused")
			String header = br.readLine();

			String line;
			while ((line = br.readLine()) != null) {

				String[] columns = line.split(delim, -1);

				if (!hasPrintedDebugInfo) {
					// System.out.println("Line being read: " + line);
					// System.out.println("Number of columns: " + columns.length);
					hasPrintedDebugInfo = true;
				}

				// String columns[] = line.split(DELIM, -1);
				System.out.println(line); // Console will print what the data looks like as its read in

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

				if (!hasPrintedAddingPrecursor) {
					// System.out.println("Adding " + peptide + " with precursor at " + targetMz + "
					// m/z, charge " + charge + ", RT " + rtCenter + " min " + rtMin/60 + " max " +
					// rtMax/60 + " isDecoy = " + isDecoy);
					hasPrintedAddingPrecursor = true;
				}

			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		return isolationWindows;

	}

	private static String getDelimiterByChar(String filePath, String header) {
		int tabs = countChar(header, '\t');
		int commas = countChar(header, ',');
		if (tabs > 0 || commas > 0) {
			return tabs >= commas ? "\t" : ",";
		}

		String lowerPath = filePath.toLowerCase();
		if (lowerPath.endsWith(".csv")) return ",";
		if (lowerPath.endsWith(".txt") || lowerPath.endsWith(".tsv")) return "\t";

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

	private static int countChar(String s, char c) {
		int n = 0;
		for (int i = 0; i < s.length(); i++) {
			if (s.charAt(i) == c) n++;
		}
		return n;
	}
}
