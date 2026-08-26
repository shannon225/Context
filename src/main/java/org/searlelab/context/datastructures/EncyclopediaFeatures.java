package org.searlelab.context.datastructures;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

public class EncyclopediaFeatures {
	private static final String DELIM = "\t";

	public static final double NEAR_CONSTANT_SIGMA = 1e-6;

	public static final int ID_INDEX = 0;
	public static final int LABEL_INDEX = 1;

	private final String[] header;
	private final ArrayList<String[]> rows;

	private EncyclopediaFeatures(String[] header, ArrayList<String[]> rows) {
		this.header = header;
		this.rows = rows;
	}

	public static EncyclopediaFeatures read(File file) throws IOException {
		ArrayList<String[]> rows = new ArrayList<>();
		String[] header;

		try (BufferedReader in = new BufferedReader(new FileReader(file))) {
			String headerLine = in.readLine();
			if (headerLine == null) {
				throw new IOException("Feature file is empty: " + file.getAbsolutePath());
			}
			header = headerLine.split(DELIM, -1);
			if (header.length < 5) {
				throw new IOException("Feature file needs at least id, Label, ScanNr, sequence and Proteins columns: "
						+ file.getAbsolutePath());
			}

			String line;
			int lineNumber = 1;
			while ((line = in.readLine()) != null) {
				lineNumber++;
				if (line.trim().isEmpty()) continue;

				String[] columns = line.split(DELIM, -1);
				if (columns.length > header.length) {
					columns = mergeTrailingColumns(columns, header.length);
				}
				if (columns.length != header.length) {
					throw new IOException(file.getName() + " line " + lineNumber + ": found " + columns.length
							+ " fields but the header declares " + header.length);
				}
				rows.add(columns);
			}
		}

		return new EncyclopediaFeatures(header, rows);
	}

	private static String[] mergeTrailingColumns(String[] columns, int width) {
		String[] merged = Arrays.copyOf(columns, width);
		StringBuilder tail = new StringBuilder(columns[width - 1]);
		for (int i = width; i < columns.length; i++) {
			tail.append(',').append(columns[i]);
		}
		merged[width - 1] = tail.toString();
		return merged;
	}

	public static String[] requireMatchingHeaders(File background, File reference) throws IOException {
		String[] backgroundHeader = readHeader(background);
		String[] referenceHeader = readHeader(reference);

		if (!Arrays.equals(backgroundHeader, referenceHeader)) {
			throw new IOException("Background and reference feature files do not declare the same columns."
					+ "\n  background: " + background.getAbsolutePath()
					+ "\n  reference:  " + reference.getAbsolutePath()
					+ "\nThey must come from the same EncyclopeDIA run.");
		}
		return backgroundHeader;
	}

	public static String[] readHeader(File file) throws IOException {
		try (BufferedReader in = new BufferedReader(new FileReader(file))) {
			String headerLine = in.readLine();
			if (headerLine == null) {
				throw new IOException("Feature file is empty: " + file.getAbsolutePath());
			}
			return headerLine.split(DELIM, -1);
		}
	}

	public String[] getHeader() {
		return header;
	}

	public ArrayList<String[]> getRows() {
		return rows;
	}

	public int size() {
		return rows.size();
	}

	public int getPeptideIndex() {
		return header.length - 2;
	}

	public int getProteinIndex() {
		return header.length - 1;
	}

	public int getScanIndex() {
		for (int i = 0; i < header.length; i++) {
			if ("ScanNr".equalsIgnoreCase(header[i]) || "scan".equalsIgnoreCase(header[i])) return i;
		}
		return -1;
	}

	public int[] getFeatureIndices() {
		return getFeatureIndices(header);
	}

	public static int[] getFeatureIndices(String[] header) {
		Set<Integer> skip = new LinkedHashSet<>();
		skip.add(ID_INDEX);
		skip.add(LABEL_INDEX);
		skip.add(header.length - 2);
		skip.add(header.length - 1);
		for (int i = 0; i < header.length; i++) {
			if ("ScanNr".equalsIgnoreCase(header[i]) || "scan".equalsIgnoreCase(header[i])) skip.add(i);
		}

		int[] indices = new int[header.length - skip.size()];
		int n = 0;
		for (int i = 0; i < header.length; i++) {
			if (!skip.contains(i)) indices[n++] = i;
		}
		return indices;
	}

	public List<String> getFeatureNames() {
		ArrayList<String> names = new ArrayList<>();
		for (int index : getFeatureIndices()) {
			names.add(header[index]);
		}
		return names;
	}

	public boolean isDecoy(String[] row) {
		return parseLabel(row) < 0;
	}

	public int parseLabel(String[] row) {
		try {
			return Integer.parseInt(row[LABEL_INDEX].trim());
		} catch (NumberFormatException nfe) {
			return 1;
		}
	}

	public double[] getFeatureValues(String[] row, int[] featureIndices) {
		double[] values = new double[featureIndices.length];
		for (int i = 0; i < featureIndices.length; i++) {
			values[i] = parseDouble(row[featureIndices[i]]);
		}
		return values;
	}

	private static double parseDouble(String value) {
		if (value == null) return 0.0;
		String trimmed = value.trim();
		if (trimmed.isEmpty()) return 0.0;
		try {
			double parsed = Double.parseDouble(trimmed);
			return Double.isFinite(parsed) ? parsed : 0.0;
		} catch (NumberFormatException nfe) {
			return 0.0;
		}
	}

	public List<String> getNearConstantFeatureNames() {
		return getNearConstantFeatureNames(NEAR_CONSTANT_SIGMA);
	}

	public List<String> getNearConstantFeatureNames(double sigmaThreshold) {
		ArrayList<String> constant = new ArrayList<>();
		if (rows.isEmpty()) return constant;

		for (int index : getFeatureIndices()) {
			double sum = 0.0;
			for (String[] row : rows) {
				sum += parseDouble(row[index]);
			}
			double mean = sum / rows.size();

			double sumSquaredDeviation = 0.0;
			for (String[] row : rows) {
				double deviation = parseDouble(row[index]) - mean;
				sumSquaredDeviation += deviation * deviation;
			}
			double sigma = Math.sqrt(sumSquaredDeviation / rows.size());

			if (sigma < sigmaThreshold) constant.add(header[index]);
		}
		return constant;
	}

	public void writeWithout(File outputFile, Set<String> columnsToDrop) throws IOException {
		boolean[] keep = new boolean[header.length];
		int kept = 0;
		for (int i = 0; i < header.length; i++) {
			keep[i] = !columnsToDrop.contains(header[i]);
			if (keep[i]) kept++;
		}
		if (kept < 5) {
			throw new IOException("Dropping " + columnsToDrop + " would leave too few columns to run Percolator.");
		}

		File parent = outputFile.getParentFile();
		if (parent != null && !parent.exists() && !parent.mkdirs()) {
			throw new IOException("Could not create directory " + parent.getAbsolutePath());
		}

		try (BufferedWriter out = new BufferedWriter(new FileWriter(outputFile))) {
			writeKeptFields(out, header, keep);
			for (String[] row : rows) {
				writeKeptFields(out, row, keep);
			}
		}
	}

	private static void writeKeptFields(BufferedWriter out, String[] fields, boolean[] keep) throws IOException {
		boolean first = true;
		for (int i = 0; i < fields.length; i++) {
			if (!keep[i]) continue;
			if (!first) out.write(DELIM);
			out.write(fields[i]);
			first = false;
		}
		out.write("\n");
	}
}
