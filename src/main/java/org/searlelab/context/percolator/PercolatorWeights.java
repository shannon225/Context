package org.searlelab.context.percolator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class PercolatorWeights {
	private static final String DELIM = "\t";
	private static final String BIAS_ROW_NAME = "__bias__";

	private final List<String> featureNames;
	private final double[] weights;
	private final double bias;
	private final int numberOfBins;

	PercolatorWeights(List<String> featureNames, double[] weights, double bias, int numberOfBins) {
		this.featureNames = Collections.unmodifiableList(new ArrayList<>(featureNames));
		this.weights = weights;
		this.bias = bias;
		this.numberOfBins = numberOfBins;
	}

	public static PercolatorWeights parse(File weightsFile) throws IOException {
		List<String> lines = readMeaningfulLines(weightsFile);
		if (lines.isEmpty()) {
			throw new IOException("Percolator wrote no weights to " + weightsFile.getAbsolutePath()
					+ ". Check the Percolator log next to it.");
		}

		List<String> featureNames = null;
		ArrayList<double[]> rawWeightsPerBin = new ArrayList<>();

		List<String> currentHeader = null;
		double[] currentRaw = null;

		for (String line : lines) {
			String[] fields = line.split(DELIM, -1);

			if (isNumeric(fields[0])) {
				if (currentHeader == null) {
					throw new IOException(weightsFile.getName()
							+ ": found a row of weights before any feature-name header.");
				}
				if (fields.length != currentHeader.size() + 1) {
					throw new IOException(weightsFile.getName() + ": a weights row has " + fields.length
							+ " values but its header names " + currentHeader.size()
							+ " features (plus a bias column).");
				}
				currentRaw = parseRow(fields, weightsFile);
			} else {
				if (currentHeader != null) {
					rawWeightsPerBin.add(requireRow(currentRaw, weightsFile));
				}
				currentHeader = namesFromHeader(fields);
				currentRaw = null;

				if (featureNames == null) {
					featureNames = currentHeader;
				} else if (!featureNames.equals(currentHeader)) {
					throw new IOException(weightsFile.getName()
							+ ": the cross-validation bins name different features.");
				}
			}
		}
		if (currentHeader != null) {
			rawWeightsPerBin.add(requireRow(currentRaw, weightsFile));
		}

		if (featureNames == null || rawWeightsPerBin.isEmpty()) {
			throw new IOException("Could not parse any weights from " + weightsFile.getAbsolutePath());
		}

		int numberOfFeatures = featureNames.size();
		double[] averagedWeights = new double[numberOfFeatures];
		double averagedBias = 0.0;
		for (double[] row : rawWeightsPerBin) {
			for (int i = 0; i < numberOfFeatures; i++) {
				averagedWeights[i] += row[i];
			}
			averagedBias += row[numberOfFeatures];
		}
		for (int i = 0; i < numberOfFeatures; i++) {
			averagedWeights[i] /= rawWeightsPerBin.size();
		}
		averagedBias /= rawWeightsPerBin.size();

		return new PercolatorWeights(featureNames, averagedWeights, averagedBias, rawWeightsPerBin.size());
	}

	private static List<String> readMeaningfulLines(File file) throws IOException {
		ArrayList<String> lines = new ArrayList<>();
		try (BufferedReader in = new BufferedReader(new FileReader(file))) {
			String line;
			while ((line = in.readLine()) != null) {
				String trimmed = line.trim();
				if (trimmed.isEmpty() || trimmed.startsWith("#")) continue;
				lines.add(line);
			}
		}
		return lines;
	}

	private static List<String> namesFromHeader(String[] fields) {
		ArrayList<String> names = new ArrayList<>(fields.length - 1);
		for (int i = 0; i < fields.length - 1; i++) {
			names.add(fields[i].trim());
		}
		return names;
	}

	private static double[] parseRow(String[] fields, File weightsFile) throws IOException {
		double[] values = new double[fields.length];
		for (int i = 0; i < fields.length; i++) {
			try {
				values[i] = Double.parseDouble(fields[i].trim());
			} catch (NumberFormatException nfe) {
				throw new IOException(weightsFile.getName() + ": could not parse weight [" + fields[i] + "]", nfe);
			}
		}
		return values;
	}

	private static double[] requireRow(double[] row, File weightsFile) throws IOException {
		if (row == null) {
			throw new IOException(weightsFile.getName()
					+ ": a cross-validation bin has a feature-name header but no weights.");
		}
		return row;
	}

	public List<String> getFeatureNames() {
		return featureNames;
	}

	public double[] getWeights() {
		return weights.clone();
	}

	public double getBias() {
		return bias;
	}

	/** How many cross-validation bins were averaged together. */
	public int getNumberOfBins() {
		return numberOfBins;
	}

	public double score(double[] featureValues) {
		if (featureValues.length != weights.length) {
			throw new IllegalArgumentException("Model has " + weights.length + " weights but was handed "
					+ featureValues.length + " feature values.");
		}
		double score = bias;
		for (int i = 0; i < weights.length; i++) {
			score += weights[i] * featureValues[i];
		}
		return score;
	}

	public void writeAveraged(File outputFile) throws IOException {
		File parent = outputFile.getParentFile();
		if (parent != null && !parent.exists() && !parent.mkdirs()) {
			throw new IOException("Could not create directory " + parent.getAbsolutePath());
		}

		try (BufferedWriter out = new BufferedWriter(new FileWriter(outputFile))) {
			out.write("feature" + DELIM + "weight");
			out.write("\n");
			for (int i = 0; i < weights.length; i++) {
				out.write(featureNames.get(i) + DELIM + weights[i]);
				out.write("\n");
			}
			out.write(BIAS_ROW_NAME + DELIM + bias);
			out.write("\n");
		}
	}

	private static boolean isNumeric(String value) {
		if (value == null) return false;
		String trimmed = value.trim();
		if (trimmed.isEmpty()) return false;
		try {
			Double.parseDouble(trimmed);
			return true;
		} catch (NumberFormatException nfe) {
			return false;
		}
	}

	@Override
	public String toString() {
		return "PercolatorWeights[" + weights.length + " features, " + numberOfBins + " CV bins, bias=" + bias + "]";
	}
}
