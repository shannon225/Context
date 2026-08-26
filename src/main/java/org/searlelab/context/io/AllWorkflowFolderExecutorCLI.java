package org.searlelab.context.io;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Stream;

import org.searlelab.context.percolator.ContextPercolator;
import org.searlelab.context.percolator.PyIsoPEPRunner;

import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

/**
 * Finds matched complete, background, and reference feature files in a folder
 * and runs AllWorkflowExecutorCLI once for each matched set.
 */
public final class AllWorkflowFolderExecutorCLI {

	private static final String COMPLETE_SUFFIX = ".features.txt";
	private static final String BACKGROUND_SUFFIX = "_background.features.txt";
	private static final String REFERENCE_SUFFIX = "_reference.features.txt";

	private AllWorkflowFolderExecutorCLI() {
	}

	public static void main(String[] args) {
		HashMap<String, String> arguments = CommandLineParser.parseArguments(args);

		if (args.length == 0 || containsHelpFlag(arguments)) {
			printHelp();
			System.exit(args.length == 0 ? 1 : 0);
		}

		try {
			File featureFolder = requiredReadableDirectory(arguments, "--features-folder");
			File fasta = requiredReadableFile(arguments, "-f");

			float fdr = Float.parseFloat(arguments.getOrDefault("-fdr", Float.toString(ContextPercolator.DEFAULT_FDR)));

			if (!(fdr > 0.0f && fdr <= 1.0f)) {
				throw new IllegalArgumentException("FDR must be greater than 0 and at most 1.");
			}

			File outputRoot = DirectoryOptions.outputDirectory(arguments, featureFolder.getAbsoluteFile());

			PyIsoPEPRunner pyIsoPEP = new PyIsoPEPRunner(arguments.get("--pyisopep"));

			HashMap<String, String> encyclopediaArguments = encyclopediaArguments(arguments);

			List<FeatureFileSet> featureSets = findFeatureFileSets(featureFolder);

			Logger.logLine("Found " + featureSets.size() + " matched feature-file sets.");

			for (int i = 0; i < featureSets.size(); i++) {
				FeatureFileSet set = featureSets.get(i);

				Logger.logLine("");
				Logger.logLine("[" + (i + 1) + "/" + featureSets.size() + "] Processing " + set.prefix);

				File setOutputDirectory = new File(outputRoot, set.prefix);

				AllWorkflowExecutorCLI.runAll(set.completeFeatures, set.backgroundFeatures, set.referenceFeatures,
						fasta, pyIsoPEP, fdr, setOutputDirectory, set.prefix, new HashMap<>(encyclopediaArguments));
			}

			Logger.logLine("");
			Logger.logLine("Finished " + featureSets.size() + " matched feature-file sets.");

		} catch (Exception e) {
			Logger.errorLine("Folder workflow execution failed: " + e.getMessage());
			e.printStackTrace();
			System.exit(2);
		}
	}

	/**
	 * Finds complete feature files and validates that each has corresponding
	 * background and reference files.
	 */
	private static List<FeatureFileSet> findFeatureFileSets(File folder) throws IOException {

		List<Path> completeFiles = new ArrayList<>();

		try (Stream<Path> paths = Files.list(folder.toPath())) {
			paths.filter(Files::isRegularFile).filter(AllWorkflowFolderExecutorCLI::isCompleteFeatureFile)
					.sorted(Comparator.comparing(path -> path.getFileName().toString())).forEach(completeFiles::add);
		}

		if (completeFiles.isEmpty()) {
			throw new IOException("No complete *.features.txt files were found in " + folder.getAbsolutePath());
		}

		List<FeatureFileSet> featureSets = new ArrayList<>();
		List<String> missingFiles = new ArrayList<>();

		for (Path completePath : completeFiles) {
			File completeFile = completePath.toFile();
			String completeName = completeFile.getName();

			String prefix = completeName.substring(0, completeName.length() - COMPLETE_SUFFIX.length());

			File backgroundFile = new File(folder, prefix + BACKGROUND_SUFFIX);

			File referenceFile = new File(folder, prefix + REFERENCE_SUFFIX);

			if (!backgroundFile.isFile() || !backgroundFile.canRead()) {
				missingFiles.add(backgroundFile.getName());
			}

			if (!referenceFile.isFile() || !referenceFile.canRead()) {
				missingFiles.add(referenceFile.getName());
			}

			if (backgroundFile.isFile() && backgroundFile.canRead() && referenceFile.isFile()
					&& referenceFile.canRead()) {

				featureSets.add(new FeatureFileSet(prefix, completeFile, backgroundFile, referenceFile));
			}
		}

		/*
		 * Validate every set before running anything. This prevents the program from
		 * completing some seeds and then failing halfway through because a later set is
		 * incomplete.
		 */
		if (!missingFiles.isEmpty()) {
			throw new IOException("Missing matched feature files: " + String.join(", ", missingFiles));
		}

		return featureSets;
	}

	/**
	 * A complete feature file ends in .features.txt but is not one of the
	 * background or reference partitions.
	 */
	private static boolean isCompleteFeatureFile(Path path) {
		String name = path.getFileName().toString();

		return name.endsWith(COMPLETE_SUFFIX) && !name.endsWith(BACKGROUND_SUFFIX) && !name.endsWith(REFERENCE_SUFFIX);
	}

	private static HashMap<String, String> encyclopediaArguments(HashMap<String, String> arguments) {

		HashMap<String, String> parameters = SearchParameterParser.getDefaultParameters();

		HashMap<String, String> remaining = new HashMap<>(arguments);

		for (String workflowFlag : new String[] { "--features-folder", "--features", "--background", "--reference", "-f",
				"-fdr", "-o", "--outdir", "--prefix", "--pyisopep", "-h", "-help", "--help" }) {
			remaining.remove(workflowFlag);
		}

		parameters.putAll(remaining);
		return parameters;
	}

	private static File requiredReadableDirectory(HashMap<String, String> arguments, String flag) throws IOException {

		String value = arguments.get(flag);

		if (value == null || value.trim().isEmpty()) {
			throw new IllegalArgumentException("Missing required argument: " + flag);
		}

		File directory = new File(value);

		if (!directory.isDirectory() || !directory.canRead()) {
			throw new IOException(flag + " is not a readable directory: " + directory.getAbsolutePath());
		}

		return directory;
	}

	private static File requiredReadableFile(HashMap<String, String> arguments, String flag) throws IOException {

		String value = arguments.get(flag);

		if (value == null || value.trim().isEmpty()) {
			throw new IllegalArgumentException("Missing required argument: " + flag);
		}

		File file = new File(value);

		if (!file.isFile() || !file.canRead()) {
			throw new IOException(flag + " is not a readable file: " + file.getAbsolutePath());
		}

		return file;
	}

	private static boolean containsHelpFlag(HashMap<String, String> arguments) {

		return arguments.containsKey("-h") || arguments.containsKey("-help") || arguments.containsKey("--help");
	}

	private static void printHelp() {
		Logger.timelessLogLine("AllWorkflowFolderExecutorCLI");
		Logger.timelessLogLine("Run all four workflows on every matched feature-file set in a folder.");
		Logger.timelessLogLine("");
		Logger.timelessLogLine("Required:");
		Logger.timelessLogLine("  --features-folder <dir>  folder containing matched feature files");
		Logger.timelessLogLine("  -f               <file> FASTA protein database");
		Logger.timelessLogLine("");
		Logger.timelessLogLine("Optional:");
		Logger.timelessLogLine("  -fdr      <float> peptide FDR threshold");
		Logger.timelessLogLine("  -o        <dir>   output root; defaults to the feature folder");
		Logger.timelessLogLine("  --pyisopep <path>  pyIsoPEP executable; defaults to PATH");
		Logger.timelessLogLine("");
		Logger.timelessLogLine("Expected naming:");
		Logger.timelessLogLine("  X.features.txt");
		Logger.timelessLogLine("  X_background.features.txt");
		Logger.timelessLogLine("  X_reference.features.txt");
	}

	private static final class FeatureFileSet {

		private final String prefix;
		private final File completeFeatures;
		private final File backgroundFeatures;
		private final File referenceFeatures;

		private FeatureFileSet(String prefix, File completeFeatures, File backgroundFeatures, File referenceFeatures) {

			this.prefix = prefix;
			this.completeFeatures = completeFeatures;
			this.backgroundFeatures = backgroundFeatures;
			this.referenceFeatures = referenceFeatures;
		}
	}
}