package org.searlelab.context.io;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.HashMap;

import org.searlelab.context.encyclopedia.MProphetReiter;
import org.searlelab.context.percolator.ContextPercolator;
import org.searlelab.context.percolator.ContextPercolatorExecutor;
import org.searlelab.context.percolator.ContextPercolatorResult;
import org.searlelab.context.percolator.PyIsoPEPRunner;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearDiscriminantAnalysis;

/**
 * Runs four confidence-estimation workflows on one matched feature-file set:
 *
 * <ol>
 *   <li>Context Percolator on the background and reference feature files</li>
 *   <li>standard Percolator on the complete feature file</li>
 *   <li>Context mProphet on the background and reference feature files</li>
 *   <li>standard mProphet on the complete feature file</li>
 * </ol>
 *
 * The complete, background, and reference files must have been produced from
 * the same acquisition. The background and reference files must partition the
 * peptides represented by the complete feature file.
 */
public final class AllWorkflowExecutorCLI {

	private static final String CONTEXT_MPROPHET_ENGINE_NAME = "context-mprophet";
	private static final String STANDARD_MPROPHET_ENGINE_NAME = "standard-mprophet";
	private static final float DEFAULT_FDR = ContextPercolator.DEFAULT_FDR;
	private static final int MPROPHET_SEED = 1;
	private static final int MPROPHET_ROUND = 1;

	private AllWorkflowExecutorCLI() {
	}

	public static void main(String[] args) {
		HashMap<String, String> arguments = CommandLineParser.parseArguments(args);

		if (args.length == 0 || containsHelpFlag(arguments)) {
			printHelp();
			System.exit(args.length == 0 ? 1 : 0);
		}

		try {
			File allFeatures = requiredReadableFile(arguments, "-features");
			File backgroundFeatures = requiredReadableFile(arguments, "-background");
			File referenceFeatures = requiredReadableFile(arguments, "-reference");
			File fasta = requiredReadableFile(arguments, "-f");

			float fdr = Float.parseFloat(
					arguments.getOrDefault("-fdr", Float.toString(DEFAULT_FDR)));
			File outputDirectory = DirectoryOptions.outputDirectory(
					arguments, allFeatures.getAbsoluteFile().getParentFile());
			String prefix = arguments.getOrDefault("-prefix", stripFeatureSuffix(allFeatures));
			PyIsoPEPRunner pyIsoPEP = new PyIsoPEPRunner(arguments.get("-pyisopep"));
			HashMap<String, String> encyclopediaArguments = encyclopediaArguments(arguments);

			runAll(allFeatures, backgroundFeatures, referenceFeatures, fasta, pyIsoPEP,
					fdr, outputDirectory, prefix, encyclopediaArguments);
		} catch (Exception e) {
			Logger.errorLine("Workflow execution failed: " + e.getMessage());
			e.printStackTrace();
			System.exit(2);
		}
	}

	/**
	 * Programmatic entry point used by the CLI and suitable for integration tests.
	 */
	public static AllWorkflowResult runAll(File allFeatures, File backgroundFeatures,
			File referenceFeatures, File fasta, PyIsoPEPRunner pyIsoPEP, float fdr,
			File outputDirectory, String prefix,
			HashMap<String, String> encyclopediaArguments) throws Exception {

		requireReadableFile(allFeatures, "Complete feature file");
		requireReadableFile(backgroundFeatures, "Background feature file");
		requireReadableFile(referenceFeatures, "Reference feature file");
		requireReadableFile(fasta, "FASTA file");

		if (!(fdr > 0.0f && fdr <= 1.0f)) {
			throw new IllegalArgumentException("FDR must be greater than 0 and at most 1.");
		}
		if (outputDirectory == null) {
			throw new IllegalArgumentException("Output directory must not be null.");
		}
		if (prefix == null || prefix.trim().isEmpty()) {
			throw new IllegalArgumentException("Output prefix must not be empty.");
		}

		HashMap<String, String> parameters = encyclopediaArguments == null
				? SearchParameterParser.getDefaultParameters()
						: new HashMap<>(encyclopediaArguments);

		Logger.logLine("[1/4] Running Context Percolator");
		ContextPercolatorResult contextPercolator = ContextPercolator.trainAndApply(
				backgroundFeatures, referenceFeatures, fasta, parameters, pyIsoPEP,
				fdr, outputDirectory, prefix);

		Logger.logLine("[2/4] Running standard Percolator");
		PercolatorExecutionData standardPercolator =
				ContextPercolatorExecutor.runStandardPercolator(
						allFeatures, fasta, pyIsoPEP, fdr, outputDirectory, prefix,
						new HashMap<>(parameters));

		SearchParameters searchParameters = parseSearchParameters(parameters, fdr);

		Logger.logLine("[3/4] Running Context mProphet");
		MProphetResult contextMProphet = runContextMProphet(
				backgroundFeatures, referenceFeatures, fasta, searchParameters,
				fdr, outputDirectory, prefix);

		Logger.logLine("[4/4] Running standard mProphet");
		MProphetResult standardMProphet = runStandardMProphet(
				allFeatures, fasta, searchParameters, fdr, outputDirectory, prefix);

		Logger.logLine("Finished all four workflows. Results are under "
				+ outputDirectory.getAbsolutePath());

		return new AllWorkflowResult(
				contextPercolator, standardPercolator, contextMProphet, standardMProphet);
	}

	private static MProphetResult runContextMProphet(File backgroundFeatures,
			File referenceFeatures, File fasta, SearchParameters parameters, float fdr,
			File outputDirectory, String prefix) throws Exception {

		File engineDirectory = DirectoryOptions.engineDirectory(
				outputDirectory, CONTEXT_MPROPHET_ENGINE_NAME);

		MProphetExecutionData backgroundData = new MProphetExecutionData(
				backgroundFeatures,
				fasta,
				new File(engineDirectory, prefix + ".peptide.background.target.txt"),
				new File(engineDirectory, prefix + ".peptide.background.decoy.txt"),
				parameters);

		MProphetExecutionData referenceData = new MProphetExecutionData(
				referenceFeatures,
				fasta,
				new File(engineDirectory, prefix + ".peptide.reference.target.txt"),
				new File(engineDirectory, prefix + ".peptide.reference.decoy.txt"),
				parameters);

		deleteMProphetOutputs(backgroundData);
		deleteMProphetOutputs(referenceData);

		MProphetResult backgroundResult = MProphetReiter.executeMProphetTSV(
				backgroundData, fdr, MPROPHET_SEED, parameters.getAAConstants(),
				MPROPHET_ROUND);
		LinearDiscriminantAnalysis backgroundModel = backgroundResult.getLDA();

		MProphetResult referenceResult = MProphetReiter.executeMProphetTSVWithModel(
				referenceData, fdr, backgroundModel, parameters.getAAConstants());

		Logger.logLine("Context mProphet found "
				+ referenceResult.getPassingPeptides().size()
				+ " reference peptides at " + (fdr * 100.0f) + "% FDR");
		return referenceResult;
	}

	private static MProphetResult runStandardMProphet(File allFeatures, File fasta,
			SearchParameters parameters, float fdr, File outputDirectory, String prefix)
					throws Exception {

		File engineDirectory = DirectoryOptions.engineDirectory(
				outputDirectory, STANDARD_MPROPHET_ENGINE_NAME);
		MProphetExecutionData data = new MProphetExecutionData(
				allFeatures,
				fasta,
				new File(engineDirectory, prefix + ".peptide.target.txt"),
				new File(engineDirectory, prefix + ".peptide.decoy.txt"),
				parameters);

		deleteMProphetOutputs(data);
		MProphetResult result = MProphetReiter.executeMProphetTSV(
				data, fdr, MPROPHET_SEED, parameters.getAAConstants(), MPROPHET_ROUND);

		Logger.logLine("Standard mProphet found " + result.getPassingPeptides().size()
				+ " peptides at " + (fdr * 100.0f) + "% FDR");
		return result;
	}

	private static SearchParameters parseSearchParameters(
			HashMap<String, String> encyclopediaArguments, float fdr) {
		HashMap<String, String> parameterMap =
				SearchParameterParser.getDefaultParametersObject().toParameterMap();
		parameterMap.putAll(encyclopediaArguments);
		parameterMap.put("-percolatorThreshold", Float.toString(fdr));
		return SearchParameterParser.parseParameters(parameterMap);
	}

	private static void deleteMProphetOutputs(MProphetExecutionData data) throws IOException {
		Files.deleteIfExists(data.getPeptideOutputFile().toPath());
		Files.deleteIfExists(data.getPeptideDecoyFile().toPath());
	}

	private static HashMap<String, String> encyclopediaArguments(
			HashMap<String, String> arguments) {
		HashMap<String, String> parameters = SearchParameterParser.getDefaultParameters();
		HashMap<String, String> remaining = new HashMap<>(arguments);

		for (String workflowFlag : new String[] {
				"-features", "-background", "-reference", "-f", "-fdr", "-o",
				"-outdir", "-prefix", "-pyisopep", "-h", "-help", "--help" }) {
			remaining.remove(workflowFlag);
		}

		parameters.putAll(remaining);
		return parameters;
	}

	private static File requiredReadableFile(HashMap<String, String> arguments,
			String flag) throws IOException {
		String value = arguments.get(flag);
		if (value == null || value.trim().isEmpty()) {
			throw new IllegalArgumentException("Missing required argument: " + flag);
		}

		File file = new File(value);
		requireReadableFile(file, flag);
		return file;
	}

	private static void requireReadableFile(File file, String description)
			throws IOException {
		if (file == null || !file.isFile() || !file.canRead()) {
			throw new IOException(description + " is not a readable file: "
					+ (file == null ? "null" : file.getAbsolutePath()));
		}
	}

	private static boolean containsHelpFlag(HashMap<String, String> arguments) {
		return arguments.containsKey("-h") || arguments.containsKey("-help")
				|| arguments.containsKey("--help");
	}

	private static String stripFeatureSuffix(File featureFile) {
		return featureFile.getName().replaceFirst("\\.features\\.txt$", "");
	}

	private static void printHelp() {
		Logger.timelessLogLine("AllWorkflowExecutorCLI");
		Logger.timelessLogLine("Run Context Percolator, standard Percolator, Context mProphet,");
		Logger.timelessLogLine("and standard mProphet on one matched feature-file set.");
		Logger.timelessLogLine("");
		Logger.timelessLogLine("Required:");
		Logger.timelessLogLine("  -features   <file> complete, unsplit feature TSV");
		Logger.timelessLogLine("  -background <file> background feature TSV from the same experiment");
		Logger.timelessLogLine("  -reference  <file> reference feature TSV from the same experiment");
		Logger.timelessLogLine("  -f          <file> FASTA protein database");
		Logger.timelessLogLine("");
		Logger.timelessLogLine("Optional:");
		Logger.timelessLogLine("  -fdr      <float> peptide FDR threshold (default: " + DEFAULT_FDR + ")");
		Logger.timelessLogLine("  -o        <dir>   output root (default: next to -features)");
		Logger.timelessLogLine("  -prefix   <name>  output filename prefix");
		Logger.timelessLogLine("  -pyisopep <path>  pyIsoPEP executable (default: PATH)");
		Logger.timelessLogLine("");
		Logger.timelessLogLine("Other flags are passed through to EncyclopeDIA parameters.");
		Logger.timelessLogLine("");
		Logger.timelessLogLine("Output directories:");
		Logger.timelessLogLine("  <o>/" + ContextPercolator.ENGINE_NAME);
		Logger.timelessLogLine("  <o>/standard-percolator");
		Logger.timelessLogLine("  <o>/" + CONTEXT_MPROPHET_ENGINE_NAME);
		Logger.timelessLogLine("  <o>/" + STANDARD_MPROPHET_ENGINE_NAME);
	}

	public static final class AllWorkflowResult {
		private final ContextPercolatorResult contextPercolator;
		private final PercolatorExecutionData standardPercolator;
		private final MProphetResult contextMProphet;
		private final MProphetResult standardMProphet;

		private AllWorkflowResult(ContextPercolatorResult contextPercolator,
				PercolatorExecutionData standardPercolator,
				MProphetResult contextMProphet, MProphetResult standardMProphet) {
			this.contextPercolator = contextPercolator;
			this.standardPercolator = standardPercolator;
			this.contextMProphet = contextMProphet;
			this.standardMProphet = standardMProphet;
		}

		public ContextPercolatorResult getContextPercolator() {
			return contextPercolator;
		}

		public PercolatorExecutionData getStandardPercolator() {
			return standardPercolator;
		}

		public MProphetResult getContextMProphet() {
			return contextMProphet;
		}

		public MProphetResult getStandardMProphet() {
			return standardMProphet;
		}
	}
}
