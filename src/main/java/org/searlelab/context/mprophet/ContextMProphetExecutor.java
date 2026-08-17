package org.searlelab.context.mprophet;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;

import org.searlelab.context.encyclopedia.MProphetReiter;
import org.searlelab.context.io.ContextOptions;
import org.searlelab.context.io.MassListDecoyGenerator;
import org.searlelab.context.io.PreparedFeatures;
import org.searlelab.context.io.RawFiles;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearDiscriminantAnalysis;

public class ContextMProphetExecutor {

	public static final String ENGINE_NAME = "mprophet";

	private static final float DEFAULT_FDR = 0.01f;
	private static final int DEFAULT_SEED = 1;

	private static final int FIRST_ROUND = 1;

	public static void main(String[] args) {
		HashMap<String, String> arguments = CommandLineParser.parseArguments(args);

		if (args.length == 0 || arguments.containsKey("-h") || arguments.containsKey("-help")
				|| arguments.containsKey("--help")) {
			printHelp();
			System.exit(args.length == 0 ? 1 : 0);
		}

		try {
			float fdr = Float.parseFloat(arguments.getOrDefault("-fdr", Float.toString(DEFAULT_FDR)));
			int seed = Integer.parseInt(arguments.getOrDefault("-seed", Integer.toString(DEFAULT_SEED)));

			if (arguments.containsKey("-background") || arguments.containsKey("-reference")) {
				File background = requiredFile(arguments, "-background");
				File reference = requiredFile(arguments, "-reference");
				File fasta = requiredFile(arguments, "-f");
				File outputDirectory = ContextOptions.outputDirectory(arguments,
						background.getAbsoluteFile().getParentFile());
				String prefix = arguments.getOrDefault("-prefix", ContextOptions.stripSplitSuffix(background));

				runTrainApply(reference, background, fasta, fdr, seed, outputDirectory, prefix);
			} else {
				File dia = requiredFile(arguments, "-i");
				File library = requiredFile(arguments, "-l");
				File fasta = requiredFile(arguments, "-f");
				File massList = requiredFile(arguments, "-massList");
				File outputDirectory = ContextOptions.outputDirectory(arguments,
						dia.getAbsoluteFile().getParentFile());
				String prefix = arguments.get("-prefix");
				boolean generateDecoys = ContextOptions.isEnabled(arguments, "-generateDecoys");

				runEndToEnd(library, fasta, dia, massList, fdr, seed, outputDirectory, prefix, generateDecoys);
			}
		} catch (Exception e) {
			Logger.errorLine("ContextMProphetExecutor failed: " + e.getMessage());
			e.printStackTrace();
			System.exit(2);
		}
	}

	public static MProphetResult runEndToEnd(File library, File fasta, File dia, File massList, float fdr, int seed,
			File outputDirectory, String prefix, boolean generateDecoys) throws Exception {

		if (!dia.exists()) throw new IOException("Input file not found: " + dia);
		if (!library.exists()) throw new IOException("Library file not found: " + library);
		if (!fasta.exists()) throw new IOException("FASTA file not found: " + fasta);
		if (!massList.exists()) throw new IOException("Mass list file not found: " + massList);
		RawFiles.requireSupported(dia);

		String baseName = RawFiles.baseName(dia);
		String resolvedPrefix = prefix != null ? prefix : new File(baseName).getName();

		File splitOn = MassListDecoyGenerator.resolveForSplit(massList, outputDirectory, resolvedPrefix,
				generateDecoys);

		Logger.logLine("Scoring " + dia.getName() + " against " + library.getName());
		ContextFeatureScorer.scoreFeatures(library, dia, fasta, baseName, splitOn.getAbsolutePath());

		File referenceFeatures = new File(baseName + "_reference.features.txt");
		File backgroundFeatures = new File(baseName + "_background.features.txt");

		return runTrainApply(referenceFeatures, backgroundFeatures, fasta, fdr, seed, outputDirectory, resolvedPrefix);
	}

	public static MProphetResult runTrainApply(File referenceFeatures, File backgroundFeatures, File fasta, float fdr,
			int seed, File outputDirectory, String prefix) throws Exception {

		File engineDirectory = ContextOptions.engineDirectory(outputDirectory, ENGINE_NAME);
		File workingDirectory = ContextOptions.subdirectory(engineDirectory, ContextOptions.WORK_DIRECTORY);
		File modelDirectory = ContextOptions.subdirectory(engineDirectory, ContextOptions.MODEL_DIRECTORY);

		PreparedFeatures prepared = PreparedFeatures.prepare(backgroundFeatures, referenceFeatures, workingDirectory,
				prefix);

		SearchParameters params = SearchParameterParser.getDefaultParametersObject();

		MProphetExecutionData backgroundData = buildExecutionData(prepared.getBackground(), fasta, engineDirectory,
				prefix + ".peptide.background", params);
		MProphetExecutionData referenceData = buildExecutionData(prepared.getReference(), fasta, engineDirectory,
				prefix + ".peptide.reference", params);

		Logger.logLine("Training mProphet LDA on the background...");
		MProphetResult backgroundResult = MProphetReiter.executeMProphetTSV(backgroundData, fdr, seed,
				params.getAAConstants(), FIRST_ROUND);
		LinearDiscriminantAnalysis backgroundLDA = backgroundResult.getLDA();

		writeModel(backgroundLDA, backgroundResult.getFeatureNames(), new File(modelDirectory, prefix + ".lda.txt"));

		Logger.logLine("Applying the background-trained LDA to the reference peptides...");
		MProphetResult referenceResult = MProphetReiter.executeMProphetTSVWithModel(referenceData, fdr, backgroundLDA,
				params.getAAConstants());

		Logger.logLine("Reference passing peptides: " + referenceResult.getPassingPeptides().size() + " of "
				+ prepared.getReferenceTargets() + " targets at " + (fdr * 100f) + "% FDR");
		Logger.logLine("Results are under " + engineDirectory.getAbsolutePath());

		return referenceResult;
	}

	private static void writeModel(LinearDiscriminantAnalysis lda, List<String> featureNames, File modelFile)
			throws IOException {
		double[] coefficients = lda.getCoefficients();
		try (BufferedWriter out = new BufferedWriter(new FileWriter(modelFile))) {
			out.write("feature\tweight\n");
			for (int i = 0; i < coefficients.length; i++) {
				String name = featureNames != null && i < featureNames.size() ? featureNames.get(i) : "feature" + i;
				out.write(name + "\t" + coefficients[i] + "\n");
			}
			out.write("m0\t" + lda.getConstant() + "\n");
		}
	}

	private static MProphetExecutionData buildExecutionData(File inputFeatures, File fasta, File outputDirectory,
			String outputPrefix, SearchParameters params) {
		File peptideOutputFile = new File(outputDirectory, outputPrefix + ".txt");
		File peptideDecoyFile = new File(outputDirectory, outputPrefix + ".decoy.txt");
		return new MProphetExecutionData(inputFeatures, fasta, peptideOutputFile, peptideDecoyFile, params);
	}

	private static File requiredFile(HashMap<String, String> arguments, String flag) {
		String value = arguments.get(flag);
		if (value == null) {
			Logger.errorLine("Missing required argument: " + flag);
			printHelp();
			System.exit(1);
		}
		return new File(value);
	}

	private static void printHelp() {
		Logger.timelessLogLine("ContextMProphetExecutor");
		Logger.timelessLogLine("Context-mode mProphet: train the LDA on background peptides, then transfer it");
		Logger.timelessLogLine("to the targeted reference peptides without retraining.");
		Logger.timelessLogLine("");
		Logger.timelessLogLine("End-to-end from a raw file:");
		Logger.timelessLogLine("  -i        <file>   acquisition: " + RawFiles.supportedExtensions());
		Logger.timelessLogLine("  -l        <file>   library (.elib preferred, .dlib accepted)");
		Logger.timelessLogLine("  -f        <file>   FASTA protein database");
		Logger.timelessLogLine("  -massList <file>   assay / mass-list .txt or .csv");
		Logger.timelessLogLine("");
		Logger.timelessLogLine("Or from feature files that were already split:");
		Logger.timelessLogLine("  -background <file> background feature TSV");
		Logger.timelessLogLine("  -reference  <file> reference feature TSV");
		Logger.timelessLogLine("  -f          <file> FASTA protein database");
		Logger.timelessLogLine("");
		Logger.timelessLogLine("Optional:");
		Logger.timelessLogLine("  -fdr      <float>  peptide FDR threshold (default: " + DEFAULT_FDR + ")");
		Logger.timelessLogLine("  -seed     <int>    random seed for LDA training (default: " + DEFAULT_SEED + ")");
		Logger.timelessLogLine("  -o        <dir>    output directory (default: next to the input)");
		Logger.timelessLogLine("  -prefix   <name>   base name for output files");
		Logger.timelessLogLine("  -generateDecoys    add entrapment decoys when the mass list has none;");
		Logger.timelessLogLine("                     a list that already has decoys is used unchanged");
		Logger.timelessLogLine("");
		Logger.timelessLogLine("Outputs, under <-o>/" + ENGINE_NAME + "/:");
		Logger.timelessLogLine("  <prefix>.peptide.reference.txt          reference target peptides");
		Logger.timelessLogLine("  <prefix>.peptide.reference.decoy.txt    reference decoy peptides");
		Logger.timelessLogLine("  <prefix>.peptide.background.txt         background target peptides");
		Logger.timelessLogLine("  <prefix>.peptide.background.decoy.txt   background decoy peptides");
		Logger.timelessLogLine("  model/<prefix>.lda.txt                  the discriminant it learned");
		Logger.timelessLogLine("  work/<prefix>.*.pin	                  the pruned tables the engine saw");
		Logger.timelessLogLine("");
		Logger.timelessLogLine("and, when -generateDecoys had work to do, in <-o> itself:");
		Logger.timelessLogLine("  <prefix>.assay.with_decoys.txt          the assay it built");
	}
}
