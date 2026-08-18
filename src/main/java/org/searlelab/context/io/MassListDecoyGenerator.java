package org.searlelab.context.io;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import org.searlelab.context.datastructures.IsolationWindow;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;

public class MassListDecoyGenerator {

	public static final String GENERATED_SUFFIX = ".assay.with_decoys.txt";

	public static void main(String[] args) throws Exception {
		HashMap<String, String> arguments = CommandLineParser.parseArguments(args);

		if (args.length == 0 || arguments.containsKey("-h") || arguments.containsKey("-help")
				|| arguments.containsKey("--help")) {
			printHelp();
			System.exit(args.length == 0 ? 1 : 0);
		}

		String massListPath = arguments.get("-massList");
		if (massListPath == null) {
			Logger.errorLine("Missing required argument: -massList");
			printHelp();
			System.exit(1);
		}
		File massList = new File(massListPath);
		if (!massList.exists()) {
			Logger.errorLine("Mass list not found: " + massList.getAbsolutePath());
			System.exit(1);
		}

		String outputValue = arguments.get("-o");
		if (outputValue == null) outputValue = arguments.get("-outdir");
		File output = outputValue == null
				? defaultOutput(massList.getAbsoluteFile().getParentFile(), baseName(massList))
				: new File(outputValue);
		if (output.isDirectory()) output = defaultOutput(output, baseName(massList));

		write(addDecoys(massList.getAbsolutePath()), output);
	}

	public static File ensureDecoys(File massList, File outputDirectory, String prefix) throws IOException {
		if (hasDecoys(massList.getAbsolutePath())) {
			Logger.logLine("Mass list " + massList.getName() + " already contains decoys; using it as it is");
			return massList;
		}

		DirectoryOptions.makeDirectory(outputDirectory);
		File output = new File(outputDirectory, prefix + GENERATED_SUFFIX);
		Logger.logLine("Mass list " + massList.getName() + " has no decoys; generating entrapment decoys");
		write(addDecoys(massList.getAbsolutePath()), output);
		return output;
	}

	public static File resolveForSplit(File massList, File outputDirectory, String prefix, boolean generateDecoys)
			throws IOException {
		if (generateDecoys) {
			return ensureDecoys(massList, outputDirectory, prefix);
		}
		if (!hasDecoys(massList.getAbsolutePath())) {
			Logger.logLine("Note: " + massList.getName() + " has no decoys. Add -generateDecoys to build them.");
		}
		return massList;
	}

	public static boolean hasDecoys(String massListPath) {
		for (IsolationWindow window : IsolationWindowReader.parseMassList(massListPath)) {
			if (window.isDecoy()) return true;
		}
		return false;
	}

	private static void write(ArrayList<IsolationWindow> windows, File output) throws IOException {
		Path parent = output.getAbsoluteFile().toPath().getParent();
		if (parent != null) DirectoryOptions.makeDirectory(parent.toFile());

		new TargetedBootstrapper().writeAssayList(windows, output.toPath());

		int decoys = 0;
		for (IsolationWindow window : windows) {
			if (window.isDecoy()) decoys++;
		}
		Logger.logLine("Wrote " + output.getAbsolutePath() + ": " + (windows.size() - decoys) + " targets, "
				+ decoys + " decoys");
	}

	private static File defaultOutput(File directory, String base) {
		return new File(directory, base + GENERATED_SUFFIX);
	}

	private static String baseName(File file) {
		return file.getName().replaceFirst("\\.[^.]+$", "");
	}

	public static ArrayList<IsolationWindow> addDecoys(String massListPath) {
		ArrayList<IsolationWindow> input = IsolationWindowReader.parseMassList(massListPath);

		for (IsolationWindow window : input) {
			if (window.isDecoy()) {
				Logger.logLine("Mass list already contains decoys; leaving it unchanged");
				return input;
			}
		}

		SearchParameters params = PecanParameterParser.getDefaultParametersObject();
		AminoAcidConstants constants = new AminoAcidConstants();

		HashSet<String> takenSequences = new HashSet<>();
		ArrayList<IsolationWindow> output = new ArrayList<>();

		for (IsolationWindow target : input) {
			output.add(target);

			String sequence = target.getCompound();
			byte charge = target.getCharge();
			takenSequences.add(sequence);

			String decoy = PeptideUtils.getSmartDecoy(sequence, charge, takenSequences, params);
			String correctedDecoyMass = PeptideUtils.getCorrectedMasses(decoy, constants);
			double decoyMz = constants.getChargedMass(correctedDecoyMass, charge);
			takenSequences.add(decoy);

			output.add(new IsolationWindow(decoy, decoyMz, charge, target.getRtMin(), target.getRtMax(), true));
		}

		return output;
	}

	private static void printHelp() {
		Logger.timelessLogLine("MassListDecoyGenerator");
		Logger.timelessLogLine("Add entrapment decoys to an assay that only lists targets, so that the");
		Logger.timelessLogLine("reference peptides have something to compete against.");
		Logger.timelessLogLine("");
		Logger.timelessLogLine("  -massList <file>   assay / mass list, .csv or .txt, 7 or 8 columns");
		Logger.timelessLogLine("  -o        <path>   file to write, or a directory to write into");
		Logger.timelessLogLine("                     (default: next to the input)");
		Logger.timelessLogLine("");
	}
}
