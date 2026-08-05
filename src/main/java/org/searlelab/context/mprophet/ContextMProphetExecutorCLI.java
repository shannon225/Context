package org.searlelab.context.mprophet;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class ContextMProphetExecutorCLI {
	
	private static final String MODE_CONTEXT = "context";
	private static final String MODE_FOLDER = "folder";
	private static final String MODE_MPROPHET = "mprophet";
	private static final Set<String> ALLOWED_ARGUMENTS = new HashSet<>();

	static {
		ALLOWED_ARGUMENTS.add("--library");
		ALLOWED_ARGUMENTS.add("--fasta");
		ALLOWED_ARGUMENTS.add("--dia");
		ALLOWED_ARGUMENTS.add("--mass-list");
		ALLOWED_ARGUMENTS.add("--mode");
//		ALLOWED_ARGUMENTS.add("--dia-folder");
	}

	public static void main(String[] args) {
		try {
			if (args.length < 3 || containsHelpFlag(args)) {
				printUsage();
				System.exit(1);
				return;
			}

			Map<String, String> arguments = parseArguments(args);

			String libraryPath = requireArguments(arguments, "--library");
			String fastaPath = requireArguments(arguments, "--fasta");
			String diaFilePath = requireArguments(arguments, "--dia");
			String mode = arguments.getOrDefault("--mode", MODE_CONTEXT).toLowerCase();

			validateReadableFile("Library", libraryPath);
			validateReadableFile("FASTA", fastaPath);

			if (MODE_CONTEXT.equals(mode)) {
				
				String massListPath = requireArguments(arguments, "--mass-list");
				validateReadableFile("Mass list", massListPath);
				validateReadableFile("DIA", diaFilePath);

				ContextMProphetExecutor.executeContextMProphet(libraryPath, fastaPath, diaFilePath, massListPath);
				
			} else if (MODE_FOLDER.equals(mode)) {
				
				String diaFolderPath = requireArguments(arguments, "--dia");
				File diaFolder = validateDirectory("DIA folder", diaFolderPath);

				ContextMProphetExecutor.executeContextMProphetOnFolder(libraryPath, fastaPath, diaFolder);
				
			} else if (MODE_MPROPHET.equals(mode)) {
				
				File diaFolder = new File(diaFilePath).getAbsoluteFile().getParentFile();
				String massListPath = requireArguments(arguments, "--mass-list");
				validateReadableFile("Mass list", massListPath);
				validateReadableFile("DIA", diaFilePath);

				ContextMProphetExecutor.executeMProphet(libraryPath, fastaPath, diaFilePath, massListPath);
		
			} else {
				throw new IllegalArgumentException(
						"Unknown mode '" + mode + "'. Expected context, folder, or mprophet.");
			}
		} catch (IllegalArgumentException e) {
			System.err.println("Error: " + e.getMessage());
			System.err.println();
			printUsage();
			System.exit(1);
		}
	}

	private static Map<String, String> parseArguments(String[] args) {
		Map<String, String> arguments = new HashMap<>();

		for (int index = 0; index < args.length; index++) {
			String flag = args[index];

			if (!flag.startsWith("--")) {
				throw new IllegalArgumentException("Expected a named flag, but found '" + flag + "'.");
			}
			if (!ALLOWED_ARGUMENTS.contains(flag)) {
				throw new IllegalArgumentException("Unknown argument: " + flag + ".");
			}
			if (index + 1 >= args.length || args[index + 1].startsWith("--")) {
				throw new IllegalArgumentException("Missing input for " + flag + ".");
			}
			if (arguments.containsKey(flag)) {
				throw new IllegalArgumentException("Argument supplied more than once: " + flag + ".");
			}

			arguments.put(flag, args[++index]);
		}

		return arguments;
	}

	private static String requireArguments(Map<String, String> arguments, String flag) {
		String value = arguments.get(flag);
		if (value == null || value.trim().isEmpty()) {
			throw new IllegalArgumentException("Required arguments are missing: " + flag + ".");
		}
		return value;
	}

	private static void validateReadableFile(String description, String path) {
		File file = new File(path);
		if (!file.isFile() || !file.canRead()) {
			throw new IllegalArgumentException(
					description + " is not a readable file: " + file.getAbsolutePath());
		}
	}

	private static File validateDirectory(String description, String path) {
		File directory = new File(path);
		if (!directory.isDirectory()) {
			throw new IllegalArgumentException(
					description + " is not a directory: " + directory.getAbsolutePath());
		}
		return directory;
	}

	private static boolean containsHelpFlag(String[] args) {
		for (String argument : args) {
			if ("-h".equals(argument) || "--help".equals(argument)) {
				return true;
			}
		}
		return false;
	}

	private static void printUsage() {
		System.out.println("Usage:");
		System.out.println("  java org.searlelab.context.mprophet.ContextMProphetCLI \\");
		System.out.println("      --library <library.elib> \\");
		System.out.println("      --fasta <proteins.fasta> \\");
		System.out.println("      --dia <input.dia> \\");
		System.out.println("      --mass-list <assay.txt> \\");
		System.out.println("      [--mode context|folder|mprophet] \\");
		System.out.println("      [--dia-folder <folder>]");
		System.out.println();
		System.out.println("Modes:");
		System.out.println("  context    Run executeContextMProphet (default).");
		System.out.println("  folder     Run executeContextMProphetOnFolder; requires --dia-folder.");
		System.out.println("  mprophet   Run executeMProphet.");
	}
}

