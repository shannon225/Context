package org.searlelab.context.io;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

public class ContextOptions {

	public static final String MODEL_DIRECTORY = "model";
	public static final String WORK_DIRECTORY = "work";

	private ContextOptions() {
	}

	public static boolean isEnabled(HashMap<String, String> arguments, String flag) {
		if (!arguments.containsKey(flag)) return false;
		String value = arguments.get(flag);
		if (value == null || value.trim().isEmpty()) return true;
		return Boolean.parseBoolean(value.trim());
	}

	public static File outputDirectory(HashMap<String, String> arguments, File fallback) {
		String value = arguments.get("-o");
		if (value == null) value = arguments.get("-outdir");
		if (value == null) return fallback != null ? fallback : new File(".");
		return new File(value);
	}

	public static File engineDirectory(File outputDirectory, String engineName) throws IOException {
		File directory = new File(outputDirectory, engineName);
		makeDirectory(directory);
		return directory;
	}

	public static File subdirectory(File engineDirectory, String name) throws IOException {
		File directory = new File(engineDirectory, name);
		makeDirectory(directory);
		return directory;
	}

	public static String stripSplitSuffix(File featureFile) {
		return featureFile.getName().replaceFirst("_(background|reference)\\.features\\.txt$", "");
	}

	public static void makeDirectory(File directory) throws IOException {
		if (!directory.exists() && !directory.mkdirs()) {
			throw new IOException("Could not create directory " + directory.getAbsolutePath());
		}
	}
}
