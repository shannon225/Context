package org.searlelab.context.io;

import java.io.File;
import java.io.IOException;
import java.util.Locale;

public class RawFiles {

	private static final String[] EXTENSIONS = { ".mzml", ".dia", ".raw", ".d" };

	private RawFiles() {
	}

	public static String stripExtension(String path) {
		if (path == null) return null;

		String lower = path.toLowerCase(Locale.ROOT);
		for (String extension : EXTENSIONS) {
			if (lower.endsWith(extension)) {
				return path.substring(0, path.length() - extension.length());
			}
		}
		return path;
	}

	public static String baseName(File acquisition) {
		return stripExtension(acquisition.getAbsolutePath());
	}

	public static boolean isRecognised(String path) {
		if (path == null) return false;
		String lower = path.toLowerCase(Locale.ROOT);
		for (String extension : EXTENSIONS) {
			if (lower.endsWith(extension)) return true;
		}
		return false;
	}

	public static String supportedExtensions() {
		if (hasThermoSupport()) {
			return ".raw (Thermo), .d (Bruker), .dia, .mzML";
		}
		return ".d (Bruker), .dia, .mzML (convert Thermo .raw with msconvert first)";
	}

	public static boolean hasThermoSupport() {
		return RawFiles.class.getClassLoader().getResource(THERMO_MARKER) != null;
	}

	private static final String THERMO_MARKER = "msraw/thermo/bin/common/CommandLine.dll";

	public static void requireSupported(File acquisition) throws IOException {
		String name = acquisition.getName();
		if (name.toLowerCase(Locale.ROOT).endsWith(".raw") && !hasThermoSupport()) {
			throw new IOException("Context cannot read Thermo .raw files: it is built against "
					+ "msrawjava-core-nolice, which omits Thermo's RawFileReader for licensing reasons.\n"
					+ "Convert " + name + " to mzML first,"
					+ "then run Context on the resulting .mzML.");
		}
	}
}
