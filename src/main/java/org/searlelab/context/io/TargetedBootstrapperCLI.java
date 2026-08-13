package org.searlelab.context.io;

public class TargetedBootstrapperCLI {
	private static final int DEFAULT_MAX_SEED = 1;
	private static final int DEFAULT_NUMBER_OF_PEPTIDES = 100;
	private static final float DEFAULT_WINDOW_WIDTH_RT = 0.50f;
	private static final double DEFAULT_HALF_WINDOW_WIDTH_MZ = 1.0;

	public static void main(String[] args) throws Throwable {
		if (args.length < 2 || args.length > 6) {
			System.err.println("Usage: " + "TargetedBootstrapperCLI "
					+ "<library file location> <.dia file location> "
					+ "\n[seed] [numberOfpeptides] [halfWindowWidthRT] [halfWindowWidthMz]");
			System.exit(1);
		}

		String libraryPath = args[0];
		String rawFilePath = args[1];

		// Default parameters
		int seed = DEFAULT_MAX_SEED;
		int numberOfPeptides = DEFAULT_NUMBER_OF_PEPTIDES; // number of Peptides per assay
		float halfWindowWidthRT = DEFAULT_WINDOW_WIDTH_RT;
		double halfWindowWidthMz = DEFAULT_HALF_WINDOW_WIDTH_MZ;

		// Parameters as input
		if (args.length >= 3) {
			seed = Integer.parseInt(args[2]);
		}

		if (args.length >= 4) {
			numberOfPeptides = Integer.parseInt(args[3]);
		}

		if (args.length >= 5) {
			halfWindowWidthRT = Float.parseFloat(args[4]);
		}

		if (args.length >= 6) {
			halfWindowWidthMz = Double.parseDouble(args[5]);
		}
		
		if (args.length >= 7) {
			System.out.println("There are too many arguments, retry with a max of 7 arguments.");
		}

		TargetedBootstrapper bootstrapper = new TargetedBootstrapper();
		bootstrapper.execute(libraryPath, rawFilePath, seed, numberOfPeptides, halfWindowWidthRT,
				halfWindowWidthMz);

	}

}
