package org.searlelab.context.datastructures;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class ContextFeatures {

	private final File background;
	private final File reference;
	private final EncyclopediaFeatures referenceTable;
	private final List<String> droppedFeatures;
	private final int backgroundTargets;
	private final int backgroundDecoys;
	private final int referenceTargets;
	private final int referenceDecoys;

	private ContextFeatures(File background, File reference, EncyclopediaFeatures referenceTable,
			List<String> droppedFeatures, int backgroundTargets, int backgroundDecoys, int referenceTargets,
			int referenceDecoys) {
		this.background = background;
		this.reference = reference;
		this.referenceTable = referenceTable;
		this.droppedFeatures = droppedFeatures;
		this.backgroundTargets = backgroundTargets;
		this.backgroundDecoys = backgroundDecoys;
		this.referenceTargets = referenceTargets;
		this.referenceDecoys = referenceDecoys;
	}

	public static ContextFeatures prepare(File backgroundFeatures, File referenceFeatures, File workingDirectory,
			String prefix) throws IOException {

		requireReadable(backgroundFeatures, "background feature file");
		requireReadable(referenceFeatures, "reference feature file");
		makeDirectory(workingDirectory);

		EncyclopediaFeatures.requireMatchingHeaders(backgroundFeatures, referenceFeatures);

		Logger.logLine("Reading background features: " + backgroundFeatures.getName());
		EncyclopediaFeatures background = EncyclopediaFeatures.read(backgroundFeatures);
		Logger.logLine("Reading reference features: " + referenceFeatures.getName());
		EncyclopediaFeatures reference = EncyclopediaFeatures.read(referenceFeatures);

		int[] backgroundCounts = countLabels(background);
		int[] referenceCounts = countLabels(reference);

		requireTrainable(backgroundCounts);
		requireCompetable(referenceCounts);

		Logger.logLine("Background: " + backgroundCounts[0] + " targets, " + backgroundCounts[1] + " decoys");
		Logger.logLine("Reference: " + referenceCounts[0] + " targets, " + referenceCounts[1] + " decoys");

		List<String> dropped = background.getNearConstantFeatureNames();
		if (!dropped.isEmpty()) {
			Logger.logLine("Dropping " + dropped.size() + " feature(s) that are constant across the background: "
					+ dropped);
		}
		Set<String> toDrop = new HashSet<>(dropped);

		File prunedBackground = new File(workingDirectory, prefix + ".background.pin");
		File prunedReference = new File(workingDirectory, prefix + ".reference.pin");
		background.writeWithout(prunedBackground, toDrop);
		reference.writeWithout(prunedReference, toDrop);

		return new ContextFeatures(prunedBackground, prunedReference, EncyclopediaFeatures.read(prunedReference),
				Collections.unmodifiableList(dropped), backgroundCounts[0], backgroundCounts[1], referenceCounts[0],
				referenceCounts[1]);
	}

	private static int[] countLabels(EncyclopediaFeatures table) {
		int targets = 0;
		int decoys = 0;
		for (String[] row : table.getRows()) {
			if (table.isDecoy(row)) {
				decoys++;
			} else {
				targets++;
			}
		}
		return new int[] { targets, decoys };
	}

	private static void requireTrainable(int[] counts) throws IOException {
		if (counts[0] == 0 || counts[1] == 0) {
			throw new IOException("The background has " + counts[0] + " targets and " + counts[1]
					+ " decoys; training needs both. Check that the mass list is not matching every peptide in the "
					+ "library.");
		}
	}

	private static void requireCompetable(int[] counts) throws IOException {
		if (counts[1] == 0) {
			throw new IOException("The reference set has " + counts[0] + " targets and no decoys.");
		}
	}

	public File getBackground() {
		return background;
	}

	public File getReference() {
		return reference;
	}

	public EncyclopediaFeatures getReferenceTable() {
		return referenceTable;
	}

	public List<String> getDroppedFeatures() {
		return droppedFeatures;
	}

	public int getBackgroundTargets() {
		return backgroundTargets;
	}

	public int getBackgroundDecoys() {
		return backgroundDecoys;
	}

	public int getReferenceTargets() {
		return referenceTargets;
	}

	public int getReferenceDecoys() {
		return referenceDecoys;
	}

	private static void makeDirectory(File directory) throws IOException {
		if (!directory.exists() && !directory.mkdirs()) {
			throw new IOException("Could not create directory " + directory.getAbsolutePath());
		}
	}

	private static void requireReadable(File file, String description) throws IOException {
		if (file == null || !file.exists() || !file.canRead()) {
			throw new IOException("Cannot read the " + description + ": "
					+ (file == null ? "null" : file.getAbsolutePath()));
		}
	}
}
