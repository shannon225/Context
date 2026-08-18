package org.searlelab.context.percolator;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

import org.searlelab.context.io.MassListDecoyGenerator;
import org.searlelab.context.io.RawFiles;
import org.searlelab.context.mprophet.ContextFeatureScorer;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class ContextPercolatorExecutor {

	public static ContextPercolatorResult runEndToEnd(File library, File fasta, File dia, File massList,
			PyIsoPEPRunner pyIsoPEP, float fdr, File outputDirectory, String prefix, boolean generateDecoys,
			HashMap<String, String> encyclopediaArgs) throws Exception {

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

		File backgroundFeatures = new File(baseName + "_background.features.txt");
		File referenceFeatures = new File(baseName + "_reference.features.txt");

		return ContextPercolator.trainAndApply(backgroundFeatures, referenceFeatures, fasta, encyclopediaArgs,
				pyIsoPEP, fdr, outputDirectory, resolvedPrefix);
	}


}
