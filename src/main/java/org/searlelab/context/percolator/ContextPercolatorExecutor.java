package org.searlelab.context.percolator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.nio.file.Files;
import java.util.ArrayList;

import org.searlelab.context.io.ContextFeatureScorer;
import org.searlelab.context.io.DirectoryOptions;
import org.searlelab.context.io.MassListDecoyGenerator;
import org.searlelab.context.io.RawFiles;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;


public class ContextPercolatorExecutor {
	private static final String STANDARD_ENGINE_NAME = "standard-percolator";
	private static final int FINAL_ROUND = 2;
	
	public static ContextPercolatorResult runContextPercolator(File library, File fasta, File dia, File massList,
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
		ContextFeatureScorer.scoreFeaturesForContext(library, dia, fasta, baseName, splitOn.getAbsolutePath());

		File backgroundFeatures = new File(baseName + "_background.features.txt");
		File referenceFeatures = new File(baseName + "_reference.features.txt");

		return ContextPercolator.trainAndApply(backgroundFeatures, referenceFeatures, fasta, encyclopediaArgs,
				pyIsoPEP, fdr, outputDirectory, resolvedPrefix);
	}
	
	public static PercolatorExecutionData runStandardPercolator(File allFeatures, File fasta, PyIsoPEPRunner pyIsoPEP, float fdr, File outputDirectory, String prefix, HashMap<String, String> encyclopediaArgs) throws IOException, InterruptedException {
		
		if (!allFeatures.exists() || !allFeatures.canRead()) {
			throw new IOException("Feature file not found or is unreadable: " + allFeatures.getAbsolutePath());
		}
		
		if (!fasta.exists() || !fasta.canRead()) {
			throw new IOException("FASTA file not found or is unreadable: " + fasta);
		}
		
		if (!(fdr > 0.0f && fdr <= 1.0f)) {
			throw new IllegalArgumentException("FDR must a value between 0 and 1.");
		}
		
		HashMap<String, String> parameterMap = SearchParameterParser.getDefaultParametersObject().toParameterMap();
		
		if (encyclopediaArgs !=null) {
			parameterMap.putAll(encyclopediaArgs);
		}
		
		parameterMap.put("-percolatorThreshold",  Float.toString(fdr));
		SearchParameters parameters = SearchParameterParser.parseParameters(parameterMap);
		
		File engineDirectory = DirectoryOptions.engineDirectory(outputDirectory, STANDARD_ENGINE_NAME);
		File peptideTargets = new File(engineDirectory, prefix + ".peptide.target.txt");
		File peptideDecoys = new File(engineDirectory, prefix + ".peptide.decoy.txt");
		File proteinTargets = new File(engineDirectory, prefix + ".protein.target.txt");
		File proteinDecoys = new File(engineDirectory, prefix + ".protein.decoy.txt");
		
		PercolatorExecutionData run = new PercolatorExecutionData(allFeatures, fasta, peptideTargets, peptideDecoys, proteinTargets, proteinDecoys, parameters);
		
		 Files.deleteIfExists(run.getModelFile().toPath());
		 Files.deleteIfExists(run.getWeightsFile(FINAL_ROUND).toPath());
		 
		 Logger.logLine("Running standard Percolator cross-validation on " + allFeatures.getName());
		 
		 Pair<ArrayList<PercolatorPeptide>, Float> result = PercolatorExecutor.executePercolatorTSV(parameters.getPercolatorVersionNumber(), run, fdr, parameters.getAAConstants(), FINAL_ROUND);
		 
		 Logger.logLine("Standard Percolator run and resulted in " + result.x.size() + " target peptides below " + (fdr*100f) + "% FDR");
		 Logger.logLine("Results are under " + engineDirectory.getAbsolutePath());
		 
		 
		return run;
		
	}
	
	private static void preparePercolatorOutputForPyIsoPEP(File source, File destination) throws IOException {

	    try (BufferedReader reader = new BufferedReader(new FileReader(source));
	            BufferedWriter writer = new BufferedWriter(new FileWriter(destination))) {
	        String line;

	        while ((line = reader.readLine()) != null) {
	            if (line.startsWith(PercolatorExecutor.PI_0_TAG)) {
	                continue;
	            }

	            if (line.trim().isEmpty()) {
	                continue;
	            }

	            writer.write(line);
	            writer.newLine();
	        }
	    }
	}
	
	private static int countPassingPyIsoPEPPeptides(PyIsoPEPRunner.Table table,float fdr) throws IOException {

	    if (table.indexOf(PyIsoPEPRunner.Q_VALUE_COLUMN) < 0) {
	        throw new IOException("pyIsoPEP output is missing the q-value column: "+ PyIsoPEPRunner.Q_VALUE_COLUMN);
	    }

	    int passing = 0;

	    for (String[] row : table.getRows()) {
	        String value = table.get(row, PyIsoPEPRunner.Q_VALUE_COLUMN);

	        try {
	            if (Double.parseDouble(value.trim()) <= fdr) {
	                passing++;
	            }
	            
	        } catch (NumberFormatException e) {
	            throw new IOException(
	                    "Could not parse pyIsoPEP q-value: " + value,
	                    e
	            );
	        }
	    }

	    return passing;
	}
}
