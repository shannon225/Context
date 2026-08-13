package org.searlelab.context.mprophet;

import java.io.File;

import org.searlelab.context.encyclopedia.MProphetReiter;
import org.searlelab.context.mprophet.ContextFeatureScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearDiscriminantAnalysis;

public class ContextMProphetExecutor {

	public static void executeContextMProphet(String libraryPath, String fastaPath, String diaFilePath, String massListPath) {
		File fasta = new File(fastaPath);
		File diaFile = new File(diaFilePath);
		File library = new File(libraryPath);

		String baseName = diaFilePath.replaceFirst("\\.dia$", "");

		SearchParameters params = SearchParameterParser.getDefaultParametersObject();


		// Score features in the .dia file against the library, split the results
		try {
			ContextFeatureScorer.scoreFeatures(library, diaFile, fasta, baseName, massListPath); // run this if the feature file hasn't been processed yet
			String featureFileName = baseName.replaceAll("\\.txt$", "");

			File backgroundFeatureFile = new File(featureFileName + "_background.features.txt");
			File referenceFeatureFile = new File(featureFileName + "_reference.features.txt");

			MProphetExecutionData backgroundData = makeMProphetExecutionData(backgroundFeatureFile, fasta, params, ".pep");
			MProphetExecutionData referenceData = makeMProphetExecutionData(referenceFeatureFile, fasta, params, ".pep");

			float peptideFDRThreshold = 0.01f;
			int seed = 1;
			int round = 1;

			MProphetResult backgroundMProphetResult = MProphetReiter.executeMProphetTSV(backgroundData, peptideFDRThreshold, seed, params.getAAConstants(), round);
			LinearDiscriminantAnalysis backgroundLDA = backgroundMProphetResult.getLDA();

			// 	Use the background LDA model on the reference feature file without retraining
			MProphetResult referenceMProphetResult = MProphetReiter.executeMProphetTSVWithModel(referenceData, peptideFDRThreshold, backgroundLDA, params.getAAConstants());

			System.out.println("The lda model has been trained on background feature. Now we'll use reference features from " + referenceFeatureFile.getAbsolutePath());

			System.out.println("Finished scoring peptides with background-trained lda model. "
					+ "\nReference passing peptides: " + referenceMProphetResult.getPassingPeptides().size());

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static void executeContextMProphetOnFolder1(String libraryPath, String fastaPath, File diaFolder) {
		File fasta = new File(fastaPath);
		File library = new File(libraryPath);
		File[] diaFilesInFolder = diaFolder.listFiles();
		
		System.out.println("Running ContextMProphetOnFolder for " + diaFolder.getAbsolutePath());
		
		SearchParameters params = SearchParameterParser.getDefaultParametersObject();

		// Score features in the .dia file against the library, split the results
		try {
			if (diaFilesInFolder != null) {
				for (File diaFile : diaFilesInFolder) {
				
					// Ignore files that do not end in .dia 
					if (!diaFile.isFile() || !diaFile.getName().endsWith(".dia")) {
						continue;
					}

					String diaName = diaFile.getName();
					String baseName = diaName.replaceFirst("\\.dia$", "");
					
					File massListFile = new File(diaFolder, baseName + ".txt");
					String massListPath = massListFile.getAbsolutePath();
	
					System.out.println("Processesing " + diaFile.getName());
					
					if (!massListFile.exists()) {
						System.out.println("Skipping " + diaFile.getName() + " because mass list was not found.");
						continue;
						
					}

					ContextFeatureScorer.scoreFeatures(library, diaFile, fasta, baseName, massListPath); // run this if the feature file hasn't been processed yet
					String featureFileName = baseName.replaceAll("\\.txt$", "");

					File backgroundFeatureFile = new File(featureFileName + "_background.features.txt");
					File referenceFeatureFile = new File(featureFileName + "_reference.features.txt");

					MProphetExecutionData backgroundData = makeMProphetExecutionData(backgroundFeatureFile, fasta, params, ".pep");
					MProphetExecutionData referenceData = makeMProphetExecutionData(referenceFeatureFile, fasta, params, ".pep");

					float peptideFDRThreshold = 0.01f;
					int seed = 1;
					int round = 1;

					MProphetResult backgroundMProphetResult = MProphetReiter.executeMProphetTSV(backgroundData, peptideFDRThreshold, seed, params.getAAConstants(), round);
					LinearDiscriminantAnalysis backgroundLDA = backgroundMProphetResult.getLDA();

					// 	Use the background LDA model on the reference feature file without retraining
					MProphetResult referenceMProphetResult = MProphetReiter.executeMProphetTSVWithModel(referenceData, peptideFDRThreshold, backgroundLDA, params.getAAConstants());

					System.out.println("The lda model has been trained on background feature. Now we'll use reference features from " + referenceFeatureFile.getAbsolutePath());
					System.out.println("Finished scoring peptides with background-trained lda model. "
							+ "\nReference passing peptides: " + referenceMProphetResult.getPassingPeptides().size());
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void executeContextMProphetOnFolder(String libraryPath, String fastaPath, File diaFolder) {
		File fasta = new File(fastaPath);
		File library = new File(libraryPath);
		File[] diaFilesInFolder = diaFolder.listFiles();
		
		System.out.println("Running ContextMProphetOnFolder for " + diaFolder.getAbsolutePath());
		
		SearchParameters params = SearchParameterParser.getDefaultParametersObject();

		// Score features in the .dia file against the library, split the results
		try {
			if (diaFilesInFolder != null) {
				for (File diaFile : diaFilesInFolder) {
				
					// Ignore files that do not end in .dia 
					if (!diaFile.isFile() || !diaFile.getName().endsWith(".dia")) {
						continue;
					}

			//		String diaName = diaFile.getName();
					String baseName = diaFile.getAbsolutePath().replaceFirst("\\.dia$", "");
					
					File massListFile = new File(diaFolder, baseName + ".txt");
					String massListPath = massListFile.getAbsolutePath();
	
					System.out.println("Processesing " + diaFile.getName());
					
					if (!massListFile.exists()) {
						System.out.println("Skipping " + diaFile.getName() + " because mass list was not found.");
						continue;
						
					}

					ContextFeatureScorer.scoreFeatures(library, diaFile, fasta, baseName, massListPath); // run this if the feature file hasn't been processed yet
					String featureFileName = baseName.replaceAll("\\.txt$", "");

					File backgroundFeatureFile = new File(featureFileName + "_background.features.txt");
					File referenceFeatureFile = new File(featureFileName + "_reference.features.txt");

					MProphetExecutionData backgroundData = makeMProphetExecutionData(backgroundFeatureFile, fasta, params, ".pep");
					MProphetExecutionData referenceData = makeMProphetExecutionData(referenceFeatureFile, fasta, params, ".pep");

					float peptideFDRThreshold = 0.01f;
					int seed = 1;
					int round = 1;

					MProphetResult backgroundMProphetResult = MProphetReiter.executeMProphetTSV(backgroundData, peptideFDRThreshold, seed, params.getAAConstants(), round);
					LinearDiscriminantAnalysis backgroundLDA = backgroundMProphetResult.getLDA();

					// 	Use the background LDA model on the reference feature file without retraining
					MProphetResult referenceMProphetResult = MProphetReiter.executeMProphetTSVWithModel(referenceData, peptideFDRThreshold, backgroundLDA, params.getAAConstants());

					System.out.println("The lda model has been trained on background feature. Now we'll use reference features from " + referenceFeatureFile.getAbsolutePath());
					System.out.println("Finished scoring peptides with background-trained lda model. "
							+ "\nReference passing peptides: " + referenceMProphetResult.getPassingPeptides().size());
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void executeMProphet(String libraryPath, String fastaPath, String diaFilePath, String massListPath) {
		
		File fasta = new File(fastaPath);
		File diaFile = new File(diaFilePath);
		File library = new File(libraryPath);

		String baseName = diaFilePath.replaceFirst("\\.dia$", "");

		SearchParameters params = SearchParameterParser.getDefaultParametersObject();

		// Score features in the .dia file against the library, split the results
		try {
			ContextFeatureScorer.scoreFeatures(library, diaFile, fasta, baseName, massListPath); // run this if the feature file hasn't been processed yet
			String featureFileName = baseName.replaceAll("\\.txt$", "");

			File featureFile = new File(featureFileName + ".features.txt");

			MProphetExecutionData featureData = makeMProphetExecutionData(featureFile, fasta, params, ".pep");


//			MProphetExecutionData referenceData = makeMProphetExecutionData(referenceFeatureFile, fasta, params, ".pep");
			float peptideFDRThreshold = 0.01f;
			int seed = 1;
			int round = 1;

			MProphetResult mprophetResult = MProphetReiter.executeMProphetTSV(featureData, peptideFDRThreshold, seed, params.getAAConstants(), round);

			System.out.println("The lda model has been trained on background feature. Now we'll use reference features from " + featureFile.getAbsolutePath());
			System.out.println("Finished scoring peptides with background-trained lda model. "
					+ "\nReference passing peptides: " + mprophetResult.getPassingPeptides().size());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}


	public static void executeMProphetOnFolder(String libraryPath, String fastaPath, File diaFolder) {

		File fasta = new File(fastaPath);
		File library = new File(libraryPath);
		File[] diaFilesInFolder = diaFolder.listFiles();

		System.out.println("Running mProphet on folder: " + diaFolder.getAbsolutePath());

		if (diaFilesInFolder == null) {
			System.out.println("Could not read DIA folder: " + diaFolder.getAbsolutePath());
			return;
		}

		SearchParameters params = SearchParameterParser.getDefaultParametersObject();

		float peptideFDRThreshold = 0.01f;
		int seed = 1;
		int round = 1;

		for (File diaFile : diaFilesInFolder) {
			
			// Ignore folders and files that do not end in .dia
			if (!diaFile.isFile() || !diaFile.getName().endsWith(".dia")) {
				continue;
			}

			String baseName = diaFile.getAbsolutePath().replaceFirst("\\.dia$", "");
			File massListFile = new File(baseName + ".txt");
			
			System.out.println("Processing " + diaFile.getName());

			if (!massListFile.isFile()) {
				System.out.println("Skipping " + diaFile.getName() + " because its mass list was not found: " + massListFile.getAbsolutePath());
				continue;
			}

			try {
				ContextFeatureScorer.scoreFeatures(library, diaFile, fasta, baseName, massListFile.getAbsolutePath());

				File featureFile = new File(baseName + ".features.txt");

				if (!featureFile.isFile()) {
					System.out.println("Skipping mProphet because the feature file " + "was not created: " + featureFile.getAbsolutePath());
					continue;
				}

				MProphetExecutionData featureData = makeMProphetExecutionData(featureFile, fasta, params, ".pep");

				MProphetResult mprophetResult = MProphetReiter.executeMProphetTSV(featureData, peptideFDRThreshold, seed,params.getAAConstants(),round);

				System.out.println("Finished mProphet analysis for "+ diaFile.getName()+ "\nPassing peptides: "+ mprophetResult.getPassingPeptides().size());

			} catch (Exception e) {
				
				System.err.println("mProphet failed for "+ diaFile.getAbsolutePath());
				e.printStackTrace();
				
			} finally {
				System.out.println("");
			}
		}
	}
	

	private static  MProphetExecutionData makeMProphetExecutionData(File inputFeatureFile, File fasta, SearchParameters params, String outputSuffix) {

		File peptideOutputFile = new File(inputFeatureFile.getAbsolutePath().replaceAll("\\.txt$", "") + outputSuffix + ".output.txt");
		File peptideDecoyFile = new File(inputFeatureFile.getAbsolutePath().replaceAll("\\.txt$", "") + outputSuffix + ".decoy.txt");

		return new MProphetExecutionData(inputFeatureFile, fasta, peptideOutputFile, peptideDecoyFile, params);
	}
}


