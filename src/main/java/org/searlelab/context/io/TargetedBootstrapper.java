package org.searlelab.context.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
//import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;

import org.searlelab.msrawjava.model.FragmentScan;
import org.searlelab.msrawjava.model.PrecursorScan;
import org.searlelab.msrawjava.model.WindowData;
import org.searlelab.msrawjava.model.Range;
import org.searlelab.msrawjava.io.encyclopedia.EncyclopeDIAFile;

import org.searlelab.context.mprophet.IsolationWindow;

public class TargetedBootstrapper {

	public void execute(String libraryPath, String rawFilePath, int maxSeed, int numberOfPeptides,
			float halfWindowWidthRT, double halfWindowWidthMz) throws Throwable {
		Path rawFile = Paths.get(rawFilePath);
		String rawFileName = rawFile.getFileName().toString();
		String baseName = rawFileName.replaceFirst("\\.dia$", "");

		AminoAcidConstants aaConstants = new AminoAcidConstants();

		for (int seed = 1; seed <= maxSeed; seed++) {
			ArrayList<IsolationWindow> isolationWindows = selectMask(numberOfPeptides, aaConstants, seed, libraryPath,
					halfWindowWidthRT, halfWindowWidthMz);

			Path maskedFileOutputPath = rawFile.getParent().resolve(baseName + "_masked" + seed + "_assay.dia");
			Path assayOutputPath = rawFile.getParent().resolve(baseName + "_masked" + seed + "_assay.txt");

			writeMaskedFile(isolationWindows, seed, rawFilePath, maskedFileOutputPath, halfWindowWidthMz);

			writeAssayList(isolationWindows, assayOutputPath);
		}
	}

	public ArrayList<IsolationWindow> selectTargets(int numberOfPeptides, AminoAcidConstants aaConstants, int i,
			String libraryPath, float halfWindowWidthRT, double halfWindowWidthMz, HashSet<String> excludedTargetSequences) throws Exception {

		ArrayList<IsolationWindow> targetWindows = new ArrayList<>();

		int randomValue = 1 + i; // Add haliburton's number to get a random number

		HashSet<Integer> simulatedAssaySet = new HashSet<>();
		HashSet<String> targetSequencesSelected = new HashSet<>();

		LibraryFile library = new LibraryFile();
		File file = new File(libraryPath);
		library.openFile(file);

		// Randomly select precursors loop
		try {
			// Load all entries
			ArrayList<LibraryEntry> entries = library.getAllEntries(false, aaConstants);

			System.out.println("There are " + entries.size() + " entries from the library.");
//			System.out.println("The highest m/z is " + Collections.max(entries).getPrecursorMZ());
//			System.out.println("The lowest m/z is " + Collections.min(entries).getPrecursorMZ());
			
			if (numberOfPeptides > entries.size()) {
				throw new IllegalArgumentException("The requested number of targets is larger than the number of library entries. Select a smaller number of peptides to put into each assay.");
			}

			// Select target precursors
			while (simulatedAssaySet.size() < numberOfPeptides) {
				randomValue = RandomGenerator.randomInt(randomValue);
				int index = Math.abs(randomValue) % entries.size();
				simulatedAssaySet.add(index);
			}

			// While loop for selecting targets and decoys
			while (targetSequencesSelected.size() < numberOfPeptides) {
				for (Integer index : simulatedAssaySet) {

					if (targetSequencesSelected.size() >= numberOfPeptides) {
						break;
					}
					
					// Retrieve the library entry at the random index
					LibraryEntry entry = entries.get(index);
					
					// Get the m/z, RT and sequence
					double targetMz = entry.getPrecursorMZ();
					float rtCenter = entry.getRetentionTimeInSec();
					String sequence = entry.getAccuratePeptideModSeq(aaConstants);
					byte charge = entry.getPrecursorCharge();

					// Calculate a RT ranges for the isolationWindows object
					float rtMin = (float) (rtCenter - (60f * (halfWindowWidthRT)));
					float rtMax = (float) (rtCenter + (60f * (halfWindowWidthRT)));

					if (!targetSequencesSelected.contains(sequence) && !excludedTargetSequences.contains(sequence)) {

						IsolationWindow window = new IsolationWindow(sequence, targetMz, charge, rtMin, rtMax, false);
						targetWindows.add(window);
						targetSequencesSelected.add(sequence);

						System.out.println("Target is " + targetMz + " for " + sequence + " at " + rtCenter / 60);
					}
				} // If there aren't enough peptides selected in the target list, select more indicies
				
				if (targetSequencesSelected.size() < numberOfPeptides) {
					
					if (simulatedAssaySet.size() >= entries.size()) {
						throw new IllegalStateException("The library does not contain enough unique peptides to bootstrap from the assay. Select a new library");
					}
					
					int previousSetSize = simulatedAssaySet.size(); 
					
					while (simulatedAssaySet.size() == previousSetSize) {
						randomValue = RandomGenerator.randomInt(randomValue);
						
						int index = Math.abs(randomValue) % entries.size();
						simulatedAssaySet.add(index);
					}
				}
			}
		} catch (Exception e) {
			System.out.println(
					"There was an error with selecting targets. Check the file path for the library. If it fails again, you may not have enough peptides in your dataset to bootstrap a synthetic assay.");
			throw e;
		} finally {
			System.out.println(targetWindows.size() + " Precursors marked for extraction.");
			library.close();
		}
		return targetWindows;
	}

	public ArrayList<IsolationWindow> selectDecoys(ArrayList<IsolationWindow> targetWindows, AminoAcidConstants aaConstants, int i, String libraryPath, float halfWindowWidthRT, double halfWindowWidthMz, ArrayList<IsolationWindow> unmatchedTargets) throws Exception {

		ArrayList<IsolationWindow> decoyWindows = new ArrayList<IsolationWindow>();
		
		HashSet<String> targetSequencesSelected = new HashSet<>();
		HashSet<String> decoySequencesSelected = new HashSet<>();
		HashMap<String, String> discardedTDPairs = new HashMap<>();
		
		LibraryFile library = new LibraryFile();
		File file = new File(libraryPath);
		library.openFile(file);

		try {
			for (IsolationWindow targetWindow : targetWindows) {
				targetSequencesSelected.add(targetWindow.getCompound());
			}
			
			for (IsolationWindow targetWindow : targetWindows) {

				// Get the m/z, RT and sequence for targets
				double targetMz = targetWindow.getTargetMz();
				String targetSequence = targetWindow.getCompound();
				
				boolean decoyFound = false;
				double mzTolerancePPM = 10;
				
				while (!decoyFound && mzTolerancePPM <= 320) {
					double mzToleranceDa = targetMz * mzTolerancePPM / 1_000_000;
					double upperMz = targetMz + mzToleranceDa;
					double lowerMz = targetMz - mzToleranceDa;

					Range libraryMzRange = new Range(lowerMz, upperMz);

					ArrayList<LibraryEntry> candidateDecoys = library.getEntries(libraryMzRange, false, aaConstants);
	
					// Loop to find entrapment decoys at a different window from target peptides
					for (LibraryEntry candidate : candidateDecoys) {
						double candidateMinMz = candidate.getPrecursorMZ() - halfWindowWidthMz;
						double candidateMaxMz = candidate.getPrecursorMZ() + halfWindowWidthMz;

						float candidateRT = candidate.getRetentionTime();
						float candidateMinRT = candidateRT - (halfWindowWidthRT * 60f);
						float candidateMaxRT = candidateRT + (halfWindowWidthRT * 60f);
						
						String decoySequence = candidate.getAccuratePeptideModSeq(aaConstants);
						boolean decoyWindowOverlapsTarget = false;
						
					for (IsolationWindow comparisonTarget : targetWindows) {
						
						// Calculate isolation window for each window
						double comparisonTargetMzMin = comparisonTarget.getTargetMz() - (10 * halfWindowWidthMz);
						double comparisonTargetMzMax = comparisonTarget.getTargetMz() + (10 * halfWindowWidthMz);

							
						// Calculate target windows + buffer range - targets will be selected outside of this window
						float minTargetRTWithBuffer = comparisonTarget.getRtMin() - (halfWindowWidthRT * 60f);
						float maxTargetRTWithBuffer = comparisonTarget.getRtMax() + (halfWindowWidthRT * 60f);
						
						boolean overlapsThisTargetRT = candidateMinRT <= maxTargetRTWithBuffer && candidateMaxRT >= minTargetRTWithBuffer;
						boolean overlapsThisTargetMz = candidateMinMz <= comparisonTargetMzMax && candidateMaxMz >= comparisonTargetMzMin;

						if (overlapsThisTargetRT && overlapsThisTargetMz) {
							decoyWindowOverlapsTarget = true;
							break;
						}
					}
						// If the candidate RT is outside of the target RT range, and the decoy sequence
						// is different from the target, then add the decoys
						if (!decoyWindowOverlapsTarget 
								&& !decoySequence.equals(targetSequence) 
								&& !targetSequencesSelected.contains(decoySequence)
								&& !decoySequencesSelected.contains(decoySequence) // has the decoy been used? 
								&& !discardedTDPairs.containsValue(targetSequence) // is the target used in target-decoy pair
								&& !discardedTDPairs.containsKey(decoySequence)) { // checks if decoy was used in target-decoy pair

							double decoyMz = candidate.getPrecursorMZ();
							byte decoyCharge = candidate.getPrecursorCharge();

							float decoyRTMin = candidateRT - (60f * (halfWindowWidthRT));
							float decoyRTMax = candidateRT + (60f * (halfWindowWidthRT));

							IsolationWindow decoyWindow = new IsolationWindow(decoySequence, decoyMz, decoyCharge, decoyRTMin, decoyRTMax, true);
							decoyWindows.add(decoyWindow);
							
							decoySequencesSelected.add(decoySequence);
							discardedTDPairs.put(decoySequence, targetSequence);
							decoyFound = true;

							// Print info to the console to see what decoys are selected
							System.out.println("Decoy candidate is " + decoyMz + " for " + decoySequence + " at "	+ candidateRT / 60);

							break; // end loop after adding one decoy per target
						}
					}
					if (!decoyFound) {
						mzTolerancePPM = mzTolerancePPM *2;
					}
				}
				
				if (!decoyFound) {
					unmatchedTargets.add(targetWindow);
				//	System.out.println("No valid decoy was found per target " + targetSequence + " after searching up to " + mzTolerancePPM + " ppm.");
					continue;
				}
			}

		} catch (Exception e) {
			throw e;
		} finally {
			System.out.println(decoyWindows.size() + " decoy windows selected.");
			library.close();
		}
		
		return decoyWindows;

	}

	// First function - Randomly Selects Precursors from a library, places on a list
	public ArrayList<IsolationWindow> selectMask(int numberOfPeptides, AminoAcidConstants aaConstants, int i, String libraryPath, float halfWindowWidthRT, double halfWindowWidthMz)
					throws IOException, SQLException, Throwable {

		// START TIMER 1
		long startTime = System.nanoTime();
		HashSet<String> rejectedTargetSequences = new HashSet<>();
		
		ArrayList<IsolationWindow> targetWindows = selectTargets(numberOfPeptides, aaConstants, i, libraryPath, halfWindowWidthRT, halfWindowWidthMz, rejectedTargetSequences);
		ArrayList<IsolationWindow> decoyWindows;
		
		int replacementRound = 0; 
		
		while(true) {
			ArrayList<IsolationWindow> unmatchedTargets = new ArrayList<>();
			
			// Select deoys
			decoyWindows = selectDecoys(targetWindows, aaConstants, i, libraryPath, halfWindowWidthRT, halfWindowWidthMz, unmatchedTargets);
		
			if (unmatchedTargets.isEmpty()) {
				break;
			}
			
			// Remove targets that could not find appropriate decoys
			for (IsolationWindow unmatchedTarget : unmatchedTargets) {
				targetWindows.remove(unmatchedTarget);
				
				rejectedTargetSequences.add(unmatchedTarget.getCompound());
			}
			
			HashSet<String> unavailableTargetSequences = new HashSet<>(rejectedTargetSequences);
			// Exclude the targets that are still in the assay so they can't be selected again 
			for (IsolationWindow existingTarget : targetWindows) {
				unavailableTargetSequences.add(existingTarget.getCompound());
			}
			
			int replacementsNeeded = numberOfPeptides - targetWindows.size();
			replacementRound++;
			
			ArrayList<IsolationWindow> replacementTargets = selectTargets(replacementsNeeded, aaConstants, i + replacementRound, libraryPath, halfWindowWidthRT, halfWindowWidthMz, unavailableTargetSequences);
			
			targetWindows.addAll(replacementTargets);
		}
		
		if (targetWindows.size() !=decoyWindows.size() ) {
			throw new IllegalStateException("The number of targets does not match the number of decoys");
		}
		
		ArrayList<IsolationWindow> isolationWindows = new ArrayList<>(targetWindows.size() + decoyWindows.size());
		
		for (int index = 0; index < targetWindows.size(); index++) {
			isolationWindows.add(targetWindows.get(index));
			isolationWindows.add(decoyWindows.get(index));
		}
		
		System.out.println(isolationWindows.size() + " total target and decoy windows selected.");
		
		long endTime = System.nanoTime();
		System.out.println("selectMask(): Time taken (ms): " + (endTime-startTime) / 1_000_000);
		// Randomly select precursors loop
		return isolationWindows;
}

// Second function - Uses the IsolationWindow List to mask the raw data
public EncyclopeDIAFile writeMaskedFile(ArrayList<IsolationWindow> isolationWindows, int i, String diaFilePath,
		Path outputPath, double halfWindowWidthMz) throws Throwable {
	// START TIMER 2
	long startTime = System.nanoTime();

	File rawFile = new File(diaFilePath);
	EncyclopeDIAFile maskedFile = new EncyclopeDIAFile();
	EncyclopeDIAFile rawLibraryFile = new EncyclopeDIAFile();
	File outputFile = outputPath.toFile();

	HashSet<Integer> addedPrecursors = new HashSet<>();
	HashSet<Integer> addedFragments = new HashSet<>();

	// System.out.println("Is the .dia file open? " + rawLibraryFile.isOpen());
	rawLibraryFile.openFile(rawFile);

	try {
		maskedFile.openFile();

		// Add Ranges
		HashMap<Range, WindowData> dutyCycleMap = new HashMap<>();
		System.out.println("Masking DIA file based on the selected precursors...");

		for (IsolationWindow window : isolationWindows) {
			boolean isDecoy = window.isDecoy();

			if (!isDecoy) {

				double windowMz = window.getTargetMz();
				float windowStartTime = window.getRtMin();
				float windowStopTime = window.getRtMax();
				boolean sqrt = false;
				double mzStart = windowMz - halfWindowWidthMz;
				double mzStop = windowMz + halfWindowWidthMz;
				Range mzRange = new Range(mzStart, mzStop);

				ArrayList<org.searlelab.msrawjava.model.FragmentScan> fragmentScansFromWindow = rawLibraryFile
						.getStripes(windowMz, windowStartTime, windowStopTime, sqrt);
				ArrayList<FragmentScan> matchingScans = new ArrayList<>();

				// Add Fragment Scans
				for (FragmentScan scan : fragmentScansFromWindow) {
					double scanMz = scan.getPrecursorMZ();
					int scanIndex = scan.getSpectrumIndex();
					if (mzRange.contains(scanMz) && !addedFragments.contains(scanIndex)) {
						matchingScans.add(scan);
						addedFragments.add(scanIndex);
					}
				}

				for (Entry<org.searlelab.msrawjava.model.Range, WindowData> entry : rawLibraryFile.getRanges()
						.entrySet()) {
					if (mzRange.contains(entry.getKey().getMiddle())) {
						dutyCycleMap.put(entry.getKey(), entry.getValue());
					}
				}

				maskedFile.setRanges(dutyCycleMap);
				maskedFile.addStripe(matchingScans);

				// Add Precursor Scans
				ArrayList<PrecursorScan> precursorScanFromWindow = rawLibraryFile.getPrecursors(windowStartTime,
						windowStopTime);

				ArrayList<PrecursorScan> matchingPrecursors = new ArrayList<>();

				for (PrecursorScan precursor : precursorScanFromWindow) {
					Range precursorRange = new Range(precursor.getIsolationWindowLower(),
							precursor.getIsolationWindowUpper());
					int spectrumIndex = precursor.getSpectrumIndex();
					if ((precursorRange.contains(mzRange) && !addedPrecursors.contains(spectrumIndex))) {
						matchingPrecursors.add(precursor);
						addedPrecursors.add(spectrumIndex);
					}
				}
				maskedFile.addPrecursor(matchingPrecursors);
			}
			maskedFile.setFileName(rawFile.getName(), null, rawFile.getAbsolutePath());
			maskedFile.addMetadata(rawLibraryFile.getMetadata());
			maskedFile.setFractionNames(rawLibraryFile.getFractionNames());
		}
	} catch (IOException e) {
		System.out.println("Unable to open raw file.");
		throw e;
	}

	// END TIMER 2
	long endTime = System.nanoTime();
	long duration = endTime - startTime;
	rawLibraryFile.close();

	System.out.println("maskDIAFileBasedOnIsolationWindows(): Time taken (ms) : " + duration / 1_000_000);

	maskedFile.saveAsFile(outputFile);
	System.out.println("Target mass list for the masked file  was written to " + outputPath
			+ "\n Number of added Precursor scans: " + addedPrecursors.size()
			+ "\n Number of added Fragment scans: " + addedFragments.size());

	return maskedFile;
}

	public void writeAssayList(ArrayList<IsolationWindow> isolationWindows, Path outputPath) throws IOException {
		try (BufferedWriter writer = Files.newBufferedWriter(outputPath)) {
			writer.write("Compound\tFormula\tAdduct\tm/z\tz\tRT Time (min)\tWindow (min)\tisDecoy");
			writer.newLine();

			for (IsolationWindow window : isolationWindows) {
				String compound = window.getCompound();
				double targetMz = window.getTargetMz();
				byte charge = window.getCharge();
				boolean isDecoy = window.isDecoy();

				float rtCenterMin = ((window.getRtMin() + window.getRtMax()) / 2.0f) / 60.0f;
				float windowMin = (window.getRtMax() - window.getRtMin()) / 60.0f;

				writer.write(compound + "\t" + "\t" + "(no adduct)" + "\t" + targetMz + "\t" + charge + "\t"
						+ rtCenterMin + "\t" + windowMin + "\t" + isDecoy);
				writer.newLine();
			}
		} catch (Exception e) {
			throw e;
		}
	}

	@SuppressWarnings("unused")
	private void writeTargetDecoyMap(HashMap<String, String> targetDecoyMap, Path outputPath) throws IOException {

		try (BufferedWriter writer = Files.newBufferedWriter(outputPath)) {

			writer.write("decoySequence\ttargetSequence");
			writer.newLine();

			for (Entry<String, String> entry : targetDecoyMap.entrySet()) {
				String decoySequence = entry.getKey();
				String targetSequence = entry.getValue();

				writer.write(decoySequence + "\t" + targetSequence);
				writer.newLine();
			}
		} catch (Exception e) {
			throw e;
		}
	}
}
