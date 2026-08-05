package org.searlelab.context.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.ArrayList;
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
			ArrayList<IsolationWindow> isolationWindows = selectMask(numberOfPeptides, aaConstants, seed, libraryPath,halfWindowWidthRT, halfWindowWidthMz);

			Path maskedFileOutputPath = rawFile.getParent().resolve(baseName + "_masked" + seed + "_assay.dia");
			Path assayOutputPath = rawFile.getParent().resolve(baseName + "_masked" + seed + "_assay.txt");

			writeMaskedFile(isolationWindows, seed, rawFilePath, maskedFileOutputPath, halfWindowWidthMz);

			writeAssayList(isolationWindows, assayOutputPath);
		}
	}

	// First function - Randomly Selects Precursors from a library and compiles them
	// into a list

	public ArrayList<IsolationWindow> selectMask(int numberOfPeptides, AminoAcidConstants aaConstants, int i,
			String libraryPath, float halfWindowWidthRT, double halfWindowWidthMz)
					throws IOException, SQLException, Throwable {

		// START TIMER 1
		long startTime = System.nanoTime();
		ArrayList<IsolationWindow> isolationWindows = new ArrayList<>();

		// For mapping targets and decoys later
	//	HashMap<String, String> targetDecoyOriginMap = new HashMap<>();

		// Set parameters before the loop
		int randomValue = 1 + i; // Add haliburton's number to get a random number

		HashSet<Integer> simulatedAssaySet = new HashSet<>();
		HashSet<String> targetSequencesSelected = new HashSet<>();
		HashSet<String> decoySequencesSelected = new HashSet<>();
		HashMap<String, String> discardedTDPairs = new HashMap<>();

		LibraryFile library = new LibraryFile();
		File file = new File(libraryPath);
		library.openFile(file);

		// Randomly select precursors loop
		try {

			// Load all entries
			ArrayList<LibraryEntry> entries = library.getAllEntries(false, aaConstants);

			// Select target precursors
			while (simulatedAssaySet.size() < numberOfPeptides) {
				randomValue = RandomGenerator.randomInt(randomValue);
				int index = Math.abs(randomValue) % entries.size();
				simulatedAssaySet.add(index);
			}
			
			
			// While loop for selecting targets and decoys

			while (targetSequencesSelected.size() < numberOfPeptides) {

				for (Integer index : simulatedAssaySet) {

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

					float windowRTMin = rtMin - 240f;
					float windowRTMax = rtMax + 240f;

					Range rtRange = new Range(windowRTMin, windowRTMax);

					// Add target sequences to the isolationWindows object
					if (!targetSequencesSelected.contains(sequence) && !decoySequencesSelected.contains(sequence)) {

						// Mark the target as not having a decoy
						boolean decoyFound = false;

						// m/z tolerance for decoys
						double mzTolerancePPM = 10;

						while (!decoyFound && mzTolerancePPM < 320) { 
							double mzToleranceDa = targetMz * mzTolerancePPM / 1_000_000;
							double upperMz = targetMz + mzToleranceDa;
							double lowerMz = targetMz - mzToleranceDa;

							Range libraryMzRange = new Range(lowerMz, upperMz);

							ArrayList<LibraryEntry> candidateDecoys = library.getEntries(libraryMzRange, false,
									aaConstants);

							// Loop to find entrapment decoys at a different window from target peptides

							for (LibraryEntry candidate : candidateDecoys) {

								float candidateRT = candidate.getRetentionTime();
								String decoySequence = candidate.getAccuratePeptideModSeq(aaConstants);

								// If the candidate RT is outside of the target RT range, and the decoy sequence
								// is different from the target, then add the decoy
								if (!rtRange.contains(candidateRT) 
										&& !decoySequence.equals(sequence)
										&& !targetSequencesSelected.contains(sequence)
										&& !decoySequencesSelected.contains(decoySequence) // checks if the decoy has been used
										&& !targetSequencesSelected.contains(decoySequence) // checks if decoy is also a target
										&& !discardedTDPairs.containsValue(sequence)  // checks if target was used in target-decoy pair 
										&& !discardedTDPairs.containsKey(decoySequence)) {  // checks if decoy was used in target-decoy pair
									
									double decoyMz = candidate.getPrecursorMZ();
									byte decoyCharge = candidate.getPrecursorCharge();

									float decoyRTMin = candidateRT - (60f * (halfWindowWidthRT));
									float decoyRTMax = candidateRT + (60f * (halfWindowWidthRT));

									IsolationWindow decoyWindow = new IsolationWindow(decoySequence, decoyMz,
											decoyCharge, decoyRTMin, decoyRTMax, true);

									decoyFound = true;

									//  Add the target once a matching decoy is found

									IsolationWindow window = new IsolationWindow(sequence, targetMz, charge, rtMin,
											rtMax, false);
									isolationWindows.add(window);
									targetSequencesSelected.add(sequence);

									System.out.println(
											"Target is " + targetMz + " for " + sequence + " at " + rtCenter / 60);

									isolationWindows.add(decoyWindow);
									decoySequencesSelected.add(decoySequence);
									discardedTDPairs.put(decoySequence, sequence);

									// Print info to the console to see what decoys are selected
									System.out.println("Decoy candidate is " + decoyMz + " for " + decoySequence
											+ " at " + candidateRT / 60);

									break; // end loop after adding one decoy per target
								}
							}

							if (!decoyFound && mzTolerancePPM < 320) {
								mzTolerancePPM = mzTolerancePPM * 2;
								continue;
							}
						} // end finding decoy loop
					} else {
						
					}
				} // end indexing through the HashSet loop 
				
				int peptidesLeftToSelect = simulatedAssaySet.size() - targetSequencesSelected.size();
				int newAssaySetSize = numberOfPeptides + peptidesLeftToSelect;

  			  while (newAssaySetSize <entries.size() && simulatedAssaySet.size() < newAssaySetSize) {
					randomValue = RandomGenerator.randomInt(randomValue);
					int index = Math.abs(randomValue) % entries.size();
					simulatedAssaySet.add(index);
				}
			}
			
			
		} catch (Exception e) {
			System.out.println("There was an error with selecting precursors. Check the file path. If it fails again, you may not have enough peptides in your dataset to bootstrap a synthetic assay.");
			throw e;
		} finally {
	//		writeTargetDecoyMap(targetDecoyOriginMap, mapOutputPath);

			System.out.println(isolationWindows.size() + " Precursors marked for extraction.");
			
			library.close();

		}

		// END TIMER 1
		long endTime = System.nanoTime();
		long duration = endTime - startTime;
		System.out.println("randomlySelectPrecursors(): Time taken (ms) : " + duration / 1_000_000);

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

					ArrayList<org.searlelab.msrawjava.model.FragmentScan> fragmentScansFromWindow = rawLibraryFile.getStripes(windowMz, windowStartTime, windowStopTime, sqrt);
					ArrayList<FragmentScan> matchingScans = new ArrayList<>();

					// Add Fragment Scans
					for (FragmentScan scan : fragmentScansFromWindow) {
						double scanMz = scan.getPrecursorMZ();
//						float scanRT = scan.getScanStartTime();
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
	//			maskedFile.addMetadata(diaFilePath, diaFilePath);
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
