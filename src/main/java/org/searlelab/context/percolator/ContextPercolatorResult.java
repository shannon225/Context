package org.searlelab.context.percolator;

import java.io.File;

public class ContextPercolatorResult {
	private final PercolatorWeights model;
	private final File nativeWeightsFile;
	private final File averagedWeightsFile;
	private final File rescoredFeaturesFile;
	private final File psmOutputFile;
	private final File peptideOutputFile;
	private final File rawPsmReportFile;
	private final File rawPeptideReportFile;
	private final int referencePSMCount;
	private final int referencePeptideCount;
	private final int passingPeptideCount;
	private final float fdrThreshold;

	ContextPercolatorResult(PercolatorWeights model, File nativeWeightsFile, File averagedWeightsFile,
			File rescoredFeaturesFile, File psmOutputFile, File peptideOutputFile, File rawPsmReportFile,
			File rawPeptideReportFile, int referencePSMCount, int referencePeptideCount, int passingPeptideCount,
			float fdrThreshold) {
		this.model = model;
		this.nativeWeightsFile = nativeWeightsFile;
		this.averagedWeightsFile = averagedWeightsFile;
		this.rescoredFeaturesFile = rescoredFeaturesFile;
		this.psmOutputFile = psmOutputFile;
		this.peptideOutputFile = peptideOutputFile;
		this.rawPsmReportFile = rawPsmReportFile;
		this.rawPeptideReportFile = rawPeptideReportFile;
		this.referencePSMCount = referencePSMCount;
		this.referencePeptideCount = referencePeptideCount;
		this.passingPeptideCount = passingPeptideCount;
		this.fdrThreshold = fdrThreshold;
	}

	public PercolatorWeights getModel() {
		return model;
	}

	public File getNativeWeightsFile() {
		return nativeWeightsFile;
	}

	public File getAveragedWeightsFile() {
		return averagedWeightsFile;
	}

	public File getRescoredFeaturesFile() {
		return rescoredFeaturesFile;
	}

	public File getPsmOutputFile() {
		return psmOutputFile;
	}

	public File getPeptideOutputFile() {
		return peptideOutputFile;
	}

	public File getRawPsmReportFile() {
		return rawPsmReportFile;
	}

	public File getRawPeptideReportFile() {
		return rawPeptideReportFile;
	}

	public int getReferencePSMCount() {
		return referencePSMCount;
	}

	public int getReferencePeptideCount() {
		return referencePeptideCount;
	}

	public int getPassingPeptideCount() {
		return passingPeptideCount;
	}

	public float getFdrThreshold() {
		return fdrThreshold;
	}

	@Override
	public String toString() {
		return "ContextPercolatorResult[" + passingPeptideCount + "/" + referencePeptideCount
				+ " reference peptides at " + (fdrThreshold * 100f) + "% FDR]";
	}
}
