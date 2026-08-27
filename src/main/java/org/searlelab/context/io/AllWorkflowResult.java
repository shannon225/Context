package org.searlelab.context.io;

import org.searlelab.context.percolator.ContextPercolatorResult;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;

public final class AllWorkflowResult {
	
		private final ContextPercolatorResult contextPercolator;
		private final PercolatorExecutionData standardPercolator;
		private final MProphetResult contextMProphet;
		private final MProphetResult standardMProphet;

		AllWorkflowResult(ContextPercolatorResult contextPercolator, PercolatorExecutionData standardPercolator,
				MProphetResult contextMProphet, MProphetResult standardMProphet) {
			this.contextPercolator = contextPercolator;
			this.standardPercolator = standardPercolator;
			this.contextMProphet = contextMProphet;
			this.standardMProphet = standardMProphet;
		}

		public ContextPercolatorResult getContextPercolator() {
			return contextPercolator;
		}

		public PercolatorExecutionData getStandardPercolator() {
			return standardPercolator;
		}

		public MProphetResult getContextMProphet() {
			return contextMProphet;
		}

		public MProphetResult getStandardMProphet() {
			return standardMProphet;
		}
	
}

