package org.searlelab.context.io;

import org.searlelab.context.percolator.ContextPercolatorResult;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;

public final class AllWorkflowResult {
	
		private final ContextPercolatorResult contextPercolator;
		private final PercolatorExecutionData standardPercolator;
		private final PercolatorExecutionData targetPercolator;
		private final MProphetResult contextMProphet;
		private final MProphetResult standardMProphet;
		private final MProphetResult targetMProphet;

		AllWorkflowResult(ContextPercolatorResult contextPercolator, PercolatorExecutionData standardPercolator, PercolatorExecutionData targetPercolator,
				MProphetResult contextMProphet, MProphetResult standardMProphet, MProphetResult targetMProphet) {
			this.contextPercolator = contextPercolator;
			this.standardPercolator = standardPercolator;
			this.targetPercolator = targetPercolator;
			this.contextMProphet = contextMProphet;
			this.standardMProphet = standardMProphet;
			this.targetMProphet = targetMProphet;
		}

		public ContextPercolatorResult getContextPercolator() {
			return contextPercolator;
		}

		public PercolatorExecutionData getStandardPercolator() {
			return standardPercolator;
		}
		
		public PercolatorExecutionData getTargetPercolator() {
			return targetPercolator;
		}

		public MProphetResult getContextMProphet() {
			return contextMProphet;
		}

		public MProphetResult getStandardMProphet() {
			return standardMProphet;
		}
		
		public MProphetResult getTargetMProphet() {
			return targetMProphet;
		}
	
}

