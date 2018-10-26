package normCuts;

import org.apache.commons.math3.linear.MatrixUtils;

import normCuts.NormalizedCuts;

/**
 * This is a Normalized Cuts stub class that we are using to write a test for
 * the Normalized Cuts algorithm.
 * 
 * Hence this is a trivially simple implementation of Normalized Cuts.
 * 
 * @author Andrew P Talbot
 */
public class NormalizedCutsTestStub extends NormalizedCuts {
	
	private double[][] W;

	/**
	 * {@inheritDoc}
	 * 
	 * @param W
	 *            is the similarity matrix that we use to test for the
	 *            Normalized Cuts method.
	 */
	public NormalizedCutsTestStub(double sigmaI, double sigmaX, double r, int l, double nCutThreshold, int nCutType,
			int nClusters, double accuracy, boolean useEigenvectors, double[][] W) {
		super(sigmaI, sigmaX, r, l, nCutThreshold, nCutType, nClusters, accuracy, useEigenvectors);
		this.W = W;
	}

	@Override
	public void setSimilarityMatrix() {
		setW(MatrixUtils.createRealMatrix(W));
	}

	@Override
	public double getFeatureFactor(int i, int j) {
		return 0;
	}

	@Override
	public double getWeight(int i, int j) {
		return 0;
	}

	@Override
	public void displayClusters() {
		// Do nothing
	}

}
