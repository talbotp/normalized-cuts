package normCuts;

/**
 * Here is a class that we use to implement Normalized Cuts for the use of MSI data.
 * 
 * @author Andrew P Talbot
 */
public class NormalizedCutsSpectrum extends NormalizedCuts {

	public NormalizedCutsSpectrum(double sigmaI, double sigmaX, double r, double l, double nCutThreshold, int nCutType,
			int nClusters, double accuracy, boolean useEigenvectors) {
		super(sigmaI, sigmaX, r, l, nCutThreshold, nCutType, nClusters, accuracy, useEigenvectors);
		
		
	}
	
	@Override
	public void setSimilarityMatrix() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double getFeatureFactor(int i, int j) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getWeight(int i, int j) {
		// TODO Auto-generated method stub
		return 0;
	}

}
