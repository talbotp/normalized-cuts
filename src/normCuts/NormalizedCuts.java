package normCuts;

import java.awt.FlowLayout;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;

import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import normCuts.LanczosSolver.EigenValueVectorPair;

/**
 * Class which defines all generic methods to be used for any Normalized Cuts
 * implementations.
 * 
 * @author Andrew P Talbot
 */
public abstract class NormalizedCuts {

	// Types of methods to choose the partitions.
	public static final int PARTITION_BY_ZERO = 0;
	public static final int PARTITION_BY_MEDIAN = 1;
	public static final int PARTITION_BY_NCUT = 2;

	// Parameters of algorithm
	private double sigmaI, sigmaX, r, l, nCutThreshold;

	// Fields specific to the image.
	private int imgWidth, imgHeight, nPixels, nCutType, nClusters;
	
	private boolean useEigenvectors;
	private double accuracy;

	// Matrices created for the algorithm.
	private RealMatrix W;

	// Array which stores the value of all of the clusters.
	private int[] clusters;
	private int currentNClusters;
	// THIS FIELD IS NOT THE CURRENTNCLUSTERS etc.
	private int tmpNClusters;

	private BufferedImage img;

	/**
	 * Constructor for the Normalized Cuts image segmentation.
	 * 
	 * REQUIRE : Please set the Similarity matrix W in the implemented
	 * constructor.
	 * 
	 * @param sigmaI
	 *            is the value which divides the features.
	 * @param sigmaX
	 *            is the value which divides the spatial value.
	 * @param r
	 *            is the length which if the distance between two nodes is
	 *            greater than or equal to, then we award a weight of zero.
	 * @param l
	 *            is the gap size used to estimate the best n cut partition.
	 * @param nCutThreshold
	 *            is the threshold that we use to consider repartitioning and
	 *            recursively calling the algorithm.
	 * @param nCutType
	 *            is the method that we use to repartition our graph.
	 * @param nClusters
	 *            is the number of clusters, is only used when the nCutType is
	 *            ZERO or MEDIAN.
	 */
	public NormalizedCuts(double sigmaI, double sigmaX, double r, double l, double nCutThreshold, int nCutType,
			int nClusters, double accuracy, boolean useEigenvectors) {
		setSigmaI(sigmaI);
		setSigmaX(sigmaX);
		setR(r);
		setL(l);
		setNCutThreshold(nCutThreshold);
		setNCutType(nCutType);
		setNClusters(nClusters);
		setAccuracy(accuracy);
		setUseEigenvectors(useEigenvectors);
		this.currentNClusters = 0;
	}

	/**
	 * This is NOT the recommended constructor to use for NCuts. Please use the
	 * only other NCuts Constructor.
	 * 
	 * Empty constructor if the user wishes to make a NCuts implementation from
	 * scratch using setters.
	 */
	public NormalizedCuts() {

	}

	/**
	 * Method we call to actually perform clustering on our Normalized Cuts
	 * object, to get the result, get the clusters integer array which contains
	 * a pixel array containing the cluster value for a given pixel.
	 */
	public void clusterAndSetW() {

		setSimilarityMatrix();

		initClusterArray();

		recursiveCut(new BiPartition(this.W, null));
	}

	/**
	 * Method that we call to cluster our data, without setting the similarity
	 * matrix W.
	 * 
	 * NB : Please Set the similarity matrix before calling this method.
	 */
	public void cluster() {

		initClusterArray();

		recursiveCut(new BiPartition(this.W, null));
	}

	/*********************************************************************/
	/* Here are methods to be implemented by extending classes. */
	/*********************************************************************/

	/**
	 * We call this method to take the input data, create a graph of nodes and
	 * edges, then create and set the Weights matrix W, calling setW().
	 */
	public abstract void setSimilarityMatrix();

	/**
	 * Returns the feature factor. As seen in the NCuts paper, it is generalized
	 * as l2Norm((F(i) - F(j))^2 / sigmaX
	 * 
	 * @param i
	 *            is the index of the pixel array for node i.
	 * @param j
	 *            is the index of the pixel array for node j.
	 * @return the featured weight between the two nodes i and j.
	 */
	public abstract double getFeatureFactor(int i, int j);

	/**
	 * Return the weight of an edge between the index i and index j nodes in the
	 * pixel array.
	 * 
	 * Typically, this method is to be called in the setSimilarityMatrix method.
	 * 
	 * @param i
	 *            is the index of the pixel array for node i.
	 * @param j
	 *            is the index of the pixel array for node j.
	 * @return the weight between two nodes with index i and j in the pixel
	 *         array.
	 */
	public abstract double getWeight(int i, int j);

	/*********************************************************************/
	/* Here are the recursive utility methods for the clustering. */
	/*********************************************************************/

	/**
	 * This is the recursive method that we call to recursively partition and
	 * cluster our image, when the settings for the nCutType = PARTITION_BY_ZERO
	 * or PARTITION_BY_MEDIAN.
	 * 
	 * @param W
	 *            is the similarity matrix for this iteration in the algorithm.
	 * @param nClusters
	 *            is the maximum number of clusters, we use this so we know when
	 *            to stop partitioning.
	 * @param originalIndices
	 *            is the mapping from the indexes of the new weight matrix to
	 *            their original indices in the first similarity matrix W, so it
	 *            can properly be mapped to the right pixel.
	 * @param originalIndices2
	 *            is the same mapping as originalIndices but for the partition
	 *            number 2.
	 */
	public void recursiveCut(BiPartition... partitions) {

		// Make the new array for the recursive call.
		ArrayList<BiPartition> bipartitions = new ArrayList<>();

		// For each leaf in the recursive call we do cluster.
		for (BiPartition leaf : partitions) {

			System.out.println("Size of bipartition is " + bipartitions.size());

			if (leaf.W.getColumnDimension() <= 2)
				continue;

			RealMatrix D = getD(leaf.W);
			RealMatrix D_1_2 = getD_minus_1_2(D);
			RealMatrix toEig = getMatrixToEigenDecompose(D, D_1_2, leaf.W);

			// IF too big then we use Lanczos. - CHANGE TO TRUE. TODO : CHANGE
			int val = 30;
			boolean lessThanVal = (leaf.W.getColumnDimension() < val) ? true : false;
			double[] eigvec2 = getSecondEigenvector(toEig, lessThanVal, this.getAccuracy(), this.getUseEigenvectors());

			double splitter = getSplitter(eigvec2, leaf.W);
			System.out.println("splitter = " + splitter);

			// Make the new maps for the new output.
			Map<Integer, Integer> p1 = new HashMap<>();
			Map<Integer, Integer> p2 = new HashMap<>();

			// Invalid splitter value
			System.out.println(splitter > this.nCutThreshold);
			if (this.nCutType == PARTITION_BY_NCUT && splitter > this.nCutThreshold) {
				System.out.println("End the ncuts partition");
				continue;
			}

			if (leaf.partition == null) { // First iteration of method.
				for (int i = 0; i < eigvec2.length; i++) {
					if (eigvec2[i] < splitter) {
						clusters[i] = this.currentNClusters + 1;
						p1.put(i, i);
					} else {
						clusters[i] = this.currentNClusters + 2;
						p2.put(i, i);
					}
				}

				this.currentNClusters = 2;
				this.tmpNClusters = 2;

			} else {

				// We only do this if type is ncut and (the type is not ncut and
				// the number of clusters is less than the max n clusters.)
				if ((getNCutType() != PARTITION_BY_NCUT && getCurrentNClusters() < getNClusters())
						|| getNCutType() == PARTITION_BY_NCUT) {

					int sizeP1 = 0;
					int sizeP2 = 0;
					Set<Integer> keySet = leaf.partition.keySet();
					int counter = 0;
					for (Integer i : keySet) {
						if (eigvec2[counter] < splitter) {
							clusters[leaf.partition.get(i)] = this.tmpNClusters + 1;
							p1.put(sizeP1, leaf.partition.get(i));
							sizeP1++;
						} else {
							clusters[leaf.partition.get(i)] = this.tmpNClusters + 2;
							p2.put(sizeP2, leaf.partition.get(i));
							sizeP2++;
						}
						counter++;
					}

					// Increment the number of clusters
					this.currentNClusters += 1;
					this.tmpNClusters += 2;

				} else {
					// TODO : ISSUE HERE
					break;
				}
			}

			// We have now created the two new partitions for this, now we
			// decide whether to add it to the new array for the recursive call
			// or not.

			// Make the image of the clusters.
			displayClusters();

			if (this.nCutType == PARTITION_BY_NCUT) {
				// if (splitter <= this.nCutThreshold) {
				if (p1.size() != 0)
					bipartitions.add(new BiPartition(getNewW(leaf.W, p1), p1));
				if (p2.size() != 0)
					bipartitions.add(new BiPartition(getNewW(leaf.W, p2), p2));
				// }
			} else {
				if (this.currentNClusters >= nClusters) {
					return;
				} else {
					if (p1.size() != 0)
						bipartitions.add(new BiPartition(getNewW(leaf.W, p1), p1));
					if (p2.size() != 0)
						bipartitions.add(new BiPartition(getNewW(leaf.W, p2), p2));
				}
			}
		}

		// Make the recursive call
		if (bipartitions.size() != 0)
			recursiveCut(bipartitions.toArray(new BiPartition[0]));
	}

	/*********************************************************************/
	/* Here are the recursive utility methods for the clustering. */
	/*********************************************************************/

	/**
	 * Get the appropriate splitter according to the partitioning method
	 * specified by the user in the field nCutType.
	 * 
	 * @param eigvec2
	 *            is the vector that we partition.
	 * @param W
	 *            is the weights matrix we use to partition.
	 * @return the value that we use to partition the eigenvector.
	 */
	public double getSplitter(double[] eigvec2, RealMatrix W) {
		switch (this.nCutType) {
		case PARTITION_BY_ZERO:
			return 0;
		case PARTITION_BY_MEDIAN:
			return getMedian(eigvec2);
		case PARTITION_BY_NCUT:
			return getBestCut(eigvec2, getL(), W);
		default:
			// This should NEVER happen.
			throw new IllegalStateException("An invalid partitioning method is being used.");
		}
	}

	/**
	 * Calculate the nth cut and the best partition.
	 * 
	 * @param eigvec2
	 *            is the second eigenvector of the weights matrix W for this
	 *            iteration of the algorithm.
	 * @param l
	 *            is the size variable of the gaps used to estimate the best cut
	 *            value.
	 * @param W
	 *            is the weight matrix for the current iteration of the
	 *            algorithm.
	 */
	public static double getBestCut(double[] eigvec2, double l, RealMatrix W) {

		// Find the min of the second eigenvector, so we know the intensity to
		// start partitioning from.
		double min = Double.MAX_VALUE;
		for (int i = 0; i < eigvec2.length; i++)
			if (eigvec2[i] < min)
				min = eigvec2[i];

		ArrayList<Double> nCutValues = new ArrayList<>();

		boolean validCutVal = true;
		double splitter = min;

		/*
		 * We fill the array p2 with all indices of pixels with an intensity of
		 * greater than the splitter. As we start from below (min), we stop the
		 * loop when the element finally has nothing in it.
		 */
		while (validCutVal) {

			ArrayList<Integer> p1 = new ArrayList<>();
			ArrayList<Integer> p2 = new ArrayList<>();

			for (int i = 0; i < eigvec2.length; i++) {
				double tmp = eigvec2[i];
				if (tmp < splitter)
					p1.add(i);
				else
					p2.add(i);
			}

			double cutVal = NCut(p1.toArray(new Integer[0]), p2.toArray(new Integer[0]), W);
			nCutValues.add(cutVal);

			/*
			 * Here, we alter the loop invariant, and increment the value of
			 * splitter by l.
			 */
			if (p2.size() == 0)
				validCutVal = false;

			splitter += l;
		}

		// Set the bestCut value
		Collections.sort(nCutValues);

		return nCutValues.get(0);
	}

	/**
	 * We use this to create a new similarity matrix based off of the new
	 * partition and the old matrix W.
	 * 
	 * @param W
	 *            is the original similarity matrix.
	 * @param partition
	 *            is the values in the new partition.
	 * @return the new similarity matrix.
	 */
	public static RealMatrix getNewW(RealMatrix W, Map<Integer, Integer> partition) {
		Integer[] p = partition.keySet().toArray(new Integer[0]);
		int sizeW = p.length;

		if (sizeW > W.getColumnDimension())
			throw new IllegalArgumentException("That is an invalid map size.");

		RealMatrix newW = new OpenMapRealMatrix(sizeW, sizeW);

		// Add the cut values to the new matrix.
		for (int i = 0; i < sizeW; i++) {
			for (int j = 0; j < sizeW; j++) {
				newW.setEntry(i, j, W.getEntry(p[i], p[j]));
			}
		}

		return newW;
	}

	/**
	 * initialize the cluster array, must be done once we get the number of
	 * pixels.
	 */
	public void initClusterArray() {
		clusters = new int[nPixels];
	}

	/**
	 * Returns the spatial factor for a given distance between two nodes in a
	 * graph.
	 * 
	 * @param spatialDistance
	 *            is the distance between the two nodes in the graph.
	 * @return the factor that is multiplied to calculate their weight.
	 */
	public double getSpatialFactor(double spatialDistance) {
		return Math.exp(-Math.pow(spatialDistance, 2) / this.sigmaX);
	}

	/**
	 * Get the euclidean distance between two vectors (arrays).
	 * 
	 * @param i
	 *            is the array of the location of pixel i.
	 * @param j
	 *            is the array of the location of pixel j.
	 * @return the euclidean distance between these two arrays.
	 */
	public static double getSpatialWeight(double[] i, double[] j) {
		return l2NormDiff(i, j);
	}

	/**
	 * Calculate the Matrix D^(-1/2).
	 * 
	 * REQUIRES : D is a diagonal matrix.
	 * 
	 * @param D
	 *            is the diagonal matrix we calculate on.
	 * @return D^(-1/2)
	 */
	public static RealMatrix getD_minus_1_2(RealMatrix D) {
		int nPixels = D.getColumnDimension();

		RealMatrix D_minus = new OpenMapRealMatrix(nPixels, nPixels);

		for (int i = 0; i < nPixels; i++)
			D_minus.setEntry(i, i, 1d / Math.sqrt(D.getEntry(i, i)));

		return D_minus;
	}

	/**
	 * Get the second eigenvector for the Matrix X. Typically this matrix is
	 * D^(-1/2) * (D - W) * D^(-1/2)
	 * 
	 * @param X
	 *            is the matrix we get the second eigenvector for.
	 * @param trivialMethod
	 *            true if we trivially calculate the second eigenvector,
	 *            requires us to calculate all eigenvectors, else we calulate it
	 *            efficiently, still need to figure this out.
	 * @return the eigenvector corresponding to the second smallest eigenvalue.
	 */
	public static double[] getSecondEigenvector(RealMatrix X, boolean trivialMethod, double accuracy, boolean useEigenvectors) {

		if (trivialMethod) {
			// Use library to calculate all eigenvectors, get the second.
			EigenDecomposition eigDec = new EigenDecomposition(X);
			double[] eigenvalues = eigDec.getRealEigenvalues();

			System.out.println(Arrays.toString(eigenvalues));

			ArrayList<EigenvalueIndexPair> eigPairs = new ArrayList<>();
			for (int i = 0; i < eigenvalues.length; i++)
				eigPairs.add(new EigenvalueIndexPair(i, eigenvalues[i]));
			Collections.sort(eigPairs);
			int index_eigval2 = eigPairs.get(1).index;
			return eigDec.getEigenvector(index_eigval2).toArray();
		} else {
			// Use Lanczos to get only the second eigenvector.

			// Here are the Lanczos variables, perhaps these should be static
			// variables at the top of this class?
			int eigIndex = 1;
			int m = X.getColumnDimension(); // TODO : CHANGE THIS TO ALTER THE
											// ACCURACY OF THE APPROXIMATION.

			LanczosSolver ls = new LanczosSolver(X, eigIndex, m, useEigenvectors, accuracy);
			EigenValueVectorPair evp = ls.solve();
			System.out.println(evp.eigenvalue);

			return evp.eigenvector;
		}
	}

	/**
	 * Method to check if a matrix is symmetric.
	 * 
	 * We use this as a asanity check that the matrix we are checking is in fact
	 * symmetric, otherwise Lanczos will simply fail.
	 * 
	 * @param X
	 *            is the matrix we check if is symmetric.
	 * @param accuracy
	 *            is the maximum difference between values for the matrix to be
	 *            considered symmetric, we do this if we wish to account for
	 *            rounding errors.
	 * @return true if the matrix is symmetric, false otherwise.
	 */
	public static boolean isSymmetric(RealMatrix X, double accuracy) {

		// Case when matrix is not square it cannot be symmetric.
		if (X.getColumnDimension() != X.getRowDimension()) {
			System.out.println("unequal size row and columns, not symmetric");
			return false;
		}

		// Work through the upper triangle of the matrix.
		for (int i = 0; i < X.getColumnDimension(); i++) {
			for (int j = i; j < X.getRowDimension(); j++) {
				if (Math.abs(X.getEntry(i, j) - X.getEntry(j, i)) > accuracy)
					return false;
			}
		}
		return true;
	}

	/**
	 * Gets the matrix X.
	 * 
	 * Method which calculates the matrix which we are going to eigen decompose.
	 * 
	 * We calculate : D^(-1/2) * (D - W) * D^(-1/2);
	 * 
	 * @param D
	 *            is the diagonal matrix with the sum of the rows of W in the
	 *            diagonal.
	 * @param D_1_2
	 *            is D^(-1/2)
	 * @param W
	 *            is the original weights matrix W.
	 */
	public static RealMatrix getMatrixToEigenDecompose(RealMatrix D, RealMatrix D_1_2, RealMatrix W) {

		RealMatrix X = D.subtract(W);

		X = X.multiply(D_1_2);

		X = D_1_2.multiply(X);

		return X;
	}

	/**
	 * Get the diagonal matrix D, which is a diagonal matrix with the sum of the
	 * rows of W as its columns.
	 * 
	 * @param W
	 *            is the weight matrix W.
	 * @return the matrix D.
	 */
	public static RealMatrix getD(RealMatrix W) {
		int nPixels = W.getColumnDimension();

		RealMatrix D = new OpenMapRealMatrix(nPixels, nPixels);

		for (int i = 0; i < nPixels; i++) {
			D.setEntry(i, i, sumRow(W, i));
		}

		return D;
	}

	/*********************************************************************/
	/* Here are methods we use to calculate NCut and Nassoc etc. */
	/*********************************************************************/

	/**
	 * Get the normalized cut value with partitions p1 and p2 on the weights
	 * matrix W.
	 * 
	 * Requires : 1) the union of p1 and p2 should contain be exactly equal to
	 * the indexes of the number of elements in W. 2) W must be a square matrix.
	 * 
	 * @param p1
	 *            is a partition on graph W
	 * @param p2
	 *            is a partition on graph W, such that p1 and p2 have a pixel
	 *            value for each pixel in the matrix W.
	 * @param W
	 *            is weight matrix W.
	 * @return NCut value.
	 */
	public static double NCut(Integer[] p1, Integer[] p2, RealMatrix W) {
		// As we get NaN in these cases, and we wish to minimise the NCut
		// values, we simply set these to be the maximum double value so that
		// these partitioning values will never be chosen.
		if (p1.length == 0 || p2.length == 0)
			return Double.MAX_VALUE;
		return 2 - Nassoc(p1, p2, W);
	}

	/**
	 * Return Nassoc of p1, p2 on weight matrix W.
	 * 
	 * REQUIRES : p1 and p2 non empty.
	 * 
	 * SEE DOCS for NCut().
	 */
	public static double Nassoc(Integer[] p1, Integer[] p2, RealMatrix W) {
		return assoc(p1, false, W) / assoc(p1, true, W) + assoc(p2, false, W) / assoc(p2, true, W);
	}

	/**
	 * Returns asso of p1 and p2 on weight matrix W.
	 * 
	 * REQUIRES : p1 and p2 non empty.
	 * 
	 * SEE DOCS for NCut().
	 */
	public static double assoc(Integer[] p1, boolean V, RealMatrix W) {
		double sumW = 0;

		if (V) {
			// Sum of all weights

			for (int i = 0; i < p1.length; i++) {
				for (int j = 0; j < W.getColumnDimension(); j++) {
					sumW += W.getEntry(p1[i], j);
				}
			}

		} else {

			// Sum of weights between p1 and W
			for (int i = 0; i < p1.length; i++) {
				for (int j = 0; j < p1.length; j++) {
					sumW += W.getEntry(p1[i], p1[j]);
				}
			}
		}

		return sumW;
	}

	/*********************************************************************/
	/* Here are some utility methods */
	/*********************************************************************/

	/**
	 * Add the values of a column of a RealMatrix.
	 * 
	 * @param A
	 *            is the matrix we add the values in column col of.
	 * @param col
	 *            is the column of A we wish to sum.
	 * @return the sum of the values in column col of matrix A.
	 */
	public static double sumColumn(RealMatrix A, int col) {
		double sumCol = 0;
		for (int i = 0; i < A.getRowDimension(); i++)
			sumCol += A.getEntry(i, col);
		return sumCol;
	}

	/**
	 * Add the values of a row of a RealMatrix.
	 * 
	 * @param A
	 *            is the matrix we add the values in row col of.
	 * @param row
	 *            is the row of A we wish to sum.
	 * @return the sum of the values in column row of matrix A.
	 */
	public static double sumRow(RealMatrix A, int row) {
		double sumRow = 0;
		for (int i = 0; i < A.getColumnDimension(); i++)
			sumRow += A.getEntry(row, i);
		return sumRow;
	}

	/**
	 * Method to visualize a BufferedImage.
	 */
	public static void show(BufferedImage img) {
		JFrame frame = new JFrame();
		frame.getContentPane().setLayout(new FlowLayout());
		frame.getContentPane().add(new JLabel(new ImageIcon(img)));
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}

	/**
	 * Method to visualize a resized BufferedImage.
	 */
	public static void show(BufferedImage img, int newWidth, int newHeight) {
		Image tmp = img.getScaledInstance(newWidth, newHeight, Image.SCALE_SMOOTH);
		BufferedImage dimg = new BufferedImage(newWidth, newHeight, BufferedImage.TYPE_INT_ARGB);
		Graphics2D g2d = dimg.createGraphics();
		g2d.drawImage(tmp, 0, 0, null);
		g2d.dispose();

		JFrame frame = new JFrame();
		frame.getContentPane().setLayout(new FlowLayout());
		frame.getContentPane().add(new JLabel(new ImageIcon(dimg)));
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}

	/**
	 * Create a BufferedImage with pixel array pixels, create a frame for it and
	 * display it.
	 * 
	 * @param pixels
	 *            is the pixel array to visualize.
	 */
	public void show(int[] pixels) {
		BufferedImage img = new BufferedImage(imgWidth, imgHeight, BufferedImage.TYPE_BYTE_GRAY);
		for (int i = 0; i < imgWidth; i++) {
			for (int j = 0; j < imgHeight; j++) {
				img.setRGB(i, j, pixels[i * imgWidth + j]);
			}
		}
		show(img);
	}

	/**
	 * Method we use to make a BufferedImage of our clusters and display it in a
	 * 500x500 pixels frame.
	 */
	public void displayClusters() {
		BufferedImage tmp_img = new BufferedImage(getImgWidth(), getImgHeight(), BufferedImage.TYPE_BYTE_GRAY);

		setImgPixelsGrayscale(tmp_img, getClusterPixel(getclusters(true), getCurrentNClusters()));
		System.out.println("Now we show the new image with " + getCurrentNClusters() + " clusters.");
		System.out.println("But there are " + this.getNClusters() + "clusters");
		show(tmp_img, 300, 300);
	}

	/**
	 * Method for setting the values of a BufferedImage to match the values of
	 * that in an RGB array.
	 * 
	 * @param img
	 *            is the images we are to set the value of pixels for.
	 * @param pixels
	 *            is the new RGB values.
	 */
	public static BufferedImage setImgPixelsGrayscale(BufferedImage img, int[] pixels) {
		for (int i = 0; i < pixels.length; i++)
			img.getRaster().getDataBuffer().setElem(i, pixels[i] - 128);
		return img;
	}

	/**
	 * Convert the pixels array to an array which can be visualized by a
	 * grayscale BufferedImage
	 * 
	 * @param clusters
	 *            is the cluster array we wish to view.
	 * @param nClusters
	 *            is the number of clusters there are, we do this so that we
	 *            know how to scale the pixel values.
	 * @return the pixel array for clusters.
	 */
	public static int[] getClusterPixel(int[] clusters, int nClusters) {
		int[] pixels = new int[clusters.length];
		int gap = (int) Math.ceil(256d / nClusters);
		for (int i = 0; i < clusters.length; i++) {
			int clusterVal = (clusters[i] - 1) * gap;
			pixels[i] = clusterVal;
		}
		return pixels;
	}

	/**
	 * Convert a BufferedImage to a grayscale BufferedImage.
	 * 
	 * @param biColor
	 *            is the color BufferedImage.
	 * @return the grayscale BufferedImage version of biColor.
	 */
	public static BufferedImage convertToGrayBI(BufferedImage biColor) {
		BufferedImage biGray = new BufferedImage(biColor.getWidth(), biColor.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
		biGray.getGraphics().drawImage(biColor, 0, 0, null);
		return biGray;
	}

	/**
	 * Convert an array from bytes to int.
	 * 
	 * @param pixels
	 *            is the byte array to convert to int.
	 * @return an int array of the byte pixels.
	 */
	public static int[] convertToIntArray(byte[] pixels) {
		int[] arr = new int[pixels.length];
		for (int i = 0; i < pixels.length; i++) {
			arr[i] = 128 + pixels[i];
		}
		return arr;
	}

	// TODO : REPLACE WITH MATHMSK VERSION.
	public static double l2Norm(double[] spec) {

		if (spec.length < 1)
			throw new IllegalArgumentException("The array requires a length of at least 1.");

		double sum = 0;
		for (double val : spec) {
			sum += Math.pow(Math.abs(val), 2);
		}
		return Math.sqrt(sum);
	}

	// TODO : REPLACE WITH MATHMSK VERSION.
	public static double l2NormDiff(double[] spec1, double[] spec2) {
		int length = spec1.length;
		double[] arr = new double[spec1.length];
		for (int i = 0; i < length; ++i) {
			arr[i] = spec1[i] - spec2[i];
		}
		return l2Norm(arr);
	}

	// TODO: Replace with MathMSK.median().
	public static double getMedian(double[] arr) {
		double[] tmp = arr.clone();
		Arrays.sort(tmp);
		double median;
		if (tmp.length % 2 == 0)
			median = ((double) tmp[tmp.length / 2] + (double) tmp[tmp.length / 2 - 1]) / 2;
		else
			median = (double) tmp[tmp.length / 2];
		return median;
	}

	/*********************************************************************/
	/* Here are getters and setters. */
	/*********************************************************************/

	public double getSigmaI() {
		return sigmaI;
	}

	public void setSigmaI(double sigmaI) {
		if (sigmaI > 0)
			this.sigmaI = sigmaI;// Math.pow(sigmaI, 2);
		else
			throw new IllegalArgumentException("Please use a positive value of sigmaI.");
	}

	public double getSigmaX() {
		return sigmaX;
	}

	public void setSigmaX(double sigmaX) {
		if (sigmaX > 0)
			this.sigmaX = sigmaX;// Math.pow(sigmaX, 2);
		else
			throw new IllegalArgumentException("Please use a positive value of sigmaX.");
	}

	public double getR() {
		return r;
	}

	public void setR(double r) {
		if (r > 0)
			this.r = r;
		else
			throw new IllegalArgumentException("Please use a positive value of r.");
	}

	public double getNCutThreshold() {
		return nCutThreshold;
	}

	public void setNCutThreshold(double nCutThreshold) {
		if (nCutThreshold > 0)
			this.nCutThreshold = nCutThreshold;
		else
			throw new IllegalArgumentException("Please use a positive value of the N Cut threshold.");
	}

	public int getImgWidth() {
		return imgWidth;
	}

	public void setImgWidth(int imgWidth) {
		this.imgWidth = imgWidth;
	}

	public int getImgHeight() {
		return imgHeight;
	}

	public void setImgHeight(int imgHeight) {
		this.imgHeight = imgHeight;
	}

	public int getnPixels() {
		return nPixels;
	}

	public void setnPixels(int nPixels) {
		this.nPixels = nPixels;
	}

	public RealMatrix getW() {
		return W;
	}

	public void setW(RealMatrix w) {
		W = w;
	}

	public double getL() {
		return this.l;
	}

	public void setL(double l2) {
		if (l2 >= 0)
			this.l = l2;
		else
			throw new IllegalArgumentException("Please ensure that l is non-negative.");
	}

	public int getNCutType() {
		return nCutType;
	}

	public void setNCutType(final int TYPE) {
		if (TYPE == PARTITION_BY_ZERO || TYPE == PARTITION_BY_MEDIAN || TYPE == PARTITION_BY_NCUT)
			nCutType = TYPE;
		else
			throw new IllegalArgumentException("Please choose a valid method of partitioning the graph.");
	}

	/**
	 * Method to clean the cluster array, because it in order to correctly get
	 * the cluster values, we
	 * 
	 * @param cleanClusters
	 * @return the clean clusters with values ranging from 0 <= val <=
	 *         getCurrentNClusters().
	 */
	public int[] getclusters(boolean cleanClusters) {
		if (cleanClusters) {
			int[] c = clusters.clone();
			int min = Integer.MAX_VALUE;
			for (int i = 0; i < c.length; i++) {
				if (c[i] < min)
					min = c[i];

			}

			for (int i = 0; i < c.length; i++) {
				c[i] -= min;
			}
			return c;
		} else
			return clusters;
	}

	public int getNClusters() {
		return nClusters;
	}

	public int getCurrentNClusters() {
		return this.currentNClusters;
	}

	public void setNClusters(int nClusters) {
		if (nClusters > 0)
			this.nClusters = nClusters;
		else
			throw new IllegalArgumentException("Please ensure the number of clusters is greater than 0.");
	}

	public BufferedImage getImg() {
		return this.img;
	}

	public void setImg(BufferedImage img) {
		this.img = img;
	}

	public double getAccuracy() {
		return accuracy;
	}
	
	public void setAccuracy(double accuracy) {
		if (accuracy < 0)
			throw new IllegalArgumentException("That is an invalid value of the accuracy.");
		else
			this.accuracy = accuracy;
	}
	
	public boolean getUseEigenvectors() {
		return useEigenvectors;
	}
	
	public void setUseEigenvectors(boolean toUseEigenvectors) {
		this.useEigenvectors = toUseEigenvectors;
	}

}
