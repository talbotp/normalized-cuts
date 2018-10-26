package normCuts;

import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import javax.imageio.ImageIO;

import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * Normalized Cuts implementation to work on RGB Images.
 * 
 * @author Andrew P Talbot
 */
public class NormalizedCutsGray extends NormalizedCuts {

	private BufferedImage img;
	private int[] pixels;

	/**
	 * {@inheritDoc}
	 * 
	 * Constructor for when we wish to perform Normalized Cuts on an image saved
	 * on our system.
	 * 
	 * REQUIRES : We set the similarity matrix in this constructor.
	 * 
	 * @param imgPath
	 *            is the path to the grayscale image that we wish to perform
	 *            NCuts on.
	 */
	public NormalizedCutsGray(double sigmaI, double sigmaX, double r, double l, double nCutThreshold, int nCutType,
			int nClusters, double accuracy, boolean useEigenvectors, String imgPath) {
		super(sigmaI, sigmaX, r, l, nCutThreshold, nCutType, nClusters, accuracy, useEigenvectors);

		// Read in the image.
		try {
			img = ImageIO.read(new File(imgPath));
		} catch (IOException e) {
			System.out.println("That is not a valid Image.");
			System.exit(-1);
		}

		// Convert to grayscale image.
		setImg(convertToGrayBI(img));

		// Set the number of pixels of the image.
		setImgWidth(img.getWidth());
		setImgHeight(img.getHeight());
		setnPixels(getImgWidth() * getImgHeight());
		pixels = convertToIntArray(((DataBufferByte) img.getRaster().getDataBuffer()).getData());

		setSimilarityMatrix();

		// System.out.println("We print the original pixel array below.");
		// System.out.println(Arrays.toString(pixels));
	}

	/**
	 * {@inheritDoc}
	 * 
	 * Constructor for when we wish to perform Normalized Cuts on a
	 * BufferedImage.
	 * 
	 * @param img
	 *            is the BufferedImage we perform Normalized Cuts on.
	 */
	public NormalizedCutsGray(double sigmaI, double sigmaX, double r, double l, double nCutThreshold, int nCutType,
			int nClusters, double accuracy, boolean useEigenvectors, BufferedImage img) {
		super(sigmaI, sigmaX, r, l, nCutThreshold, nCutType, nClusters, accuracy, useEigenvectors);

		// Convert to grayscale image.
		this.img = convertToGrayBI(img);

		// Set the number of pixels of the image.
		setImgWidth(this.img.getWidth());
		setImgHeight(this.img.getHeight());
		setnPixels(getImgWidth() * getImgHeight());

		pixels = convertToIntArray(((DataBufferByte) img.getRaster().getDataBuffer()).getData());

		System.out.println("We print the original pixel array below.");
		System.out.println(Arrays.toString(pixels));

		setSimilarityMatrix();
	}

	/**
	 * Empty Constructor to use if we wish to create the Normalized Cuts object
	 * uisng the setters instead of just one constructor.
	 */
	public NormalizedCutsGray() {

	}

	/**
	 * {@inheritDoc}}
	 * 
	 * We set it to decimal format so that the integrity is maintained.
	 */
	@Override
	public void setSimilarityMatrix() {
		int nPixels = getnPixels();

		RealMatrix W = new OpenMapRealMatrix(nPixels, nPixels);
		
		// Loop through the array to set the symmetric matrix W.
		for (int i = 0; i < nPixels; i++)
			for (int j = i; j < nPixels; j++) {
				
				double w = getWeight(i, j);
				W.setEntry(i, j, w);
				W.setEntry(j, i, w);
				
			}

		setW(W);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public double getFeatureFactor(int i, int j) {
		return Math.exp(-Math.abs(i - j) / getSigmaI());
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public double getWeight(int i, int j) {

		// Get the pixel coordinates for point i.
		int i_x = i % getImgWidth();
		int i_y = i / getImgWidth();

		// Get the pixel coordinates for point j.
		int j_x = j % getImgWidth();
		int j_y = j / getImgWidth();

		// Work out the spatial value X(i) - X(j).
		double spatialDistance = getSpatialWeight(new double[] { i_x, i_y }, new double[] { j_x, j_y });

		return (spatialDistance < getR())
				? getFeatureFactor(getPixels()[i], getPixels()[j]) * getSpatialFactor(spatialDistance) : 0;
	}

	public BufferedImage getImg() {
		return img;
	}

	public void setImg(BufferedImage img) {
		this.img = img;
	}

	public int[] getPixels() {
		return pixels;
	}

	public void setPixels(int[] pixels) {
		this.pixels = pixels;
	}

	public static void main(String[] args) {
		System.out.println("BEGIN");
		
		double sigmaI = 15;
		double sigmaX = 16;
		double r = 2.2;
		double accuracy = 0.0001;
		boolean useEigenvectors = true;
		double l = 0.02;
		double nCutThreshold = 0.05;
		int nCutType = PARTITION_BY_ZERO;
		int nClusters = 3;

		// Big image to test on
		String imgPath = "images/test_30x30.png";

		// Small image to test on.
		int[] pixels = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
		BufferedImage img = new BufferedImage(3, 3, BufferedImage.TYPE_BYTE_GRAY);
		for (int i = 0; i < pixels.length; i++)
			img.getRaster().getDataBuffer().setElem(i, pixels[i]);

		NormalizedCutsGray nc = new NormalizedCutsGray(sigmaI, sigmaX, r, l, nCutThreshold, nCutType, nClusters,
				accuracy, useEigenvectors, imgPath);
		show(nc.getImg(), 300, 300);

		nc.clusterAndSetW();
		System.out.println("END");
	}

}
