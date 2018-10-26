package normCuts;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * Here is a Normalized Cuts implementation which can be used to perform the
 * algorithm on RGB images.
 * 
 * @author Andrew P Talbot
 *
 */
public class NormalizedCutsRGB extends NormalizedCuts {

	private BufferedImage img;
	private int[] pixels;
	private boolean hasAlpha;

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
	public NormalizedCutsRGB(double sigmaI, double sigmaX, double r, double l, double nCutThreshold, int nCutType,
			int nClusters, double accuracy, boolean useEigenvectors, String imgPath) {
		super(sigmaI, sigmaX, r, l, nCutThreshold, nCutType, nClusters, accuracy, useEigenvectors);

		// Read in the image.
		try {
			img = ImageIO.read(new File(imgPath));
		} catch (IOException e) {
			System.out.println("That is not a valid Image.");
			System.exit(-1);
		}

		// Set the number of pixels of the image.
		setImgWidth(img.getWidth());
		setImgHeight(img.getHeight());
		setnPixels(getImgWidth() * getImgHeight());
		setHasAlpha(this.img.getColorModel().hasAlpha());

		this.pixels = img.getRGB(0, 0, getImgWidth(), getImgHeight(), null, 0, getImgWidth());

		setSimilarityMatrix();
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
	public NormalizedCutsRGB(double sigmaI, double sigmaX, double r, double l, double nCutThreshold, int nCutType,
			int nClusters, double accuracy, boolean useEigenvectors, BufferedImage img) {
		super(sigmaI, sigmaX, r, l, nCutThreshold, nCutType, nClusters, accuracy, useEigenvectors);

		this.img = img;

		// Set the number of pixels of the image.
		setImgWidth(this.img.getWidth());
		setImgHeight(this.img.getHeight());
		setnPixels(getImgWidth() * getImgHeight());
		setHasAlpha(this.img.getColorModel().hasAlpha());

		pixels = img.getRGB(0, 0, getImgWidth(), getImgHeight(), null, 0, getImgWidth());

		setSimilarityMatrix();
	}

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

	@Override
	public double getFeatureFactor(int i, int j) {

		Color c_i = new Color(i, this.hasAlpha);
		Color c_j = new Color(j, this.hasAlpha);

		double[] c_i_arr = (this.hasAlpha)
				? new double[] { c_i.getRed(), c_i.getGreen(), c_i.getBlue(), c_i.getAlpha() }
				: new double[] { c_i.getRed(), c_i.getGreen(), c_i.getBlue() };

		double[] c_j_arr = (this.hasAlpha)
				? new double[] { c_j.getRed(), c_j.getGreen(), c_j.getBlue(), c_j.getAlpha() }
				: new double[] { c_j.getRed(), c_j.getGreen(), c_j.getBlue() };

		return Math.exp(-l2NormDiff(c_i_arr, c_j_arr) / getSigmaI());
	}

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

	public boolean getHasAlpha() {
		return this.hasAlpha;
	}

	public void setHasAlpha(boolean hasAlpha) {
		this.hasAlpha = hasAlpha;
	}

	public static void main(String[] args) {
		double sigmaI = 0.07;
		double sigmaX = 15;
		double r = 4;
		double l = 0.02;
		double accuracy = 0.0001;
		boolean useEigenvectors = true;
		double nCutThreshold = 0.0002;
		int nCutType = PARTITION_BY_ZERO;
		int nClusters = 4;

		// Big image to test on
		String imgPath = "images/color_50x50.png";
		NormalizedCutsRGB ncrgb = new NormalizedCutsRGB(sigmaI, sigmaX, r, l, nCutThreshold, nCutType, nClusters,
				accuracy, useEigenvectors, imgPath);
		show(ncrgb.getImg(), 500, 500);
		ncrgb.clusterAndSetW();

	}

}
