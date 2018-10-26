package normCuts;

import static org.junit.Assert.assertEquals;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import normCuts.NormalizedCuts;
import normCuts.NormalizedCutsGray;

/**
 * Test cases for NormalizedCutsGrey, NB We extend the abstract class
 * NormalizedCuts and hence we MUST also test the methods in the class.
 * 
 * NB : To test some of the methods, we create an empty NormalizedCuts object,
 * and then set the imgPixels object and use that, it seems to be the easiest
 * way to set the imgPixels array to something we can use sutiably for testing.
 * 
 * @author Andrew P Talbot
 */
public class NormalizedCutsGreyTest {

	// Here are simple similarity matrices to be used to test the Normalized
	// Cuts Methods.
	public static final int[] DATA_1 = { 1 };
	public static final int[] DATA_2 = { 1, 2, 3, 4 };
	public static final int[] DATA_3 = { 1, 0, 0, 2 };
	public static final int[] DATA_4 = { 1, 0, 1, 1, 1, 2, 0, 0, 1 };

	private static double sigmaI = 2;
	private static double sigmaX = 5;
	private static double r = 8;
	private static int l = 2;
	private static double nCutThreshold = 2;
	private static int nClusters = 1;
	private static double accuracy = 0.001;
	private static boolean useEigenvectors = true;
	private static int nCutType = NormalizedCuts.PARTITION_BY_ZERO;
	private static String imgPath = "images/cat_grayscale_30x30.gif";

	private static NormalizedCutsGray nc, nc1, nc2, nc3, nc4;

	@Rule
	public final ExpectedException exception = ExpectedException.none();

	/*
	 * Construct the Normalized Cuts objects that we test on.
	 */
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		nc = new NormalizedCutsGray(sigmaI, sigmaX, r, l, nCutThreshold, nCutType, nClusters, accuracy, useEigenvectors,
				imgPath);

		nc1 = new NormalizedCutsGray();
		nc1.setPixels(DATA_1);
		nc1.setImgWidth(1);
		nc1.setImgHeight(1);
		nc1.setnPixels(1);
		nc1.setSigmaI(sigmaI);
		nc1.setSigmaX(sigmaX);
		nc1.setR(r);

		nc2 = new NormalizedCutsGray();
		nc2.setPixels(DATA_2);
		nc2.setImgWidth(2);
		nc2.setImgHeight(2);
		nc2.setSigmaI(sigmaI);
		nc2.setSigmaX(sigmaX);
		nc2.setR(r);

		nc3 = new NormalizedCutsGray();
		nc3.setPixels(DATA_3);
		nc3.setImgWidth(2);
		nc3.setImgHeight(2);
		nc3.setnPixels(4);
		nc3.setSigmaI(sigmaI);
		nc3.setSigmaX(sigmaX);
		nc3.setR(r);

		nc4 = new NormalizedCutsGray();
		nc4.setPixels(DATA_4);
		nc4.setImgWidth(3);
		nc4.setImgHeight(3);
		nc4.setnPixels(9);
		nc4.setSigmaI(sigmaI);
		nc4.setSigmaX(sigmaX);
		nc4.setR(r);
	}

	// Test constructor getters, setters.

	@Test
	public void test_1_1() {
		double expectedI = sigmaI;
		double actualI = nc.getSigmaI();
		assertEquals(expectedI, actualI, 0);

		double expectedX = sigmaX;
		double actualX = nc.getSigmaX();
		assertEquals(expectedX, actualX, 0);

		double expectedR = r;
		double actualR = nc.getR();
		assertEquals(expectedR, actualR, 0);

		int expectedL = l;
		double actualL = nc.getL();
		assertEquals(expectedL, actualL, 0);

		double expectedThreshold = nCutThreshold;
		double actualThreshold = nc.getNCutThreshold();
		assertEquals(expectedThreshold, actualThreshold, 0);

		int expectedNClusters = nClusters;
		int actualNClusters = nc.getNClusters();
		assertEquals(expectedNClusters, actualNClusters);

		int expectedCutType = nCutType;
		int actualCutType = nc.getNCutType();
		assertEquals(expectedCutType, actualCutType);
	}

	// Test the method getFeatureFactor(). This method is what creates the
	// feature factor when calculating the weight between two nodes in the
	// graph. It should always return e ^ (- abs(i-j)).

	@Test
	public void test_2_1() {
		int i = 0;
		int j = 0;

		double expected = 1;
		double actual = nc1.getFeatureFactor(i, j);
		assertEquals(expected, actual, 0);
	}

	@Test
	public void test_2_2() {
		int i = 1;
		int j = 0;

		double expected = Math.exp(-1 / sigmaI);
		double actual = nc1.getFeatureFactor(i, j);
		assertEquals(expected, actual, 0);

		double actualSwap = nc1.getFeatureFactor(j, i);
		assertEquals(expected, actualSwap, 0);
	}

	@Test
	public void test_2_3() {
		int i = 1;
		int j = 3;

		double expected = Math.exp(-2 / sigmaI);
		double actual = nc1.getFeatureFactor(i, j);
		assertEquals(expected, actual, 0);

		double actualSwap = nc1.getFeatureFactor(j, i);
		assertEquals(expected, actualSwap, 0);
	}

	// Test method getWeight().
	// This method takes on input pixels indices and gets their weight in the
	// similarity matrix W.

	@Test
	public void test_3_1() {
		int i = 0;
		int j = 0;

		double expected = 1;
		double actual = nc1.getWeight(i, j);
		assertEquals(expected, actual, 0);

		double actualSwap = nc1.getWeight(j, i);
		assertEquals(expected, actualSwap, 0);
	}

	@Test
	public void test_3_2() {
		int i = 0;
		int j = 0;

		double expected = 1;
		double actual = nc2.getWeight(i, j);
		assertEquals(expected, actual, 0);

		double actualSwap = nc2.getWeight(j, i);
		assertEquals(expected, actualSwap, 0);
	}

	@Test
	public void test_3_3() {
		int i = 1;
		int j = 0;

		double expected = Math.exp(-Math.abs(1 - 2) / sigmaI) * Math.exp(-Math.abs(1 / sigmaX));
		double actual = nc2.getWeight(i, j);
		assertEquals(expected, actual, 0);

		double actualSwap = nc2.getWeight(j, i);
		assertEquals(expected, actualSwap, 0);
	}

	@Test
	public void test_3_4() {
		int i = 0;
		int j = 1;

		double expected = Math.exp(-Math.abs(0 - 1) / sigmaI) * Math.exp(-Math.abs(1 / sigmaX));
		double actual = nc3.getWeight(i, j);
		assertEquals(expected, actual, 0);

		double actualSwap = nc3.getWeight(j, i);
		assertEquals(expected, actualSwap, 0);
	}

	@Test
	public void test_3_5() {
		int i = 5;
		int j = 1;

		double expected = Math.exp(-Math.abs(0 - 2) / sigmaI) * Math.exp(-Math.pow(Math.sqrt(2), 2) / sigmaX);
		double actual = nc4.getWeight(i, j);
		assertEquals(expected, actual, 0);

		double actualSwap = nc4.getWeight(j, i);
		assertEquals(expected, actualSwap, 0);
	}

	@Test
	public void test_3_6() {
		int i = 3;
		int j = 8;

		double expected = Math.exp(-Math.abs(1 - 1) / sigmaI) * Math.exp(-Math.pow(Math.sqrt(5), 2) / sigmaX);
		double actual = nc4.getWeight(i, j);
		assertEquals(expected, actual, 0);

		double actualSwap = nc4.getWeight(j, i);
		assertEquals(expected, actualSwap, 0);
	}

	// Test method setSimilarityMatrix()

	@Test
	public void test_4_1() {
		nc1.setSimilarityMatrix();

		double[][] expectedW = { { 1 } };
		double[][] actualW = nc1.getW().getData();
		Assert.assertArrayEquals(expectedW, actualW);
	}

}
