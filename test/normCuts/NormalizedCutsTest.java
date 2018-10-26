package normCuts;

import static org.junit.Assert.assertEquals;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import normCuts.NormalizedCuts;

/**
 * Test cases for Normalized Cuts Methods.
 *
 * We have created a basic extension of the Normalized Cuts named
 * NormalizedeCutsTestStub which is a very basic extension of The Abstract class
 * NormalizedCuts so we can safely use it to test its methods using the basic
 * Similarity Matrix W's defined below.
 * 
 * @author Andrew P Talbot
 */
public class NormalizedCutsTest {

	private final double TEST_ACCURACY = 0.000001;

	// Here are simple similarity matrices to be used to test the Normalized
	// Cuts Methods.
	public static final double[][] DATA_1 = { { 1 } };
	public static final double[][] DATA_2 = { { 1, 2 }, { 3, 4 } };
	public static final double[][] DATA_3 = { { 1, 0 }, { 0, 2 } };
	public static final double[][] DATA_4 = { { 1, 0, 1 }, { 1, 1, 2 }, { 0, 0, 1 } };

	private static double sigmaI = 2;
	private static double sigmaX = 5;
	private static double r = 8;
	private static int l = 2;
	private static double accuracy = 0.001;
	private static boolean useEigenvectors = true;
	private static double nCutThreshold = 2;
	private static int nClusters = 1;
	private static int nCutType = NormalizedCuts.PARTITION_BY_ZERO;

	private static NormalizedCuts nc1;

	@Rule
	public final ExpectedException exception = ExpectedException.none();

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		double[][] W = DATA_1;
		nc1 = new NormalizedCutsTestStub(sigmaI, sigmaX, r, l, nCutThreshold, nCutType, nClusters, accuracy, useEigenvectors, W);
	}

	// Test the Constructor, Getters and Setters.

	@Test
	public void test_1_1() {
		double expectedI = sigmaI;
		double actualI = nc1.getSigmaI();
		assertEquals(expectedI, actualI, 0);

		double expectedX = sigmaX;
		double actualX = nc1.getSigmaX();
		assertEquals(expectedX, actualX, 0);

		double expectedR = r;
		double actualR = nc1.getR();
		assertEquals(expectedR, actualR, 0);

		int expectedL = l;
		double actualL = nc1.getL();
		assertEquals(expectedL, actualL, 0);

		double expectedThreshold = nCutThreshold;
		double actualThreshold = nc1.getNCutThreshold();
		assertEquals(expectedThreshold, actualThreshold, 0);

		int expectedNClusters = nClusters;
		int actualNClusters = nc1.getNClusters();
		assertEquals(expectedNClusters, actualNClusters);

		int expectedCutType = nCutType;
		int actualCutType = nc1.getNCutType();
		assertEquals(expectedCutType, actualCutType);
	}

	// Test Method NormalizedCuts.getSpatialWeight().

	@Test
	public void test_2_1() {
		double[] arr1 = { 0 };
		double[] arr2 = { 0 };

		double expected = 0;
		double actual = NormalizedCuts.getSpatialWeight(arr1, arr2);
		assertEquals(expected, actual, 0);
	}

	@Test
	public void test_2_2() {
		double[] arr1 = { 0 };
		double[] arr2 = { 1 };

		double expected = 1;
		double actual = NormalizedCuts.getSpatialWeight(arr1, arr2);
		assertEquals(expected, actual, 0);

		double expectedSwap = 1;
		double actualSwap = NormalizedCuts.getSpatialWeight(arr2, arr1);
		assertEquals(expectedSwap, actualSwap, 0);
	}

	@Test
	public void test_2_3() {
		double[] arr1 = { 2, 2 };
		double[] arr2 = { 0, 0 };

		double expected = Math.sqrt(8);
		double actual = NormalizedCuts.getSpatialWeight(arr1, arr2);
		assertEquals(expected, actual, 0);
	}

	@Test
	public void test_2_4() {
		double[] arr1 = { 3, 3 };
		double[] arr2 = { 2, 2 };

		double expected = Math.sqrt(2);
		double actual = NormalizedCuts.getSpatialWeight(arr1, arr2);
		assertEquals(expected, actual, 0);
	}

	// Test methods assoc()

	@Test
	public void test_3_2() {
		double[][] arr = { { 0 } };
		RealMatrix W = MatrixUtils.createRealMatrix(arr);

		Integer[] p1 = { 0 };

		double expected = 0;
		double actual = NormalizedCuts.assoc(p1, true, W);
		assertEquals(expected, actual, 0);

		double actual2 = NormalizedCuts.assoc(p1, false, W);
		assertEquals(expected, actual2, 0);
	}

	@Test
	public void test_3_3() {
		double[][] arr = { { 1, 2 }, { 3, 4 } };
		RealMatrix W = MatrixUtils.createRealMatrix(arr);

		Integer[] p1 = {};

		double expected1 = 0;
		double actual1 = NormalizedCuts.assoc(p1, true, W);
		assertEquals(expected1, actual1, 0);

		double expected2 = 0;
		double actual2 = NormalizedCuts.assoc(p1, false, W);
		assertEquals(expected2, actual2, 0);

		Integer[] p2 = { 0 };

		double expected3 = 3;
		double actual3 = NormalizedCuts.assoc(p2, true, W);
		assertEquals(expected3, actual3, 0);

		double expected4 = 1;
		double actual4 = NormalizedCuts.assoc(p2, false, W);
		assertEquals(expected4, actual4, 0);

		Integer[] p3 = { 1 };

		double expected5 = 7;
		double actual5 = NormalizedCuts.assoc(p3, true, W);
		assertEquals(expected5, actual5, 0);

		double expected6 = 4;
		double actual6 = NormalizedCuts.assoc(p3, false, W);
		assertEquals(expected6, actual6, 0);

		Integer[] p4 = { 0, 1 };
		double expected7 = 10;
		double actual7 = NormalizedCuts.assoc(p4, true, W);
		assertEquals(expected7, actual7, 0);

		double expected8 = 10;
		double actual8 = NormalizedCuts.assoc(p4, false, W);
		assertEquals(expected8, actual8, 0);
	}

	// Test Nassoc()

	@Test
	public void test_3_7() {
		double[][] arr = { { 1, 2 }, { 3, 4 } };
		RealMatrix W = MatrixUtils.createRealMatrix(arr);

		Integer[] p1 = { 0 };
		Integer[] p2 = { 1 };

		double expected1 = 1d / 3 + 4d / 7;
		double actual1 = NormalizedCuts.Nassoc(p1, p2, W);
		assertEquals(expected1, actual1, 0);

		double expected2 = 1d / 3 + 4d / 7;
		double actual2 = NormalizedCuts.Nassoc(p2, p1, W);
		assertEquals(expected2, actual2, 0);
	}

	// Test NCut()

	@Test
	public void test_4_1() {
		double[][] arr = { { 1 } };
		RealMatrix W = MatrixUtils.createRealMatrix(arr);

		Integer[] p1 = {};
		Integer[] p2 = {};

		double expected = Double.MAX_VALUE;
		double actual = NormalizedCuts.NCut(p1, p2, W);
		assertEquals(expected, actual, 0);
	}

	// Test method getD()

	@Test
	public void test_5_1() {
		RealMatrix W = new OpenMapRealMatrix(1, 1);
		W.setEntry(0, 0, 1);

		RealMatrix expected = W;
		RealMatrix actual = NormalizedCuts.getD(W);
		assertEquals(expected, actual);
	}

	@Test
	public void test_5_2() {
		RealMatrix W = new OpenMapRealMatrix(2, 2);
		W.setRow(0, new double[] { 0, 0 });
		W.setRow(1, new double[] { 0, 0 });

		RealMatrix expected = W;
		RealMatrix actual = NormalizedCuts.getD(W);
		assertEquals(expected, actual);
	}

	@Test
	public void test_5_3() {
		RealMatrix W = new OpenMapRealMatrix(2, 2);
		W.setRow(0, new double[] { 1, 0 });
		W.setRow(1, new double[] { 3, 2 });

		RealMatrix expected = new OpenMapRealMatrix(2, 2);
		expected.setRow(0, new double[] { 1, 0 });
		expected.setRow(1, new double[] { 0, 5 });

		RealMatrix actual = NormalizedCuts.getD(W);
		assertEquals(expected, actual);
	}

	@Test
	public void test_5_4() {
		RealMatrix W = new OpenMapRealMatrix(3, 3);
		W.setRow(0, new double[] { 1, 2, 3 });
		W.setRow(1, new double[] { 4, 5, 6 });
		W.setRow(2, new double[] { 7, 8, 9 });

		RealMatrix expected = new OpenMapRealMatrix(3, 3);
		expected.setRow(0, new double[] { 6, 0, 0 });
		expected.setRow(1, new double[] { 0, 15, 0 });
		expected.setRow(2, new double[] { 0, 0, 24 });

		RealMatrix actual = NormalizedCuts.getD(W);
		assertEquals(expected, actual);
	}

	// Test method getD_minus_1_2()

	@Test
	public void test_6_1() {
		RealMatrix W = new OpenMapRealMatrix(1, 1);
		W.setEntry(0, 0, 1);

		RealMatrix expected = W;
		RealMatrix actual = NormalizedCuts.getD_minus_1_2(W);
		assertEquals(expected, actual);
	}

	@Test
	public void test_6_2() {
		RealMatrix W = new OpenMapRealMatrix(1, 1);
		W.setEntry(0, 0, 2);

		RealMatrix expected = new OpenMapRealMatrix(1, 1);
		expected.setEntry(0, 0, 1 / Math.sqrt(2));

		RealMatrix actual = NormalizedCuts.getD_minus_1_2(W);
		assertEquals(expected, actual);
	}

	@Test
	public void test_6_3() {
		RealMatrix W = new OpenMapRealMatrix(2, 2);
		W.setEntry(0, 0, 4);
		W.setEntry(1, 1, 9);

		RealMatrix expected = new OpenMapRealMatrix(2, 2);
		expected.setEntry(0, 0, 1d / 2);
		expected.setEntry(1, 1, 1d / 3);

		RealMatrix actual = NormalizedCuts.getD_minus_1_2(W);
		assertEquals(expected, actual);
	}

	// Test method getMatrixToEigenDecompose()

	@Test
	public void test_7_1() {
		RealMatrix W = new OpenMapRealMatrix(1, 1);
		W.setEntry(0, 0, 1);

		RealMatrix expected = new OpenMapRealMatrix(1, 1);
		expected.setEntry(0, 0, 0);

		RealMatrix D = NormalizedCuts.getD(W);
		RealMatrix actual = NormalizedCuts.getMatrixToEigenDecompose(D, NormalizedCuts.getD_minus_1_2(D), W);
		assertEquals(expected, actual);
	}

	@Test
	public void test_7_2() {
		RealMatrix W = new OpenMapRealMatrix(1, 1);
		W.setEntry(0, 0, 5);

		RealMatrix expected = new OpenMapRealMatrix(1, 1);
		expected.setEntry(0, 0, 0);

		RealMatrix D = NormalizedCuts.getD(W);
		RealMatrix actual = NormalizedCuts.getMatrixToEigenDecompose(D, NormalizedCuts.getD_minus_1_2(D), W);

		assertEquals(expected, actual);
	}

	@Test
	public void test_7_3() {
		RealMatrix W = new OpenMapRealMatrix(2, 2);
		W.setRow(0, new double[] { 1, 0 });
		W.setRow(1, new double[] { 1, 1 });

		RealMatrix expected = new OpenMapRealMatrix(2, 2);
		expected.setRow(0, new double[] { 0, 0 });
		expected.setRow(1, new double[] { -1d / Math.sqrt(2), 1d / 2 });

		RealMatrix D = NormalizedCuts.getD(W);
		RealMatrix actual = NormalizedCuts.getMatrixToEigenDecompose(D, NormalizedCuts.getD_minus_1_2(D), W);

		// Slight inaccuracues so we dont test the objects, we test their
		// columns.
		for (int i = 0; i < W.getColumnDimension(); i++) {
			System.out.println();
			Assert.assertArrayEquals(expected.getColumn(i), actual.getColumn(i), TEST_ACCURACY);
		}
	}

	@Test
	public void test_7_4() {
		RealMatrix W = new OpenMapRealMatrix(2, 2);
		W.setRow(0, new double[] { 1, 2 });
		W.setRow(1, new double[] { 2, 1 });

		RealMatrix expected = new OpenMapRealMatrix(2, 2);
		expected.setRow(0, new double[] { 2d / 3d, -2d / 3d });
		expected.setRow(1, new double[] { -2d / 3d, 2d / 3d });

		RealMatrix D = NormalizedCuts.getD(W);
		RealMatrix actual = NormalizedCuts.getMatrixToEigenDecompose(D, NormalizedCuts.getD_minus_1_2(D), W);

		// Slight inaccuracues so we dont test the objects, we test their
		// columns.
		for (int i = 0; i < W.getColumnDimension(); i++) {
			Assert.assertArrayEquals(expected.getColumn(i), actual.getColumn(i), TEST_ACCURACY);
		}
	}

	// Test method NormalizedCuts.getNewW()

	@Test
	public void test_8_1() {
		RealMatrix W = new OpenMapRealMatrix(1, 1);
		W.setEntry(0, 0, 1);

		Map<Integer, Integer> partition = new HashMap<>();
		partition.put(0, 0);

		RealMatrix expected = new OpenMapRealMatrix(1, 1);
		expected.setEntry(0, 0, 1);

		RealMatrix actual = NormalizedCuts.getNewW(W, partition);

		assertEquals(expected, actual);
	}

	@Test
	public void test_8_2() {
		RealMatrix W1 = new OpenMapRealMatrix(2, 2);
		W1.setRow(0, new double[] { 1, 3 });
		W1.setRow(1, new double[] { 3, 2 });

		Map<Integer, Integer> partition1 = new HashMap<>();
		partition1.put(0, 154);
		RealMatrix expected1 = new OpenMapRealMatrix(1, 1);
		expected1.setEntry(0, 0, 1);
		RealMatrix actual1 = NormalizedCuts.getNewW(W1, partition1);
		assertEquals(expected1, actual1);

		Map<Integer, Integer> partition2 = new HashMap<>();
		partition2.put(1, 154);
		RealMatrix expected2 = new OpenMapRealMatrix(1, 1);
		expected2.setEntry(0, 0, 2);
		RealMatrix actual2 = NormalizedCuts.getNewW(W1, partition2);
		assertEquals(expected2, actual2);
	}

	@Test
	public void test_8_3() {
		RealMatrix W1 = new OpenMapRealMatrix(2, 2);
		W1.setRow(0, new double[] { 1, 3 });
		W1.setRow(1, new double[] { 3, 2 });

		Map<Integer, Integer> partition1 = new HashMap<>();
		partition1.put(0, 154);
		partition1.put(1, 213214);
		RealMatrix expected1 = W1;
		RealMatrix actual1 = NormalizedCuts.getNewW(W1, partition1);
		assertEquals(expected1, actual1);
	}

	@Test
	public void test_8_4() {
		RealMatrix W1 = new OpenMapRealMatrix(3, 3);
		W1.setRow(0, new double[] { 1, 3, 5 });
		W1.setRow(1, new double[] { 3, 2, 11 });
		W1.setRow(2, new double[] { 55, 3, 41 });

		Map<Integer, Integer> partition1 = new HashMap<>();
		partition1.put(1, 154);
		RealMatrix expected1 = new OpenMapRealMatrix(1, 1);
		expected1.setEntry(0, 0, 2);
		RealMatrix actual1 = NormalizedCuts.getNewW(W1, partition1);
		assertEquals(expected1, actual1);

		Map<Integer, Integer> partition2 = new HashMap<>();
		partition2.put(0, 154);
		partition2.put(2, 234);
		RealMatrix expected2 = new OpenMapRealMatrix(2, 2);
		expected2.setRow(0, new double[] { 1, 5 });
		expected2.setRow(1, new double[] { 55, 41 });
		RealMatrix actual2 = NormalizedCuts.getNewW(W1, partition2);
		assertEquals(expected2, actual2);

		Map<Integer, Integer> partition3 = new HashMap<>();
		partition3.put(0, 154);
		partition3.put(2, 234);
		partition3.put(1, 123);
		RealMatrix expected3 = W1;
		RealMatrix actual3 = NormalizedCuts.getNewW(W1, partition3);
		assertEquals(expected3, actual3);
	}

	@Test
	public void test_8_5() {
		RealMatrix W1 = new OpenMapRealMatrix(3, 3);
		W1.setRow(0, new double[] { 1, 3, 5 });
		W1.setRow(1, new double[] { 3, 2, 11 });
		W1.setRow(2, new double[] { 55, 3, 41 });

		Map<Integer, Integer> partition1 = new HashMap<>();
		partition1.put(1, 154);
		partition1.put(2, 4352);
		partition1.put(4234, 2342);
		partition1.put(3, 3);

		exception.expect(IllegalArgumentException.class);
		NormalizedCuts.getNewW(W1, partition1);
	}

}
