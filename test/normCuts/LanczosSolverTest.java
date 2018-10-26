package normCuts;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Assert;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import normCuts.LanczosSolver;
import normCuts.LanczosSolver.EigenValueVectorPair;

/**
 * Unit test cases for class src.normCuts.LanczosSolver
 * 
 * @author P Talbot
 */
public class LanczosSolverTest {

	private final double EIGVAL_ACCURACY = 0.0001;

	@Rule
	public final ExpectedException exception = ExpectedException.none();

	// Test constructors, getters setters etc.

	@Test
	public void test_1_1() {
		RealMatrix A = new OpenMapRealMatrix(1, 2);
		int eigenIndex = 1;
		int maxIterations = 10000;
		boolean compareEigenvectors = true;
		double accuracy = 0;

		// Unsymmetric matrix
		exception.expect(IllegalArgumentException.class);

		@SuppressWarnings("unused")
		LanczosSolver ls = new LanczosSolver(A, eigenIndex, maxIterations, compareEigenvectors, accuracy);
	}

	@Test
	public void test_1_2() {
		RealMatrix A = new OpenMapRealMatrix(2, 2);
		int eigenIndex = 1;
		int maxIterations = 10000;
		boolean compareEigenvectors = true;
		double accuracy = 0;

		LanczosSolver ls = new LanczosSolver(A, eigenIndex, maxIterations, compareEigenvectors, accuracy);

		int expectedIndex = 1;
		int actualIndex = ls.getEigenIndex();
		assertEquals(expectedIndex, actualIndex);

		int expectedIterations = 2;
		int actualIterations = ls.getMaxIterations();
		assertEquals(expectedIterations, actualIterations);

		boolean expCompare = true;
		boolean actualCmp = ls.isCompareEigenvectors();
		assertEquals(expCompare, actualCmp);

		double expAcc = 0;
		double actualAcc = ls.getAccuracy();
		assertEquals(expAcc, actualAcc, 0);
	}

	@Test
	public void test_1_3() {
		RealMatrix A = new OpenMapRealMatrix(2, 2);
		int eigenIndex = 1;
		int maxIterations = 10000;
		boolean compareEigenvectors = true;
		double accuracy = -1;

		// Negative accuracy
		exception.expect(IllegalArgumentException.class);
		@SuppressWarnings("unused")
		LanczosSolver ls = new LanczosSolver(A, eigenIndex, maxIterations, compareEigenvectors, accuracy);
	}

	@Test
	public void test_1_4() {
		RealMatrix A = new OpenMapRealMatrix(2, 2);
		int eigenIndex = -1;
		int maxIterations = 10000;
		boolean compareEigenvectors = true;
		double accuracy = -1;

		// Invalid range of eigen index.
		exception.expect(IllegalStateException.class);
		@SuppressWarnings("unused")
		LanczosSolver ls = new LanczosSolver(A, eigenIndex, maxIterations, compareEigenvectors, accuracy);
	}

	@Test
	public void test_1_5() {
		RealMatrix A = new OpenMapRealMatrix(2, 2);
		int eigenIndex = 3;
		int maxIterations = 10000;
		boolean compareEigenvectors = true;
		double accuracy = -1;

		// Invalid range of eigen index.
		exception.expect(IllegalStateException.class);
		@SuppressWarnings("unused")
		LanczosSolver ls = new LanczosSolver(A, eigenIndex, maxIterations, compareEigenvectors, accuracy);
	}

	@Test
	public void test_1_6() {
		RealMatrix A = new OpenMapRealMatrix(1, 1);
		int eigenIndex = 3;
		int maxIterations = 10000;
		boolean compareEigenvectors = true;
		double accuracy = -1;

		// Matrix too small.
		exception.expect(IllegalArgumentException.class);
		@SuppressWarnings("unused")
		LanczosSolver ls = new LanczosSolver(A, eigenIndex, maxIterations, compareEigenvectors, accuracy);
	}

	// Test method makeV()

	@Test
	public void test_2_1() {
		int rowSize = 1;
		ArrayList<double[]> V = new ArrayList<>();
		V.add(new double[] { 1 });
		V.add(new double[] { 2 });
		V.add(new double[] { 3 });

		RealMatrix expected = new Array2DRowRealMatrix(rowSize, V.size());
		expected.setRow(0, new double[] { 1, 2, 3 });
		RealMatrix actual = LanczosSolver.makeV(rowSize, V);
		assertEquals(expected, actual);
	}

	@Test
	public void test_2_2() {
		int rowSize = 2;
		ArrayList<double[]> V = new ArrayList<>();
		V.add(new double[] { 1, 2 });
		V.add(new double[] { 2, 3 });
		V.add(new double[] { 3, 7 });

		RealMatrix expected = new Array2DRowRealMatrix(rowSize, V.size());
		expected.setRow(0, new double[] { 1, 2, 3 });
		expected.setRow(1, new double[] { 2, 3, 7 });
		RealMatrix actual = LanczosSolver.makeV(rowSize, V);
		assertEquals(expected, actual);
	}

	@Test
	public void test_2_3() {
		int rowSize = 3;
		ArrayList<double[]> V = new ArrayList<>();
		V.add(new double[] { 1, 2, -1 });
		V.add(new double[] { 2, 3, 2 });
		V.add(new double[] { 3, 7, 55 });

		RealMatrix expected = new Array2DRowRealMatrix(rowSize, V.size());
		expected.setRow(0, new double[] { 1, 2, 3 });
		expected.setRow(1, new double[] { 2, 3, 7 });
		expected.setRow(2, new double[] { -1, 2, 55 });
		RealMatrix actual = LanczosSolver.makeV(rowSize, V);
		assertEquals(expected, actual);
	}

	// Test toPrimitive

	@Test
	public void test_3_1() {
		ArrayList<Double> list = new ArrayList<>();

		double[] expArr = new double[0];
		double[] actArr = LanczosSolver.toPrimitive(list);

		Assert.assertArrayEquals(expArr, actArr, 0);
	}

	@Test
	public void test_3_2() {
		ArrayList<Double> list = new ArrayList<>(Arrays.asList(1d));

		double[] expArr = new double[] { 1 };
		double[] actArr = LanczosSolver.toPrimitive(list);

		Assert.assertArrayEquals(expArr, actArr, 0);
	}

	@Test
	public void test_3_3() {
		ArrayList<Double> list = new ArrayList<>(Arrays.asList(1d, 2d, 3d, 11d));

		double[] expArr = new double[] { 1, 2, 3, 11 };
		double[] actArr = LanczosSolver.toPrimitive(list);

		Assert.assertArrayEquals(expArr, actArr, 0);
	}

	// Test solve

	@Test
	public void test_4_1() {
		RealMatrix A = new OpenMapRealMatrix(2, 2);
		A.setRow(0, new double[] { 1, 2 });
		A.setRow(1, new double[] { 2, 1 });

		int eigenIndex = 0;
		int maxIterations = 100;
		boolean compareEigenvectors = true;
		double accuracy = 0;

		LanczosSolver ls = new LanczosSolver(A, eigenIndex, maxIterations, compareEigenvectors, accuracy);
		EigenValueVectorPair pair = ls.solve();

		double expected_eigval = -1;
		double actual_eigval = pair.eigenvalue;
		assertEquals(expected_eigval, actual_eigval, EIGVAL_ACCURACY);

		ls.setEigenIndex(1);
		pair = ls.solve();
		expected_eigval = 3;
		actual_eigval = pair.eigenvalue;
		assertEquals(expected_eigval, actual_eigval, EIGVAL_ACCURACY);
	}

	@Test
	public void test_4_2() {
		RealMatrix A = new OpenMapRealMatrix(3,3);
		A.setRow(0, new double[] { 1, 2, 0 });
		A.setRow(1, new double[] { 2, 1, 1 });
		A.setRow(2, new double[] { 0, 1, 0 });

		int eigenIndex = 0;
		int maxIterations = 100;
		boolean compareEigenvectors = true;
		double accuracy = 0;

		LanczosSolver ls = new LanczosSolver(A, eigenIndex, maxIterations, compareEigenvectors, accuracy);
		EigenValueVectorPair pair = ls.solve();

		double expected_eigval = -1.39138;
		double actual_eigval = pair.eigenvalue;
		assertEquals(expected_eigval, actual_eigval, EIGVAL_ACCURACY);

		ls.setEigenIndex(1);
		pair = ls.solve();
		expected_eigval = 0.227134;
		actual_eigval = pair.eigenvalue;
		assertEquals(expected_eigval, actual_eigval, EIGVAL_ACCURACY);
		
		ls.setEigenIndex(2);
		pair = ls.solve();
		expected_eigval = 3.16425;
		actual_eigval = pair.eigenvalue;
		assertEquals(expected_eigval, actual_eigval, EIGVAL_ACCURACY);
	}

}
