package util;

import static org.junit.Assert.*;

import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Test;

/**
 * Unit test cases for util.MathNC.
 * 
 * @author P Talbot
 */
public class MathNCTest {

	// Test method sumRow()

	@Test
	public void test_1_1() {
		RealMatrix W = new OpenMapRealMatrix(1, 1);
		W.setEntry(0, 0, 2);

		double expected = 2d;
		double actual = util.MathNC.sumRow(W, 0);
		assertEquals(expected, actual, 0);
	}

	@Test
	public void test_1_2() {
		RealMatrix W = new OpenMapRealMatrix(2, 2);
		W.setRow(0, new double[] { 1, 3 });
		W.setRow(1, new double[] { 5, 11 });

		double expected_row0 = 4d;
		double actual_row0 = util.MathNC.sumRow(W, 0);
		assertEquals(expected_row0, actual_row0, 0);

		double expected_row1 = 16d;
		double actual_row1 = util.MathNC.sumRow(W, 1);
		assertEquals(expected_row1, actual_row1, 0);
	}

	// Test method sumColumn()

	@Test
	public void test_2_1() {
		RealMatrix W = new OpenMapRealMatrix(1, 1);
		W.setEntry(0, 0, 2);

		double expected = 2d;
		double actual = util.MathNC.sumColumn(W, 0);
		assertEquals(expected, actual, 0);
	}

	@Test
	public void test_2_2() {
		RealMatrix W = new OpenMapRealMatrix(2, 2);
		W.setColumn(0, new double[] { 1, 3 });
		W.setColumn(1, new double[] { 5, 11 });

		double expected_row0 = 4d;
		double actual_row0 = util.MathNC.sumColumn(W, 0);
		assertEquals(expected_row0, actual_row0, 0);

		double expected_row1 = 16d;
		double actual_row1 = util.MathNC.sumColumn(W, 1);
		assertEquals(expected_row1, actual_row1, 0);
	}

	// Test method isSymmetric()

	@Test
	public void test_9_1() {
		RealMatrix A = new OpenMapRealMatrix(1, 1);
		A.setEntry(0, 0, 1);

		double accuracy = 0;

		boolean expected = true;
		boolean actual = util.MathNC.isSymmetric(A, accuracy);
		assertEquals(expected, actual);
	}

	@Test
	public void test_9_2() {
		RealMatrix A = new OpenMapRealMatrix(1, 1);
		A.setEntry(0, 0, -1);

		double accuracy = 0;

		boolean expected = true;
		boolean actual = util.MathNC.isSymmetric(A, accuracy);
		assertEquals(expected, actual);
	}

	@Test
	public void test_9_3() {
		RealMatrix A = new OpenMapRealMatrix(2, 2);
		A.setRow(0, new double[] { 1, 0 });
		A.setRow(1, new double[] { 0, 1 });

		double accuracy = 0;

		boolean expected = true;
		boolean actual = util.MathNC.isSymmetric(A, accuracy);
		assertEquals(expected, actual);
	}

	@Test
	public void test_9_4() {
		RealMatrix A = new OpenMapRealMatrix(2, 2);
		A.setRow(0, new double[] { 1, 1 });
		A.setRow(1, new double[] { 0, 1 });

		double accuracy = 0;

		boolean expected = false;
		boolean actual = util.MathNC.isSymmetric(A, accuracy);
		assertEquals(expected, actual);
	}

	@Test
	public void test_9_5() {
		RealMatrix A = new OpenMapRealMatrix(2, 2);
		A.setRow(0, new double[] { 1, 0.2 });
		A.setRow(1, new double[] { 0, 1 });

		double accuracy = 0.2;

		boolean expected = true;
		boolean actual = util.MathNC.isSymmetric(A, accuracy);
		assertEquals(expected, actual);
	}

	@Test
	public void test_9_6() {
		RealMatrix A = new OpenMapRealMatrix(2, 2);
		A.setRow(0, new double[] { 1, 0.2 });
		A.setRow(1, new double[] { 0, 1 });

		double accuracy = 0.1;

		boolean expected = false;
		boolean actual = util.MathNC.isSymmetric(A, accuracy);
		assertEquals(expected, actual);
	}

	@Test
	public void test_9_7() {
		RealMatrix A = new OpenMapRealMatrix(3, 3);
		A.setRow(0, new double[] { 1, 0, 0 });
		A.setRow(1, new double[] { 0, 1, 0 });
		A.setRow(2, new double[] { 0, 0, 1 });

		double accuracy = 0;

		boolean expected = true;
		boolean actual = util.MathNC.isSymmetric(A, accuracy);
		assertEquals(expected, actual);
	}

	@Test
	public void test_9_8() {
		RealMatrix A = new OpenMapRealMatrix(3, 3);
		A.setRow(0, new double[] { 1, 2, 6 });
		A.setRow(1, new double[] { 2, 1, -1 });
		A.setRow(2, new double[] { 6, -1, 1 });

		double accuracy = 0;

		boolean expected = true;
		boolean actual = util.MathNC.isSymmetric(A, accuracy);
		assertEquals(expected, actual);
	}

	@Test
	public void test_9_9() {
		RealMatrix A = new OpenMapRealMatrix(3, 3);
		A.setRow(0, new double[] { 1, 2, 0.1 });
		A.setRow(1, new double[] { 2, 1, -1 });
		A.setRow(2, new double[] { 6, -1, 1 });

		double accuracy = 0;

		boolean expected = false;
		boolean actual = util.MathNC.isSymmetric(A, accuracy);
		assertEquals(expected, actual);
	}

	@Test
	public void test_9_10() {
		RealMatrix A = new OpenMapRealMatrix(3, 3);
		A.setRow(0, new double[] { 1, 2, 6 });
		A.setRow(1, new double[] { -1, 15, -1 });
		A.setRow(2, new double[] { 6, -1, 1 });

		double accuracy = 0;

		boolean expected = false;
		boolean actual = util.MathNC.isSymmetric(A, accuracy);
		assertEquals(expected, actual);
	}

}
