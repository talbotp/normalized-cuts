package normCuts;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import normCuts.EigenvalueIndexPair;

/**
 * Test cases for src.normCuts.EigenvalueIndexPair
 * 
 * @author P Talbot
 */
public class EigenvalueIndexPairTest {

	// Test constructor, getter toString etc

	@Test
	public void test_1_1() {
		int index = 1;
		double evalue = 17.5;
		EigenvalueIndexPair e = new EigenvalueIndexPair(index, evalue);

		int expectedIndex = 1;
		int actualIndex = e.index;
		assertEquals(expectedIndex, actualIndex);

		double expectedVal = 17.5;
		double actualVal = e.eigenvalue;
		assertEquals(expectedVal, actualVal, 0);

		String expectedStr = 17.5 + "";
		String actualStr = e.toString();
		System.out.println(actualStr);
		assertEquals(expectedStr, actualStr);
	}

	@Test
	public void test_1_2() {
		int index = -1;
		double evalue = -123.123;
		EigenvalueIndexPair e = new EigenvalueIndexPair(index, evalue);

		int expectedIndex = -1;
		int actualIndex = e.index;
		assertEquals(expectedIndex, actualIndex);

		double expectedVal = -123.123;
		double actualVal = e.eigenvalue;
		assertEquals(expectedVal, actualVal, 0);

		String expectedStr = -123.123 + "";
		String actualStr = e.toString();
		System.out.println(actualStr);
		assertEquals(expectedStr, actualStr);
	}

	// Test compareTo

	@Test
	public void test_2_1() {
		int index1 = 0;
		double evalue1 = 1;
		EigenvalueIndexPair e1 = new EigenvalueIndexPair(index1, evalue1);

		int index2 = 0;
		double evalue2 = 1;
		EigenvalueIndexPair e2 = new EigenvalueIndexPair(index2, evalue2);

		double expected = 0;
		double actual = e1.compareTo(e2);

		assertEquals(expected, actual, 0);
	}

	@Test
	public void test_2_2() {
		int index1 = 0;
		double evalue1 = 1;
		EigenvalueIndexPair e1 = new EigenvalueIndexPair(index1, evalue1);

		int index2 = 432432234;
		double evalue2 = 1;
		EigenvalueIndexPair e2 = new EigenvalueIndexPair(index2, evalue2);

		double expected = 0;
		double actual = e1.compareTo(e2);

		assertEquals(expected, actual, 0);
	}

	@Test
	public void test_2_3() {
		int index1 = 0;
		double evalue1 = 1;
		EigenvalueIndexPair e1 = new EigenvalueIndexPair(index1, evalue1);

		int index2 = 0;
		double evalue2 = 3;
		EigenvalueIndexPair e2 = new EigenvalueIndexPair(index2, evalue2);

		double expected = -1;
		double actual = e1.compareTo(e2);

		assertEquals(expected, actual, 0);
	}

	@Test
	public void test_2_4() {
		int index1 = 0;
		double evalue1 = 1;
		EigenvalueIndexPair e1 = new EigenvalueIndexPair(index1, evalue1);

		int index2 = 3242;
		double evalue2 = 3;
		EigenvalueIndexPair e2 = new EigenvalueIndexPair(index2, evalue2);

		double expected = -1;
		double actual = e1.compareTo(e2);

		assertEquals(expected, actual, 0);
	}

	@Test
	public void test_2_5() {
		int index1 = 0;
		double evalue1 = 22;
		EigenvalueIndexPair e1 = new EigenvalueIndexPair(index1, evalue1);

		int index2 = 3242;
		double evalue2 = 3;
		EigenvalueIndexPair e2 = new EigenvalueIndexPair(index2, evalue2);

		double expected = 1;
		double actual = e1.compareTo(e2);

		assertEquals(expected, actual, 0);
	}

	@Test
	public void test_2_6() {
		int index1 = 3142;
		double evalue1 = 22;
		EigenvalueIndexPair e1 = new EigenvalueIndexPair(index1, evalue1);

		int index2 = 3242;
		double evalue2 = 21.9999;
		EigenvalueIndexPair e2 = new EigenvalueIndexPair(index2, evalue2);

		double expected = 1;
		double actual = e1.compareTo(e2);

		assertEquals(expected, actual, 0);
	}

}
