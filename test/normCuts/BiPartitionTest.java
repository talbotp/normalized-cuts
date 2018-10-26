package normCuts;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Unit test cases for trivial class src.normCuts.BiPartition.
 *
 * @author P Talbot
 */
public class BiPartitionTest {

	private static double[][] W_1;
	private static Map<Integer, Integer> partition_1, partition_2, partition_3;
	
	@BeforeClass
	public static void setup() {
		W_1 = new double[][] { { 1, 0 }, { 0, 1 } };

		partition_1 = new HashMap<>();
		partition_1.put(1, 0);
		partition_1.put(1, 0);

		// Same size as partition_1
		partition_2 = new HashMap<>();
		partition_2.put(1, 1);
		partition_2.put(1, 1);

		// Larger partition.
		partition_3 = new HashMap<>();
		partition_3.put(2, 2);
		partition_3.put(3, 3);
		partition_3.put(4, 4);
	}

	/*
	 * In these next tests we test results of the compare method.
	 */

	// Test equal size.
	@Test
	public void test_1_1() {
		RealMatrix W = new Array2DRowRealMatrix(W_1);

		BiPartition bp_1 = new BiPartition(W, partition_1);
		BiPartition bp_2 = new BiPartition(W, partition_2);

		int expected = 0;
		int actual = bp_1.compare(bp_1, bp_2);
		Assert.assertEquals(expected, actual);
	}
	
	// Test smaller
	@Test
	public void test_1_2() {
		RealMatrix W = new Array2DRowRealMatrix(W_1);

		BiPartition bp_1 = new BiPartition(W, partition_1);
		BiPartition bp_3 = new BiPartition(W, partition_3);

		int expected = 2;
		int actual = bp_1.compare(bp_1, bp_3);
		Assert.assertEquals(expected, actual);
	}

	// Test larger
	@Test
	public void test_1_3() {
		RealMatrix W = new Array2DRowRealMatrix(W_1);

		BiPartition bp_3 = new BiPartition(W, partition_3);
		BiPartition bp_1 = new BiPartition(W, partition_2);

		int expected = -2;
		int actual = bp_3.compare(bp_3, bp_1);
		Assert.assertEquals(expected, actual);
	}
	
}
