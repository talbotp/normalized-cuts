package util;

import java.util.Arrays;

import org.apache.commons.math3.linear.RealMatrix;

/**
 * Here are some useful static Math related methods that either have been used
 * in the development of NormalizedCuts or are currently in use in Normalized
 * Cuts.
 * 
 * @author P Talbot
 *
 */
public class MathNC {

	/**
	 * Method to check if a matrix is symmetric.
	 * 
	 * We use this as a a sanity check that the matrix we are checking is in
	 * fact symmetric, otherwise Lanczos will simply fail.
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
	 * Calculate the l2 norm, i.e the Euclidean Norm.
	 * 
	 * @param arr
	 *            is the array we calculate the l2 norm of.
	 * @return the Euclidean norm of array arr.
	 */
	public static double l2Norm(double[] arr) {

		if (arr.length < 1)
			throw new IllegalArgumentException(
					"The array requires a length of at least 1.");

		double sum = 0;
		for (double val : arr) {
			sum += Math.pow(Math.abs(val), 2);
		}
		return Math.sqrt(sum);
	}

	/**
	 * Calculate the l2 norm of the difference between two arrays.
	 * 
	 * @param arr1
	 *            is the first array
	 * @param arr2
	 *            is the second array.
	 * @return ||arr1 - arr2||
	 */
	public static double l2NormDiff(double[] arr1, double[] arr2) {
		int length = arr1.length;
		double[] arr = new double[arr1.length];
		for (int i = 0; i < length; ++i) {
			arr[i] = arr1[i] - arr2[i];
		}
		return l2Norm(arr);
	}

	/**
	 * Return the median value in an array.
	 * 
	 * @param arr
	 *            is an array.
	 * @return the median value in arr.
	 */
	public static double getMedian(double[] arr) {
		double[] tmp = arr.clone();
		Arrays.sort(tmp);
		double median;
		if (tmp.length % 2 == 0)
			median = ((double) tmp[tmp.length / 2]
					+ (double) tmp[tmp.length / 2 - 1]) / 2;
		else
			median = (double) tmp[tmp.length / 2];
		return median;
	}

}
