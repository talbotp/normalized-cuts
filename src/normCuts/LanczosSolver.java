package normCuts;

import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.random.UnitSphereRandomVectorGenerator;

import util.MathNC;

/**
 * Logic behind the Lanczos Algorithm for finding the k most useful extreme
 * eigenvectors of a Hermittian Matrix.
 * 
 * In our case a Hermittian matrix is relaxed to a symmetric matrix. By the k
 * most useful eigenvectors, we mean the k most extreme eigenvectors. For
 * example, given a 100x100 symmetric real matrix, we can use this algorithm to
 * calculate the smallest and largest eigenvectors. (eigenvectors 0-4 and 94-99
 * might be appropriate for this algorithm.
 * 
 * @author P Talbot
 * @version 13/08/2018
 */
public class LanczosSolver {

	private RealMatrix A;
	private int maxIterations;
	private int eigenIndex;
	private boolean compareEigenvectors;
	private double accuracy;

	/**
	 * Create a LanczosSolver object in order to find the eigenvectors that we
	 * need.
	 * 
	 * @param A
	 *            is the matrix we solve for.
	 * @param eigenIndex
	 *            refers to the eigenIndex'th smallest eigenvalue.
	 * @param maxIterations
	 *            refers to the maximum number of Lanczos iterations that we
	 *            will allow before ending the algorithm if convergence doesn't
	 *            already occur.
	 * @param compareEigenvectors
	 *            is true if we wish to compare eigenvectors to see if
	 *            convergence has occured, false if we compare eigenvalues to
	 *            see if convergence has occured.
	 * @param accuracy
	 *            is the parameter which defines how close we require values to
	 *            be for convergence to be considered true.
	 */
	public LanczosSolver(RealMatrix A, int eigenIndex, int maxIterations,
			boolean compareEigenvectors, double accuracy) {
		setA(A);
		setEigenIndex(eigenIndex);
		setMaxIterations(maxIterations);
		setAccuracy(accuracy);
		setCompareEigenvectors(compareEigenvectors);
	}

	/**
	 * Create a LanczosSolver object in order to find the eigenvectors that we
	 * need, when we wish to set the max number of iterations as equal to the
	 * size of our matrix, ie m = n.
	 * 
	 * @param A
	 *            is the matrix we solve for.
	 * @param compareEigenvectors
	 *            is true if we wish to compare eigenvectors to see if
	 *            convergence has occured, false if we compare eigenvalues to
	 *            see if convergence has occured.
	 * @param accuracy
	 *            is the parameter which defines how close we require values to
	 *            be for convergence to be considered true.
	 */
	public LanczosSolver(RealMatrix A, boolean compareEigenvectors,
			double accuracy) {
		setA(A);
		setEigenIndex(eigenIndex);
		setMaxIterations(A.getColumnDimension());
		setAccuracy(accuracy);
		setCompareEigenvectors(compareEigenvectors);
	}

	/**
	 * Solve to get a specific eigenvalue, eigenvector pair from the given field
	 * variables.
	 * 
	 * @return the EigenValueVectorPair corresponding to the index given as a
	 *         field.
	 */
	public EigenValueVectorPair solve() {

		// Init the storage for required variables.
		ArrayList<Double> alphas = new ArrayList<>();
		ArrayList<Double> betas = new ArrayList<>();
		ArrayList<double[]> V = new ArrayList<>();

		// Step 1 : Let v be a random unit vector.
		int sizeA = this.A.getColumnDimension();
		UnitSphereRandomVectorGenerator randVecGen = new UnitSphereRandomVectorGenerator(
				sizeA);
		RealVector v_i = new ArrayRealVector(randVecGen.nextVector());
		V.add(v_i.toArray());

		// Step 2 : Initial Step (reuse these in the loop.)
		RealVector w_dash_i = A.operate(v_i);
		alphas.add(w_dash_i.dotProduct(v_i));
		RealVector w_i = w_dash_i.subtract(v_i.mapMultiply(alphas.get(0)));

		// This is variable which tell us about the convergence.
		boolean iterateAgain = true;
		double eigval = 0;
		double[] eigvec = null;
		double previousEigval = 0;
		double[] previousEigvec = null;
		ArrayList<EigenvalueIndexPair> eigvalpairs = null;

		// Step 3 : Main loop
		for (int i = 1; i < this.maxIterations; i++) {

			betas.add(w_i.getNorm()); // betas.size() == i.

			if (betas.get(i - 1) != 0) {
				v_i = w_i.mapDivide(betas.get(i - 1));
				V.add(v_i.toArray());
			} else {
				System.out.println(
						"beta_" + i + " = 0, we calculate a new vector.");
				throw new IllegalStateException(
						"Please rerun the algorithm, rng has failed.");
			}

			// Can reuse old variables for W as they are not needed.
			w_dash_i = this.A.operate(v_i);

			alphas.add(w_dash_i.dotProduct(v_i));

			RealVector tmp_v = new ArrayRealVector(V.get(i - 1));
			w_i = (w_dash_i.subtract(v_i.mapMultiply(alphas.get(i))))
					.subtract(tmp_v.mapMultiply(betas.get(i - 1)));

			// Make matrix V containing v_i as columns.
			RealMatrix V_matrix = makeV(this.A.getRowDimension(), V);

			// Eigendecompose matrix T.
			EigenDecomposition eigdec = new EigenDecomposition(
					toPrimitive(alphas), toPrimitive(betas));

			// Now we sort the eigenvalues.
			double[] eigvals = eigdec.getRealEigenvalues();
			eigvalpairs = new ArrayList<>();
			for (int j = 0; j < eigvals.length; j++) {
				eigvalpairs.add(new EigenvalueIndexPair(j, eigvals[j]));
			}

			Collections.sort(eigvalpairs);

			// Get the appropriate index for all sized current matrices.
			int required_index = (this.eigenIndex > i)
					? eigvalpairs.get(i - 1).index
					: eigvalpairs.get(this.eigenIndex).index;

			// Get the eigenvalue and corresponding eigenvector.
			eigval = eigdec.getRealEigenvalue(required_index);
			eigvec = V_matrix.operate(eigdec.getEigenvector(required_index))
					.toArray();

			// Skip on first iter.
			if (i != 1) {

				if (this.compareEigenvectors) {

					// Case where we compare components of eigenvectors.
					iterateAgain = false;
					for (int k = 0; k < eigvec.length; k++) {
						if (Math.abs(Math.abs(previousEigvec[k])
								- Math.abs(eigvec[k])) > this.accuracy) {
							iterateAgain = true;
							break;
						}
					}

				} else {

					// Case where we check convergence based on the eigenvalues.
					if (Math.abs(previousEigval - eigval) < accuracy) {
						iterateAgain = false;
					} else
						iterateAgain = true;

				}

				// Check if we iterate again :
				if (iterateAgain) {
					// If we iterate again, set the variables so we can check
					// similarity later.
					previousEigval = eigval;
					previousEigvec = eigvec;
				} else {
					// If we dont, then we return the most recent values.
					break;
				}

			} else {

				// In first iteration we set previous eigval to this.
				previousEigval = eigval;
				previousEigvec = eigvec;
			}

		}

		return new EigenValueVectorPair(eigval, eigvec);
	}

	/**
	 * Compute an arbitrary unit vector which is orthogonal to every vector in
	 * the ArrayList V.
	 * 
	 * @param V
	 *            is the set of vectors we wish to find a vector orthogonal to.
	 * @return a vector orthogonal to all vectors in ArrayList V.
	 */
	public static double[] getOrthogonalUnitVectorOfV(ArrayList<double[]> V) {

		if (V.isEmpty())
			throw new IllegalArgumentException(
					"V is empt, we cannot create a vector for this.");

		// Create the matrix of vectors in V as rows.
		RealMatrix V_rows = new Array2DRowRealMatrix(V.size(), V.get(0).length);
		for (int i = 0; i < V.size(); i++) {
			V_rows.setRow(i, V.get(i));
		}

		// Now we solve V_rows x = 0
		DecompositionSolver solver = new QRDecomposition(V_rows).getSolver();
		RealVector zeros = new ArrayRealVector(new double[V.size()]);

		return solver.solve(zeros).toArray();
	}

	/**
	 * Return a unit vector of the vector provided.
	 * 
	 * @param vec
	 *            is the vector we wish to transform into a linear vector.
	 * @return the unit vector.
	 */
	public double[] toUnitVector(double[] vec) {
		double euc_dist = MathNC.l2Norm(vec);
		for (int i = 0; i < vec.length; i++)
			vec[i] /= euc_dist;
		return vec;
	}

	/**
	 * Creates a matrix V of row dimensions rowSize, where the V_i are the
	 * columns of the matrix.
	 * 
	 * REQUIRES rowSize > 0 && V.size() > 0
	 * 
	 * @param rowSize
	 *            is the row dimensions of A.
	 * @param V
	 *            is the ArrayList containing v_i arrays.
	 * @return a matrix with v_i as its columns.
	 */
	public static RealMatrix makeV(int rowSize, ArrayList<double[]> V) {
		RealMatrix V_matrix = new Array2DRowRealMatrix(rowSize, V.size());
		for (int col = 0; col < V.size(); col++)
			V_matrix.setColumn(col, V.get(col));
		return V_matrix;
	}

	/**
	 * Convert a Double[] to double
	 * 
	 * @return the primitive double array of list.
	 */
	public static double[] toPrimitive(ArrayList<Double> list) {
		double[] arr = new double[list.size()];
		for (int i = 0; i < arr.length; i++)
			arr[i] = list.get(i);
		return arr;
	}

	/***********************************************************************/
	/* Small class here which stores Eigenvalue and Eigenvector */
	/***********************************************************************/

	/**
	 * This is a basic class that we use to store an Eigenvalue and an
	 * Eigenvector together as a pair.
	 */
	public static class EigenValueVectorPair {

		public double eigenvalue;
		public double[] eigenvector;

		/**
		 * Stores an eigenvalue and eigenvector pair.
		 * 
		 * @param eigenvalue
		 *            is the eigenvalue.
		 * @param eigenvector
		 *            is the corresponding eigenvector.
		 */
		public EigenValueVectorPair(double eigenvalue, double[] eigenvector) {
			this.eigenvalue = eigenvalue;
			this.eigenvector = eigenvector;
		}

	}

	/***********************************************************************/
	/* Getters and Setters below here */
	/***********************************************************************/

	/**
	 * Get matrix A that we perform Lanczos on.
	 * 
	 * @return A.
	 */
	public RealMatrix getA() {
		return A;
	}

	/**
	 * Set the matrix A
	 * 
	 * @throws IllegalArgumentException
	 *             if the matrix is not square.
	 * @param A
	 *            is the matrix we perform lanczos on.
	 */
	public void setA(RealMatrix A) {
		if (A == null)
			throw new IllegalArgumentException(
					"The matrix cannot be set as null.");
		if (!A.isSquare())
			throw new IllegalArgumentException("The matrix must be square");
		if (A.getColumnDimension() < 2)
			throw new IllegalArgumentException(
					"The matrix must have a size greater than 1.");
		this.A = A;
	}

	/**
	 * Get the max number of Lanczos iterations
	 * 
	 * @return m, the max number of Lanczos iterations.
	 */
	public int getMaxIterations() {
		return maxIterations;
	}

	/**
	 * Set the max number of Lanczos iterations.
	 * 
	 * @throws IllegalStateException
	 *             if A is null.
	 * @param maxIterations
	 *            of the Lanczos algorithm.
	 */
	public void setMaxIterations(int maxIterations) {
		if (A == null)
			throw new IllegalStateException(
					"A is equal to null, we must set A before setting anything else.");
		if (maxIterations <= 1)
			throw new IllegalArgumentException(
					"Please use a value of maxIterations which is greater than 1.");
		if (maxIterations > this.A.getColumnDimension())
			maxIterations = this.A.getColumnDimension();
		this.maxIterations = maxIterations;
	}

	/**
	 * @return the index of the eigenvector/eigenvalue that we are getting.
	 */
	public int getEigenIndex() {
		return eigenIndex;
	}

	/**
	 * Sets the index of the eigenvalue/eigenvector that we are getting.
	 * 
	 * @param eigenIndex
	 *            is the index of the eigenvalue we are retrieving.
	 */
	public void setEigenIndex(int eigenIndex) {
		if (eigenIndex < 0 || eigenIndex >= A.getColumnDimension())
			throw new IllegalStateException(
					"That index is outside the range of the matrix.");
		this.eigenIndex = eigenIndex;
	}

	/**
	 * Get the boolean which describes how we compare values between Lanczos
	 * iterations.
	 * 
	 * @return boolean compareEigenvectors.
	 */
	public boolean isCompareEigenvectors() {
		return compareEigenvectors;
	}

	/**
	 * Set the field compareEigenvectos.
	 * 
	 * @param compareEigenvectors
	 *            is true if we compare eigenvectors, false if we compare
	 *            eigenvalues.
	 */
	public void setCompareEigenvectors(boolean compareEigenvectors) {
		this.compareEigenvectors = compareEigenvectors;
	}

	/**
	 * Get the accuracy which once reached we stop the algorithm.
	 * 
	 * @return the accuracy which halts the algorithm.
	 */
	public double getAccuracy() {
		return accuracy;
	}

	/**
	 * Set the accuracy which stops the algorithm.
	 * 
	 * @param accuracy
	 *            the new accuracy which stops the algorithm.
	 */
	public void setAccuracy(double accuracy) {
		if (accuracy < 0)
			throw new IllegalArgumentException(
					"Please use a non negative accuracy value.");
		this.accuracy = accuracy;
	}

}
