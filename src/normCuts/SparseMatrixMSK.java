package normCuts;

/**
 * Class which allows us to store SparseMatrix, specifically catered for our
 * implementation of NormalizedCuts.
 * 
 * Notably, in our implementation of NormalizedCuts, EVERY row and column has
 * one or more entry, we need it to be quick for traversing elements in the
 * matrix more than anything.
 * 
 * TODO : Implement this in the future!
 * 
 * @author Andrew P Talbot
 *
 */
public class SparseMatrixMSK {
	
	protected Header head;

	/**
	 * Construct a SparseMatrixMSK object,
	 */
	public SparseMatrixMSK() {

	}

	/**
	 * Header class which stores information about its associated row/column.
	 * 
	 * @author Andrew P Talbot
	 */
	private static class Header extends MatrixEntity {
		
		public long row_col;
		
		public Header(MatrixEntity right, MatrixEntity down, long row_col) {
			super(right, down);
			this.row_col = row_col;
		}
		
	}

	/**
	 * Stores information about an Element in the SparseMatrix, holds pointer to
	 * right and down Elements.
	 * 
	 * @author Andrew P Talbot
	 */
	private static class Element extends MatrixEntity {

		public double value;
		public long row;
		public long col;

		public Element(MatrixEntity right, MatrixEntity down, double value, long row, long col) {
			super(right, down);
			this.value = value;
			this.row = row;
			this.col = col;
		}

	}

	/**
	 * Superclass which we write so that Header and Element classes can point to
	 * either one or the other.
	 * 
	 * @author Andrew P Talbot
	 */
	private static class MatrixEntity {

		public MatrixEntity right;
		public MatrixEntity down;
		
		public MatrixEntity(MatrixEntity right, MatrixEntity down) {
			this.right = right;
			this.down = down;
		}

	}

}
