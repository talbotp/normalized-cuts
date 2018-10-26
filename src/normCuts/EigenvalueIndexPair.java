package normCuts;

/**
 * Trivial class that we create to store the index of a eigenvalue so it can be
 * properly sorted using Collections.
 * 
 * @author P Talbot
 */
public class EigenvalueIndexPair implements Comparable<EigenvalueIndexPair> {

	public int index;
	public double eigenvalue;

	/**
	 * Object to store an index and its eigenvalue so that we properly sort it.
	 * 
	 * @param index
	 *            is the index of the eigenvalue in the array.
	 * @param eigenvalue
	 *            is the eigenvalue at index given.
	 */
	public EigenvalueIndexPair(int index, double eigenvalue) {
		this.index = index;
		this.eigenvalue = eigenvalue;
	}

	/**
	 * 
	 * @param other
	 *            is the EigenvalueIndexPair that we compare this to.
	 * @return 0 if this eigenvalue is equal to other eigenvalue, -1 if this
	 *         eigenvalue is less than other eigenvalue, 1 otherwise.
	 */
	@Override
	public int compareTo(EigenvalueIndexPair other) {
		double diff = this.eigenvalue - other.eigenvalue;
		if (diff == 0)
			return 0;
		return (diff < 0) ? -1 : 1;
	}

	/**
	 * @return the eigenvalue in string form.
	 */
	@Override
	public String toString() {
		return "" + eigenvalue;
	}

}
