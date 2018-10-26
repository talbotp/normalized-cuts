package normCuts;

import org.apache.commons.math3.linear.RealMatrix;

import java.util.Comparator;
import java.util.Map;

/**
 * Trivial Class to store a matrix with its original indices Mapping.
 *
 * We implement Comparator as we wish to sort an array of BiPartition and this
 * is NOT a natural ordering, as we sort so that the biggest matrix will be
 * first.
 *
 * @author P Talbot
 */
public class BiPartition implements Comparator<BiPartition> {

    public RealMatrix W;
    public Map<Integer, Integer> partition;

    /**
     * Create a Bipartition consisting of a similarity matrix W. and a map
     * containing partition values.
     *
     * @param W         is the similarity matrix to be partitioned.
     * @param partition is the map from points in W to partition A or B.
     */
    public BiPartition(RealMatrix W, Map<Integer, Integer> partition) {
        this.W = W;
        this.partition = partition;
    }

    /**
     * Not a natural ordering of the object.
     * Sort such that matrix with the largest size is first.
     *
     * @param o1 is Bipartition to be compared against o2.
     * @param o2 is Bipartition to be compared against o1.
     * @return integer x. x < 0 if o1 is bigger than 02, x == 0 if o1 and o2 are
     *         equal size, x > 0 otherwise.
     */
    @Override
    public int compare(BiPartition o1, BiPartition o2) {
        return o2.partition.size() - o1.partition.size();
    }

}
