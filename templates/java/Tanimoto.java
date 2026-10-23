import java.util.BitSet;
class Tanimoto 
{
    public static double calculate(byte[] v1, byte[] v2) 
    {
        // Note, we don't have to check the args for null,
        // because this UDF is declared with RETURNS NULL ON NULL INPUT,
        // so it won't be called if either of both v1 and v1 are null's.
        if (v1.length != v2.length)
        {
            String msg = String.format("Vectors representing bitsets must be of equal size. v1.length: %d, v2.length: %d.", 
                                v1.length, v2.length);
            throw new IllegalArgumentException(msg);
        }
        BitSet bitset1 = BitSet.valueOf(v1); 
        BitSet bitset2 = BitSet.valueOf(v2);
        // Create a copy of bitset1 to find the intersection
        BitSet intersection = (BitSet) bitset1.clone();
        intersection.and(bitset2);
    
        int nIntersection = intersection.cardinality(); // Number of "on" bits in intersection
        int nA = bitset1.cardinality(); // Number of "on" bits in bitset1
        int nB = bitset2.cardinality(); // Number of "on" bits in bitset2
    
        if (nA + nB - nIntersection == 0) 
        {
            return 0.0; // Avoid division by zero if both sets are empty.
        }
    
        return (double) nIntersection / (nA + nB - nIntersection);
    }
};
