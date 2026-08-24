# NOTE: due to problems with the bitarray package, this function is currently not used.
# ChemCoreApi.tpl.sql now defines two identical functions, Tanimoto and Tanimoto_J, 
# both using the Java implementation.
from typing import Optional
from bitarray import bitarray
def tanimoto(b1: Optional[bytes], b2: Optional[bytes])->Optional[float]:
    if not b1 or not b2:
        return None
    b1b = bitarray()
    b2b = bitarray()
    b1b.frombytes(b1)
    b2b.frombytes(b2)
    if len(b1b) != len(b2b):
        raise ValueError('Bitsets are of different lengths')
    ba_and = b1b & b2b  #  intersection (common bits)
    and_count = ba_and.count(1)
    union_cnt = b1b.count(1) + b2b.count(1) - and_count  # could also use (b1b | b2b).count(1), but it would be slower
    if union_cnt == 0:
        return 1.0  # can only happen if both b1 and b2 have no bits set to 1
    return and_count / float(union_cnt)
