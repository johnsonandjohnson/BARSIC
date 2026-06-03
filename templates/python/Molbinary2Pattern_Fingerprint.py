from rdkit import Chem, DataStructs
from typing import Optional
import threading

try:
    from Util import *
except ImportError:
    # module is inlined, ignore the import error
    pass


_lck = threading.RLock()
_prev_bin = None
_prev_fp = None


@safe_call_decorator
def molbinary2pattern_fingerprint(molbinary: Optional[bytes]) -> Optional[bytes]:
    global _prev_bin
    global _prev_fp
    if not molbinary:
        return None
    with _lck:
        if molbinary == _prev_bin:
            return _prev_fp
        _prev_bin = molbinary
        m = Chem.Mol(molbinary)
        if not m or m.GetNumAtoms() == 0:
            _prev_fp = None
            return None
        fp = Chem.PatternFingerprint(m)
        _prev_fp = DataStructs.BitVectToBinaryText(fp)
        return _prev_fp
