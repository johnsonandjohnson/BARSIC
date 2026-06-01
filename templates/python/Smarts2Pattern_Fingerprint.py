from rdkit import Chem, DataStructs
from typing import Optional
import threading

try:
    from Util import *
except ImportError:
    # module is inlined, ignore the import error
    pass


_lck = threading.RLock()
_prev_smarts = None
_prev_fp = None

def smarts2pattern_fingerprint(smarts: Optional[str]) -> Optional[bytes]:
    global _prev_smarts
    global _prev_fp
    if not smarts:
        return None
    with _lck:
        if _prev_smarts == smarts:
            return _prev_fp
        m = Chem.MolFromSmarts(smarts)
        if not m:
            raise ValueError(f'Error parsing SMARTS: {smarts}')
        _prev_smarts = smarts
        if m.GetNumAtoms() == 0:
            _prev_fp = None
            return None
        fp = Chem.PatternFingerprint(m)        
        _prev_fp = DataStructs.BitVectToBinaryText(fp)
        return _prev_fp
