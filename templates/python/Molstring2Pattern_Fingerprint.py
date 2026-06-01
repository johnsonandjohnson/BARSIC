from rdkit import Chem, DataStructs
from typing import Optional
import threading

try:
    from Util import *
except ImportError:
    # module is inlined, ignore the import error
    pass


_lck = threading.RLock()
_prev_molstring = None
_prev_fp = None


def is_molfile(molstring: Optional[str]) -> bool:
    if not molstring:
        return False
    return '\n' in molstring


def molstring2pattern_fingerprint(molstring: Optional[str]) -> Optional[bytes]:
    global _prev_molstring
    global _prev_fp
    if not molstring:
        return None
    with _lck:
        if molstring == _prev_molstring:
            return _prev_fp
        _prev_molstring = molstring
        m = Chem.MolFromMolBlock(molstring) if is_molfile(molstring) else Chem.MolFromSmiles(molstring)
        if not m or m.GetNumAtoms() == 0:
            _prev_fp = None
            return None
        fp = Chem.PatternFingerprint(m)
        _prev_fp = DataStructs.BitVectToBinaryText(fp)
        return _prev_fp
