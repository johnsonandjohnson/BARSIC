from rdkit import Chem
from typing import Optional
from functools import lru_cache

try:
    from Util import *
except ImportError:
    # module is inlined, ignore the import error
    pass


def molstring_matches_smarts(molstring: Optional[str], smarts: Optional[str], screen_pass: Optional[bool]) -> bool:
    if not screen_pass:
        return False
    if molstring is None or smarts is None:
        return False
    p = get_pattern_mol(smarts)
    if not p:
        return False
    m = Chem.MolFromMolBlock(molstring) if is_molfile(molstring) else Chem.MolFromSmiles(molstring)
    if not m:
        return False
    return m.HasSubstructMatch(p)
