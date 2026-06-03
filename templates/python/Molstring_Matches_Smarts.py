from rdkit import Chem
from typing import Optional
from functools import lru_cache

try:
    from Util import *
except ImportError:
    # module is inlined, ignore the import error
    pass


@lru_cache(maxsize=128)
def get_pattern_mol(smarts: str) -> Optional[Chem.Mol]:
    m = Chem.MolFromSmarts(smarts)
    if not m:
       raise ValueError(f'Error parsing SMARTS: {smarts}')
    if m.GetNumAtoms() == 0:
        return None
    return m

    
def is_molfile(molstring: Optional[str]) -> bool:
    if not molstring:
        return False
    return '\n' in molstring

    
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
