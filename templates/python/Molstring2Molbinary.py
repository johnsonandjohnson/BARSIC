from rdkit import Chem
from typing import Optional

try:
    from Util import *
except ImportError:
    # module is inlined, ignore the import error
    pass


def is_molfile(molstring: Optional[str]) -> bool:
    if not molstring:
        return False
    return '\n' in molstring

def molstring_to_binary(molstring: Optional[str])->Optional[bytes]:
    if molstring is None:
        return None
    m = Chem.MolFromMolBlock(molstring) if is_molfile(molstring) else Chem.MolFromSmiles(molstring)
    if not m:
        return None
    
    if m.GetNumAtoms() == 0: # don't need empty mols
        return None

    return m.ToBinary()
