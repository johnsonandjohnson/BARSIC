from functools import lru_cache
from typing import Callable, Optional

from rdkit import Chem
from rdkit.Chem.SaltRemover import InputFormat
from rdkit.Chem.rdchem import MolSanitizeException
from rdkit.Chem import SaltRemover


def safe_call_decorator(func: Callable):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except MolSanitizeException:
            return None
        except RuntimeError:
            return None
    return wrapper


# todo: auto-detect not only molfile/sdf and SMILES, but also other
# encodings (InChi, etc.).
def is_molfile(molstring: Optional[str]) -> bool:
    if not molstring:
        return False
    return '\n' in molstring


@lru_cache(128)
def getmol(molstring: Optional[str], molbinary: Optional[bytes]):
    if not molstring and not molbinary:
        return None
    if molstring and molbinary:
        raise ValueError('Either molstring or molbinary can be not NULL, but not both')
    if molbinary:
        m = Chem.Mol(molbinary)
    else:
        m = Chem.MolFromMolBlock(molstring) if is_molfile(molstring) else Chem.MolFromSmiles(molstring)
    if m and m.GetNumAtoms() == 0:
        return None
    return m


@lru_cache(maxsize=128)
def get_salt_remover(additional_patterns: str | None):
    if not additional_patterns:
        return SaltRemover.SaltRemover()
    sr0 = SaltRemover.SaltRemover()
    additional_patterns = additional_patterns.replace('|', '\n')
    sr1 = SaltRemover.SaltRemover(defnData=additional_patterns, defnFormat=InputFormat.SMARTS)
    # hack
    sr0.salts = sr0.salts + sr1.salts
    return sr0


@lru_cache(maxsize=128)
def get_pattern_mol(smarts: Optional[str]) -> Optional[Chem.Mol]:
    if not smarts:
        return None
    m = Chem.MolFromSmarts(smarts)
    if not m:
       raise ValueError(f'Error parsing SMARTS: {smarts}')
    if m.GetNumAtoms() == 0:
        return None
    return m

