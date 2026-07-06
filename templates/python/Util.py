from functools import lru_cache
from typing import Callable

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

