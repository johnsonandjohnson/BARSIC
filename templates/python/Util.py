from typing import Callable
from rdkit.Chem.rdchem import MolSanitizeException

def safe_call_decorator(func: Callable):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except MolSanitizeException:
            return None
        except RuntimeError:
            return None
    return wrapper

