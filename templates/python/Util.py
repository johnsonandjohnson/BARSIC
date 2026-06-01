from typing import Callable
from rdkit.Chem.rdchem import MolSanitizeException

def safe_call_decorator(func: Callable, exception_type=MolSanitizeException):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except exception_type:
            return None
    return wrapper
