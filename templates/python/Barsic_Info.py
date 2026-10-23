import rdkit
import sys

def barsic_info() -> str:
    return f'BARSIC 3.0, Python: {sys.version}, RDKit: {rdkit.__version__}'

