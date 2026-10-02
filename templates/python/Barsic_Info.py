import rdkit
import sys

def barsic_info() -> str:
    return f'BARSIC 2.0, Python: {sys.version}, RDKit: {rdkit.__version__}'

