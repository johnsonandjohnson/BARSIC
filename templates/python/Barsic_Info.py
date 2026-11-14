import rdkit
import sys
import snowflake.snowpark

def barsic_info() -> str:
    return f'BARSIC 3.1, Python: {sys.version}, RDKit: {rdkit.__version__}, snowflake-snowpark-python: {snowflake.snowpark.__version__}'

