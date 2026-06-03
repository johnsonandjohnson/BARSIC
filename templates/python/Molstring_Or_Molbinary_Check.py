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

class MolChecker:
    @staticmethod
    def check_str(molstring: str, raise_exception: bool) -> tuple[bool|None, str|None, str|None]:
        if is_molfile(molstring):
            fmt = 'MOLBLOCK'
            m = Chem.MolFromMolBlock(molstring, sanitize=False)
        else:
            fmt = 'SMILES'
            m = Chem.MolFromSmiles(molstring, sanitize=False)

        if not m:
            if raise_exception:
                raise Exception(f'{fmt} parser error')
            return False, fmt, 'Parser error'

        try:
            Chem.SanitizeMol(m)
        except Exception as ex:
            if raise_exception:
                raise
            return False, fmt, str(ex)

        return True, fmt, None

    @staticmethod
    def check_bin(molbinary: Optional[bytes], raise_exception: bool) -> tuple[bool|None, str|None, str|None]:
        fmt = 'BINARY'
        try:
            m = Chem.Mol(molbinary)
            Chem.SanitizeMol(m)
        except Exception as ex:
            if raise_exception:
                raise
            return False, fmt, str(ex)
        return True, fmt, None

    @staticmethod
    def check(molstring: Optional[str], molbinary: Optional[bytes], raise_exception: bool) -> tuple[bool|None, str|None, str|None]:
        if not molstring and not molbinary:
            return None, None, None
        if molstring and molbinary:
            raise ValueError('Either molstring or molbinary can be not NULL, but not both')

        if molstring:
            return MolChecker.check_str(molstring, raise_exception)
        return MolChecker.check_bin(molbinary, raise_exception)


    def process(self, molstring: Optional[str], molbinary: Optional[bytes], raise_exception: bool):
        yield self.check(molstring, molbinary, raise_exception)
