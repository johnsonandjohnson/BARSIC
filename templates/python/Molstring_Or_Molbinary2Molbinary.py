import sys
import hashlib
from dataclasses import dataclass
from functools import lru_cache
import argparse
from rdkit import Chem
from rdkit.Chem import SaltRemover
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


@dataclass
class M2MOptions:
    desalt: bool
    remove_stereo: bool
    # options to be added as needed

class ArgumentParserX(argparse.ArgumentParser):
    def __init__(self, prog=None):
        super().__init__(prog, add_help=True, exit_on_error=False)

    def error(self, message):
        sys.tracebacklimit = 0
        raise ValueError(message) from None

    def exit(self, status = 0, message = None):
        sys.tracebacklimit = 0
        raise ValueError(message) from None

    def print_usage(self, file=None):
        sys.tracebacklimit = 0
        raise ValueError('Usage: ' + self.format_usage()) from  None

    def print_help(self, file=None):
        sys.tracebacklimit = 0
        raise ValueError('Help: ' + self.format_help()) from None


@lru_cache(maxsize=128)
def parse_options(option_str: str) -> M2MOptions:
    argp = ArgumentParserX('Convert molstring (SMILES or MOLBLOCK (molfile, sdf) to the RDKit binary encoding '
                          'format, with additional options')

    group1 = argp.add_argument_group('Transform')
    # future work: standardize tautomers w/options.
    group1.add_argument('--desalt', action='store_true', default=False, help='Remove (strip) salt')
    group1.add_argument('--remove_stereo', action='store_true', default=False, help='Remove stereo')


    args = argp.parse_args(option_str.split())

    return M2MOptions(desalt=args.desalt, remove_stereo=args.remove_stereo)

_salt_remover = SaltRemover.SaltRemover()

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


@safe_call_decorator
def molstring_or_molbinary_to_molbinary(molstring: Optional[str], molbinary: Optional[bytes], option_str: str):
    # parse options first and show usage help if option_str has --help or -h flags
    opt = parse_options(option_str)
    m = getmol(molstring, molbinary)
    if not m:
        return None
    is_clone = False

    # we need to clone molecules before we change them so we won't change molecules in the getmol lru cache
    def clone_if_needed():
        nonlocal is_clone
        nonlocal m
        if not is_clone:
            # Create a deep copy of the original molecule to avoid modifying the cached
            # instance returned by getmol. the overhead of copying is insignificant
            m = Chem.Mol(m, quickCopy=True)
            is_clone = True

    if opt.desalt:
        # note that StripMol returns a new Mol instance and does not change the molecule
        # passed to the function, so we don't need to clone it
        m = _salt_remover.StripMol(m, dontRemoveEverything=True)
    if opt.remove_stereo:
        clone_if_needed()
        Chem.RemoveStereochemistry(m)
        Chem.ClearMolSubstanceGroups(m)
        
    if m.GetNumAtoms() == 0: # don't need empty mols
        return None
        
    return m.ToBinary()

