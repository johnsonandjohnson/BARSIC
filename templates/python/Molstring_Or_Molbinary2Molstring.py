import sys
import re
import hashlib
from dataclasses import dataclass
from functools import lru_cache
import argparse
from enum import Enum, auto
from rdkit import Chem
from rdkit.Chem import SaltRemover
from typing import Optional

# future work: auto-detect not only molfile/sdf and SMILES, but also other
# encodings (InChi, etc.)
def is_molfile(molstring: Optional[str]) -> bool:
    if not molstring:
        return False
    return '\n' in molstring

class MolEnc(Enum):
    MOLBLOCK = auto()
    SMILES = auto()
    MOLHASH = auto()
    MOLHASHFULL = auto()
    INCHI = auto()
    INCHIKEY = auto()


@dataclass
class M2MOptions:
    out_enc: MolEnc
    desalt: bool
    remove_stereo: bool
    smiles_kekule: bool
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
    argp = ArgumentParserX('Convert molstring (SMILES or MOLBLOCK (molfile, sdf) to the specified text-based '
                          'format, with additional options')
    group0 = argp.add_argument_group('Output Encoding')
    group0.add_argument('--out', required = True, help='The output molecule encoding '
                                    '((smiles|smi)|(molfile|molblock)|molhash|molhashfull|inchi|inchikey)')

    group1 = argp.add_argument_group('Transform')
    group1.add_argument('--desalt', action='store_true', default=False, help='Remove (strip) salt')
    group1.add_argument('--remove_stereo', action='store_true', default=False, help='Remove stereo')

    group2 = argp.add_argument_group('SMILES encoder options')
    group2.add_argument('--smiles_kekule', action='store_true', default=False,
                        help='Kekulize the molecule before generating the SMILES and output single/double '
                             'rather than aromatic bonds')

    args = argp.parse_args(option_str.split())

    match args.out:
        case 'molfile' | 'molblock':
            enc = MolEnc.MOLBLOCK
        case 'smiles' | 'smi':
            enc = MolEnc.SMILES
        case 'molhash':
            enc = MolEnc.MOLHASH
        case 'molhashfull':
            enc = MolEnc.MOLHASHFULL            
        case 'inchi':
            enc = MolEnc.INCHI
        case 'inchikey':
            enc = MolEnc.INCHIKEY
        case _:
            sys.tracebacklimit = 0
            raise ValueError('Invalid/unknown --out parameter, must be one of '
                             '((smiles|smi)|(molfile|molblock)|molhash|molhashfull|inchi|inchikey)')
    return M2MOptions(out_enc=enc, desalt=args.desalt, remove_stereo=args.remove_stereo,
                      smiles_kekule=args.smiles_kekule)

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

def get_hashstring(s: str) -> str:
    h = hashlib.sha1()
    h.update(s.encode())
    return h.hexdigest()

# Removing data sgroups from extended SMILES strings is much easier and faster than removing
# those groups from molecules before encoding them into extended SMILES.
# See: https://docs.chemaxon.com/display/docs/formats_chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.md
def remove_data_sgroups(smiles: Optional[str]) -> Optional[str]:
    if not smiles:
        return smiles
    # Don't touch stereolabels. In some known use cases, stereolabels are used to resolve ambiguities
    # in the extended relative stereo notation.
    data_sgroup_rx = r'SgD:(\d|,)*:(?!stereolabel)[^:]*:[^:]*:[^:]*:[^:]*:[^:]*:[^,|:]*,?'
    # If we ended up with trailing ' ||', remove 3 last characters.
    s = re.sub(data_sgroup_rx, '', smiles)
    if s.endswith(' ||'):
        return s[:-3]
    return s

def molstring_or_molbinary_to_molstring_internal(molstring: Optional[str], molbinary: Optional[bytes], option_str: str):
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

    match opt.out_enc:
        case MolEnc.SMILES:
            p = Chem.SmilesWriteParams()
            p.canonical = True
            p.doKekule = opt.smiles_kekule
            f = Chem.rdmolfiles.CXSmilesFields.CX_ALL_BUT_COORDS
            # Note that we are not removing data sgroups from extended canonical SMILES,
            # but we do remove them when computing the hash (see below).
            return Chem.MolToCXSmiles(m, p, f)
        case MolEnc.MOLBLOCK:
            return Chem.MolToMolBlock(m)
        case MolEnc.MOLHASH:
            # Note: RegistrationHash.GetMolLayers(m) and RegistrationHash.GetMolHash(...)
            # are unnecessarily complex, just use the canonical SMILES.
            # Add tautomer hash (stereo/non-stereo) if requested later...
            p = Chem.SmilesWriteParams()
            p.canonical = True
            f = Chem.rdmolfiles.CXSmilesFields.CX_ALL_BUT_COORDS
            return get_hashstring(remove_data_sgroups(Chem.MolToCXSmiles(m, p, f)))
        case MolEnc.MOLHASHFULL:  # take everything into account, including the coordinates and data sgroups
            p = Chem.SmilesWriteParams()
            p.canonical = True
            f = Chem.rdmolfiles.CXSmilesFields.CX_ALL
            return get_hashstring(Chem.MolToCXSmiles(m, p, f))            
        case MolEnc.INCHI:
            return Chem.MolToInchi(m)
        case MolEnc.INCHIKEY:
            return Chem.MolToInchiKey(m)
        case _:
            raise ValueError('Invalid/unknown output encoding')
    # return is not needed, but keeps Sonar happy
    return None

def molstring_or_molbinary_to_molstring(molstring: Optional[str], molbinary: Optional[bytes], option_str: str):
    try:
        return molstring_or_molbinary_to_molstring_internal(molstring, molbinary, option_str)
    except:
        return None
