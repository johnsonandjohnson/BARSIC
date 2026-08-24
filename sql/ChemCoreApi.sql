-- Set db/schema in which you want these UDF's to be defined.
-- It is recommended that the schema be called chem_api.
--USE <database_name_here>;
CREATE SCHEMA IF NOT EXISTS chem_api;
USE schema chem_api;


-- Chemical structure conversion UDF's ---------------------------

CREATE OR REPLACE FUNCTION Molstring_Or_Molbinary2Molstring(molstring VARCHAR DEFAULT NULL, molbinary VARBINARY DEFAULT NULL, options VARCHAR DEFAULT '--out smiles')
     RETURNS VARCHAR 
     LANGUAGE PYTHON 
     IMMUTABLE    
     RUNTIME_VERSION = '3.11' 
     PACKAGES = ('rdkit')
     HANDLER = 'molstring_or_molbinary_to_molstring'
     COMMENT='Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) or RDKit binary molecule encoding to a variety of different string-based formats, optionally transforming the input structure. Only one of molstring or molbinary arguments can be non-NULL. All available options can be listed by by running the following SQL: select Molstring_Or_Molbinary2Molstring(NULL, NULL, ''-h''); usage info will be returned as part of the error message. If the options argument value is not specified, computes canonical ChemAxon-compatible extended SMILES.'
     AS
$$
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

import enum
import sys
import re
import hashlib
from dataclasses import dataclass
from functools import lru_cache
import argparse
from enum import Enum, auto
from rdkit import Chem
from rdkit.Chem import SaltRemover
from rdkit.Chem import RegistrationHash as rh
from typing import Optional

try:
    from Util import *
except ImportError:
    # module is inlined, ignore the import error
    pass


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
    TAUTOHASH = auto()
    INCHI = auto()
    INCHIKEY = auto()


@dataclass
class M2MOptions:
    out_enc: MolEnc
    desalt: bool
    additional_desalt_patterns: str | None
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
                                    '((smiles|smi)|(molfile|molblock)|molhash|molhashfull|tautohash|inchi|inchikey)')

    group1 = argp.add_argument_group('Transform')
    group1.add_argument('--desalt', action='store_true', default=False, help='Remove (strip) salt')
    group1.add_argument('--desalt_smarts_list', required=False, default=None, help='Additional desalt SMARTS patterns separated with |')
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
        case 'tautohash':
            enc = MolEnc.TAUTOHASH
        case 'inchi':
            enc = MolEnc.INCHI
        case 'inchikey':
            enc = MolEnc.INCHIKEY
        case _:
            sys.tracebacklimit = 0
            raise ValueError('Invalid/unknown --out parameter, must be one of '
                             '((smiles|smi)|(molfile|molblock)|molhash|molhashfull|tautohash|inchi|inchikey)')
    return M2MOptions(out_enc=enc, desalt=args.desalt,
                      additional_desalt_patterns=args.desalt_smarts_list,
                      remove_stereo=args.remove_stereo,
                      smiles_kekule=args.smiles_kekule)


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

_DATA_SGROUP_RX = re.compile(r'SgD:(\d|,)*:(?!stereolabel)[^:]*:[^:]*:[^:]*:[^:]*:[^:]*:[^,|:]*,?')

# Removing data sgroups from extended SMILES strings is much easier and faster than removing
# those groups from molecules before encoding them into extended SMILES.
# See: https://docs.chemaxon.com/display/docs/formats_chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.md
def remove_data_sgroups(smiles: Optional[str]) -> Optional[str]:
    if not smiles:
        return smiles
    # Don't touch stereolabels. In some known use cases, stereolabels are used to resolve ambiguities
    # in the extended relative stereo notation.
    # If we ended up with trailing ' ||', remove 3 last characters.
    s = _DATA_SGROUP_RX.sub('', smiles)
    if s.endswith(' ||'):
        return s[:-3]
    return s


# for some reason, the 'original' HashScheme does not have these options:
@enum.unique
class HashSchemeX(enum.Enum):
    EXACT_LAYERS = (
        rh.HashLayer.CANONICAL_SMILES, # trailing commas are important! This must be a tuple
    )

    STEREO_INSENSITIVE_EXACT_LAYERS = (
        rh.HashLayer.NO_STEREO_SMILES,
    )

    STEREO_INSENSITIVE_TAUTOMER_INSENSITIVE_LAYERS = (
        rh.HashLayer.NO_STEREO_TAUTOMER_HASH,
    )

    TAUTOMER_INSENSITIVE_LAYERS = (
        rh.HashLayer.TAUTOMER_HASH,
    )


@lru_cache(128)
@safe_call_decorator
def mol_to_reg_layers(m: Optional[Chem.Mol]) -> Optional[dict]:
    if m is None:
        return None
    return rh.GetMolLayers(m, enable_tautomer_hash_v2=True)  # todo: find out how it is implemented internally


@lru_cache(128)
@safe_call_decorator
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
        m = get_salt_remover(opt.additional_desalt_patterns).StripMol(m, dontRemoveEverything=True)

    # we don't need to remove stereo for the tautomer hash,
    # because rh.GetMolLayers takes care of that internally
    if opt.remove_stereo and opt.out_enc != MolEnc.TAUTOHASH:
        clone_if_needed()
        Chem.RemoveStereochemistry(m)
        Chem.ClearMolSubstanceGroups(m)

    # TODO: see if we can avoid calling GetMolLayers, because we need it only for the tautomeric hashes
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
        case MolEnc.TAUTOHASH:
            layers = mol_to_reg_layers(m)
            if not layers:
                return None
            scheme = HashSchemeX.STEREO_INSENSITIVE_TAUTOMER_INSENSITIVE_LAYERS if opt.remove_stereo \
                     else HashSchemeX.TAUTOMER_INSENSITIVE_LAYERS
            # ignore wrong type warning here
            return rh.GetMolHash(layers, hash_scheme=scheme)
        case MolEnc.INCHI:
            return Chem.MolToInchi(m)
        case MolEnc.INCHIKEY:
            return Chem.MolToInchiKey(m)
        case _:
            raise ValueError('Invalid/unknown output encoding')
    # return is not needed, but keeps Sonar happy
    return None


# need this extra layer because of the @lru_cache(128) and @safe_call_decorator used on the handler result in
# Python Interpreter Error: AttributeError: 'functools._lru_cache_wrapper' object has no attribute '__code__' error
def molstring_or_molbinary_to_molstring(molstring: Optional[str], molbinary: Optional[bytes], option_str: str):
    return molstring_or_molbinary_to_molstring_internal(molstring, molbinary, option_str)

$$
;


CREATE OR REPLACE FUNCTION Molstring_Or_Molbinary2Molbinary(molstring VARCHAR DEFAULT NULL, molbinary VARBINARY DEFAULT NULL, options VARCHAR DEFAULT '')
     RETURNS VARBINARY 
     LANGUAGE PYTHON 
     IMMUTABLE    
     RUNTIME_VERSION = '3.11' 
     PACKAGES = ('rdkit')
     HANDLER = 'molstring_or_molbinary_to_molbinary'
     COMMENT='Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) or RDKit binary molecule encoding to the RDKit binary molecule encoding, optionally transforming the input structure. Only one of molstring or molbinary arguments can be non-NULL. All available options can be listed by by running the following SQL: select Molstring_Or_Molbinary2Molbinary(NULL, NULL, ''-h''); usage info will be returned as part of the error message. To convert molstrings to RDKit binary molecule encoding w/o applying any transforms, use the Molstring2Molbinary function, which has fewer arguments and is faster.'
     AS
$$
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

import sys
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
    additional_desalt_patterns: str | None
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
    group1.add_argument('--desalt_smarts_list', required=False, default=None, help='Additional desalt SMARTS patterns separated with |')
    group1.add_argument('--remove_stereo', action='store_true', default=False, help='Remove stereo')


    args = argp.parse_args(option_str.split())

    return M2MOptions(desalt=args.desalt,
                      additional_desalt_patterns=args.desalt_smarts_list,
                      remove_stereo=args.remove_stereo)


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
        m = get_salt_remover(opt.additional_desalt_patterns).StripMol(m, dontRemoveEverything=True)
    if opt.remove_stereo:
        clone_if_needed()
        Chem.RemoveStereochemistry(m)
        Chem.ClearMolSubstanceGroups(m)
        
    if m.GetNumAtoms() == 0: # don't need empty mols
        return None
        
    return m.ToBinary()

$$
;


CREATE OR REPLACE FUNCTION Molstring2Molbinary(molstring VARCHAR)
     RETURNS VARBINARY 
     LANGUAGE PYTHON 
     RETURNS NULL ON NULL INPUT
     IMMUTABLE    
     RUNTIME_VERSION = '3.11' 
     PACKAGES = ('rdkit')
     HANDLER = 'molstring_to_binary'
     COMMENT='Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) to RDKit binary molecule encoding. Returns NULL if molstring is NULL, empty, invalid, or represents an empty molecule with 0 atoms and 0 bonds.'
     AS
$$
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

def molstring_to_binary(molstring: Optional[str])->Optional[bytes]:
    if molstring is None:
        return None
    m = Chem.MolFromMolBlock(molstring) if is_molfile(molstring) else Chem.MolFromSmiles(molstring)
    if not m:
        return None
    
    if m.GetNumAtoms() == 0: # don't need empty mols
        return None

    return m.ToBinary()
$$
;

-- Fingerprinting UDF's ---------------------------

-- Fingerprints for substructure screening --------

CREATE OR REPLACE FUNCTION Molstring2Pattern_Fingerprint(molstring VARCHAR)
     RETURNS VARBINARY 
     LANGUAGE PYTHON 
     RETURNS NULL ON NULL INPUT
     IMMUTABLE
     RUNTIME_VERSION = '3.11' 
     PACKAGES = ('rdkit')
     HANDLER = 'molstring2pattern_fingerprint'
     COMMENT='Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) to RDKit substructure pattern fingerprint commonly used for fingerprint-based screening to speed up substructure searches. Returns NULL if molstring is NULL, empty, invalid, or represents an empty molecule with 0 atoms and 0 bonds.'
     AS
$$
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

from rdkit import Chem, DataStructs
from typing import Optional
import threading

try:
    from Util import *
except ImportError:
    # module is inlined, ignore the import error
    pass


_lck = threading.RLock()
_prev_molstring = None
_prev_fp = None


def is_molfile(molstring: Optional[str]) -> bool:
    if not molstring:
        return False
    return '\n' in molstring


def molstring2pattern_fingerprint(molstring: Optional[str]) -> Optional[bytes]:
    global _prev_molstring
    global _prev_fp
    if not molstring:
        return None
    with _lck:
        if molstring == _prev_molstring:
            return _prev_fp
        _prev_molstring = molstring
        m = Chem.MolFromMolBlock(molstring) if is_molfile(molstring) else Chem.MolFromSmiles(molstring)
        if not m or m.GetNumAtoms() == 0:
            _prev_fp = None
            return None
        fp = Chem.PatternFingerprint(m)
        _prev_fp = DataStructs.BitVectToBinaryText(fp)
        return _prev_fp
$$
;


CREATE OR REPLACE FUNCTION Molbinary2Pattern_Fingerprint(molbinary VARBINARY)
     RETURNS VARBINARY 
     LANGUAGE PYTHON 
     RETURNS NULL ON NULL INPUT
     IMMUTABLE
     RUNTIME_VERSION = '3.11' 
     PACKAGES = ('rdkit')
     HANDLER = 'molbinary2pattern_fingerprint'
     COMMENT='Converts RDKit binary-encoded molecule to RDKit substructure pattern fingerprint commonly used for fingerprint-based screening to speed up substructure searches. Returns NULL if molbinary is NULL, empty or represents an empty molecule with 0 atoms and 0 bonds. Returns error if molbinary is not a valid RDKit binary-encoded molecule.'     
     AS
$$
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

from rdkit import Chem, DataStructs
from typing import Optional
import threading

try:
    from Util import *
except ImportError:
    # module is inlined, ignore the import error
    pass


_lck = threading.RLock()
_prev_bin = None
_prev_fp = None


@safe_call_decorator
def molbinary2pattern_fingerprint(molbinary: Optional[bytes]) -> Optional[bytes]:
    global _prev_bin
    global _prev_fp
    if not molbinary:
        return None
    with _lck:
        if molbinary == _prev_bin:
            return _prev_fp
        _prev_bin = molbinary
        m = Chem.Mol(molbinary)
        if not m or m.GetNumAtoms() == 0:
            _prev_fp = None
            return None
        fp = Chem.PatternFingerprint(m)
        _prev_fp = DataStructs.BitVectToBinaryText(fp)
        return _prev_fp
$$
;


CREATE OR REPLACE FUNCTION Smarts2Pattern_Fingerprint(smarts VARCHAR)
     RETURNS VARBINARY 
     LANGUAGE PYTHON 
     RETURNS NULL ON NULL INPUT
     IMMUTABLE
     RUNTIME_VERSION = '3.11' 
     PACKAGES = ('rdkit')
     HANDLER = 'smarts2pattern_fingerprint'
     COMMENT='Converts SMARTS substructure pattern to RDKit substructure pattern fingerprint commonly used for fingerprint-based screening to speed up substructure searches. Returns NULL if smarts is NULL, empty, or represents an empty molecule with 0 atoms and 0 bonds. Returns error if smarts is invalid and cannot be parsed.'
     AS
$$
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

from rdkit import Chem, DataStructs
from typing import Optional
import threading

try:
    from Util import *
except ImportError:
    # module is inlined, ignore the import error
    pass


_lck = threading.RLock()
_prev_smarts = None
_prev_fp = None

def smarts2pattern_fingerprint(smarts: Optional[str]) -> Optional[bytes]:
    global _prev_smarts
    global _prev_fp
    if not smarts:
        return None
    with _lck:
        if _prev_smarts == smarts:
            return _prev_fp
        m = Chem.MolFromSmarts(smarts)
        if not m:
            raise ValueError(f'Error parsing SMARTS: {smarts}')
        _prev_smarts = smarts
        if m.GetNumAtoms() == 0:
            _prev_fp = None
            return None
        fp = Chem.PatternFingerprint(m)        
        _prev_fp = DataStructs.BitVectToBinaryText(fp)
        return _prev_fp
$$
;


-- Fingerprints for similarity search

CREATE OR REPLACE FUNCTION Molstring2Morgan_Fingerprint(molstring VARCHAR)
     RETURNS VARBINARY 
     LANGUAGE PYTHON 
     RETURNS NULL ON NULL INPUT
     IMMUTABLE
     RUNTIME_VERSION = '3.11' 
     PACKAGES = ('rdkit')
     HANDLER = 'molstring2morgan_fingerprint'
     COMMENT='Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) to RDKit Morgan fingerprint commonly used for fingerprint-based similarity search. Returns NULL if molstring is NULL, empty, invalid, or represents an empty molecule with 0 atoms and 0 bonds. Generator options: radius=2, fpSize=2048, atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen()'
     AS
$$
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

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from typing import Optional
import threading

try:
    from Util import *
except ImportError:
    # module is inlined, ignore the import error
    pass


_lck = threading.RLock()
_prev_molstring = None
_prev_fp = None
_fpg = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048,
                                                 atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen())

    
def is_molfile(molstring: Optional[str]) -> bool:
    if not molstring:
        return False
    return '\n' in molstring


def molstring2morgan_fingerprint(molstring: Optional[str]) -> Optional[bytes]:
    global _prev_molstring
    global _prev_fp
    if not molstring:
        return None
    with _lck:
        if molstring == _prev_molstring:
            return _prev_fp
        _prev_molstring = molstring
        m = Chem.MolFromMolBlock(molstring) if is_molfile(molstring) else Chem.MolFromSmiles(molstring)
        if not m or m.GetNumAtoms() == 0:
            _prev_fp = None
            return None
        fp = _fpg.GetFingerprint(m)
        _prev_fp = DataStructs.BitVectToBinaryText(fp)
        return _prev_fp
$$
;


CREATE OR REPLACE FUNCTION Molbinary2Morgan_Fingerprint(molbinary VARBINARY)
     RETURNS VARBINARY 
     LANGUAGE PYTHON 
     RETURNS NULL ON NULL INPUT
     IMMUTABLE
     RUNTIME_VERSION = '3.11' 
     PACKAGES = ('rdkit')
     HANDLER = 'molbinary2morgan_fingerprint'
     COMMENT='Converts RDKit binary-encoded molecule to RDKit Morgan fingerprint commonly used for fingerprint-based similarity search. Returns NULL if molbinary is NULL or empty, or represents an empty molecule with 0 atoms and 0 bonds. Generator options: radius=2, fpSize=2048, atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen(). Returns error if molbinary is not a valid RDKit binary-encoded molecule.'          
     AS
$$
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

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from typing import Optional
import threading

try:
    from Util import *
except ImportError:
    # module is inlined, ignore the import error
    pass

_lck = threading.RLock()
_prev_bin = None
_prev_fp = None
_fpg = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen())


@safe_call_decorator
def molbinary2morgan_fingerprint(molbinary: Optional[bytes]) -> Optional[bytes]:
    global _prev_bin
    global _prev_fp
    if not molbinary:
        return None
    with _lck:
        if molbinary == _prev_bin:
            return _prev_fp
        _prev_bin = molbinary
        m = Chem.Mol(molbinary)
        if not m:
            _prev_fp = None
            return None
        fp = _fpg.GetFingerprint(m)
        _prev_fp = DataStructs.BitVectToBinaryText(fp)
        return _prev_fp
$$
;


-- Substructure matching UDF's ---------------------------

CREATE OR REPLACE FUNCTION Molbinary_Matches_Smarts(molbinary VARBINARY, smarts VARCHAR, screen_pass BOOLEAN)
     RETURNS BOOLEAN 
     LANGUAGE PYTHON 
     RETURNS NULL ON NULL INPUT
     IMMUTABLE    
     RUNTIME_VERSION = '3.11' 
     PACKAGES = ('rdkit')
     HANDLER = 'molbinary_matches_smarts'
     COMMENT='Tests whether RDKit binary-encoded molecule matches the specified SMARTS pattern and returns TRUE iff it does and the screen_pass argument value is TRUE. The extra screen_pass argument is used for substructure query optimization based on fingerprint screening (see examples in the documentation and example workbooks). Returns NULL if any of the args are NULL. Returns error if molbinary is not a valid RDKit binary-encoded molecule of if smarts is invalid and cannot be parsed. Note: an empty SMARTS will not match any molecule, even an empty one. This seems to be illogical, since, in theory, a subgraph with 0 nodes and 0 edges must match any graph (or, at least, an empty one), but, in practice, this approach leads to fewer problems than the theoretically correct one. In RDKit itself, a molecule representing an empty pattern does not match anything either.'
     AS
$$
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

from rdkit import Chem
from typing import Optional
from functools import lru_cache

try:
    from Util import *
except ImportError:
    # module is inlined, ignore the import error
    pass


@lru_cache(maxsize=128)
def get_pattern_mol(smarts: str) -> Optional[Chem.Mol]:
    m = Chem.MolFromSmarts(smarts)
    if not m:
       raise ValueError(f'Error parsing SMARTS: {smarts}')
    if m.GetNumAtoms() == 0:
        return None
    return m

def molbinary_matches_smarts(molbinary: Optional[bytes], smarts: Optional[str], screen_pass: Optional[bool]) -> bool:
    if not screen_pass:
        return False
    if molbinary is None or smarts is None:
        return False
    p = get_pattern_mol(smarts)
    if not p:
        return False
    m = Chem.Mol(molbinary)
    if not m:
        return False
    return m.HasSubstructMatch(p)
$$
;


CREATE OR REPLACE FUNCTION Molstring_Matches_Smarts(molstring VARCHAR, smarts VARCHAR, screen_pass BOOLEAN)
     RETURNS BOOLEAN 
     LANGUAGE PYTHON 
     RETURNS NULL ON NULL INPUT
     IMMUTABLE    
     RUNTIME_VERSION = '3.11' 
     PACKAGES = ('rdkit')
     HANDLER = 'molstring_matches_smarts'
     COMMENT='Tests whether the molecule encoded as molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) matches the specified SMARTS pattern and returns TRUE iff it does and the screen_pass argument value is TRUE. The extra screen_pass argument is used for substructure query optimization based on fingerprint screening (see examples in the documentation and example workbooks). Returns NULL if any of the args are NULL. Returns FALSE if molstring is not a valid SMILES or molblock. Returns error if smarts is invalid and cannot be parsed. Note: an empty SMARTS will not match any molecule, even an empty one. This seems to be illogical, since, in theory, a subgraph with 0 nodes and 0 edges must match any graph (or, at least, an empty one), but, in practice, this approach leads to fewer problems than the theoretically correct one. In RDKit itself, a molecule representing an empty pattern does not match anything either.'     
     AS
$$
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

from rdkit import Chem
from typing import Optional
from functools import lru_cache

try:
    from Util import *
except ImportError:
    # module is inlined, ignore the import error
    pass


@lru_cache(maxsize=128)
def get_pattern_mol(smarts: str) -> Optional[Chem.Mol]:
    m = Chem.MolFromSmarts(smarts)
    if not m:
       raise ValueError(f'Error parsing SMARTS: {smarts}')
    if m.GetNumAtoms() == 0:
        return None
    return m

    
def is_molfile(molstring: Optional[str]) -> bool:
    if not molstring:
        return False
    return '\n' in molstring

    
def molstring_matches_smarts(molstring: Optional[str], smarts: Optional[str], screen_pass: Optional[bool]) -> bool:
    if not screen_pass:
        return False
    if molstring is None or smarts is None:
        return False
    p = get_pattern_mol(smarts)
    if not p:
        return False
    m = Chem.MolFromMolBlock(molstring) if is_molfile(molstring) else Chem.MolFromSmiles(molstring)
    if not m:
        return False
    return m.HasSubstructMatch(p)
$$
;


-- Similarity search UDF's ---------------------------

-- NOTE: due to problems with the bitarray package, this function is currently not used.
-- ChemCoreApi.tpl.sql now defines two identical functions, Tanimoto and Tanimoto_J,
-- both using the Java implementation.

--CREATE OR REPLACE FUNCTION Tanimoto(v1 VARBINARY, v2 VARBINARY)
--     RETURNS FLOAT
--     LANGUAGE PYTHON
--     RETURNS NULL ON NULL INPUT
--     IMMUTABLE
--     RUNTIME_VERSION = '3.11'
--     PACKAGES = ('bitarray==2.5.1')  -- note: the latest version in Snowflake (3.4.2 as of now) seems to be broken!
--     HANDLER = 'tanimoto'
--     COMMENT='Computes Tanimoto similarity between two bitsets represented by binary vectors. Returns BITCOUNT(v1 bitand v2) / BITCOUNT(v1 bitor v2) as float. If both v1 and v2 have no bits set to 1, returns 1.0, that is, treats two empty bitsets or two bitsets filled with only zeroes as being equal to each other. Returns error if two bitsets are of different lengths. Returns NULL if one of both arguments are NULLs'
--     AS
--$$
-- -- @donotinclude python/Tanimoto.py
-- $$
--;


-- A version of Tanimoto implemented in Java. Can be faster compared to the Python version.
CREATE STAGE IF NOT EXISTS java_handlers; --must be created once

CREATE OR REPLACE FUNCTION Tanimoto_J(v1 VARBINARY, v2 VARBINARY)
     RETURNS FLOAT
     LANGUAGE JAVA
     RETURNS NULL ON NULL INPUT
     IMMUTABLE
     HANDLER = 'Tanimoto.calculate'
     -- the Java code will be pre-compiled and stored in the jar file:
     TARGET_PATH = '@java_handlers/tanimoto.jar'
     COMMENT='A version of TANIMOTO implemented in Java. Computes Tanimoto similarity between two bitsets represented by binary vectors. Returns BITCOUNT(v1 bitand v2) / BITCOUNT(v1 bitor v2) as float. If both v1 and v2 have no bits set to 1, returns 1.0, that is, treats two empty bitsets or two bitsets filled with only zeroes as being equal to each other. Returns error if two bitsets are of different lengths. Returns NULL if one of both arguments are NULLs'
     AS
$$
import java.util.BitSet;
class Tanimoto 
{
    public static double calculate(byte[] v1, byte[] v2) 
    {
        // Note, we don't have to check the args for null,
        // because this UDF is declared with RETURNS NULL ON NULL INPUT,
        // so it won't be called if either of both v1 and v1 are null's.
        if (v1.length != v2.length)
        {
            String msg = String.format("Vectors representing bitsets must be of equal size. v1.length: %d, v2.length: %d.", 
                                v1.length, v2.length);
            throw new IllegalArgumentException(msg);
        }
        BitSet bitset1 = BitSet.valueOf(v1); 
        BitSet bitset2 = BitSet.valueOf(v2);
        // Create a copy of bitset1 to find the intersection
        BitSet intersection = (BitSet) bitset1.clone();
        intersection.and(bitset2);
    
        int nIntersection = intersection.cardinality(); // Number of "on" bits in intersection
        int nA = bitset1.cardinality(); // Number of "on" bits in bitset1
        int nB = bitset2.cardinality(); // Number of "on" bits in bitset2
    
        if (nA + nB - nIntersection == 0) 
        {
            return 1.0; // Avoid division by zero if both sets are empty. Consider them identical.
        }
    
        return (double) nIntersection / (nA + nB - nIntersection);
    }
};
$$
;


-- Defined identically to the Tanimoto_J above. Need both for backward compatibility.
-- See the comments above.
CREATE OR REPLACE FUNCTION Tanimoto(v1 VARBINARY, v2 VARBINARY)
     RETURNS FLOAT
     LANGUAGE JAVA
     RETURNS NULL ON NULL INPUT
     IMMUTABLE
     HANDLER = 'Tanimoto.calculate'
     -- The Java code will be pre-compiled and stored in the jar file.
     -- Note that the path must be different from the one in Tanimoto_J above.
     TARGET_PATH = '@java_handlers/tanimoto_1.jar'
     COMMENT='A version of TANIMOTO implemented in Java. Computes Tanimoto similarity between two bitsets represented by binary vectors. Returns BITCOUNT(v1 bitand v2) / BITCOUNT(v1 bitor v2) as float. If both v1 and v2 have no bits set to 1, returns 1.0, that is, treats two empty bitsets or two bitsets filled with only zeroes as being equal to each other. Returns error if two bitsets are of different lengths. Returns NULL if one of both arguments are NULLs'
     AS
$$
import java.util.BitSet;
class Tanimoto 
{
    public static double calculate(byte[] v1, byte[] v2) 
    {
        // Note, we don't have to check the args for null,
        // because this UDF is declared with RETURNS NULL ON NULL INPUT,
        // so it won't be called if either of both v1 and v1 are null's.
        if (v1.length != v2.length)
        {
            String msg = String.format("Vectors representing bitsets must be of equal size. v1.length: %d, v2.length: %d.", 
                                v1.length, v2.length);
            throw new IllegalArgumentException(msg);
        }
        BitSet bitset1 = BitSet.valueOf(v1); 
        BitSet bitset2 = BitSet.valueOf(v2);
        // Create a copy of bitset1 to find the intersection
        BitSet intersection = (BitSet) bitset1.clone();
        intersection.and(bitset2);
    
        int nIntersection = intersection.cardinality(); // Number of "on" bits in intersection
        int nA = bitset1.cardinality(); // Number of "on" bits in bitset1
        int nB = bitset2.cardinality(); // Number of "on" bits in bitset2
    
        if (nA + nB - nIntersection == 0) 
        {
            return 1.0; // Avoid division by zero if both sets are empty. Consider them identical.
        }
    
        return (double) nIntersection / (nA + nB - nIntersection);
    }
};
$$
;


-- Molecular properties and descriptor UDF's ---------------------------

CREATE OR REPLACE FUNCTION Molstring_Or_Molbinary2Medchem_Descriptors(molstring VARCHAR DEFAULT NULL, molbinary VARBINARY DEFAULT NULL)
     RETURNS TABLE (Q_Estim_Drug_Likeness float, Mol_Wt float, Heavy_Atom_Mol_Wt float, Exact_Mol_Wt float, Num_Valence_Electrons int, Num_Radical_Electrons int, Max_Partial_Charge float, Min_Partial_Charge float, Max_Abs_Partial_Charge float, Min_Abs_Partial_Charge float, TPSA float, Fraction_CSP3 float, Heavy_Atom_Count int, Num_Aliphatic_Carbocycles int, Num_Aliphatic_Heterocycles int, Num_Aliphatic_Rings int, Num_Aromatic_Carbocycles int, Num_Aromatic_Heterocycles int, Num_Aromatic_Rings int, Num_HAcceptors int, Num_HDonors int, Num_Heteroatoms int, Num_Rotatable_Bonds int, Num_Saturated_Carbocycles int, Num_Saturated_Heterocycles int, Num_Saturated_Rings int, Ring_Count int, Mol_Log_P float, Mol_MR float)
     LANGUAGE PYTHON
     IMMUTABLE
     RUNTIME_VERSION = '3.11'
     PACKAGES = ('rdkit')
     HANDLER = 'DGen'
     COMMENT='Computes commonly used medchem properties (descriptors) for a molecule encoded as molstring (standard or ChemAxon SMILES or molfile/molblock) or as RDKit binary-encoded molecule. Only one of molstring or molbinary arguments can be non-NULL, otherwise, an error will be returned. Returns a table with one row and multiple columns corresponding to the computed descriptors. If molstring and molbinary are both NULLs, or molstring is invalid and cannot be parsed, returns a table with one row filled with NULLs. Returns an error if molbinary is not a valid RDKit binary-encoded molecule.'
     AS
$$
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

from functools import lru_cache
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
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

@lru_cache(128)
def getmol(molstring: Optional[str], molbinary: Optional[bytes]):
    if not molstring and not molbinary:  # does it make sense to return molfile or SMILES strings representing empty molecules?
        return None
    if molstring and molbinary:
        raise ValueError('Either molstring or molbinary can be not NULL, but not both')
    if molbinary:
        m = Chem.Mol(molbinary)
    else:
        m = Chem.MolFromMolBlock(molstring) if is_molfile(molstring) else Chem.MolFromSmiles(molstring)
    return m


_rdkCommonMedchem = [
    'qed',
    'MolWt',
    'HeavyAtomMolWt',
    'ExactMolWt',
    'NumValenceElectrons',
    'NumRadicalElectrons',
    'MaxPartialCharge',
    'MinPartialCharge',
    'MaxAbsPartialCharge',
    'MinAbsPartialCharge',
    'TPSA',
    'FractionCSP3',
    'HeavyAtomCount',
    'NumAliphaticCarbocycles',
    'NumAliphaticHeterocycles',
    'NumAliphaticRings',
    'NumAromaticCarbocycles',
    'NumAromaticHeterocycles',
    'NumAromaticRings',
    'NumHAcceptors',
    'NumHDonors',
    'NumHeteroatoms',
    'NumRotatableBonds',
    'NumSaturatedCarbocycles',
    'NumSaturatedHeterocycles',
    'NumSaturatedRings',
    'RingCount',
    'MolLogP',
    'MolMR'
]
_dc = MoleculeDescriptors.MolecularDescriptorCalculator(_rdkCommonMedchem)
_null_desc = (None,) * len(_rdkCommonMedchem)

class DGen:
    def process(self, molstring: Optional[str], molbinary: Optional[bytes]):
        m = getmol(molstring, molbinary)
        if not m:
            yield _null_desc
        else:
            yield _dc.CalcDescriptors(m)
$$
;


CREATE OR REPLACE FUNCTION Molstring_Or_Molbinary_Check(molstring VARCHAR DEFAULT NULL, molbinary VARBINARY DEFAULT NULL, raise_exception BOOLEAN DEFAULT FALSE)
     RETURNS TABLE (is_ok boolean, encoding varchar, error_msg varchar) 
     LANGUAGE PYTHON 
     IMMUTABLE    
     RUNTIME_VERSION = '3.11' 
     PACKAGES = ('rdkit')
     HANDLER = 'MolChecker'
     COMMENT='Checks a molecule encoded as molstring (standard or ChemAxon SMILES or molfile/molblock) or as RDKit binary-encoded molecule. Only one of molstring or molbinary arguments can be non-NULL, otherwise, an error will be returned. Returns a table with one row and three columns: is_ok boolean, encoding varchar, and error_msg varchar. If both molstring and molbinary are NULL, the entire result row will be filled with NULLs. Otherwise, is_ok will contain True iff the input can be parsed into a molecule with no errors, encoding will contain a string representation of the encoding (MOLBLOCK, SMILES, or BINARY), and the error_msg will contain a description of the error or NULL. If raise_exception parameter is TRUE (it is FALSE by default) and the input cannot be parsed into a valid molecule, the method will raise an exception and quit instead of returning.'
     AS
$$
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
$$
;

-- HELM support (Proof of Concept) ---------------------------------------------------------
CREATE OR REPLACE PROCEDURE helm2molbinary(monomer_table VARCHAR, helm_table VARCHAR)
RETURNS table(ID INT, MOLBINARY VARBINARY)
LANGUAGE PYTHON
RUNTIME_VERSION = 3.12
PACKAGES = ('snowflake-snowpark-python', 'rdkit')
HANDLER = 'main'
COMMENT=$$Converts HELM strings in the specified helm_table to full molecular representation and returns structures
encoded in the RDKit binary format as a tabile with the ID INT and BINARY_MOL VARBINARY columns. The IDs correspond
to those in the table specified via helm_table arg, which must contain at least two columns, ID INT and HELM VARCHAR.
Reads monomers from the table specified via monomer_table argument.
The table specified via the monomer_table argument must contain the MOLFILE VARCHAR column populated with
SDF records with data fields as documented in https://github.com/adaliaramon/helmkit
The current implementation relies on re-factored code from the above repo (see additional comments in HelmMolecule.py).
Example:
call helm2molbinary(TABLE(monomer), TABLE(peptide)) ->>
SELECT p.id, p.helm, p.original_name_in_source_literature
FROM $1 as h
join peptide p on p.id = h.id
where molbinary_matches_smarts(molbinary, 'c1cc2ccccc2n1', TRUE);
$$
AS
$$

# This code is based upon https://github.com/adaliaramon/helmkit
# Revision: 4401bccdfbce14a56e594740a2e2783967d66641
# Author: Ramon Adàlia <adaliaramon@gmail.com>
# Date: 6/15/2026 3:57:23 AM
# MIT license

# The original implementation has been significantly refactored, but the original logic has not changed (yet).
# See my comments below.

import bisect
import re
import warnings
from collections import defaultdict
from collections.abc import Callable
from collections.abc import Sequence
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
from typing import TypeVar, Iterable

from rdkit import Chem
from rdkit import rdBase


class SequenceConstants:
    max_rgroups: int = 4


def get_molecule_property(
    molecule: Chem.Mol, property_name: str, default: str | None = None
) -> str | None:
    return (
        molecule.GetProp(property_name) if molecule.HasProp(property_name) else default
    )


T = TypeVar("T")


def parse_comma_separated_property(
    molecule: Chem.Mol,
    property_name: str,
    convert_func: Callable[[str], T] | None = None,
) -> list[str | T | None]:
    property_value = get_molecule_property(molecule, property_name)
    if not property_value:
        return []

    values = property_value.split(",")
    if convert_func:
        values = [convert_func(v) if v != "None" else None for v in values]
    else:
        values = [None if v == "None" else v for v in values]

    return values


def infer_attachment_points(
    molecule: Chem.Mol, rgroup_indices: Sequence[int | None]
) -> list[int]:
    """Infer attachment points by finding atoms bonded to R-group atoms."""
    attachment_points = []

    for r_idx in rgroup_indices:
        if r_idx is None:
            attachment_points.append(None)
            continue

        atom = molecule.GetAtomWithIdx(r_idx)

        for bond in atom.GetBonds():
            other_idx = bond.GetOtherAtomIdx(r_idx)
            attachment_points.append(other_idx)
            break
        else:
            attachment_points.append(None)
            warnings.warn(
                f"R-group atom {r_idx} has no bonds to determine attachment point"
            )

    return attachment_points

# D.R.: much easier to use and more type-safe than the original typed dictionary
@dataclass
class MonomerData:
    m_romol: Chem.Mol
    m_Rgroups: list[str | None]
    m_RgroupIdx: list[int]
    m_attachmentPointIdx: list[int]
    m_type: str
    m_subtype: str
    m_abbr: str
    m_name: str | None = None
    m_chainID: str | None = None
    m_resID: int | None = None


MonomerLibrary = dict[str, dict[str, MonomerData]]

# D.R.: Note, this function can also incrementally add monomers to an existing lib
def load_monomers(molecules: Iterable[Chem.Mol], monomers_dict: MonomerLibrary = None) -> MonomerLibrary:
    if monomers_dict is None:
        monomers_dict = defaultdict(dict)
    for mol in molecules:
        if mol is None:
            continue

        symbol = get_molecule_property(mol, "symbol")
        if not symbol:
            continue

        m_type = get_molecule_property(mol, "m_type", "")

        rgroups = parse_comma_separated_property(mol, "m_Rgroups")
        rgroup_idx = parse_comma_separated_property(mol, "m_RgroupIdx", int)
        attachment_point_idx = infer_attachment_points(mol, rgroup_idx)

        monomers_dict[m_type][symbol] = MonomerData(mol, rgroups, rgroup_idx, attachment_point_idx, m_type,
                                                    get_molecule_property(mol, "m_subtype", ""),
                                                    get_molecule_property(mol, "m_abbr", ""))

    return monomers_dict

# D.R.: don't use @lru_cache here! See comments below (Molecule)
def load_monomer_library(sdf_library_path: str | None = None) -> MonomerLibrary:
    """Load and prepare monomer data from SDF file.
    If sdf_library_path is None or empty, loads the default monomer lib from monomers.sdf
    supplied as part of the package.
    """
    if not sdf_library_path:
        sdf_library_path = str(files("helmkit.data") / "monomers.sdf")
    supplier = Chem.SDMolSupplier(sdf_library_path, removeHs=False)

    return load_monomers(supplier)

# D.R.: that's where we really need the cache
@lru_cache
def _mol_from_smarts(smarts: str) -> Chem.Mol:
    m = Chem.MolFromSmarts(smarts)
    if not m:
        raise ValueError(f'Smarts {smarts} is None, empty or invalid')
    return m

def _create_missing_monomer(monomer_name: str, m_type: str = "aa") -> MonomerData:
    mol = Chem.MolFromSmiles(monomer_name, sanitize=False)
    if mol is None:
        if monomer_name.endswith("|") and not monomer_name.endswith("$|"):
            return _create_missing_monomer(monomer_name[:-1] + "$|", m_type)
        raise ValueError(
            f"Monomer {monomer_name} not in monomer library and is not a valid SMILES string"
        )
    with rdBase.BlockLogs():
        error = Chem.SanitizeMol(mol, catchErrors=True)
    if error == Chem.SanitizeFlags.SANITIZE_PROPERTIES:
        mol = Chem.RWMol(mol)
        pattern = _mol_from_smarts("O[CX4]=O")
        matches = mol.GetSubstructMatches(pattern)
        for drop_idx, *_ in matches:
            mol.RemoveAtom(drop_idx)
        error = Chem.SanitizeMol(mol, catchErrors=True)
    if error:
        raise ValueError(
            f"Monomer {monomer_name} not in monomer library and is not a valid SMILES string"
        )

    r_group_map = {}
    main_atoms = []

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        label = atom.GetProp("atomLabel") if atom.HasProp("atomLabel") else ""

        if label.startswith("_R"):
            try:
                r_num = int(label[2:])
                atom.SetProp("dummyLabel", f"R{r_num}")
                atom.SetIntProp("_MolFileRLabel", r_num)
                atom.SetProp("molFileValue", "*")
                r_group_map[r_num] = idx
            except ValueError:
                continue
        else:
            main_atoms.append(idx)

    sorted_r = sorted(r_group_map.items())
    r_group_idx = [idx for _, idx in sorted_r]
    mol = Chem.RenumberAtoms(mol, main_atoms + r_group_idx)

    rgroup_idx_full: list[int | None] = [None] * SequenceConstants.max_rgroups
    for i, (r_num, _) in enumerate(sorted_r):
        if 1 <= r_num <= SequenceConstants.max_rgroups:
            rgroup_idx_full[r_num - 1] = len(main_atoms) + i

    attachment_points = infer_attachment_points(mol, rgroup_idx_full)
    rgroup_vals = [None] * SequenceConstants.max_rgroups

    if m_type == "aa" and "_R1" not in monomer_name:
        matches = {
            idx
            for _, idx, _ in mol.GetSubstructMatches(
                _mol_from_smarts("[#6][NX3H][#6]")
            )
        }
        if len(matches) == 0:
            matches = {
                idx
                for idx, _ in mol.GetSubstructMatches(_mol_from_smarts("[NX3H2][#6]"))
            }
        if len(matches) == 1:
            attachment_id = matches.pop()

            mol = Chem.RWMol(mol)
            new_idx = mol.AddAtom(Chem.Atom(0))
            mol.AddBond(attachment_id, new_idx, Chem.BondType.SINGLE)
            rgroup_idx_full[0] = new_idx
            attachment_points[0] = attachment_id

    if m_type == "aa" and "_R2" not in monomer_name:
        aldehyde = _mol_from_smarts("[CX3H1]=O")
        matches = mol.GetSubstructMatches(aldehyde)
        if len(matches) == 0:
            matches = mol.GetSubstructMatches(_mol_from_smarts("[CX3](=O)[OH]"))
        if len(matches) == 1:
            attachment_id, *_ = matches[0]

            mol = Chem.RWMol(mol)
            new_idx = mol.AddAtom(Chem.Atom(0))
            mol.AddBond(attachment_id, new_idx, Chem.BondType.SINGLE)
            rgroup_idx_full[1] = new_idx
            attachment_points[1] = attachment_id

    mol.SetProp("m_name", monomer_name)

    mol.SetProp("symbol", monomer_name)
    mol.SetProp("m_abbr", monomer_name)
    mol.SetProp("m_type", m_type)
    mol.SetProp("m_subtype", "non-natural" if m_type == "aa" else "")
    mol.SetProp("m_RgroupIdx", ",".join(map(str, rgroup_idx_full)))
    mol.SetProp("m_Rgroups", ",".join(map(str, rgroup_vals)))
    mol.SetProp("m_attachmentPointIdx", ",".join(map(str, attachment_points)))
    mol.SetProp("natAnalog", "")
    return MonomerData(mol, rgroup_vals, rgroup_idx_full, attachment_points, m_type,
                       "non-natural" if m_type == "aa" else "", monomer_name)


# D.R.: it does not make sense to initialize monomer lib in the __init__ method.
# Using the same monomer lib for all instances of Molecule is also a bad idea.
# A much cleaner and less bug-prone approach is to always pass an instance of monomer lib as arg to __init__
class Molecule:
    """Single class for HELM to RDKit Mol conversion."""

    _bracket_re = re.compile(r"{(.*?)}")
    _pipe_outside_brackets = re.compile(r"\|(?![^\[]*\])")
    _dollar_outside_brackets = re.compile(r"\$(?![^\[]*\])")

    def __init__(self, helm: str, monomer_lib: MonomerLibrary):
        """Initialize a Molecule object from a HELM string."""
        self.mol = None
        self.offset = []
        self.bondlist = []
        self.monomers: list[MonomerData] = []
        self.chain_offset = {}
        self.residue_reps = defaultdict(list)
        self.has_ambiguous_monomers = False
        self.hydrogen_bonds = []
        assert monomer_lib
        self.monomer_lib = monomer_lib
        self._parse_helm_string(helm)
        self._build_molecule()

        if not isinstance(self.mol, Chem.rdchem.Mol):
            raise TypeError("Failed to initialize RDKit Mol object")

    def _parse_helm_string(self, helm: str) -> None:
        """Parse a HELM string into molecular components."""
        helm_parts = self._split_helm_sections(helm)

        if len(helm_parts) < 5:
            warnings.warn(f"Problem with HELM string - not enough sections: {helm}")
            return

        polymer_sections, connection_sections, hydrogen_bonds_sections = (
            helm_parts[0],
            helm_parts[1],
            helm_parts[2],
        )

        if not polymer_sections:
            warnings.warn(f"No simple polymers in HELM string {helm}")
            return

        self._process_polymers(polymer_sections)
        self._process_connections(connection_sections)
        self._process_hydrogen_bonds(hydrogen_bonds_sections)

    def _split_helm_sections(self, helm: str) -> list[str | list[str]]:
        """Split a HELM string into its components."""
        parts = self._dollar_outside_brackets.split(helm, 4)
        parts.extend([""] * (5 - len(parts)))

        parts[0] = (
            self._pipe_outside_brackets.split(parts[0])
            if "|" in parts[0]
            else [parts[0]]
        )

        if parts[1]:
            parts[1] = parts[1].split("|") if "|" in parts[1] else [parts[1]]
        else:
            parts[1] = []

        if parts[2]:
            parts[2] = parts[2].split("|") if "|" in parts[2] else [parts[2]]
        else:
            parts[2] = []

        return parts

    @staticmethod
    def _split_sequence_with_brackets(sequence: str) -> list[str]:
        """Split a sequence into individual monomers, respecting brackets."""
        result = []
        current = ""
        bracket_depth = 0

        for char in sequence:
            if char in "[(":
                bracket_depth += 1
                current += char
            elif char in "])":
                bracket_depth -= 1
                current += char
            elif char == "." and bracket_depth == 0:
                result.append(current)
                current = ""
            else:
                current += char

        if current:
            result.append(current)

        return result

    def _extract_chain_id(self, chain_str: str) -> tuple[str | None, bool, str | None]:
        """Extract chain ID and validate chain type."""
        match = re.match(r"([A-Z]+)(\d+)", chain_str)
        if not match:
            warnings.warn(f"Invalid chain format: {chain_str}")
            return None, False, None

        polymer_type = match.group(1)
        if polymer_type not in ("PEPTIDE", "RNA", "CHEM"):
            warnings.warn(f"Unsupported polymer type: {polymer_type}")
            return None, False, None
        else:
            return chain_str, True, polymer_type

    def _process_monomer(
        self, monomer_name: str, chain_id: str, residue_idx: int, polymer_type: str
    ) -> MonomerData | None:
        """Process a single monomer."""
        monomer_name = (
            monomer_name[1:-1]
            if monomer_name.startswith("[") and monomer_name.endswith("]")
            else monomer_name
        )
        if monomer_name == "":
            raise ValueError(f"Monomer {residue_idx + 1} has no name. Check HELM.")

        # Check for (a,[b]) pattern
        match = re.fullmatch(r"\([^,]+,\[([^\]]+)\]\)", monomer_name)
        if match:
            # Extract the 'b' from (a,[b]) and recurse
            self.has_ambiguous_monomers = True
            return self._process_monomer(
                match.group(1), chain_id, residue_idx, polymer_type
            )

        if polymer_type == "PEPTIDE":
            m_type = "aa"
        elif polymer_type == "RNA":
            m_type = "rna"
        elif polymer_type == "CHEM":
            m_type = "chem"
        else:
            m_type = "aa"

        if m_type in self.monomer_lib and monomer_name in self.monomer_lib[m_type]:
            monomer_info = self.monomer_lib[m_type][monomer_name]
        else:
            try:
                monomer_info = _create_missing_monomer(monomer_name, m_type)
                if m_type not in self.monomer_lib:
                    self.monomer_lib[m_type] = {}
                self.monomer_lib[m_type][monomer_name] = monomer_info
            except ValueError as e:
                warnings.warn(str(e))
                return None

        return MonomerData(monomer_info.m_romol, monomer_info.m_Rgroups[:], monomer_info.m_RgroupIdx,
                           monomer_info.m_attachmentPointIdx, monomer_info.m_type, monomer_info.m_subtype,
                           monomer_info.m_abbr, monomer_name, chain_id, residue_idx)

    @staticmethod
    def _parse_rna_string(sequence: str) -> list[str]:
        result = []
        current = ""
        bracket_depth = 0

        for char in sequence:
            if char in "[(":
                bracket_depth += 1
                current += char
            elif char in "])":
                bracket_depth -= 1
                current += char
            else:
                current += char
            if bracket_depth == 0:
                result.append(current)
                current = ""

        if current:
            result.append(current)

        return [r[1:-1] if r.startswith("[") and r.endswith("]") else r for r in result]

    def _process_polymers(self, polymers: list[str]) -> None:
        """Process polymer chains from HELM, creating backbone bonds on the fly."""
        monomer_idx = 0

        for chain in polymers:
            chain = chain.strip()
            match = self._bracket_re.search(chain)
            if not match:
                warnings.warn(f"No sequence in polymer: {chain}")
                continue

            id_chain = chain[: match.start()]
            chain_id, valid, polymer_type = self._extract_chain_id(id_chain)
            if not valid:
                continue

            if chain_id in self.chain_offset:
                raise ValueError(f"Duplicate chain ID: {chain_id}")

            sequence = match.group(1)
            if not sequence:
                warnings.warn(f"Empty polymer: {chain}")
                continue

            residues = self._split_sequence_with_brackets(sequence)
            self.chain_offset[chain_id] = monomer_idx

            if polymer_type == "PEPTIDE":
                for residue_idx, monomer_name in enumerate(residues):
                    monomer = self._process_monomer(
                        monomer_name, chain_id, residue_idx, polymer_type
                    )
                    if not monomer:
                        continue

                    self.monomers.append(monomer)
                    self.residue_reps[chain_id].append(monomer_idx)

                    if residue_idx > 0:
                        monomer1 = self.monomers[monomer_idx - 1]
                        monomer2 = monomer

                        attachment_point1 = monomer1.m_attachmentPointIdx[1]
                        if attachment_point1 is None:
                            raise ValueError(
                                f"R-group 2 is not present in monomer {monomer_idx} ({monomer1.m_name}). Check monomers."
                            )
                        attachment_point2 = monomer2.m_attachmentPointIdx[0]
                        if attachment_point2 is None:
                            raise ValueError(
                                f"R-group 1 is not present in monomer {monomer_idx + 1} ({monomer2.m_name}). Check monomers."
                            )

                        self.bondlist.append([
                            monomer_idx - 1,
                            attachment_point1,
                            monomer_idx,
                            attachment_point2,
                        ])
                        self._mark_used_rgroup(monomer_idx - 1, 1)
                        self._mark_used_rgroup(monomer_idx, 0)

                    monomer_idx += 1
            elif polymer_type == "RNA":
                prev_monomer = None
                for residue_idx, residue in enumerate(residues):
                    split_residue = self._parse_rna_string(residue)
                    for subresidue in split_residue:
                        is_base = subresidue.startswith("(") and subresidue.endswith(
                            ")"
                        )
                        monomer_name = subresidue[1:-1] if is_base else subresidue
                        monomer_name = (
                            monomer_name[1:-1]
                            if monomer_name.startswith("[")
                            and monomer_name.endswith("]")
                            else monomer_name
                        )
                        monomer = self._process_monomer(
                            monomer_name, chain_id, residue_idx, polymer_type
                        )

                        self.monomers.append(monomer)
                        self.residue_reps[chain_id].append(monomer_idx)

                        if prev_monomer is not None:
                            monomer1 = self.monomers[prev_monomer]
                            monomer2 = monomer

                            # Attach to R3 to R1 if the monomer is a base, R2 to R1 otherwise
                            r_index = 2 if is_base else 1
                            try:
                                attachment_point1 = monomer1.m_attachmentPointIdx[r_index]
                            except IndexError:
                                attachment_point1 = None
                            if attachment_point1 is None:
                                raise ValueError(
                                    f"R-group {r_index + 1} is not present in monomer {prev_monomer} ({monomer1.m_name}). Check monomers."
                                )
                            attachment_point2 = monomer2.m_attachmentPointIdx[0]
                            if attachment_point2 is None:
                                raise ValueError(
                                    f"R-group 1 is not present in monomer {monomer_idx} ({monomer2.m_name}). Check monomers."
                                )

                            self.bondlist.append([
                                prev_monomer,
                                attachment_point1,
                                monomer_idx,
                                attachment_point2,
                            ])
                            self._mark_used_rgroup(prev_monomer, r_index)
                            self._mark_used_rgroup(monomer_idx, 0)

                        # Only set prev_monomer if the monomer is not a base
                        if not is_base:
                            prev_monomer = monomer_idx

                        monomer_idx += 1
            elif polymer_type == "CHEM":
                if len(residues) != 1:
                    raise ValueError("CHEM polymers must have exactly one residue")
                monomer_name = residues[0]
                residue_idx = 0
                monomer_name = (
                    monomer_name[1:-1]
                    if monomer_name.startswith("[") and monomer_name.endswith("]")
                    else monomer_name
                )
                monomer = self._process_monomer(
                    monomer_name, chain_id, residue_idx, polymer_type
                )
                self.monomers.append(monomer)
                self.residue_reps[chain_id].append(monomer_idx)
                monomer_idx += 1

    def _parse_connection(
        self, connection_str: str
    ) -> tuple[str, int, int, str, int, int] | None:
        """Parse a single connection string."""
        parts = connection_str.split(",")
        if len(parts) != 3:
            warnings.warn(f"Invalid connection format: {connection_str}")
            return None

        chain_id1, chain_id2, bond_spec = parts

        try:
            bond_parts = re.split(r"[-:]", bond_spec)
            if len(bond_parts) != 4:
                warnings.warn(f"Invalid bond format: {bond_spec}")
                return None

            residue1, rgroup1, residue2, rgroup2 = bond_parts

            residue1 = int(residue1) - 1
            residue2 = int(residue2) - 1
            rgroup1 = int(rgroup1.replace("R", ""))
            rgroup2 = int(rgroup2.replace("R", ""))
        except (ValueError, IndexError) as e:
            warnings.warn(f"Error parsing connection {connection_str}: {e}")
            return None
        else:
            return chain_id1, residue1, rgroup1, chain_id2, residue2, rgroup2

    def _process_connections(self, connections: list[str]) -> None:
        """Process connections between chains."""
        if not connections:
            return

        for connection_str in connections:
            parsed = self._parse_connection(connection_str)
            if not parsed:
                continue

            chain_id1, residue1, rgroup1, chain_id2, residue2, rgroup2 = parsed
            rgroup1 -= 1
            rgroup2 -= 1

            monomer_idx1 = self.residue_reps[chain_id1][residue1]
            monomer_idx2 = self.residue_reps[chain_id2][residue2]

            monomer1 = self.monomers[monomer_idx1]
            monomer2 = self.monomers[monomer_idx2]

            attachment_idx1 = monomer1.m_attachmentPointIdx[rgroup1]
            if attachment_idx1 is None:
                raise ValueError(
                    f"R-group {rgroup1} is not present in monomer {monomer_idx1 + 1} ({monomer1.m_name}). Check connections."
                )
            attachment_idx2 = monomer2.m_attachmentPointIdx[rgroup2]
            if attachment_idx2 is None:
                raise ValueError(
                    f"R-group {rgroup2} is not present in monomer {monomer_idx2 + 1} ({monomer2.m_name}). Check connections."
                )

            self.bondlist.append([
                monomer_idx1,
                attachment_idx1,
                monomer_idx2,
                attachment_idx2,
            ])

            self._mark_used_rgroup(monomer_idx1, rgroup1)
            self._mark_used_rgroup(monomer_idx2, rgroup2)

    def _process_hydrogen_bonds(self, connections: list[str]) -> None:
        """Process hydrogen bonds."""
        if not connections:
            return

        for connection_str in connections:
            parts = connection_str.split(",")
            if len(parts) != 3:
                warnings.warn(f"Invalid hydrogen bond format: {connection_str}")
                continue
            chain_id1, chain_id2, bond_spec = parts

            bond_parts = re.split(r"[-:]", bond_spec)
            if len(bond_parts) != 4:
                warnings.warn(f"Invalid hydrogen bond format: {bond_spec}")
                continue

            residue1, _, residue2, _ = bond_parts
            residue1 = int(residue1) - 1
            residue2 = int(residue2) - 1
            self.hydrogen_bonds.append([chain_id1, residue1, chain_id2, residue2])

    def _mark_used_rgroup(self, monomer_idx: int, rgroup: int) -> None:
        """Mark an R-group as used based on its attachment point index."""
        monomer = self.monomers[monomer_idx]
        monomer.m_Rgroups[rgroup] = None

    def _build_molecule(self) -> None:
        """Build the RDKit molecule from parsed monomer and bond data."""
        if not self.monomers:
            self.mol = Chem.RWMol()
            return

        monomer = self.monomers[0]
        self.mol = Chem.RWMol(monomer.m_romol)

        rgroups = monomer.m_Rgroups
        rgroup_idx = monomer.m_RgroupIdx
        for i in range(min(len(rgroups), SequenceConstants.max_rgroups)):
            if rgroups[i] is not None:
                self._replace_rgroup(self.mol, 0, rgroup_idx[i], rgroups[i])

        current_offset = self.mol.GetNumAtoms()
        self.offset = [0, current_offset]

        for monomer in self.monomers[1:]:
            self.mol.InsertMol(monomer.m_romol)

            rgroups = monomer.m_Rgroups
            rgroup_idx = monomer.m_RgroupIdx
            for i in range(min(len(rgroups), SequenceConstants.max_rgroups)):
                if rgroups[i] is not None:
                    self._replace_rgroup(
                        self.mol, current_offset, rgroup_idx[i], rgroups[i]
                    )

            atom_count = monomer.m_romol.GetNumAtoms()
            current_offset += atom_count
            self.offset.append(current_offset)

        self._add_bonds()
        self._sanitize()

    def _add_bonds(self) -> None:
        """Add bonds between monomers based on bond list."""
        for monomer1_idx, atom1_idx, monomer2_idx, atom2_idx in self.bondlist:
            absolute_atom1_idx = self.offset[monomer1_idx] + atom1_idx
            absolute_atom2_idx = self.offset[monomer2_idx] + atom2_idx

            self.mol.AddBond(
                absolute_atom1_idx, absolute_atom2_idx, Chem.BondType.SINGLE
            )

    def _replace_rgroup(
        self, rdkit_mol: Chem.RWMol, atom_offset: int, atom_idx: int, atom_type: str
    ) -> None:
        """Replace an R-group with the appropriate atom type."""
        absolute_idx = atom_offset + atom_idx

        if atom_type == "OH":
            try:
                oxygen_atom = Chem.Atom(8)  # Oxygen
                rdkit_mol.ReplaceAtom(absolute_idx, oxygen_atom)
            except (RuntimeError, OverflowError) as e:
                warnings.warn(f"Failed to replace R-group with OH: {e}")
        elif atom_type != "H":
            warnings.warn(f"Unrecognized R-group type: {atom_type}")

    def _sanitize(self) -> None:
        """Clean up the molecule by removing dummy atoms."""
        pattern = _mol_from_smarts("[#0]")
        matches = self.mol.GetSubstructMatches(pattern)
        atoms_to_delete = sorted({idx for match in matches for idx in match})
        self.mol = Chem.DeleteSubstructs(self.mol, pattern)

        def correction(offset: int, idx: int) -> int:
            return bisect.bisect_left(atoms_to_delete, idx) - bisect.bisect_left(
                atoms_to_delete, offset
            )

        for i, (m1, a1, m2, a2) in enumerate(self.bondlist):
            offset1 = self.offset[m1]
            offset2 = self.offset[m2]
            self.bondlist[i][1] -= correction(offset1, offset1 + a1)
            self.bondlist[i][3] -= correction(offset2, offset2 + a2)

        self.offset = [
            offset - sum(d < offset for d in atoms_to_delete) for offset in self.offset
        ]

    @property
    def bond_indices(self) -> list[int]:
        return [
            self.mol.GetBondBetweenAtoms(
                self.offset[monomer1_idx] + atom1_idx,
                self.offset[monomer2_idx] + atom2_idx,
            ).GetIdx()
            for monomer1_idx, atom1_idx, monomer2_idx, atom2_idx in self.bondlist
        ]

    @property
    def monomer_indices(self) -> list[int]:
        return [
            bisect.bisect_right(self.offset, i) - 1
            for i in range(self.mol.GetNumAtoms())
        ]

# D.R.: Removed load_in_parallel, it is unsafe and can lead to various hard-to-diagnose errors
# when used in runtimes incompatible with process pools.
# Users can add their own parallel loading implementations if they really need those, but in most practical use cases
# the gain in the overall performance when process pools are used will be negligible or even negative.

from snowflake.snowpark.types import *
from snowflake.snowpark import Session, DataFrame

try:
    # this is for local tests
    from HelmMolecule import *
except ImportError:
    # module is inlined, ignore the import error
    pass

def mol_from_molfile_str(ms: str) -> Chem.Mol | None:
    if not ms:
        return None
    # No idea why such a weird method needs to be used to create an instance of
    # RDKit molecule with all data fields from an SDF record. Chem.MolFromMolBlock will read just the ctab
    # and no data fields (see https://github.com/rdkit/rdkit/discussions/5747)
    try:
        s = Chem.SDMolSupplier()
        s.SetData(ms)
        return next(s)
    except Exception:
        # todo: better error handling
        return None

def main(session: Session, monomer_table: str, helm_table: str) -> DataFrame:
    df = session.table(monomer_table)
    # Obviously, it is assumed that monomer_table has a column called MOLFILE
    # and it contains SDF strings with the correct data fields
    molfiles = (row['MOLFILE'] for row in df.select('MOLFILE').collect())
    mols = (mol_from_molfile_str(mf) for mf in molfiles)
    # no None's!
    mols = (m for m in mols if m is not None)
    monomer_lib = load_monomers(mols)
    df = session.table(helm_table)
    # helm_table must contain ID (int) and HELM (varchar) columns
    helm_records = ((row['ID'], row['HELM']) for row in df.collect())

    def helm_to_binary(h: str | None) -> bytes | None:
        if not h:
            return None
        try:
            molecule = Molecule(h, monomer_lib)
            if molecule.mol.GetNumAtoms() == 0:  # don't need empty mols
                return None
            return molecule.mol.ToBinary()
        except Exception:
            # todo: better error handling (perhaps, add 'status' column to the result df and populate it with
            # success, error, warning, etc.)
            return None
    # int() is not really necessary but is used just in case the ID column contains integers, but the column type
    # is NUMBER with decimal points (non-zero scale), which is a bad idea in general, but we can handle it.
    data = [(int(hr[0]), helm_to_binary(hr[1])) for hr in helm_records]

    schema = StructType([StructField("ID", IntegerType()),
                     StructField("BINARY_MOL", BinaryType())])

    return session.create_dataframe(data, schema)

$$;
