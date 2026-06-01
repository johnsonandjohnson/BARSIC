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
from typing import Callable
from rdkit.Chem.rdchem import MolSanitizeException

def safe_call_decorator(func: Callable, exception_type=MolSanitizeException):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except exception_type:
            return None
    return wrapper
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

@safe_call_decorator
def molstring_or_molbinary_to_molstring(molstring: Optional[str], molbinary: Optional[bytes], option_str: str):
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
from typing import Callable
from rdkit.Chem.rdchem import MolSanitizeException

def safe_call_decorator(func: Callable, exception_type=MolSanitizeException):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except exception_type:
            return None
    return wrapper
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
from typing import Callable
from rdkit.Chem.rdchem import MolSanitizeException

def safe_call_decorator(func: Callable, exception_type=MolSanitizeException):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except exception_type:
            return None
    return wrapper
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
from typing import Callable
from rdkit.Chem.rdchem import MolSanitizeException

def safe_call_decorator(func: Callable, exception_type=MolSanitizeException):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except exception_type:
            return None
    return wrapper
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
from typing import Callable
from rdkit.Chem.rdchem import MolSanitizeException

def safe_call_decorator(func: Callable, exception_type=MolSanitizeException):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except exception_type:
            return None
    return wrapper
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
from typing import Callable
from rdkit.Chem.rdchem import MolSanitizeException

def safe_call_decorator(func: Callable, exception_type=MolSanitizeException):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except exception_type:
            return None
    return wrapper
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
from typing import Callable
from rdkit.Chem.rdchem import MolSanitizeException

def safe_call_decorator(func: Callable, exception_type=MolSanitizeException):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except exception_type:
            return None
    return wrapper
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
from typing import Callable
from rdkit.Chem.rdchem import MolSanitizeException

def safe_call_decorator(func: Callable, exception_type=MolSanitizeException):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except exception_type:
            return None
    return wrapper
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
from typing import Callable
from rdkit.Chem.rdchem import MolSanitizeException

def safe_call_decorator(func: Callable, exception_type=MolSanitizeException):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except exception_type:
            return None
    return wrapper
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
from typing import Callable
from rdkit.Chem.rdchem import MolSanitizeException

def safe_call_decorator(func: Callable, exception_type=MolSanitizeException):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except exception_type:
            return None
    return wrapper
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

CREATE OR REPLACE FUNCTION Tanimoto(v1 VARBINARY, v2 VARBINARY)
     RETURNS FLOAT 
     LANGUAGE PYTHON 
     RETURNS NULL ON NULL INPUT
     IMMUTABLE    
     RUNTIME_VERSION = '3.11' 
     PACKAGES = ('bitarray==2.5.1')  -- note: the latest version in Snowflake (3.4.2 as of now) seems to be broken!
     HANDLER = 'tanimoto'
     COMMENT='Computes Tanimoto similarity between two bitsets represented by binary vectors. Returns BITCOUNT(v1 bitand v2) / BITCOUNT(v1 bitor v2) as float. If both v1 and v2 have no bits set to 1, returns 1.0, that is, treats two empty bitsets or two bitsets filled with only zeroes as being equal to each other. Returns error if two bitsets are of different lengths. Returns NULL if one of both arguments are NULLs'
     AS
$$
from typing import Optional
from bitarray import bitarray
def tanimoto(b1: Optional[bytes], b2: Optional[bytes])->Optional[float]:
    if not b1 or not b2:
        return None
    b1b = bitarray()
    b2b = bitarray()
    b1b.frombytes(b1)
    b2b.frombytes(b2)
    if len(b1b) != len(b2b):
        raise ValueError('Bitsets are of different lengths')
    ba_and = b1b & b2b  #  intersection (common bits)
    and_count = ba_and.count(1)
    union_cnt = b1b.count(1) + b2b.count(1) - and_count  # could also use (b1b | b2b).count(1), but it would be slower
    if union_cnt == 0:
        return 1.0  # can only happen if both b1 and b2 have no bits set to 1
    return and_count / float(union_cnt)
$$
;


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
from typing import Callable
from rdkit.Chem.rdchem import MolSanitizeException

def safe_call_decorator(func: Callable, exception_type=MolSanitizeException):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except exception_type:
            return None
    return wrapper
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
from typing import Callable
from rdkit.Chem.rdchem import MolSanitizeException

def safe_call_decorator(func: Callable, exception_type=MolSanitizeException):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except exception_type:
            return None
    return wrapper
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
