import base64

from rdkit import Chem
from typing import Optional
from functools import lru_cache

from rdkit.Chem.Draw import rdMolDraw2D

try:
    from Util import *
except ImportError:
    # module is inlined, ignore the import error
    pass


@lru_cache(maxsize=128)
def get_pattern_mol(smarts: Optional[str]) -> Optional[Chem.Mol]:
    if not smarts:
        return None
    m = Chem.MolFromSmarts(smarts)
    if not m:
       raise ValueError(f'Error parsing SMARTS: {smarts}')
    if m.GetNumAtoms() == 0:
        return None
    return m


# todo: auto-detect not only molfile/sdf and SMILES, but also other
# encodings (InChi, etc.). Also, move all these duplicate defs to Util.py
def is_molfile(molstring: Optional[str]) -> bool:
    if not molstring:
        return False
    return '\n' in molstring


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


def mol_to_svg(mol: Chem.Mol | None, atoms_to_highlight: list | None) -> str | None:
    if not mol:
        return None
    data = rdMolDraw2D.MolToSVG(mol, highlightAtoms=atoms_to_highlight)
    # data:image/svg+xml;utf8,<svg xmlns=...</svg>
    # remove the first line, which looks like this: <?xml version='1.0' encoding='iso-8859-1'?>
    index = data.find('\n')
    data = data[index + 1:] if index != -1 else ''
    b64 = base64.b64encode(data.encode('utf-8')).decode('utf-8')
    # not sure why, but the plain utf8 won't work
    return f'data:image/svg+xml;base64,{b64}'


def molstring_or_molbinary_to_svg(molstring: Optional[str], molbinary: Optional[bytes],
                                  highlight_smarts: Optional[str]) -> Optional[str]:
    m = getmol(molstring, molbinary)
    if not m:
        return None
    highlight_mol = get_pattern_mol(highlight_smarts)
    atoms_to_highlight = None
    if highlight_mol:
        matches = m.GetSubstructMatches(highlight_mol)
        atoms_to_highlight = [idx for match in matches for idx in match]
    return mol_to_svg(m, atoms_to_highlight)
