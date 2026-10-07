import base64

import rdkit
from rdkit import Chem
from typing import Optional

from rdkit.Chem.Draw import rdMolDraw2D

try:
    from Util import *
except ImportError:
    # module is inlined, ignore the import error
    pass


def mol_to_svg(mol: Chem.Mol | None, atoms_to_highlight: list | None, output_raw_xml: bool) -> str | None:
    if not mol:
        return None
    data = rdkit.Chem.Draw.MolToSVG(mol, kekulize=False, wedgeBonds=True, highlightAtoms=atoms_to_highlight)
    if output_raw_xml:
        return data
    # data:image/svg+xml;utf8,<svg xmlns=...</svg>
    # remove the first line, which looks like this: <?xml version='1.0' encoding='iso-8859-1'?>
    index = data.find('\n')
    data = data[index + 1:] if index != -1 else ''
    b64 = base64.b64encode(data.encode('utf-8')).decode('utf-8')
    # not sure why, but the plain utf8 won't work
    return f'data:image/svg+xml;base64,{b64}'


def molstring_or_molbinary_to_svg(molstring: Optional[str], molbinary: Optional[bytes],
                                  highlight_smarts: Optional[str],
                                  draw_options: Optional[str], output_raw_xml: bool) -> Optional[str]:
    m = getmol(molstring, molbinary)
    if not m:
        return None
    highlight_mol = get_pattern_mol(highlight_smarts)
    atoms_to_highlight = None
    if highlight_mol:
        matches = m.GetSubstructMatches(highlight_mol)
        atoms_to_highlight = [idx for match in matches for idx in match]
    return mol_to_svg(m, atoms_to_highlight, output_raw_xml)
