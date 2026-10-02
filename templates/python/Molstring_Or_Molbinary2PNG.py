import io
from rdkit import Chem

try:
    from Util import *
except ImportError:
    # module is inlined, ignore the import error
    pass

def mol_to_png(mol: Chem.Mol | None, width_px: int, height_px: int, atoms_to_highlight: list | None) -> bytes | None:
    if not mol:
        return None
    pil_image = Chem.Draw.MolToImage(mol, size=(width_px, height_px), highlightAtoms=atoms_to_highlight)
    with io.BytesIO() as buffer:
        pil_image.save(buffer, 'png')
        return buffer.getvalue()


def molstring_or_molbinary_to_png(width_px: int, height_px: int,
                                  molstring: Optional[str], molbinary: Optional[bytes],
                                  highlight_smarts: Optional[str],
                                  draw_options: Optional[str]) -> Optional[bytes]:
    if width_px <= 0 or height_px <= 0:
        raise ValueError('width_px and height_px must be > 0')
    m = getmol(molstring, molbinary)
    if not m:
        return None
    highlight_mol = get_pattern_mol(highlight_smarts)
    atoms_to_highlight = None
    if highlight_mol:
        matches = m.GetSubstructMatches(highlight_mol)
        atoms_to_highlight = [idx for match in matches for idx in match]
    return mol_to_png(m, width_px, height_px, atoms_to_highlight)