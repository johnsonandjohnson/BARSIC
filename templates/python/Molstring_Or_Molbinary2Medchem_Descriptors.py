from functools import lru_cache
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
from typing import Optional

try:
    from Util import *
except ImportError:
    # module is inlined, ignore the import error
    pass


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
