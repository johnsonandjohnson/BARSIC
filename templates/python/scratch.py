import unittest
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
import io
import pandas as pd
from rdkit.Chem.SaltRemover import InputFormat

import Molstring_Or_Molbinary2Molstring as mbs
import Molstring_Or_Molbinary2Molbinary as mbb

# charged + imine
enamine_mol1 = R"""
  Mt7.00  02222222222D

 22 24  0  0  1  0  0  0  0  0999 V2000
   -2.0015    1.2872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5015    2.1532    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5015    2.1532    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0015    1.2872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5015    0.4212    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5015    0.4212    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0015    1.2872    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4137    0.4782    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4137    2.0962    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7227   -0.4729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5374    0.7872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5374    1.7872    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   -1.7009   -0.6808    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0536   -1.2160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3626   -2.1671    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3407   -2.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0099   -1.6318    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3464    2.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4034    0.2872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2694    0.7872    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.1354    0.2872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0015    0.7872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  6  1  1  0  0  0  0
  1  7  1  0  0  0  0
  7  8  1  0  0  0  0
  7  9  1  0  0  0  0
  8 10  1  0  0  0  0
  8 11  1  0  0  0  0
 11 12  1  0  0  0  0
  9 12  1  0  0  0  0
 13 10  2  0  0  0  0
 10 14  1  0  0  0  0
 14 15  2  0  0  0  0
 15 16  1  0  0  0  0
 16 17  2  0  0  0  0
 17 13  1  0  0  0  0
 12 18  1  6  0  0  0
 11 19  1  0  0  0  0
 19 20  1  0  0  0  0
 20 21  2  0  0  0  0
 21 22  1  0  0  0  0
M  END
$$$$"""

# note the y-coordinate difference of the first atom
enamine_mol2 = R"""
  Mt7.00  02222222222D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 22 24 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.673929 2.287201 0.000000 0
M  V30 2 C -1.173929 2.153226 0.000000 0
M  V30 3 C -2.173929 2.153226 0.000000 0
M  V30 4 C -2.673929 1.287201 0.000000 0
M  V30 5 C -2.173929 0.421175 0.000000 0
M  V30 6 C -1.173929 0.421175 0.000000 0
M  V30 7 N 0.326071 1.287201 0.000000 0
M  V30 8 C 0.913856 0.478184 0.000000 0
M  V30 9 C 0.913856 2.096218 0.000000 0
M  V30 10 C 0.604839 -0.472873 0.000000 0
M  V30 11 C 1.864912 0.787201 0.000000 0
M  V30 12 C 1.864912 1.787201 0.000000 0 CFG=3
M  V30 13 C -0.373309 -0.680785 0.000000 0
M  V30 14 C 1.273969 -1.216018 0.000000 0
M  V30 15 C 0.964952 -2.167074 0.000000 0
M  V30 16 C -0.013195 -2.374986 0.000000 0
M  V30 17 C -0.682326 -1.631841 0.000000 0
M  V30 18 C 2.673929 2.374986 0.000000 0
M  V30 19 C 2.730938 0.287201 0.000000 0
M  V30 20 N 3.596963 0.787201 0.000000 0
M  V30 21 C 4.462988 0.287201 0.000000 0
M  V30 22 C 5.329014 0.787201 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 1
M  V30 7 1 1 7
M  V30 8 1 7 8
M  V30 9 1 7 9
M  V30 10 1 8 10
M  V30 11 1 8 11
M  V30 12 1 11 12
M  V30 13 1 9 12
M  V30 14 2 13 10
M  V30 15 1 10 14
M  V30 16 2 14 15
M  V30 17 1 15 16
M  V30 18 2 16 17
M  V30 19 1 17 13
M  V30 20 1 12 18 CFG=3
M  V30 21 1 11 19
M  V30 22 1 19 20
M  V30 23 1 20 21
M  V30 24 2 21 22
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STERAC1 ATOMS=(1 12)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
$$$$"""

testCsv = R"""SmilesText
"CCN/C(=N/CCO)NCCCC[C@H](N)C(=O)O"
"CC/N=C(\NCCO)NCCCC[C@H](N)C(=O)O"
"CCNC(=NCCC(=O)O)NCCCC[C@H](N)C(=O)O"
"CCN/C(=N\CCC(=O)O)NCCCC[C@H](N)C(=O)O"
"CCN/C(=N/CCC(=O)O)NCCCC[C@H](N)C(=O)O"
"CC/N=C(/NCCCC[C@H](N)C(=O)O)NCCN(C)C"
"CCNC(=NCCN(C)C)NCCCCC[C@H](N)C(=O)O"
CC1=CC=NC(=O)N1
CC1=NC(=O)NC=C1
CC1=NC(O)=NC=C1
CNC1=NC(C)=CC=C1
CN=C1NC(C)=CC=C1
CNN=C1CC=CC=C1
CN=NC1CC=CC=C1
OC1=CC=NC=C1
O=C1C=CNC=C1
C1C=CN=C1
N1C=CC=C1
"""

keto_enol_smiles1 = r'C=C(O)[C@@H](C)CC'
keto_enol_smiles2 = r'CC[C@@H](C)C(C)=O'

enamine_smiles1 = r'C=C(NC)[C@@H](C)CC'
enamine_smiles2 = r'CC[C@@H](C)/C(C)=N\C'

s1 = r'CC/N=C(\NCCO)NCCCC[C@H](N)C(=O)O'
s2 = r'CCNC(=N/CCO)\NCCCC[C@@H](N)C(=O)O'

ss1 = r'CCNC(=NCCC(=O)O)NCCCC[C@H](N)C(=O)O'
ss2 = r'CCN/C(=N\CCC(=O)O)NCCCC[C@H](N)C(=O)O'
ss3 = r'CCN/C(=N/CCC(=O)O)NCCCC[C@H](N)C(=O)O'

@enum.unique
class HashSchemeX(enum.Enum):
    EXACT_LAYERS = (
        rh.HashLayer.CANONICAL_SMILES,
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

def comp_hash(smiles: str):
    m = Chem.MolFromSmiles(smiles)
    layers_dict = rh.GetMolLayers(m, enable_tautomer_hash_v2=False)
    tih = rh.GetMolHash(layers_dict, hash_scheme=HashSchemeX.TAUTOMER_INSENSITIVE_LAYERS)
    sitih = rh.GetMolHash(layers_dict, hash_scheme=HashSchemeX.STEREO_INSENSITIVE_TAUTOMER_INSENSITIVE_LAYERS)
    return (tih, sitih)
    pass

def comp_hash_df(row):
    smiles = row['SmilesText']
    return comp_hash(smiles)

def print_df_as_is(df):
    print(df.to_string())


class Tests(unittest.TestCase):

    def test1(self):
        df = pd.read_csv(io.StringIO(testCsv))
        df[['TAUTO_HASH', 'NO_STEREO_TAUTO_HASH']] = df.apply(comp_hash_df, axis=1, result_type='expand')

        df_copy = df.copy(deep=True)
        df.drop(columns=['NO_STEREO_TAUTO_HASH'], inplace=True)
        df_copy.drop(columns=['TAUTO_HASH'], inplace=True)

        print('**Tautomeric hash***************************')
        r = df.groupby(['TAUTO_HASH'], as_index=False).agg({'SmilesText': ', '.join})
        print_df_as_is(r)

        print('**Stereo-insensitive tautomeric hash********')
        r = df_copy.groupby(['NO_STEREO_TAUTO_HASH'], as_index=False).agg({'SmilesText': ', '.join})
        print_df_as_is(r)


        print('Done')
        pass


    def test2(self):
        m1 = Chem.MolFromSmiles(keto_enol_smiles1)
        layers_dict1 = rh.GetMolLayers(m1, enable_tautomer_hash_v2=True)
        m2 = Chem.MolFromSmiles(keto_enol_smiles2)
        layers_dict2 = rh.GetMolLayers(m2, enable_tautomer_hash_v2=True)

        for k in layers_dict1:
            print(f'{k}: layers1: {layers_dict1[k]}, layers2: {layers_dict2[k]}, match: {layers_dict1[k] == layers_dict2[k]}')

        pass

    def test3(self):
        m1 = Chem.MolFromSmiles('CNCN')
        m2 = Chem.MolFromSmiles('CNCN')
        d = dict()
        d[m1] = 1
        d[m2] = 2
        print(d[m1])

    def test4(self):
        options = ['--out molhash', '--out molhashfull', '--out tautohash',
                   '--remove_stereo --out molhash', '--remove_stereo --out molhashfull', '--remove_stereo --out tautohash']
        def p(s1: str, s2: str):
            for opt in options:
                h1 = mbs.molstring_or_molbinary_to_molstring(s1, None, opt)
                h2 = mbs.molstring_or_molbinary_to_molstring(s2, None, opt)
                print(f'Opt: {opt}, Same hash? {h1 == h2}')

        print('Keto-enol ***********')
        p(keto_enol_smiles1, keto_enol_smiles2)
        print('Enamine-imine ***********')
        p(enamine_smiles1, enamine_smiles2)

    def test5(self):
        m1 = mbs.getmol(keto_enol_smiles2, None)
        m2 = mbs.getmol(keto_enol_smiles2, None)

        mx = mbs.getmol(keto_enol_smiles2, None)
        h = hash(mx)
        mx_clone = Chem.Mol(mx, quickCopy=True)
        h_clone = hash(mx_clone)

        lst = [m1, m2]
        lst = lst * 10

        for m in lst:
            res = mbs.mol_to_reg_layers(m)
        cache_info = mbs.mol_to_reg_layers.cache_info()
        #self.assertEquals(cache_info.)
        print(mbs.mol_to_reg_layers.cache_info())
        print(mbs.getmol.cache_info())


    def test6(self):
        s = 'CNCCN(C)CCCC(C(C)C)N1CC(C)(CCc2ccccc2O)C1.CC(=O)O.O=S(=O)(O)c1ccccc1.Cl'
        mol = Chem.MolFromSmiles(s)
        patterns = 'O=S(=O)(O)c1ccccc1|O=C(O)C1CCCC1'
        patterns = patterns.replace('|', '\n')
        default_salt_remover = SaltRemover.SaltRemover()
        salt_remover = SaltRemover.SaltRemover(defnData=patterns, defnFormat=InputFormat.SMARTS)
        salt_remover.salts = salt_remover.salts + default_salt_remover.salts;
        m1 = salt_remover.StripMol(mol)
        s1 = Chem.MolToSmiles(m1)
        print(s1)

    def test7(self):
        s = 'CNCCN(C)CCCC(C(C)C)N1CC(C)(CCc2ccccc2O)C1.CC(=O)O.O=S(=O)(O)c1ccccc1.Cl'
        expected_default = 'CNCCN(C)CCCC(C(C)C)N1CC(C)(CCc2ccccc2O)C1.O=S(=O)(O)c1ccccc1'
        expected_with_additional_patterns = 'CNCCN(C)CCCC(C(C)C)N1CC(C)(CCc2ccccc2O)C1'
        patterns = 'O=S(=O)(O)c1ccccc1|O=C(O)C1CCCC1'
        b0 = mbb.molstring_or_molbinary_to_molbinary(s, None, '--desalt')
        b1 = mbb.molstring_or_molbinary_to_molbinary(s, None, f'--desalt --desalt_smarts_list {patterns}')

        sb0 = mbs.molstring_or_molbinary_to_molstring(None, b0, '--out smi')
        sb1 = mbs.molstring_or_molbinary_to_molstring(None, b1, '--out smi')

        s0 = mbs.molstring_or_molbinary_to_molstring(s, None, '--out smi --desalt')
        s1 = mbs.molstring_or_molbinary_to_molstring(s, None, f'--out smi --desalt --desalt_smarts_list {patterns}')

        self.assertEquals(expected_default, sb0)
        self.assertEquals(expected_default, s0)

        self.assertEquals(expected_with_additional_patterns, sb1)
        self.assertEquals(expected_with_additional_patterns, s1)


        print('Done')
        pass

if __name__ == '__main__':
    unittest.main()
