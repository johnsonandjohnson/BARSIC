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
import Molstring_Or_Molbinary2Molstring as mbs

# charged + imine
enamine1 = R"""
  Mt7.00  02222222222D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 22 24 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.673929 1.287201 0.000000 0
M  V30 2 C -1.173929 2.153226 0.000000 0 CHG=-1
M  V30 3 C -2.173929 2.153226 0.000000 0
M  V30 4 C -2.673929 1.287201 0.000000 0
M  V30 5 C -2.173929 0.421175 0.000000 0
M  V30 6 C -1.173929 0.421175 0.000000 0
M  V30 7 N 0.326071 1.287201 0.000000 0 CHG=1
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
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 1
M  V30 7 2 1 7
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
M  V30 23 2 20 21
M  V30 24 1 21 22
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STERAC1 ATOMS=(1 12)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
$$$$"""

enamine2 = R"""
  Mt7.00  02222222222D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 22 24 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.673929 1.287201 0.000000 0
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




@enum.unique
class HashSchemeX(enum.Enum):
    STEREO_INSENSITIVE_TAUTOMER_INSENSITIVE_LAYERS = (
        rh.HashLayer.FORMULA,
        rh.HashLayer.NO_STEREO_TAUTOMER_HASH
    )

    TAUTOMER_INSENSITIVE_LAYERS = (
        rh.HashLayer.FORMULA,
        rh.HashLayer.TAUTOMER_HASH
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


def test1():
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


def test2():
    m1 = Chem.MolFromMolBlock(enamine1)
    layers_dict1 = rh.GetMolLayers(m1, enable_tautomer_hash_v2=False)
    m2 = Chem.MolFromMolBlock(enamine2)
    layers_dict2 = rh.GetMolLayers(m2, enable_tautomer_hash_v2=False)

    for k in layers_dict1:
        print(f'{k}: layers1: {layers_dict1[k]}, layers2: {layers_dict2[k]}, match: {layers_dict1[k] == layers_dict2[k]}')

    pass

if __name__ == '__main__':
    test2()
