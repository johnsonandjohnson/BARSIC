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


if __name__ == '__main__':
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

