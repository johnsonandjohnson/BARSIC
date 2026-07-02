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
