-- Set db/schema in which you want these UDF's to be defined.
-- It is recommended that the schema be called chem_api.
--USE <database_name_here>;
CREATE SCHEMA IF NOT EXISTS chem_api;
USE schema chem_api;

-- Category: Info
CREATE OR REPLACE FUNCTION Barsic_Info()
     RETURNS VARCHAR
     LANGUAGE PYTHON
     IMMUTABLE
-- @include python/Runtime_Repo.inc
     PACKAGES = (
-- @include python/Snowpark_Pkg.inc
     ,
-- @include python/Rdkit_Pkg.inc
     )
     HANDLER = 'barsic_info'
     COMMENT='Returns toolkit info'
     AS
$$
-- @include python/Barsic_Info.py
$$
;

-- Category: Chemical structure depiction

CREATE OR REPLACE FUNCTION Draw_Molstring_Or_Molbinary2SVG(molstring VARCHAR DEFAULT NULL, molbinary VARBINARY DEFAULT NULL, highlight_smarts VARCHAR DEFAULT NULL, draw_options VARCHAR DEFAULT NULL, output_raw_xml BOOLEAN DEFAULT FALSE)
     RETURNS VARCHAR
     LANGUAGE PYTHON
     IMMUTABLE
-- @include python/Runtime_Repo.inc
     PACKAGES = (
-- @include python/Rdkit_Pkg.inc
     )
     HANDLER = 'molstring_or_molbinary_to_svg'
     COMMENT='Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) or RDKit binary molecule encoding (only one of molstring and molbinary can be not NULL) to its graphical representation in the SVG format. If highlight_smarts substructure pattern is specified, highlights matching atoms. If output_raw_xml is FALSE, returns picture as varchar with the following MIME type: data:image/svg+xml;base64, otherwise, returns plain unencoded xml. draw_options arg is not used yet and is reserved for future use.'
     AS
$$
-- @include python/Util.py
-- @include python/Molstring_Or_Molbinary2SVG.py
$$
;


CREATE OR REPLACE FUNCTION Draw_Molstring_Or_Molbinary2PNG(width_px INT, height_px INT, molstring VARCHAR DEFAULT NULL, molbinary VARBINARY DEFAULT NULL, highlight_smarts VARCHAR DEFAULT NULL, draw_options VARCHAR DEFAULT NULL)
     RETURNS VARBINARY
     LANGUAGE PYTHON
     IMMUTABLE
-- @include python/Runtime_Repo.inc
     PACKAGES = (
-- @include python/Rdkit_Pkg.inc
     )
     HANDLER = 'molstring_or_molbinary_to_png'
     COMMENT='Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) or RDKit binary molecule encoding (only one of molstring and molbinary can be not NULL) to its graphical representation in the PNG format (width_px x height_px is the image size in pixels). If highlight_smarts substructure pattern is specified, highlights matching atoms. draw_options arg is not used yet and is reserved for future use.'
     AS
$$
-- @include python/Util.py
-- @include python/Molstring_Or_Molbinary2PNG.py
$$
;

-- Category: Chemical structure conversion

CREATE OR REPLACE FUNCTION Molstring_Or_Molbinary2Molstring(molstring VARCHAR DEFAULT NULL, molbinary VARBINARY DEFAULT NULL, options VARCHAR DEFAULT '--out smiles')
     RETURNS VARCHAR 
     LANGUAGE PYTHON 
     IMMUTABLE    
-- @include python/Runtime_Repo.inc
     PACKAGES = (
-- @include python/Rdkit_Pkg.inc
     )
     HANDLER = 'molstring_or_molbinary_to_molstring'
     COMMENT='Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) or RDKit binary molecule encoding to a variety of different string-based formats, optionally transforming the input structure. Only one of molstring or molbinary arguments can be non-NULL. All available options can be listed by by running the following SQL: select Molstring_Or_Molbinary2Molstring(NULL, NULL, ''-h''); usage info will be returned as part of the error message. If the options argument value is not specified, computes canonical ChemAxon-compatible extended SMILES.

Options:
   -h, --help:  Show this help message and exit
Output Encoding:
   --out OUT:  The output molecule encoding ((smiles|smi)|(molfile|mo lblock)|molhash|molhashfull|tautohash|inchi|inchikey)
Transform:
   --desalt:  Remove (strip) salt
   --desalt_smarts_list DESALT_SMARTS_LIST:  Additional desalt SMARTS patterns separated with |
   --remove_stereo:  Remove stereo
SMILES encoder options:
   --smiles_kekule: Kekulize the molecule before generating the SMILES and output single/double rather than aromatic bonds'
     AS
$$
-- @include python/Util.py
-- @include python/Molstring_Or_Molbinary2Molstring.py
$$
;


CREATE OR REPLACE FUNCTION Molstring_Or_Molbinary2Molbinary(molstring VARCHAR DEFAULT NULL, molbinary VARBINARY DEFAULT NULL, options VARCHAR DEFAULT '')
     RETURNS VARBINARY 
     LANGUAGE PYTHON 
     IMMUTABLE    
-- @include python/Runtime_Repo.inc
     PACKAGES = (
-- @include python/Rdkit_Pkg.inc
     )
     HANDLER = 'molstring_or_molbinary_to_molbinary'
     COMMENT='Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) or RDKit binary molecule encoding to the RDKit binary molecule encoding, optionally transforming the input structure. Only one of molstring or molbinary arguments can be non-NULL. All available options can be listed by by running the following SQL: select Molstring_Or_Molbinary2Molbinary(NULL, NULL, ''-h''); usage info will be returned as part of the error message. To convert molstrings to RDKit binary molecule encoding w/o applying any transforms, use the Molstring2Molbinary function, which has fewer arguments and is faster.

Options: -h, --help show this help message and exit
Transform:
   --desalt:  Remove (strip) salt
   --desalt_smarts_list DESALT_SMARTS_LIST:  Additional desalt SMARTS patterns separated with |
   --remove_stereo:  Remove stereo'
     AS
$$
-- @include python/Util.py
-- @include python/Molstring_Or_Molbinary2Molbinary.py
$$
;


CREATE OR REPLACE FUNCTION Molstring2Molbinary(molstring VARCHAR)
     RETURNS VARBINARY 
     LANGUAGE PYTHON 
     RETURNS NULL ON NULL INPUT
     IMMUTABLE    
-- @include python/Runtime_Repo.inc
     PACKAGES = (
-- @include python/Rdkit_Pkg.inc
     )
     HANDLER = 'molstring_to_binary'
     COMMENT='Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) to RDKit binary molecule encoding. Returns NULL if molstring is NULL, empty, invalid, or represents an empty molecule with 0 atoms and 0 bonds.'
     AS
$$
-- @include python/Util.py
-- @include python/Molstring2Molbinary.py
$$
;

-- Category: Fingerprinting

CREATE OR REPLACE FUNCTION Molstring2Pattern_Fingerprint(molstring VARCHAR)
     RETURNS VARBINARY 
     LANGUAGE PYTHON 
     RETURNS NULL ON NULL INPUT
     IMMUTABLE
-- @include python/Runtime_Repo.inc
     PACKAGES = (
-- @include python/Rdkit_Pkg.inc
     )
     HANDLER = 'molstring2pattern_fingerprint'
     COMMENT='Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) to RDKit substructure pattern fingerprint commonly used for fingerprint-based screening to speed up substructure searches. Returns NULL if molstring is NULL, empty, invalid, or represents an empty molecule with 0 atoms and 0 bonds.'
     AS
$$
-- @include python/Util.py
-- @include python/Molstring2Pattern_Fingerprint.py
$$
;


CREATE OR REPLACE FUNCTION Molbinary2Pattern_Fingerprint(molbinary VARBINARY)
     RETURNS VARBINARY 
     LANGUAGE PYTHON 
     RETURNS NULL ON NULL INPUT
     IMMUTABLE
-- @include python/Runtime_Repo.inc
     PACKAGES = (
-- @include python/Rdkit_Pkg.inc
     )
     HANDLER = 'molbinary2pattern_fingerprint'
     COMMENT='Converts RDKit binary-encoded molecule to RDKit substructure pattern fingerprint commonly used for fingerprint-based screening to speed up substructure searches. Returns NULL if molbinary is NULL, empty or represents an empty molecule with 0 atoms and 0 bonds. Returns error if molbinary is not a valid RDKit binary-encoded molecule.'     
     AS
$$
-- @include python/Util.py
-- @include python/Molbinary2Pattern_Fingerprint.py
$$
;


CREATE OR REPLACE FUNCTION Smarts2Pattern_Fingerprint(smarts VARCHAR)
     RETURNS VARBINARY 
     LANGUAGE PYTHON 
     RETURNS NULL ON NULL INPUT
     IMMUTABLE
-- @include python/Runtime_Repo.inc
     PACKAGES = (
-- @include python/Rdkit_Pkg.inc
     )
     HANDLER = 'smarts2pattern_fingerprint'
     COMMENT='Converts SMARTS substructure pattern to RDKit substructure pattern fingerprint commonly used for fingerprint-based screening to speed up substructure searches. Returns NULL if smarts is NULL, empty, or represents an empty molecule with 0 atoms and 0 bonds. Returns error if smarts is invalid and cannot be parsed.'
     AS
$$
-- @include python/Util.py
-- @include python/Smarts2Pattern_Fingerprint.py
$$
;

CREATE OR REPLACE FUNCTION Molstring2Morgan_Fingerprint(molstring VARCHAR)
     RETURNS VARBINARY 
     LANGUAGE PYTHON 
     RETURNS NULL ON NULL INPUT
     IMMUTABLE
-- @include python/Runtime_Repo.inc
     PACKAGES = (
-- @include python/Rdkit_Pkg.inc
     )
     HANDLER = 'molstring2morgan_fingerprint'
     COMMENT='Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) to RDKit Morgan fingerprint commonly used for fingerprint-based similarity search. Returns NULL if molstring is NULL, empty, invalid, or represents an empty molecule with 0 atoms and 0 bonds. Generator options: radius=2, fpSize=2048, atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen()'
     AS
$$
-- @include python/Util.py
-- @include python/Molstring2Morgan_Fingerprint.py
$$
;


CREATE OR REPLACE FUNCTION Molbinary2Morgan_Fingerprint(molbinary VARBINARY)
     RETURNS VARBINARY 
     LANGUAGE PYTHON 
     RETURNS NULL ON NULL INPUT
     IMMUTABLE
-- @include python/Runtime_Repo.inc
     PACKAGES = (
-- @include python/Rdkit_Pkg.inc
     )
     HANDLER = 'molbinary2morgan_fingerprint'
     COMMENT='Converts RDKit binary-encoded molecule to RDKit Morgan fingerprint commonly used for fingerprint-based similarity search. Returns NULL if molbinary is NULL or empty, or represents an empty molecule with 0 atoms and 0 bonds. Generator options: radius=2, fpSize=2048, atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen(). Returns error if molbinary is not a valid RDKit binary-encoded molecule.'          
     AS
$$
-- @include python/Util.py
-- @include python/Molbinary2Morgan_Fingerprint.py
$$
;


-- Category: Substructure matching

CREATE OR REPLACE FUNCTION Molbinary_Matches_Smarts(molbinary VARBINARY, smarts VARCHAR, screen_pass BOOLEAN)
     RETURNS BOOLEAN 
     LANGUAGE PYTHON 
     RETURNS NULL ON NULL INPUT
     IMMUTABLE    
-- @include python/Runtime_Repo.inc
     PACKAGES = (
-- @include python/Rdkit_Pkg.inc
     )
     HANDLER = 'molbinary_matches_smarts'
     COMMENT='Tests whether RDKit binary-encoded molecule matches the specified SMARTS pattern and returns TRUE iff it does and the screen_pass argument value is TRUE. The extra screen_pass argument is used for substructure query optimization based on fingerprint screening (see examples in the documentation and example workbooks). Returns NULL if any of the args are NULL. Returns error if molbinary is not a valid RDKit binary-encoded molecule of if smarts is invalid and cannot be parsed. Note: an empty SMARTS will not match any molecule, even an empty one. This seems to be illogical, since, in theory, a subgraph with 0 nodes and 0 edges must match any graph (or, at least, an empty one), but, in practice, this approach leads to fewer problems than the theoretically correct one. In RDKit itself, a molecule representing an empty pattern does not match anything either.'
     AS
$$
-- @include python/Util.py
-- @include python/Molbinary_Matches_Smarts.py
$$
;


CREATE OR REPLACE FUNCTION Molstring_Matches_Smarts(molstring VARCHAR, smarts VARCHAR, screen_pass BOOLEAN)
     RETURNS BOOLEAN 
     LANGUAGE PYTHON 
     RETURNS NULL ON NULL INPUT
     IMMUTABLE    
-- @include python/Runtime_Repo.inc
     PACKAGES = (
-- @include python/Rdkit_Pkg.inc
     )
     HANDLER = 'molstring_matches_smarts'
     COMMENT='Tests whether the molecule encoded as molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) matches the specified SMARTS pattern and returns TRUE iff it does and the screen_pass argument value is TRUE. The extra screen_pass argument is used for substructure query optimization based on fingerprint screening (see examples in the documentation and example workbooks). Returns NULL if any of the args are NULL. Returns FALSE if molstring is not a valid SMILES or molblock. Returns error if smarts is invalid and cannot be parsed. Note: an empty SMARTS will not match any molecule, even an empty one. This seems to be illogical, since, in theory, a subgraph with 0 nodes and 0 edges must match any graph (or, at least, an empty one), but, in practice, this approach leads to fewer problems than the theoretically correct one. In RDKit itself, a molecule representing an empty pattern does not match anything either.'     
     AS
$$
-- @include python/Util.py
-- @include python/Molstring_Matches_Smarts.py
$$
;


-- Category: Structure similarity

CREATE OR REPLACE FUNCTION Tanimoto_P(v1 VARBINARY, v2 VARBINARY)
     RETURNS FLOAT
     LANGUAGE PYTHON
     RETURNS NULL ON NULL INPUT
     IMMUTABLE
-- @include python/Runtime_Repo.inc
     PACKAGES = ('bitarray')
     HANDLER = 'tanimoto'
     COMMENT='A version of TANIMOTO implemented in Python. Computes Tanimoto similarity between two bitsets represented by binary vectors. Returns BITCOUNT(v1 bitand v2) / BITCOUNT(v1 bitor v2) as float. If both v1 and v2 have no bits set to 1, returns 1.0, that is, treats two empty bitsets or two bitsets filled with only zeroes as being equal to each other. Returns error if two bitsets are of different lengths. Returns NULL if one of both arguments are NULLs'
     AS
$$
-- @include python/Tanimoto.py
$$
;


-- A version of Tanimoto implemented in Java.
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
-- @include java/Tanimoto.java
$$
;

-- This implementation relies on the version of the built-in BITCOUNT function
-- that can correctly handle VARBINARY arguments and is available in Snowflake update bundle 2026_07
CREATE OR REPLACE FUNCTION Tanimoto(v1 VARBINARY, v2 VARBINARY)
     RETURNS FLOAT
     RETURNS NULL ON NULL INPUT
     COMMENT='Computes Tanimoto similarity between two bitsets represented by binary vectors. Returns BITCOUNT(v1 bitand v2) / BITCOUNT(v1 bitor v2) as float. If both v1 and v2 have no bits set to 1, returns 0.0, that is, treats two 0-length bitsets or two bitsets filled with only zeroes as being different from each other. The lengths of v1 and v2 must be the same, otheerwise, an exception will be raised. Returns NULL if one of both arguments are NULLs'
     AS
$$
    DIV0(BITCOUNT(BITAND(v1, v2)), (BITCOUNT(v1) + BITCOUNT(v2) - BITCOUNT(BITAND(v1, v2))))::float
$$
;


-- Category: Molecular properties and descriptor calculation

CREATE OR REPLACE FUNCTION Molstring_Or_Molbinary2Medchem_Descriptors(molstring VARCHAR DEFAULT NULL, molbinary VARBINARY DEFAULT NULL)
     RETURNS TABLE (Q_Estim_Drug_Likeness float, Mol_Wt float, Heavy_Atom_Mol_Wt float, Exact_Mol_Wt float, Num_Valence_Electrons int, Num_Radical_Electrons int, Max_Partial_Charge float, Min_Partial_Charge float, Max_Abs_Partial_Charge float, Min_Abs_Partial_Charge float, TPSA float, Fraction_CSP3 float, Heavy_Atom_Count int, Num_Aliphatic_Carbocycles int, Num_Aliphatic_Heterocycles int, Num_Aliphatic_Rings int, Num_Aromatic_Carbocycles int, Num_Aromatic_Heterocycles int, Num_Aromatic_Rings int, Num_HAcceptors int, Num_HDonors int, Num_Heteroatoms int, Num_Rotatable_Bonds int, Num_Saturated_Carbocycles int, Num_Saturated_Heterocycles int, Num_Saturated_Rings int, Ring_Count int, Mol_Log_P float, Mol_MR float)
     LANGUAGE PYTHON
     IMMUTABLE
-- @include python/Runtime_Repo.inc
     PACKAGES = (
-- @include python/Rdkit_Pkg.inc
     )
     HANDLER = 'DGen'
     COMMENT='Computes commonly used medchem properties (descriptors) for a molecule encoded as molstring (standard or ChemAxon SMILES or molfile/molblock) or as RDKit binary-encoded molecule. Only one of molstring or molbinary arguments can be non-NULL, otherwise, an error will be returned. Returns a table with one row and multiple columns corresponding to the computed descriptors. If molstring and molbinary are both NULLs, or molstring is invalid and cannot be parsed, returns a table with one row filled with NULLs. Returns an error if molbinary is not a valid RDKit binary-encoded molecule.'
     AS
$$
-- @include python/Util.py
-- @include python/Molstring_Or_Molbinary2Medchem_Descriptors.py
$$
;

-- Category: Structure encoding and molecular graph validity check

CREATE OR REPLACE FUNCTION Molstring_Or_Molbinary_Check(molstring VARCHAR DEFAULT NULL, molbinary VARBINARY DEFAULT NULL, raise_exception BOOLEAN DEFAULT FALSE)
     RETURNS TABLE (is_ok boolean, encoding varchar, error_msg varchar) 
     LANGUAGE PYTHON 
     IMMUTABLE    
-- @include python/Runtime_Repo.inc
     PACKAGES = (
-- @include python/Rdkit_Pkg.inc
     )
     HANDLER = 'MolChecker'
     COMMENT='Checks a molecule encoded as molstring (standard or ChemAxon SMILES or molfile/molblock) or as RDKit binary-encoded molecule. Only one of molstring or molbinary arguments can be non-NULL, otherwise, an error will be returned. Returns a table with one row and three columns: is_ok boolean, encoding varchar, and error_msg varchar. If both molstring and molbinary are NULL, the entire result row will be filled with NULLs. Otherwise, is_ok will contain True iff the input can be parsed into a molecule with no errors, encoding will contain a string representation of the encoding (MOLBLOCK, SMILES, or BINARY), and the error_msg will contain a description of the error or NULL. If raise_exception parameter is TRUE (it is FALSE by default) and the input cannot be parsed into a valid molecule, the method will raise an exception and quit instead of returning.'
     AS
$$
-- @include python/Util.py
-- @include python/Molstring_Or_Molbinary_Check.py
$$
;

-- Category: HELM support (experimental)

CREATE OR REPLACE PROCEDURE helm2molbinary(monomer_table VARCHAR, helm_table VARCHAR)
     RETURNS table(ID INT, MOLBINARY VARBINARY)
     LANGUAGE PYTHON
-- RUNTIME_VERSION = 3.12
-- @include python/Runtime_Repo.inc
     PACKAGES = (
-- @include python/Snowpark_Pkg.inc
     ,
-- @include python/Rdkit_Pkg.inc
     )
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

-- @include python/HelmMolecule.py

-- @include python/Helm2Molbinary.py

$$;


-- Category: Single-input convenience wrappers for Molstring_Or_Molbinary2* UDFs

CREATE OR REPLACE FUNCTION Conv_Mol2Molstring(molstring VARCHAR, options VARCHAR DEFAULT '--out smiles')
     RETURNS VARCHAR
     LANGUAGE SQL
     IMMUTABLE
     COMMENT='Wrapper for Molstring_Or_Molbinary2Molstring(molstring, NULL, options).'
     AS
$$
    Molstring_Or_Molbinary2Molstring(molstring, NULL::VARBINARY, options)
$$;

CREATE OR REPLACE FUNCTION Conv_Mol2Molstring(molbinary VARBINARY, options VARCHAR DEFAULT '--out smiles')
     RETURNS VARCHAR
     LANGUAGE SQL
     IMMUTABLE
     COMMENT='Wrapper for Molstring_Or_Molbinary2Molstring(NULL, molbinary, options).'
     AS
$$
    Molstring_Or_Molbinary2Molstring(NULL::VARCHAR, molbinary, options)
$$;

CREATE OR REPLACE FUNCTION Conv_Mol2Molbinary(molstring VARCHAR, options VARCHAR DEFAULT '')
     RETURNS VARBINARY
     LANGUAGE SQL
     IMMUTABLE
     COMMENT='Wrapper for Molstring_Or_Molbinary2Molbinary(molstring, NULL, options).'
     AS
$$
    Molstring_Or_Molbinary2Molbinary(molstring, NULL::VARBINARY, options)
$$;

CREATE OR REPLACE FUNCTION Conv_Mol2Molbinary(molbinary VARBINARY, options VARCHAR DEFAULT '')
     RETURNS VARBINARY
     LANGUAGE SQL
     IMMUTABLE
     COMMENT='Wrapper for Molstring_Or_Molbinary2Molbinary(NULL, molbinary, options).'
     AS
$$
    Molstring_Or_Molbinary2Molbinary(NULL::VARCHAR, molbinary, options)
$$;


CREATE OR REPLACE FUNCTION Calc_Mol_Medchem_Descriptors(molstring VARCHAR)
     RETURNS TABLE (Q_Estim_Drug_Likeness float, Mol_Wt float, Heavy_Atom_Mol_Wt float, Exact_Mol_Wt float, Num_Valence_Electrons int, Num_Radical_Electrons int, Max_Partial_Charge float, Min_Partial_Charge float, Max_Abs_Partial_Charge float, Min_Abs_Partial_Charge float, TPSA float, Fraction_CSP3 float, Heavy_Atom_Count int, Num_Aliphatic_Carbocycles int, Num_Aliphatic_Heterocycles int, Num_Aliphatic_Rings int, Num_Aromatic_Carbocycles int, Num_Aromatic_Heterocycles int, Num_Aromatic_Rings int, Num_HAcceptors int, Num_HDonors int, Num_Heteroatoms int, Num_Rotatable_Bonds int, Num_Saturated_Carbocycles int, Num_Saturated_Heterocycles int, Num_Saturated_Rings int, Ring_Count int, Mol_Log_P float, Mol_MR float)
     LANGUAGE SQL
     IMMUTABLE
     COMMENT='Wrapper for Molstring_Or_Molbinary2Medchem_Descriptors(molstring, NULL).'
     AS
$$
    SELECT * FROM TABLE(Molstring_Or_Molbinary2Medchem_Descriptors(molstring, NULL::VARBINARY))
$$;

CREATE OR REPLACE FUNCTION Calc_Mol_Medchem_Descriptors(molbinary VARBINARY)
     RETURNS TABLE (Q_Estim_Drug_Likeness float, Mol_Wt float, Heavy_Atom_Mol_Wt float, Exact_Mol_Wt float, Num_Valence_Electrons int, Num_Radical_Electrons int, Max_Partial_Charge float, Min_Partial_Charge float, Max_Abs_Partial_Charge float, Min_Abs_Partial_Charge float, TPSA float, Fraction_CSP3 float, Heavy_Atom_Count int, Num_Aliphatic_Carbocycles int, Num_Aliphatic_Heterocycles int, Num_Aliphatic_Rings int, Num_Aromatic_Carbocycles int, Num_Aromatic_Heterocycles int, Num_Aromatic_Rings int, Num_HAcceptors int, Num_HDonors int, Num_Heteroatoms int, Num_Rotatable_Bonds int, Num_Saturated_Carbocycles int, Num_Saturated_Heterocycles int, Num_Saturated_Rings int, Ring_Count int, Mol_Log_P float, Mol_MR float)
     LANGUAGE SQL
     IMMUTABLE
     COMMENT='Wrapper for Molstring_Or_Molbinary2Medchem_Descriptors(NULL, molbinary).'
     AS
$$
    SELECT * FROM TABLE(Molstring_Or_Molbinary2Medchem_Descriptors(NULL::VARCHAR, molbinary))
$$;

-- Category: Single-input convenience wrappers for Draw_Molstring_Or_Molbinary2* UDFs

CREATE OR REPLACE FUNCTION Draw_Mol2SVG(molstring VARCHAR, highlight_smarts VARCHAR DEFAULT NULL, draw_options VARCHAR DEFAULT NULL, output_raw_xml BOOLEAN DEFAULT FALSE)
     RETURNS VARCHAR
     LANGUAGE SQL
     IMMUTABLE
     COMMENT='Wrapper for Draw_Molstring_Or_Molbinary2SVG(molstring, NULL, highlight_smarts, draw_options, output_raw_xml).'
     AS
$$
    Draw_Molstring_Or_Molbinary2SVG(molstring, NULL::VARBINARY, highlight_smarts, draw_options, output_raw_xml)
$$;

CREATE OR REPLACE FUNCTION Draw_Mol2SVG(molbinary VARBINARY, highlight_smarts VARCHAR DEFAULT NULL, draw_options VARCHAR DEFAULT NULL, output_raw_xml BOOLEAN DEFAULT FALSE)
     RETURNS VARCHAR
     LANGUAGE SQL
     IMMUTABLE
     COMMENT='Wrapper for Draw_Molstring_Or_Molbinary2SVG(NULL, molbinary, highlight_smarts, draw_options, output_raw_xml).'
     AS
$$
    Draw_Molstring_Or_Molbinary2SVG(NULL::VARCHAR, molbinary, highlight_smarts, draw_options, output_raw_xml)
$$;

CREATE OR REPLACE FUNCTION Draw_Mol2PNG(width_px INT, height_px INT, molstring VARCHAR, highlight_smarts VARCHAR DEFAULT NULL, draw_options VARCHAR DEFAULT NULL)
     RETURNS VARBINARY
     LANGUAGE SQL
     IMMUTABLE
     COMMENT='Wrapper for Draw_Molstring_Or_Molbinary2PNG(width_px, height_px, molstring, NULL, highlight_smarts, draw_options).'
     AS
$$
    Draw_Molstring_Or_Molbinary2PNG(width_px, height_px, molstring, NULL::VARBINARY, highlight_smarts, draw_options)
$$;

CREATE OR REPLACE FUNCTION Draw_Mol2PNG(width_px INT, height_px INT, molbinary VARBINARY, highlight_smarts VARCHAR DEFAULT NULL, draw_options VARCHAR DEFAULT NULL)
     RETURNS VARBINARY
     LANGUAGE SQL
     IMMUTABLE
     COMMENT='Wrapper for Draw_Molstring_Or_Molbinary2PNG(width_px, height_px, NULL, molbinary, highlight_smarts, draw_options).'
     AS
$$
    Draw_Molstring_Or_Molbinary2PNG(width_px, height_px, NULL::VARCHAR, molbinary, highlight_smarts, draw_options)
$$;
