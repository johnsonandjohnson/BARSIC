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
-- @include python/Util.py
-- @include python/Molstring_Or_Molbinary2Molstring.py
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
-- @include python/Util.py
-- @include python/Molstring_Or_Molbinary2Molbinary.py
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
-- @include python/Util.py
-- @include python/Molstring2Molbinary.py
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
-- @include python/Util.py
-- @include python/Molstring2Pattern_Fingerprint.py
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
-- @include python/Util.py
-- @include python/Molbinary2Pattern_Fingerprint.py
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
-- @include python/Util.py
-- @include python/Smarts2Pattern_Fingerprint.py
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
-- @include python/Util.py
-- @include python/Molstring2Morgan_Fingerprint.py
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
-- @include python/Util.py
-- @include python/Molbinary2Morgan_Fingerprint.py
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
-- @include python/Util.py
-- @include python/Molbinary_Matches_Smarts.py
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
-- @include python/Util.py
-- @include python/Molstring_Matches_Smarts.py
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
-- @include python/Tanimoto.py
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
-- @include java/Tanimoto.java
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
-- @include python/Util.py
-- @include python/Molstring_Or_Molbinary2Medchem_Descriptors.py
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
-- @include python/Util.py
-- @include python/Molstring_Or_Molbinary_Check.py
$$
;
