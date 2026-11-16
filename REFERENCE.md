# Chem UDF Reference

## Quick Reference

| # | Function | Arguments | Returns | Purpose |
| --- | --- | --- | --- | --- |
| 1 | [Barsic_Info](#1-barsic_info) | — | VARCHAR | Toolkit and library versions |
| 2 | [Draw_Molstring_Or_Molbinary2SVG](#2-draw_molstring_or_molbinary2svg) | `[molstring]`, `[molbinary]`, `[highlight_smarts]`, `[draw_options]`, `[output_raw_xml]` | VARCHAR | Draw molecule as SVG (string or binary input) |
| 3 | [Draw_Molstring_Or_Molbinary2PNG](#3-draw_molstring_or_molbinary2png) | `width_px`, `height_px`, `[molstring]`, `[molbinary]`, `[highlight_smarts]`, `[draw_options]` | VARBINARY | Draw molecule as PNG (string or binary input) |
| 4 | [Molstring_Or_Molbinary2Molstring](#4-molstring_or_molbinary2molstring) | `[molstring]`, `[molbinary]`, `[options]` | VARCHAR | Convert to SMILES/molfile/InChI/InChIKey/hashes, with optional transforms |
| 5 | [Molstring_Or_Molbinary2Molbinary](#5-molstring_or_molbinary2molbinary) | `[molstring]`, `[molbinary]`, `[options]` | VARBINARY | Convert to RDKit binary, with optional transforms |
| 6 | [Molstring2Molbinary](#6-molstring2molbinary) | `molstring` | VARBINARY | Fast molstring → RDKit binary (no transforms) |
| 7 | [Molstring2Pattern_Fingerprint](#7-molstring2pattern_fingerprint) | `molstring` | VARBINARY | Substructure-screening fingerprint from molstring |
| 8 | [Molbinary2Pattern_Fingerprint](#8-molbinary2pattern_fingerprint) | `molbinary` | VARBINARY | Substructure-screening fingerprint from binary |
| 9 | [Smarts2Pattern_Fingerprint](#9-smarts2pattern_fingerprint) | `smarts` | VARBINARY | Substructure-screening fingerprint from SMARTS query |
| 10 | [Molstring2Morgan_Fingerprint](#10-molstring2morgan_fingerprint) | `molstring` | VARBINARY | Morgan similarity fingerprint from molstring |
| 11 | [Molbinary2Morgan_Fingerprint](#11-molbinary2morgan_fingerprint) | `molbinary` | VARBINARY | Morgan similarity fingerprint from binary |
| 12 | [Molbinary_Matches_Smarts](#12-molbinary_matches_smarts) | `molbinary`, `smarts`, `screen_pass` | BOOLEAN | SMARTS substructure match (binary input) |
| 13 | [Molstring_Matches_Smarts](#13-molstring_matches_smarts) | `molstring`, `smarts`, `screen_pass` | BOOLEAN | SMARTS substructure match (molstring input) |
| 14 | [Tanimoto_P](#14-tanimoto_p) | `v1`, `v2` | FLOAT | Tanimoto similarity (Python) |
| 15 | [Tanimoto_J](#15-tanimoto_j) | `v1`, `v2` | FLOAT | Tanimoto similarity (Java) |
| 16 | [Tanimoto](#16-tanimoto) | `v1`, `v2` | FLOAT | Tanimoto similarity (SQL) |
| 17 | [Molstring_Or_Molbinary2Medchem_Descriptors](#17-molstring_or_molbinary2medchem_descriptors) | `[molstring]`, `[molbinary]` | TABLE | 29 medchem descriptors (MW, TPSA, logP, …) |
| 18 | [Molstring_Or_Molbinary_Check](#18-molstring_or_molbinary_check) | `[molstring]`, `[molbinary]`, `[raise_exception]` | TABLE | Validate molecule and detect encoding |
| 19 | [helm2molbinary](#19-helm2molbinary) | `monomer_table`, `helm_table` | TABLE | HELM → RDKit binary (experimental) |
| 20 | [Conv_Mol2Molstring](#20-conv_mol2molstring) | `molstring`, `[options]` | VARCHAR | Molstring → molstring |
| 21 | [Conv_Mol2Molstring](#21-conv_mol2molstring) | `molbinary`, `[options]` | VARCHAR | Binary → molstring |
| 22 | [Conv_Mol2Molbinary](#22-conv_mol2molbinary) | `molstring`, `[options]` | VARBINARY | Molstring → binary |
| 23 | [Conv_Mol2Molbinary](#23-conv_mol2molbinary) | `molbinary`, `[options]` | VARBINARY | Binary → binary |
| 24 | [Calc_Mol_Medchem_Descriptors](#24-calc_mol_medchem_descriptors) | `molstring` | TABLE | Medchem descriptors from molstring |
| 25 | [Calc_Mol_Medchem_Descriptors](#25-calc_mol_medchem_descriptors) | `molbinary` | TABLE | Medchem descriptors from binary |
| 26 | [Draw_Mol2SVG](#26-draw_mol2svg) | `molstring`, `[highlight_smarts]`, `[draw_options]`, `[output_raw_xml]` | VARCHAR | SVG from molstring |
| 27 | [Draw_Mol2SVG](#27-draw_mol2svg) | `molbinary`, `[highlight_smarts]`, `[draw_options]`, `[output_raw_xml]` | VARCHAR | SVG from binary |
| 28 | [Draw_Mol2PNG](#28-draw_mol2png) | `width_px`, `height_px`, `molstring`, `[highlight_smarts]`, `[draw_options]` | VARBINARY | PNG from molstring |
| 29 | [Draw_Mol2PNG](#29-draw_mol2png) | `width_px`, `height_px`, `molbinary`, `[highlight_smarts]`, `[draw_options]` | VARBINARY | PNG from binary |

`[arg]` = optional (has a default).

## Index

### Info

| # | Name | Kind | Language | Returns |
| --- | --- | --- | --- | --- |
| 1 | [Barsic_Info](#1-barsic_info) | Function | PYTHON | `VARCHAR` |

### Chemical structure depiction

| # | Name | Kind | Language | Returns |
| --- | --- | --- | --- | --- |
| 2 | [Draw_Molstring_Or_Molbinary2SVG](#2-draw_molstring_or_molbinary2svg) | Function | PYTHON | `VARCHAR` |
| 3 | [Draw_Molstring_Or_Molbinary2PNG](#3-draw_molstring_or_molbinary2png) | Function | PYTHON | `VARBINARY` |

### Chemical structure conversion

| # | Name | Kind | Language | Returns |
| --- | --- | --- | --- | --- |
| 4 | [Molstring_Or_Molbinary2Molstring](#4-molstring_or_molbinary2molstring) | Function | PYTHON | `VARCHAR` |
| 5 | [Molstring_Or_Molbinary2Molbinary](#5-molstring_or_molbinary2molbinary) | Function | PYTHON | `VARBINARY` |
| 6 | [Molstring2Molbinary](#6-molstring2molbinary) | Function | PYTHON | `VARBINARY` |

### Fingerprinting

| # | Name | Kind | Language | Returns |
| --- | --- | --- | --- | --- |
| 7 | [Molstring2Pattern_Fingerprint](#7-molstring2pattern_fingerprint) | Function | PYTHON | `VARBINARY` |
| 8 | [Molbinary2Pattern_Fingerprint](#8-molbinary2pattern_fingerprint) | Function | PYTHON | `VARBINARY` |
| 9 | [Smarts2Pattern_Fingerprint](#9-smarts2pattern_fingerprint) | Function | PYTHON | `VARBINARY` |
| 10 | [Molstring2Morgan_Fingerprint](#10-molstring2morgan_fingerprint) | Function | PYTHON | `VARBINARY` |
| 11 | [Molbinary2Morgan_Fingerprint](#11-molbinary2morgan_fingerprint) | Function | PYTHON | `VARBINARY` |

### Substructure matching

| # | Name | Kind | Language | Returns |
| --- | --- | --- | --- | --- |
| 12 | [Molbinary_Matches_Smarts](#12-molbinary_matches_smarts) | Function | PYTHON | `BOOLEAN` |
| 13 | [Molstring_Matches_Smarts](#13-molstring_matches_smarts) | Function | PYTHON | `BOOLEAN` |

### Structure similarity

| # | Name | Kind | Language | Returns |
| --- | --- | --- | --- | --- |
| 14 | [Tanimoto_P](#14-tanimoto_p) | Function | PYTHON | `FLOAT` |
| 15 | [Tanimoto_J](#15-tanimoto_j) | Function | JAVA | `FLOAT` |
| 16 | [Tanimoto](#16-tanimoto) | Function | SQL | `FLOAT` |

### Molecular properties and descriptor calculation

| # | Name | Kind | Language | Returns |
| --- | --- | --- | --- | --- |
| 17 | [Molstring_Or_Molbinary2Medchem_Descriptors](#17-molstring_or_molbinary2medchem_descriptors) | Function | PYTHON | `TABLE(...)` |

### Structure encoding and molecular graph validity check

| # | Name | Kind | Language | Returns |
| --- | --- | --- | --- | --- |
| 18 | [Molstring_Or_Molbinary_Check](#18-molstring_or_molbinary_check) | Function | PYTHON | `TABLE(...)` |

### HELM support (experimental)

| # | Name | Kind | Language | Returns |
| --- | --- | --- | --- | --- |
| 19 | [helm2molbinary](#19-helm2molbinary) | Procedure | PYTHON | `TABLE(...)` |

### Single-input convenience wrappers for Molstring_Or_Molbinary2* UDFs

| # | Name | Kind | Language | Returns |
| --- | --- | --- | --- | --- |
| 20 | [Conv_Mol2Molstring](#20-conv_mol2molstring) | Function | SQL | `VARCHAR` |
| 21 | [Conv_Mol2Molstring](#21-conv_mol2molstring) | Function | SQL | `VARCHAR` |
| 22 | [Conv_Mol2Molbinary](#22-conv_mol2molbinary) | Function | SQL | `VARBINARY` |
| 23 | [Conv_Mol2Molbinary](#23-conv_mol2molbinary) | Function | SQL | `VARBINARY` |
| 24 | [Calc_Mol_Medchem_Descriptors](#24-calc_mol_medchem_descriptors) | Function | SQL | `TABLE(...)` |
| 25 | [Calc_Mol_Medchem_Descriptors](#25-calc_mol_medchem_descriptors) | Function | SQL | `TABLE(...)` |

### Single-input convenience wrappers for Draw_Molstring_Or_Molbinary2* UDFs

| # | Name | Kind | Language | Returns |
| --- | --- | --- | --- | --- |
| 26 | [Draw_Mol2SVG](#26-draw_mol2svg) | Function | SQL | `VARCHAR` |
| 27 | [Draw_Mol2SVG](#27-draw_mol2svg) | Function | SQL | `VARCHAR` |
| 28 | [Draw_Mol2PNG](#28-draw_mol2png) | Function | SQL | `VARBINARY` |
| 29 | [Draw_Mol2PNG](#29-draw_mol2png) | Function | SQL | `VARBINARY` |

## Info

### 1. Barsic_Info

**Kind:** Function | **Language:** PYTHON

```sql
Barsic_Info()
  RETURNS VARCHAR
```

**Parameters:**

_None._

**Returns:** VARCHAR — Toolkit, Python, RDKit and Snowpark versions.

**Description:**

Returns toolkit info

**Example:**

```sql
SELECT Barsic_Info();
-- BARSIC 3.1, Python: 3.12.13 ..., RDKit: 2026.09.2, snowflake-snowpark-python: 1.57.0
```

## Chemical structure depiction

### 2. Draw_Molstring_Or_Molbinary2SVG

**Kind:** Function | **Language:** PYTHON

```sql
Draw_Molstring_Or_Molbinary2SVG(molstring VARCHAR DEFAULT NULL, molbinary VARBINARY DEFAULT NULL, highlight_smarts VARCHAR DEFAULT NULL, draw_options VARCHAR DEFAULT NULL, output_raw_xml BOOLEAN DEFAULT FALSE)
  RETURNS VARCHAR
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molstring` | VARCHAR | `NULL` | Molecule as a string: standard or ChemAxon-extended SMILES, or molblock/molfile (auto-detected). |
| `molbinary` | VARBINARY | `NULL` | Molecule in RDKit binary encoding (e.g. from `Molstring2Molbinary`). |
| `highlight_smarts` | VARCHAR | `NULL` | Optional SMARTS pattern; atoms matching it are highlighted. |
| `draw_options` | VARCHAR | `NULL` | Reserved for future use; currently ignored. |
| `output_raw_xml` | BOOLEAN | `FALSE` | If TRUE, returns plain SVG XML; if FALSE, a `data:image/svg+xml;base64,` data URI. |

**Returns:** VARCHAR — SVG image (data URI or raw XML).

**Description:**

Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) or RDKit binary molecule encoding (only one of molstring and molbinary can be not NULL) to its graphical representation in the SVG format. If highlight_smarts substructure pattern is specified, highlights matching atoms. If output_raw_xml is FALSE, returns picture as varchar with the following MIME type: data:image/svg+xml;base64, otherwise, returns plain unencoded xml. draw_options arg is not used yet and is reserved for future use.

**Example:**

```sql
SELECT Draw_Molstring_Or_Molbinary2SVG(molstring => 'CC(=O)Oc1ccccc1C(=O)O', highlight_smarts => 'c1ccccc1');
-- data:image/svg+xml;base64,PHN2ZyB2ZXJzaW...
```

### 3. Draw_Molstring_Or_Molbinary2PNG

**Kind:** Function | **Language:** PYTHON

```sql
Draw_Molstring_Or_Molbinary2PNG(width_px INT, height_px INT, molstring VARCHAR DEFAULT NULL, molbinary VARBINARY DEFAULT NULL, highlight_smarts VARCHAR DEFAULT NULL, draw_options VARCHAR DEFAULT NULL)
  RETURNS VARBINARY
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `width_px` | INT | _required_ | Image width in pixels. |
| `height_px` | INT | _required_ | Image height in pixels. |
| `molstring` | VARCHAR | `NULL` | Molecule as a string: standard or ChemAxon-extended SMILES, or molblock/molfile (auto-detected). |
| `molbinary` | VARBINARY | `NULL` | Molecule in RDKit binary encoding (e.g. from `Molstring2Molbinary`). |
| `highlight_smarts` | VARCHAR | `NULL` | Optional SMARTS pattern; atoms matching it are highlighted. |
| `draw_options` | VARCHAR | `NULL` | Reserved for future use; currently ignored. |

**Returns:** VARBINARY — PNG image bytes.

**Description:**

Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) or RDKit binary molecule encoding (only one of molstring and molbinary can be not NULL) to its graphical representation in the PNG format (width_px x height_px is the image size in pixels). If highlight_smarts substructure pattern is specified, highlights matching atoms. draw_options arg is not used yet and is reserved for future use.

**Example:**

```sql
SELECT Draw_Molstring_Or_Molbinary2PNG(300, 200, molstring => 'CC(=O)Oc1ccccc1C(=O)O');
```

## Chemical structure conversion

### 4. Molstring_Or_Molbinary2Molstring

**Kind:** Function | **Language:** PYTHON

```sql
Molstring_Or_Molbinary2Molstring(molstring VARCHAR DEFAULT NULL, molbinary VARBINARY DEFAULT NULL, options VARCHAR DEFAULT '--out smiles')
  RETURNS VARCHAR
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molstring` | VARCHAR | `NULL` | Molecule as a string: standard or ChemAxon-extended SMILES, or molblock/molfile (auto-detected). |
| `molbinary` | VARBINARY | `NULL` | Molecule in RDKit binary encoding (e.g. from `Molstring2Molbinary`). |
| `options` | VARCHAR | `'--out smiles'` | Command-line style option string (e.g. `'--out inchikey --remove_stereo'`). Pass `'-h'` to list all options. |

**Returns:** VARCHAR — Molecule in the format selected by `--out`.

**Description:**

```text
Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) or RDKit binary molecule encoding to a variety of different string-based formats, optionally transforming the input structure. Only one of molstring or molbinary arguments can be non-NULL. All available options can be listed by by running the following SQL: select Molstring_Or_Molbinary2Molstring(NULL, NULL, '-h'); usage info will be returned as part of the error message. If the options argument value is not specified, computes canonical ChemAxon-compatible extended SMILES.

Options:
   -h, --help:  Show this help message and exit
Output Encoding:
   --out OUT:  The output molecule encoding ((smiles|smi)|(molfile|mo lblock)|molhash|molhashfull|tautohash|inchi|inchikey)
Transform:
   --desalt:  Remove (strip) salt
   --desalt_smarts_list DESALT_SMARTS_LIST:  Additional desalt SMARTS patterns separated with |
   --remove_stereo:  Remove stereo
SMILES encoder options:
   --smiles_kekule: Kekulize the molecule before generating the SMILES and output single/double rather than aromatic bonds
```

**Example:**

```sql
SELECT Molstring_Or_Molbinary2Molstring('CC(=O)Oc1ccccc1C(=O)O')                                  AS cx_smiles,
       Molstring_Or_Molbinary2Molstring('CC(=O)Oc1ccccc1C(=O)O', NULL, '--out inchikey')            AS inchikey,
       Molstring_Or_Molbinary2Molstring('CC(=O)[O-].[Na+]', NULL, '--out smiles --desalt') AS desalted  -- options must include --out;
-- List all options (returned in the error message):
-- SELECT Molstring_Or_Molbinary2Molstring(NULL, NULL, '-h');
-- CC(=O)Oc1ccccc1C(=O)O | BSYNRYMUTXBXSQ-UHFFFAOYSA-N | CC(=O)[O-]
```

### 5. Molstring_Or_Molbinary2Molbinary

**Kind:** Function | **Language:** PYTHON

```sql
Molstring_Or_Molbinary2Molbinary(molstring VARCHAR DEFAULT NULL, molbinary VARBINARY DEFAULT NULL, options VARCHAR DEFAULT '')
  RETURNS VARBINARY
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molstring` | VARCHAR | `NULL` | Molecule as a string: standard or ChemAxon-extended SMILES, or molblock/molfile (auto-detected). |
| `molbinary` | VARBINARY | `NULL` | Molecule in RDKit binary encoding (e.g. from `Molstring2Molbinary`). |
| `options` | VARCHAR | `''` | Command-line style option string (e.g. `'--out inchikey --remove_stereo'`). Pass `'-h'` to list all options. |

**Returns:** VARBINARY — RDKit binary-encoded molecule.

**Description:**

```text
Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) or RDKit binary molecule encoding to the RDKit binary molecule encoding, optionally transforming the input structure. Only one of molstring or molbinary arguments can be non-NULL. All available options can be listed by by running the following SQL: select Molstring_Or_Molbinary2Molbinary(NULL, NULL, '-h'); usage info will be returned as part of the error message. To convert molstrings to RDKit binary molecule encoding w/o applying any transforms, use the Molstring2Molbinary function, which has fewer arguments and is faster.

Options: -h, --help show this help message and exit
Transform:
   --desalt:  Remove (strip) salt
   --desalt_smarts_list DESALT_SMARTS_LIST:  Additional desalt SMARTS patterns separated with |
   --remove_stereo:  Remove stereo
```

**Example:**

```sql
SELECT Molstring_Or_Molbinary2Molbinary('CC(=O)[O-].[Na+]', NULL, '--desalt') AS molbinary;
-- round-trip: Molstring_Or_Molbinary2Molstring(NULL, <result>) -> CC(=O)[O-]
```

### 6. Molstring2Molbinary

**Kind:** Function | **Language:** PYTHON

```sql
Molstring2Molbinary(molstring VARCHAR)
  RETURNS VARBINARY
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molstring` | VARCHAR | _required_ | Molecule as a string: standard or ChemAxon-extended SMILES, or molblock/molfile (auto-detected). |

**Returns:** VARBINARY — RDKit binary-encoded molecule; NULL for NULL/empty/invalid input.

**Description:**

Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) to RDKit binary molecule encoding. Returns NULL if molstring is NULL, empty, invalid, or represents an empty molecule with 0 atoms and 0 bonds.

**Example:**

```sql
SELECT Molstring2Molbinary('CC(=O)Oc1ccccc1C(=O)O') AS molbinary;
```

## Fingerprinting

### 7. Molstring2Pattern_Fingerprint

**Kind:** Function | **Language:** PYTHON

```sql
Molstring2Pattern_Fingerprint(molstring VARCHAR)
  RETURNS VARBINARY
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molstring` | VARCHAR | _required_ | Molecule as a string: standard or ChemAxon-extended SMILES, or molblock/molfile (auto-detected). |

**Returns:** VARBINARY — Pattern fingerprint bitset.

**Description:**

Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) to RDKit substructure pattern fingerprint commonly used for fingerprint-based screening to speed up substructure searches. Returns NULL if molstring is NULL, empty, invalid, or represents an empty molecule with 0 atoms and 0 bonds.

**Example:**

```sql
SELECT Molstring2Pattern_Fingerprint('CC(=O)Oc1ccccc1C(=O)O') AS pattern_fp;
```

### 8. Molbinary2Pattern_Fingerprint

**Kind:** Function | **Language:** PYTHON

```sql
Molbinary2Pattern_Fingerprint(molbinary VARBINARY)
  RETURNS VARBINARY
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molbinary` | VARBINARY | _required_ | Molecule in RDKit binary encoding (e.g. from `Molstring2Molbinary`). |

**Returns:** VARBINARY — Pattern fingerprint bitset.

**Description:**

Converts RDKit binary-encoded molecule to RDKit substructure pattern fingerprint commonly used for fingerprint-based screening to speed up substructure searches. Returns NULL if molbinary is NULL, empty or represents an empty molecule with 0 atoms and 0 bonds. Returns error if molbinary is not a valid RDKit binary-encoded molecule.

**Example:**

```sql
SELECT Molbinary2Pattern_Fingerprint(Molstring2Molbinary('CC(=O)Oc1ccccc1C(=O)O')) AS pattern_fp;
```

### 9. Smarts2Pattern_Fingerprint

**Kind:** Function | **Language:** PYTHON

```sql
Smarts2Pattern_Fingerprint(smarts VARCHAR)
  RETURNS VARBINARY
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `smarts` | VARCHAR | _required_ | SMARTS substructure pattern. |

**Returns:** VARBINARY — Pattern fingerprint bitset.

**Description:**

Converts SMARTS substructure pattern to RDKit substructure pattern fingerprint commonly used for fingerprint-based screening to speed up substructure searches. Returns NULL if smarts is NULL, empty, or represents an empty molecule with 0 atoms and 0 bonds. Returns error if smarts is invalid and cannot be parsed.

**Example:**

```sql
SELECT Smarts2Pattern_Fingerprint('C(=O)O') AS pattern_fp;
-- 256-byte (2048-bit) fingerprint
```

### 10. Molstring2Morgan_Fingerprint

**Kind:** Function | **Language:** PYTHON

```sql
Molstring2Morgan_Fingerprint(molstring VARCHAR)
  RETURNS VARBINARY
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molstring` | VARCHAR | _required_ | Molecule as a string: standard or ChemAxon-extended SMILES, or molblock/molfile (auto-detected). |

**Returns:** VARBINARY — 2048-bit Morgan fingerprint (256 bytes).

**Description:**

Converts molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) to RDKit Morgan fingerprint commonly used for fingerprint-based similarity search. Returns NULL if molstring is NULL, empty, invalid, or represents an empty molecule with 0 atoms and 0 bonds. Generator options: radius=2, fpSize=2048, atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen()

**Example:**

```sql
SELECT Molstring2Morgan_Fingerprint('CC(=O)Oc1ccccc1C(=O)O') AS morgan_fp;
```

### 11. Molbinary2Morgan_Fingerprint

**Kind:** Function | **Language:** PYTHON

```sql
Molbinary2Morgan_Fingerprint(molbinary VARBINARY)
  RETURNS VARBINARY
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molbinary` | VARBINARY | _required_ | Molecule in RDKit binary encoding (e.g. from `Molstring2Molbinary`). |

**Returns:** VARBINARY — 2048-bit Morgan fingerprint (256 bytes).

**Description:**

Converts RDKit binary-encoded molecule to RDKit Morgan fingerprint commonly used for fingerprint-based similarity search. Returns NULL if molbinary is NULL or empty, or represents an empty molecule with 0 atoms and 0 bonds. Generator options: radius=2, fpSize=2048, atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen(). Returns error if molbinary is not a valid RDKit binary-encoded molecule.

**Example:**

```sql
SELECT Molbinary2Morgan_Fingerprint(Molstring2Molbinary('CC(=O)Oc1ccccc1C(=O)O')) AS morgan_fp;
```

## Substructure matching

### 12. Molbinary_Matches_Smarts

**Kind:** Function | **Language:** PYTHON

```sql
Molbinary_Matches_Smarts(molbinary VARBINARY, smarts VARCHAR, screen_pass BOOLEAN)
  RETURNS BOOLEAN
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molbinary` | VARBINARY | _required_ | Molecule in RDKit binary encoding (e.g. from `Molstring2Molbinary`). |
| `smarts` | VARCHAR | _required_ | SMARTS substructure pattern. |
| `screen_pass` | BOOLEAN | _required_ | Result of fingerprint pre-screening; pass `IFF(BITAND(fp, query_fp) = query_fp, TRUE, NULL)` to skip the expensive match for rows that cannot match, or TRUE to always match. |

**Returns:** BOOLEAN — TRUE iff the molecule matches and `screen_pass` is TRUE.

**Description:**

Tests whether RDKit binary-encoded molecule matches the specified SMARTS pattern and returns TRUE iff it does and the screen_pass argument value is TRUE. The extra screen_pass argument is used for substructure query optimization based on fingerprint screening (see examples in the documentation and example workbooks). Returns NULL if any of the args are NULL. Returns error if molbinary is not a valid RDKit binary-encoded molecule of if smarts is invalid and cannot be parsed. Note: an empty SMARTS will not match any molecule, even an empty one. This seems to be illogical, since, in theory, a subgraph with 0 nodes and 0 edges must match any graph (or, at least, an empty one), but, in practice, this approach leads to fewer problems than the theoretically correct one. In RDKit itself, a molecule representing an empty pattern does not match anything either.

**Example:**

```sql
-- Fingerprint-screened substructure search: the expensive match runs only for rows passing the screen
WITH mols AS (
    SELECT Molstring2Molbinary(column1) AS mb FROM VALUES ('CC(=O)Oc1ccccc1C(=O)O'), ('c1ccccc1'), ('CCO')
), q AS (SELECT Smarts2Pattern_Fingerprint('C(=O)O') AS qfp)
SELECT Molstring_Or_Molbinary2Molstring(NULL, mb) AS smiles
FROM mols, q
WHERE Molbinary_Matches_Smarts(mb, 'C(=O)O',
        IFF(BITAND(Molbinary2Pattern_Fingerprint(mb), qfp) = qfp, TRUE, NULL));
-- CC(=O)Oc1ccccc1C(=O)O
```

### 13. Molstring_Matches_Smarts

**Kind:** Function | **Language:** PYTHON

```sql
Molstring_Matches_Smarts(molstring VARCHAR, smarts VARCHAR, screen_pass BOOLEAN)
  RETURNS BOOLEAN
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molstring` | VARCHAR | _required_ | Molecule as a string: standard or ChemAxon-extended SMILES, or molblock/molfile (auto-detected). |
| `smarts` | VARCHAR | _required_ | SMARTS substructure pattern. |
| `screen_pass` | BOOLEAN | _required_ | Result of fingerprint pre-screening; pass `IFF(BITAND(fp, query_fp) = query_fp, TRUE, NULL)` to skip the expensive match for rows that cannot match, or TRUE to always match. |

**Returns:** BOOLEAN — TRUE iff the molecule matches and `screen_pass` is TRUE.

**Description:**

Tests whether the molecule encoded as molstring (standard or ChemAxon-extended SMILES or molblock/molfile, auto-detected) matches the specified SMARTS pattern and returns TRUE iff it does and the screen_pass argument value is TRUE. The extra screen_pass argument is used for substructure query optimization based on fingerprint screening (see examples in the documentation and example workbooks). Returns NULL if any of the args are NULL. Returns FALSE if molstring is not a valid SMILES or molblock. Returns error if smarts is invalid and cannot be parsed. Note: an empty SMARTS will not match any molecule, even an empty one. This seems to be illogical, since, in theory, a subgraph with 0 nodes and 0 edges must match any graph (or, at least, an empty one), but, in practice, this approach leads to fewer problems than the theoretically correct one. In RDKit itself, a molecule representing an empty pattern does not match anything either.

**Example:**

```sql
SELECT Molstring_Matches_Smarts('CC(=O)Oc1ccccc1C(=O)O', 'C(=O)O', TRUE) AS has_carbonyl_o;
-- TRUE
```

## Structure similarity

### 14. Tanimoto_P

**Kind:** Function | **Language:** PYTHON

```sql
Tanimoto_P(v1 VARBINARY, v2 VARBINARY)
  RETURNS FLOAT
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `v1` | VARBINARY | _required_ | First bitset (e.g. Morgan fingerprint). |
| `v2` | VARBINARY | _required_ | Second bitset; must be the same length as `v1`. |

**Returns:** FLOAT — Similarity in [0, 1]; 1.0 for two all-zero bitsets.

**Description:**

A version of TANIMOTO implemented in Python. Computes Tanimoto similarity between two bitsets represented by binary vectors. Returns BITCOUNT(v1 bitand v2) / BITCOUNT(v1 bitor v2) as float. If both v1 and v2 have no bits set to 1, returns 1.0, that is, treats two empty bitsets or two bitsets filled with only zeroes as being equal to each other. Returns error if two bitsets are of different lengths. Returns NULL if one of both arguments are NULLs

**Example:**

```sql
SELECT Tanimoto_P(Molstring2Morgan_Fingerprint('CC(=O)Oc1ccccc1C(=O)O'), Molstring2Morgan_Fingerprint('OC(=O)c1ccccc1O')) AS similarity;
-- 0.4285714286
```

### 15. Tanimoto_J

**Kind:** Function | **Language:** JAVA

```sql
Tanimoto_J(v1 VARBINARY, v2 VARBINARY)
  RETURNS FLOAT
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `v1` | VARBINARY | _required_ | First bitset (e.g. Morgan fingerprint). |
| `v2` | VARBINARY | _required_ | Second bitset; must be the same length as `v1`. |

**Returns:** FLOAT — Similarity in [0, 1]; 1.0 for two all-zero bitsets.

**Description:**

A version of TANIMOTO implemented in Java. Computes Tanimoto similarity between two bitsets represented by binary vectors. Returns BITCOUNT(v1 bitand v2) / BITCOUNT(v1 bitor v2) as float. If both v1 and v2 have no bits set to 1, returns 1.0, that is, treats two empty bitsets or two bitsets filled with only zeroes as being equal to each other. Returns error if two bitsets are of different lengths. Returns NULL if one of both arguments are NULLs

**Example:**

```sql
SELECT Tanimoto_J(Molstring2Morgan_Fingerprint('CC(=O)Oc1ccccc1C(=O)O'), Molstring2Morgan_Fingerprint('OC(=O)c1ccccc1O')) AS similarity;
-- 0.4285714286
```

### 16. Tanimoto

**Kind:** Function | **Language:** SQL

```sql
Tanimoto(v1 VARBINARY, v2 VARBINARY)
  RETURNS FLOAT
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `v1` | VARBINARY | _required_ | First bitset (e.g. Morgan fingerprint). |
| `v2` | VARBINARY | _required_ | Second bitset; must be the same length as `v1`. |

**Returns:** FLOAT — Similarity in [0, 1]; 0.0 for two all-zero bitsets.

**Description:**

Computes Tanimoto similarity between two bitsets represented by binary vectors. Returns BITCOUNT(v1 bitand v2) / BITCOUNT(v1 bitor v2) as float. If both v1 and v2 have no bits set to 1, returns 0.0, that is, treats two 0-length bitsets or two bitsets filled with only zeroes as being different from each other. The lengths of v1 and v2 must be the same, otheerwise, an exception will be raised. Returns NULL if one of both arguments are NULLs

**Example:**

```sql
SELECT Tanimoto(Molstring2Morgan_Fingerprint('CC(=O)Oc1ccccc1C(=O)O'), Molstring2Morgan_Fingerprint('OC(=O)c1ccccc1O')) AS similarity;
-- 0.428571
```

## Molecular properties and descriptor calculation

### 17. Molstring_Or_Molbinary2Medchem_Descriptors

**Kind:** Function | **Language:** PYTHON

```sql
Molstring_Or_Molbinary2Medchem_Descriptors(molstring VARCHAR DEFAULT NULL, molbinary VARBINARY DEFAULT NULL)
  RETURNS TABLE (Q_Estim_Drug_Likeness float, Mol_Wt float, Heavy_Atom_Mol_Wt float, Exact_Mol_Wt float, Num_Valence_Electrons int, Num_Radical_Electrons int, Max_Partial_Charge float, Min_Partial_Charge float, Max_Abs_Partial_Charge float, Min_Abs_Partial_Charge float, TPSA float, Fraction_CSP3 float, Heavy_Atom_Count int, Num_Aliphatic_Carbocycles int, Num_Aliphatic_Heterocycles int, Num_Aliphatic_Rings int, Num_Aromatic_Carbocycles int, Num_Aromatic_Heterocycles int, Num_Aromatic_Rings int, Num_HAcceptors int, Num_HDonors int, Num_Heteroatoms int, Num_Rotatable_Bonds int, Num_Saturated_Carbocycles int, Num_Saturated_Heterocycles int, Num_Saturated_Rings int, Ring_Count int, Mol_Log_P float, Mol_MR float)
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molstring` | VARCHAR | `NULL` | Molecule as a string: standard or ChemAxon-extended SMILES, or molblock/molfile (auto-detected). |
| `molbinary` | VARBINARY | `NULL` | Molecule in RDKit binary encoding (e.g. from `Molstring2Molbinary`). |

**Returns:** TABLE — One row of descriptors (columns below).

| Column | Type |
| --- | --- |
| `Q_Estim_Drug_Likeness` | FLOAT |
| `Mol_Wt` | FLOAT |
| `Heavy_Atom_Mol_Wt` | FLOAT |
| `Exact_Mol_Wt` | FLOAT |
| `Num_Valence_Electrons` | INT |
| `Num_Radical_Electrons` | INT |
| `Max_Partial_Charge` | FLOAT |
| `Min_Partial_Charge` | FLOAT |
| `Max_Abs_Partial_Charge` | FLOAT |
| `Min_Abs_Partial_Charge` | FLOAT |
| `TPSA` | FLOAT |
| `Fraction_CSP3` | FLOAT |
| `Heavy_Atom_Count` | INT |
| `Num_Aliphatic_Carbocycles` | INT |
| `Num_Aliphatic_Heterocycles` | INT |
| `Num_Aliphatic_Rings` | INT |
| `Num_Aromatic_Carbocycles` | INT |
| `Num_Aromatic_Heterocycles` | INT |
| `Num_Aromatic_Rings` | INT |
| `Num_HAcceptors` | INT |
| `Num_HDonors` | INT |
| `Num_Heteroatoms` | INT |
| `Num_Rotatable_Bonds` | INT |
| `Num_Saturated_Carbocycles` | INT |
| `Num_Saturated_Heterocycles` | INT |
| `Num_Saturated_Rings` | INT |
| `Ring_Count` | INT |
| `Mol_Log_P` | FLOAT |
| `Mol_MR` | FLOAT |

**Description:**

Computes commonly used medchem properties (descriptors) for a molecule encoded as molstring (standard or ChemAxon SMILES or molfile/molblock) or as RDKit binary-encoded molecule. Only one of molstring or molbinary arguments can be non-NULL, otherwise, an error will be returned. Returns a table with one row and multiple columns corresponding to the computed descriptors. If molstring and molbinary are both NULLs, or molstring is invalid and cannot be parsed, returns a table with one row filled with NULLs. Returns an error if molbinary is not a valid RDKit binary-encoded molecule.

**Example:**

```sql
SELECT Mol_Wt, TPSA, Mol_Log_P, Num_HDonors, Num_HAcceptors
FROM TABLE(Molstring_Or_Molbinary2Medchem_Descriptors(molstring => 'CC(=O)Oc1ccccc1C(=O)O'));
-- 180.159 | 63.6 | 1.3101 | 1 | 3
```

## Structure encoding and molecular graph validity check

### 18. Molstring_Or_Molbinary_Check

**Kind:** Function | **Language:** PYTHON

```sql
Molstring_Or_Molbinary_Check(molstring VARCHAR DEFAULT NULL, molbinary VARBINARY DEFAULT NULL, raise_exception BOOLEAN DEFAULT FALSE)
  RETURNS TABLE (is_ok boolean, encoding varchar, error_msg varchar)
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molstring` | VARCHAR | `NULL` | Molecule as a string: standard or ChemAxon-extended SMILES, or molblock/molfile (auto-detected). |
| `molbinary` | VARBINARY | `NULL` | Molecule in RDKit binary encoding (e.g. from `Molstring2Molbinary`). |
| `raise_exception` | BOOLEAN | `FALSE` | If TRUE, raises an error instead of returning when the input is not a valid molecule. |

**Returns:** TABLE — One row with validity flag, detected encoding and error message.

| Column | Type |
| --- | --- |
| `is_ok` | BOOLEAN |
| `encoding` | VARCHAR |
| `error_msg` | VARCHAR |

**Description:**

Checks a molecule encoded as molstring (standard or ChemAxon SMILES or molfile/molblock) or as RDKit binary-encoded molecule. Only one of molstring or molbinary arguments can be non-NULL, otherwise, an error will be returned. Returns a table with one row and three columns: is_ok boolean, encoding varchar, and error_msg varchar. If both molstring and molbinary are NULL, the entire result row will be filled with NULLs. Otherwise, is_ok will contain True iff the input can be parsed into a molecule with no errors, encoding will contain a string representation of the encoding (MOLBLOCK, SMILES, or BINARY), and the error_msg will contain a description of the error or NULL. If raise_exception parameter is TRUE (it is FALSE by default) and the input cannot be parsed into a valid molecule, the method will raise an exception and quit instead of returning.

**Example:**

```sql
SELECT * FROM TABLE(Molstring_Or_Molbinary_Check(molstring => 'C1CC'));  -- unclosed ring
-- FALSE | SMILES | Parser error
```

## HELM support (experimental)

### 19. helm2molbinary

**Kind:** Procedure | **Language:** PYTHON

```sql
helm2molbinary(monomer_table VARCHAR, helm_table VARCHAR)
  RETURNS table(ID INT, MOLBINARY VARBINARY)
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `monomer_table` | VARCHAR | _required_ | Table of monomers with a `MOLFILE` column (helmkit-format SDF records); pass as `TABLE(name)`. |
| `helm_table` | VARCHAR | _required_ | Table with `ID INT` and `HELM VARCHAR` columns; pass as `TABLE(name)`. |

**Returns:** TABLE — One row per input HELM record.

| Column | Type |
| --- | --- |
| `ID` | INT |
| `MOLBINARY` | VARBINARY |

**Description:**

```text
Converts HELM strings in the specified helm_table to full molecular representation and returns structures
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
```

**Example:**

```sql
-- monomer: table with a MOLFILE column (helmkit-format SDF records); peptide: table with ID and HELM columns
CALL helm2molbinary(TABLE(monomer), TABLE(peptide)) ->>
SELECT p.id, p.helm
FROM $1 AS h
JOIN peptide p ON p.id = h.id
WHERE Molbinary_Matches_Smarts(h.molbinary, 'c1cc2ccccc2n1', TRUE);
```

## Single-input convenience wrappers for Molstring_Or_Molbinary2* UDFs

### 20. Conv_Mol2Molstring

**Kind:** Function | **Language:** SQL

```sql
Conv_Mol2Molstring(molstring VARCHAR, options VARCHAR DEFAULT '--out smiles')
  RETURNS VARCHAR
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molstring` | VARCHAR | _required_ | Molecule as a string: standard or ChemAxon-extended SMILES, or molblock/molfile (auto-detected). |
| `options` | VARCHAR | `'--out smiles'` | Command-line style option string (e.g. `'--out inchikey --remove_stereo'`). Pass `'-h'` to list all options. |

**Returns:** VARCHAR — Molecule in the format selected by `--out`.

**Description:**

Wrapper for [Molstring_Or_Molbinary2Molstring](#4-molstring_or_molbinary2molstring)(molstring, NULL, options).

**Example:**

```sql
SELECT Conv_Mol2Molstring('CC(=O)Oc1ccccc1C(=O)O', '--out inchi') AS inchi;
-- InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)
```

### 21. Conv_Mol2Molstring

**Kind:** Function | **Language:** SQL

```sql
Conv_Mol2Molstring(molbinary VARBINARY, options VARCHAR DEFAULT '--out smiles')
  RETURNS VARCHAR
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molbinary` | VARBINARY | _required_ | Molecule in RDKit binary encoding (e.g. from `Molstring2Molbinary`). |
| `options` | VARCHAR | `'--out smiles'` | Command-line style option string (e.g. `'--out inchikey --remove_stereo'`). Pass `'-h'` to list all options. |

**Returns:** VARCHAR — Molecule in the format selected by `--out`.

**Description:**

Wrapper for [Molstring_Or_Molbinary2Molstring](#4-molstring_or_molbinary2molstring)(NULL, molbinary, options).

**Example:**

```sql
SELECT Conv_Mol2Molstring(Molstring2Molbinary('CC(=O)Oc1ccccc1C(=O)O')) AS cx_smiles;
-- CC(=O)Oc1ccccc1C(=O)O
```

### 22. Conv_Mol2Molbinary

**Kind:** Function | **Language:** SQL

```sql
Conv_Mol2Molbinary(molstring VARCHAR, options VARCHAR DEFAULT '')
  RETURNS VARBINARY
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molstring` | VARCHAR | _required_ | Molecule as a string: standard or ChemAxon-extended SMILES, or molblock/molfile (auto-detected). |
| `options` | VARCHAR | `''` | Command-line style option string (e.g. `'--out inchikey --remove_stereo'`). Pass `'-h'` to list all options. |

**Returns:** VARBINARY — RDKit binary-encoded molecule.

**Description:**

Wrapper for [Molstring_Or_Molbinary2Molbinary](#5-molstring_or_molbinary2molbinary)(molstring, NULL, options).

**Example:**

```sql
SELECT Conv_Mol2Molbinary('CC(=O)[O-].[Na+]', '--desalt') AS molbinary;
```

### 23. Conv_Mol2Molbinary

**Kind:** Function | **Language:** SQL

```sql
Conv_Mol2Molbinary(molbinary VARBINARY, options VARCHAR DEFAULT '')
  RETURNS VARBINARY
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molbinary` | VARBINARY | _required_ | Molecule in RDKit binary encoding (e.g. from `Molstring2Molbinary`). |
| `options` | VARCHAR | `''` | Command-line style option string (e.g. `'--out inchikey --remove_stereo'`). Pass `'-h'` to list all options. |

**Returns:** VARBINARY — RDKit binary-encoded molecule.

**Description:**

Wrapper for [Molstring_Or_Molbinary2Molbinary](#5-molstring_or_molbinary2molbinary)(NULL, molbinary, options).

**Example:**

```sql
SELECT Conv_Mol2Molbinary(Molstring2Molbinary('CC(=O)Oc1ccccc1C(=O)O'), '--remove_stereo') AS molbinary;
```

### 24. Calc_Mol_Medchem_Descriptors

**Kind:** Function | **Language:** SQL

```sql
Calc_Mol_Medchem_Descriptors(molstring VARCHAR)
  RETURNS TABLE (Q_Estim_Drug_Likeness float, Mol_Wt float, Heavy_Atom_Mol_Wt float, Exact_Mol_Wt float, Num_Valence_Electrons int, Num_Radical_Electrons int, Max_Partial_Charge float, Min_Partial_Charge float, Max_Abs_Partial_Charge float, Min_Abs_Partial_Charge float, TPSA float, Fraction_CSP3 float, Heavy_Atom_Count int, Num_Aliphatic_Carbocycles int, Num_Aliphatic_Heterocycles int, Num_Aliphatic_Rings int, Num_Aromatic_Carbocycles int, Num_Aromatic_Heterocycles int, Num_Aromatic_Rings int, Num_HAcceptors int, Num_HDonors int, Num_Heteroatoms int, Num_Rotatable_Bonds int, Num_Saturated_Carbocycles int, Num_Saturated_Heterocycles int, Num_Saturated_Rings int, Ring_Count int, Mol_Log_P float, Mol_MR float)
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molstring` | VARCHAR | _required_ | Molecule as a string: standard or ChemAxon-extended SMILES, or molblock/molfile (auto-detected). |

**Returns:** TABLE — One row of descriptors (columns below).

| Column | Type |
| --- | --- |
| `Q_Estim_Drug_Likeness` | FLOAT |
| `Mol_Wt` | FLOAT |
| `Heavy_Atom_Mol_Wt` | FLOAT |
| `Exact_Mol_Wt` | FLOAT |
| `Num_Valence_Electrons` | INT |
| `Num_Radical_Electrons` | INT |
| `Max_Partial_Charge` | FLOAT |
| `Min_Partial_Charge` | FLOAT |
| `Max_Abs_Partial_Charge` | FLOAT |
| `Min_Abs_Partial_Charge` | FLOAT |
| `TPSA` | FLOAT |
| `Fraction_CSP3` | FLOAT |
| `Heavy_Atom_Count` | INT |
| `Num_Aliphatic_Carbocycles` | INT |
| `Num_Aliphatic_Heterocycles` | INT |
| `Num_Aliphatic_Rings` | INT |
| `Num_Aromatic_Carbocycles` | INT |
| `Num_Aromatic_Heterocycles` | INT |
| `Num_Aromatic_Rings` | INT |
| `Num_HAcceptors` | INT |
| `Num_HDonors` | INT |
| `Num_Heteroatoms` | INT |
| `Num_Rotatable_Bonds` | INT |
| `Num_Saturated_Carbocycles` | INT |
| `Num_Saturated_Heterocycles` | INT |
| `Num_Saturated_Rings` | INT |
| `Ring_Count` | INT |
| `Mol_Log_P` | FLOAT |
| `Mol_MR` | FLOAT |

**Description:**

Wrapper for [Molstring_Or_Molbinary2Medchem_Descriptors](#17-molstring_or_molbinary2medchem_descriptors)(molstring, NULL).

**Example:**

```sql
SELECT Mol_Wt, TPSA FROM TABLE(Calc_Mol_Medchem_Descriptors('CC(=O)Oc1ccccc1C(=O)O'));
-- 180.159 | 63.6
```

### 25. Calc_Mol_Medchem_Descriptors

**Kind:** Function | **Language:** SQL

```sql
Calc_Mol_Medchem_Descriptors(molbinary VARBINARY)
  RETURNS TABLE (Q_Estim_Drug_Likeness float, Mol_Wt float, Heavy_Atom_Mol_Wt float, Exact_Mol_Wt float, Num_Valence_Electrons int, Num_Radical_Electrons int, Max_Partial_Charge float, Min_Partial_Charge float, Max_Abs_Partial_Charge float, Min_Abs_Partial_Charge float, TPSA float, Fraction_CSP3 float, Heavy_Atom_Count int, Num_Aliphatic_Carbocycles int, Num_Aliphatic_Heterocycles int, Num_Aliphatic_Rings int, Num_Aromatic_Carbocycles int, Num_Aromatic_Heterocycles int, Num_Aromatic_Rings int, Num_HAcceptors int, Num_HDonors int, Num_Heteroatoms int, Num_Rotatable_Bonds int, Num_Saturated_Carbocycles int, Num_Saturated_Heterocycles int, Num_Saturated_Rings int, Ring_Count int, Mol_Log_P float, Mol_MR float)
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molbinary` | VARBINARY | _required_ | Molecule in RDKit binary encoding (e.g. from `Molstring2Molbinary`). |

**Returns:** TABLE — One row of descriptors (columns below).

| Column | Type |
| --- | --- |
| `Q_Estim_Drug_Likeness` | FLOAT |
| `Mol_Wt` | FLOAT |
| `Heavy_Atom_Mol_Wt` | FLOAT |
| `Exact_Mol_Wt` | FLOAT |
| `Num_Valence_Electrons` | INT |
| `Num_Radical_Electrons` | INT |
| `Max_Partial_Charge` | FLOAT |
| `Min_Partial_Charge` | FLOAT |
| `Max_Abs_Partial_Charge` | FLOAT |
| `Min_Abs_Partial_Charge` | FLOAT |
| `TPSA` | FLOAT |
| `Fraction_CSP3` | FLOAT |
| `Heavy_Atom_Count` | INT |
| `Num_Aliphatic_Carbocycles` | INT |
| `Num_Aliphatic_Heterocycles` | INT |
| `Num_Aliphatic_Rings` | INT |
| `Num_Aromatic_Carbocycles` | INT |
| `Num_Aromatic_Heterocycles` | INT |
| `Num_Aromatic_Rings` | INT |
| `Num_HAcceptors` | INT |
| `Num_HDonors` | INT |
| `Num_Heteroatoms` | INT |
| `Num_Rotatable_Bonds` | INT |
| `Num_Saturated_Carbocycles` | INT |
| `Num_Saturated_Heterocycles` | INT |
| `Num_Saturated_Rings` | INT |
| `Ring_Count` | INT |
| `Mol_Log_P` | FLOAT |
| `Mol_MR` | FLOAT |

**Description:**

Wrapper for [Molstring_Or_Molbinary2Medchem_Descriptors](#17-molstring_or_molbinary2medchem_descriptors)(NULL, molbinary).

**Example:**

```sql
SELECT Mol_Wt, TPSA FROM TABLE(Calc_Mol_Medchem_Descriptors(Molstring2Molbinary('CC(=O)Oc1ccccc1C(=O)O')));
-- 180.159 | 63.6
```

## Single-input convenience wrappers for Draw_Molstring_Or_Molbinary2* UDFs

### 26. Draw_Mol2SVG

**Kind:** Function | **Language:** SQL

```sql
Draw_Mol2SVG(molstring VARCHAR, highlight_smarts VARCHAR DEFAULT NULL, draw_options VARCHAR DEFAULT NULL, output_raw_xml BOOLEAN DEFAULT FALSE)
  RETURNS VARCHAR
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molstring` | VARCHAR | _required_ | Molecule as a string: standard or ChemAxon-extended SMILES, or molblock/molfile (auto-detected). |
| `highlight_smarts` | VARCHAR | `NULL` | Optional SMARTS pattern; atoms matching it are highlighted. |
| `draw_options` | VARCHAR | `NULL` | Reserved for future use; currently ignored. |
| `output_raw_xml` | BOOLEAN | `FALSE` | If TRUE, returns plain SVG XML; if FALSE, a `data:image/svg+xml;base64,` data URI. |

**Returns:** VARCHAR — SVG image (data URI or raw XML).

**Description:**

Wrapper for [Draw_Molstring_Or_Molbinary2SVG](#2-draw_molstring_or_molbinary2svg)(molstring, NULL, highlight_smarts, draw_options, output_raw_xml).

**Example:**

```sql
SELECT Draw_Mol2SVG('CC(=O)Oc1ccccc1C(=O)O', 'C(=O)O');
```

### 27. Draw_Mol2SVG

**Kind:** Function | **Language:** SQL

```sql
Draw_Mol2SVG(molbinary VARBINARY, highlight_smarts VARCHAR DEFAULT NULL, draw_options VARCHAR DEFAULT NULL, output_raw_xml BOOLEAN DEFAULT FALSE)
  RETURNS VARCHAR
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `molbinary` | VARBINARY | _required_ | Molecule in RDKit binary encoding (e.g. from `Molstring2Molbinary`). |
| `highlight_smarts` | VARCHAR | `NULL` | Optional SMARTS pattern; atoms matching it are highlighted. |
| `draw_options` | VARCHAR | `NULL` | Reserved for future use; currently ignored. |
| `output_raw_xml` | BOOLEAN | `FALSE` | If TRUE, returns plain SVG XML; if FALSE, a `data:image/svg+xml;base64,` data URI. |

**Returns:** VARCHAR — SVG image (data URI or raw XML).

**Description:**

Wrapper for [Draw_Molstring_Or_Molbinary2SVG](#2-draw_molstring_or_molbinary2svg)(NULL, molbinary, highlight_smarts, draw_options, output_raw_xml).

**Example:**

```sql
SELECT Draw_Mol2SVG(Molstring2Molbinary('CC(=O)Oc1ccccc1C(=O)O'), NULL, NULL, TRUE  -- positional: named and positional args can't be mixed);
-- <?xml version='1.0' encoding=... (raw SVG XML)
```

### 28. Draw_Mol2PNG

**Kind:** Function | **Language:** SQL

```sql
Draw_Mol2PNG(width_px INT, height_px INT, molstring VARCHAR, highlight_smarts VARCHAR DEFAULT NULL, draw_options VARCHAR DEFAULT NULL)
  RETURNS VARBINARY
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `width_px` | INT | _required_ | Image width in pixels. |
| `height_px` | INT | _required_ | Image height in pixels. |
| `molstring` | VARCHAR | _required_ | Molecule as a string: standard or ChemAxon-extended SMILES, or molblock/molfile (auto-detected). |
| `highlight_smarts` | VARCHAR | `NULL` | Optional SMARTS pattern; atoms matching it are highlighted. |
| `draw_options` | VARCHAR | `NULL` | Reserved for future use; currently ignored. |

**Returns:** VARBINARY — PNG image bytes.

**Description:**

Wrapper for [Draw_Molstring_Or_Molbinary2PNG](#3-draw_molstring_or_molbinary2png)(width_px, height_px, molstring, NULL, highlight_smarts, draw_options).

**Example:**

```sql
SELECT Draw_Mol2PNG(300, 200, 'CC(=O)Oc1ccccc1C(=O)O');
```

### 29. Draw_Mol2PNG

**Kind:** Function | **Language:** SQL

```sql
Draw_Mol2PNG(width_px INT, height_px INT, molbinary VARBINARY, highlight_smarts VARCHAR DEFAULT NULL, draw_options VARCHAR DEFAULT NULL)
  RETURNS VARBINARY
```

**Parameters:**

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `width_px` | INT | _required_ | Image width in pixels. |
| `height_px` | INT | _required_ | Image height in pixels. |
| `molbinary` | VARBINARY | _required_ | Molecule in RDKit binary encoding (e.g. from `Molstring2Molbinary`). |
| `highlight_smarts` | VARCHAR | `NULL` | Optional SMARTS pattern; atoms matching it are highlighted. |
| `draw_options` | VARCHAR | `NULL` | Reserved for future use; currently ignored. |

**Returns:** VARBINARY — PNG image bytes.

**Description:**

Wrapper for [Draw_Molstring_Or_Molbinary2PNG](#3-draw_molstring_or_molbinary2png)(width_px, height_px, NULL, molbinary, highlight_smarts, draw_options).

**Example:**

```sql
SELECT Draw_Mol2PNG(300, 200, Molstring2Molbinary('CC(=O)Oc1ccccc1C(=O)O'), 'c1ccccc1');
```

