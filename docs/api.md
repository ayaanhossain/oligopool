# Oligopool Calculator - API Reference

Complete parameter reference for all oligopool modules.

**Related documentation**:
- [User Guide](docs.md) - Tutorials, examples, and workflows
- [AI Agent Guide](agent-skills.md) - Decision trees, best practices, and gotchas
- [Docker Guide](docker-notes.md) - Run oligopool in a container

---

## Table of Contents

**Design Mode**
- [barcode](#barcode) - Hamming-distance barcodes
- [primer](#primer) - Thermodynamic primers
- [motif](#motif) - Sequence motifs & anchors
- [spacer](#spacer) - Neutral spacers
- [background](#background) - K-mer screening database
- [split](#split) - Fragment long oligos
- [pad](#pad) - Assembly-ready padding
- [merge](#merge) - Collapse columns
- [revcomp](#revcomp) - Reverse complement
- [lenstat](#lenstat) - Length statistics
- [verify](#verify) - QC check
- [final](#final) - Finalize for synthesis

**Analysis Mode**
- [index](#index) - Build barcode index
- [pack](#pack) - Preprocess FastQ
- [acount](#acount) - Association counting
- [xcount](#xcount) - Combinatorial counting

**Advanced**
- [vectorDB](#vectordb) - K-mer database
- [Scry](#scry) - Barcode classifier

**Appendix**
- [IUPAC Codes](#iupac-codes)
- [Common Restriction Sites](#common-restriction-sites)
- [Type IIS Enzymes](#type-iis-enzymes)
- [File Formats](#file-formats)
- [CLI Parameter Mapping](#cli-parameter-mapping)

---

## Design Mode

### barcode

**Purpose**: Generate Hamming-distance separated barcodes for unique variant identification.

**Signature**:
```python
df, stats = op.barcode(
    # Required
    input_data,                    # str | pd.DataFrame
    oligo_length_limit,            # int
    barcode_length,                # int
    minimum_hamming_distance,      # int
    maximum_repeat_length,         # int
    barcode_column,                # str

    # Optional
    output_file=None,              # str | None
    barcode_type=0,                # int | str
    left_context_column=None,      # str | None
    right_context_column=None,     # str | None
    patch_mode=False,              # bool
    cross_barcode_columns=None,    # str | list[str] | None
    minimum_cross_distance=None,   # int | None
    excluded_motifs=None,          # list | str | pd.DataFrame | None
    verbose=True,                  # bool
    random_seed=None,              # int | None
)
```

**Required Parameters**

| Parameter | Type | Constraints | Description |
|-----------|------|-------------|-------------|
| `input_data` | str \| pd.DataFrame | - | CSV path or DataFrame with `ID` + DNA sequence columns |
| `oligo_length_limit` | int | ≥4 | Maximum allowed oligo length (bp) |
| `barcode_length` | int | ≥4 | Length of designed barcodes (bp) |
| `minimum_hamming_distance` | int | ≥1 | Minimum pairwise Hamming distance within newly designed set |
| `maximum_repeat_length` | int | ≥4 | Maximum shared repeat length between barcode and context/oligos |
| `barcode_column` | str | - | Column name to create/overwrite (or patch-fill when `patch_mode=True`) |

**Optional Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `output_file` | str \| None | None | Output CSV path (required for CLI) |
| `barcode_type` | int \| str | 0 | `0` or `'terminus'`=terminus optimized (fast), `1` or `'spectrum'`=spectrum optimized (thorough). Also accepts: `'term'`, `'fast'`, `'spec'`, `'slow'` |
| `left_context_column` | str \| None | None | Column for left DNA context (at least one context required) |
| `right_context_column` | str \| None | None | Column for right DNA context (at least one context required) |
| `patch_mode` | bool | False | Fill only missing values (`None`/`NaN`/empty/`'-'`); existing values must already be valid ATGC of length `barcode_length` |
| `cross_barcode_columns` | str \| list[str] \| None | None | Existing barcode column(s) for cross-set separation |
| `minimum_cross_distance` | int \| None | None | Min Hamming distance to cross set (requires `cross_barcode_columns`) |
| `excluded_motifs` | list \| str \| pd.DataFrame \| None | None | Motifs to exclude (list, CSV, DataFrame with `Exmotif`, or FASTA) |
| `verbose` | bool | True | Print progress output |
| `random_seed` | int \| None | None | RNG seed for reproducibility |

**Design Order**: After primers/motifs, before spacers.

**Returns**: `(DataFrame, stats_dict)`

**CLI Equivalent**:
```bash
op barcode \
    --input-data input.csv \
    --oligo-length-limit 200 \
    --barcode-length 12 \
    --minimum-hamming-distance 3 \
    --maximum-repeat-length 8 \
    --barcode-column BC1 \
    --left-context-column Variant \
    --output-file output
```

[↑ Back to TOC](#table-of-contents)

---

### primer

**Purpose**: Design thermodynamically optimal primers for amplification with Tm constraints and optional background screening.

**Signature**:
```python
df, stats = op.primer(
    # Required
    input_data,                    # str | pd.DataFrame
    oligo_length_limit,            # int
    primer_sequence_constraint,    # str
    primer_type,                   # int | str
    minimum_melting_temperature,   # float
    maximum_melting_temperature,   # float
    maximum_repeat_length,         # int
    primer_column,                 # str

    # Optional
    output_file=None,              # str | None
    left_context_column=None,      # str | None
    right_context_column=None,     # str | None
    patch_mode=False,              # bool
    oligo_sets=None,               # list | str | pd.DataFrame | None
    paired_primer_column=None,     # str | None
    excluded_motifs=None,          # list | str | pd.DataFrame | None
    background_directory=None,     # str | None
    verbose=True,                  # bool
    random_seed=None,              # int | None
)
```

**Required Parameters**

| Parameter | Type | Constraints | Description |
|-----------|------|-------------|-------------|
| `input_data` | str \| pd.DataFrame | - | CSV path or DataFrame with `ID` + DNA sequence columns |
| `oligo_length_limit` | int | ≥4 | Maximum allowed oligo length (bp) |
| `primer_sequence_constraint` | str | - | IUPAC constraint string (e.g., `'SS' + 'N'*18` for GC clamp) |
| `primer_type` | int \| str | 0 or 1 | `0` or `'forward'`=forward, `1` or `'reverse'`=reverse primer design. Also accepts: `'fwd'`, `'f'`, `'rev'`, `'r'` |
| `minimum_melting_temperature` | float | ≥25 | Minimum primer Tm (°C) |
| `maximum_melting_temperature` | float | ≤95 | Maximum primer Tm (°C) |
| `maximum_repeat_length` | int | ≥6 | Maximum shared repeat length between primer and oligos/background |
| `primer_column` | str | - | Output column name |

**Optional Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `output_file` | str \| None | None | Output CSV path (required for CLI) |
| `left_context_column` | str \| None | None | Column for left DNA context (at least one context required) |
| `right_context_column` | str \| None | None | Column for right DNA context (at least one context required) |
| `patch_mode` | bool | False | Fill only missing values |
| `oligo_sets` | list \| str \| pd.DataFrame \| None | None | Per-row grouping labels for set-specific primers (list, CSV, or DataFrame with `ID` + `OligoSet`) |
| `paired_primer_column` | str \| None | None | Column of paired primer for Tm matching (within 1°C) |
| `excluded_motifs` | list \| str \| pd.DataFrame \| None | None | Motifs to exclude |
| `background_directory` | str \| None | None | Background k-mer DB from `background()` for off-target screening |
| `verbose` | bool | True | Print progress output |
| `random_seed` | int \| None | None | RNG seed for reproducibility |

**Design Order**: Design primers early. For paired primers: design inner primer first, then outer with `paired_primer_column`.

**Returns**: `(DataFrame, stats_dict)`

**CLI Equivalent**:
```bash
op primer \
    --input-data input.csv \
    --oligo-length-limit 200 \
    --primer-sequence-constraint "'N'*20" \
    --primer-type 0 \
    --minimum-melting-temperature 52 \
    --maximum-melting-temperature 58 \
    --maximum-repeat-length 10 \
    --primer-column FwdPrimer \
    --right-context-column Variant \
    --output-file output
```

[↑ Back to TOC](#table-of-contents)

---

### motif

**Purpose**: Insert sequence motifs (per-variant or constant anchors) with constraint satisfaction.

**Signature**:
```python
df, stats = op.motif(
    # Required
    input_data,                    # str | pd.DataFrame
    oligo_length_limit,            # int
    motif_sequence_constraint,     # str
    maximum_repeat_length,         # int
    motif_column,                  # str

    # Optional
    output_file=None,              # str | None
    motif_type=0,                  # int | str
    left_context_column=None,      # str | None
    right_context_column=None,     # str | None
    patch_mode=False,              # bool
    excluded_motifs=None,          # list | str | pd.DataFrame | None
    verbose=True,                  # bool
    random_seed=None,              # int | None
)
```

**Required Parameters**

| Parameter | Type | Constraints | Description |
|-----------|------|-------------|-------------|
| `input_data` | str \| pd.DataFrame | - | CSV path or DataFrame with `ID` + DNA sequence columns |
| `oligo_length_limit` | int | ≥4 | Maximum allowed oligo length (bp) |
| `motif_sequence_constraint` | str | - | IUPAC constraint string or constant sequence |
| `maximum_repeat_length` | int | ≥4 | Maximum shared repeat length between motif/anchor and oligos |
| `motif_column` | str | - | Output column name |

**Optional Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `output_file` | str \| None | None | Output CSV path (required for CLI) |
| `motif_type` | int \| str | 0 | `0` or `'variable'`=per-variant motifs, `1` or `'constant'`=constant motif/anchor for all rows. Also accepts: `'var'`, `'per-variant'`, `'const'`, `'anchor'` |
| `left_context_column` | str \| None | None | Column for left DNA context (at least one context required) |
| `right_context_column` | str \| None | None | Column for right DNA context (at least one context required) |
| `patch_mode` | bool | False | Fill only missing values; for `motif_type=1`, existing anchor is reused |
| `excluded_motifs` | list \| str \| pd.DataFrame \| None | None | Motifs to exclude |
| `verbose` | bool | True | Print progress output |
| `random_seed` | int \| None | None | RNG seed for reproducibility |

**Design Order**: Before barcodes if designing anchors for indexing.

**Returns**: `(DataFrame, stats_dict)`

**CLI Equivalent**:
```bash
op motif \
    --input-data input.csv \
    --oligo-length-limit 200 \
    --motif-sequence-constraint NNNNNNNN \
    --maximum-repeat-length 6 \
    --motif-column Anchor \
    --motif-type 1 \
    --left-context-column Primer \
    --output-file output
```

[↑ Back to TOC](#table-of-contents)

---

### spacer

**Purpose**: Insert neutral DNA spacers to meet length requirements.

**Signature**:
```python
df, stats = op.spacer(
    # Required
    input_data,                    # str | pd.DataFrame
    oligo_length_limit,            # int
    maximum_repeat_length,         # int
    spacer_column,                 # str

    # Optional
    output_file=None,              # str | None
    spacer_length=None,            # int | list | str | pd.DataFrame | None
    left_context_column=None,      # str | None
    right_context_column=None,     # str | None
    patch_mode=False,              # bool
    excluded_motifs=None,          # list | str | pd.DataFrame | None
    verbose=True,                  # bool
    random_seed=None,              # int | None
)
```

**Required Parameters**

| Parameter | Type | Constraints | Description |
|-----------|------|-------------|-------------|
| `input_data` | str \| pd.DataFrame | - | CSV path or DataFrame with `ID` + DNA sequence columns |
| `oligo_length_limit` | int | ≥4 | Maximum allowed oligo length (bp) |
| `maximum_repeat_length` | int | ≥4 | Maximum shared repeat length between spacer and oligos |
| `spacer_column` | str | - | Output column name |

**Optional Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `output_file` | str \| None | None | Output CSV path (required for CLI) |
| `spacer_length` | int \| list \| str \| pd.DataFrame \| None | None | `None`=auto-fill to `oligo_length_limit`, `int`=fixed length, `list`=per-row lengths, DataFrame/CSV with `ID` + `Length` |
| `left_context_column` | str \| None | None | Column for left DNA context (at least one context required) |
| `right_context_column` | str \| None | None | Column for right DNA context (at least one context required) |
| `patch_mode` | bool | False | Fill only missing values |
| `excluded_motifs` | list \| str \| pd.DataFrame \| None | None | Motifs to exclude |
| `verbose` | bool | True | Print progress output |
| `random_seed` | int \| None | None | RNG seed for reproducibility |

**Design Order**: Last, after all other elements.

**Returns**: `(DataFrame, stats_dict)`

**CLI Equivalent**:
```bash
op spacer \
    --input-data input.csv \
    --oligo-length-limit 200 \
    --maximum-repeat-length 8 \
    --spacer-column Spacer \
    --left-context-column Primer2 \
    --output-file output
```

[↑ Back to TOC](#table-of-contents)

---

### background

**Purpose**: Build a k-mer database for screening primers against off-target sequences.

**Signature**:
```python
stats = op.background(
    # Required
    input_data,                    # list | str | pd.DataFrame
    maximum_repeat_length,         # int
    output_directory,              # str

    # Optional
    verbose=True,                  # bool
)
```

**Required Parameters**

| Parameter | Type | Constraints | Description |
|-----------|------|-------------|-------------|
| `input_data` | list \| str \| pd.DataFrame | - | Background sequences: list of DNA strings, CSV/DataFrame with `Sequence` column, or FASTA path |
| `maximum_repeat_length` | int | 6-20 | Maximum repeat length to screen (k-mer size = `maximum_repeat_length + 1`) |
| `output_directory` | str | - | Output directory basename (writes `<name>.oligopool.background`) |

**Optional Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `verbose` | bool | True | Print progress output |

**Returns**: `stats_dict` (stats-only, no DataFrame)

**CLI Equivalent**:
```bash
op background \
    --input-data genome.fasta \
    --maximum-repeat-length 15 \
    --output-directory ecoli_bg
```

[↑ Back to TOC](#table-of-contents)

---

### split

**Purpose**: Break long oligos into overlapping fragments for overlap-based assembly (Gibson, overlap-extension PCR, etc.).

**Signature**:
```python
df, stats = op.split(
    # Required
    input_data,                    # str | pd.DataFrame
    split_length_limit,            # int
    minimum_melting_temperature,   # float
    minimum_hamming_distance,      # int
    minimum_overlap_length,        # int
    maximum_overlap_length,        # int

    # Optional
    output_file=None,              # str | None
    random_seed=None,              # int | None
    separate_outputs=False,        # bool
    verbose=True,                  # bool
)
```

**Required Parameters**

| Parameter | Type | Constraints | Description |
|-----------|------|-------------|-------------|
| `input_data` | str \| pd.DataFrame | - | CSV path or DataFrame with `ID` + DNA sequence columns |
| `split_length_limit` | int | - | Maximum fragment length (bp) |
| `minimum_melting_temperature` | float | - | Minimum overlap Tm (°C) |
| `minimum_hamming_distance` | int | - | Minimum pairwise Hamming distance between overlaps |
| `minimum_overlap_length` | int | - | Minimum overlap length (bp) |
| `maximum_overlap_length` | int | - | Maximum overlap length (bp) |

**Optional Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `output_file` | str \| None | None | Output CSV path (required for CLI) |
| `random_seed` | int \| None | None | RNG seed for reproducibility |
| `separate_outputs` | bool | False | When enabled, return a list of per-split DataFrames; if `output_file` is set, write `{base}.SplitN.oligopool.split.csv` files |
| `verbose` | bool | True | Print progress output |

**Returns**:
- `(DataFrame, stats_dict)` when `separate_outputs` is disabled (default) - output contains `Split1`, `Split2`, ... columns.
- `([DataFrame, ...], stats_dict)` when `separate_outputs` is enabled - one DataFrame per `SplitN` column.

**Notes**:
- Number of fragments varies per oligo; even-numbered splits (`Split2`, `Split4`, ...) are reverse-complemented.
- Original annotation columns are not preserved in output.
- Treat each `SplitN` column as its own synthesis pool; run `pad` per `SplitN`, then `final` each padded DataFrame (don’t run `final()` on the raw multi-column `split` output).

**CLI Equivalent**:
```bash
op split \
    --input-data library.csv \
    --split-length-limit 150 \
    --minimum-melting-temperature 55 \
    --minimum-hamming-distance 3 \
    --minimum-overlap-length 20 \
    --maximum-overlap-length 30 \
    --output-file split_library
```
Tip: CLI defaults to separate files. Use `--no-separate-outputs` to write a single combined `.oligopool.split.csv`.

[↑ Back to TOC](#table-of-contents)

---

### pad

**Purpose**: Add primers and Type IIS restriction sites to split fragments for scarless assembly.

**Signature**:
```python
df, stats = op.pad(
    # Required
    input_data,                    # str | pd.DataFrame
    oligo_length_limit,            # int
    split_column,                  # str
    typeIIS_system,                # str
    minimum_melting_temperature,   # float
    maximum_melting_temperature,   # float
    maximum_repeat_length,         # int

    # Optional
    output_file=None,              # str | None
    random_seed=None,              # int | None
    verbose=True,                  # bool
)
```

**Required Parameters**

| Parameter | Type | Constraints | Description |
|-----------|------|-------------|-------------|
| `input_data` | str \| pd.DataFrame | - | Output from `split()` or any DataFrame with fragment column |
| `oligo_length_limit` | int | ≥60 | Maximum padded fragment length (bp) |
| `split_column` | str | - | Which fragment column to pad (e.g., `Split1`) |
| `typeIIS_system` | str | - | Type IIS enzyme name (see [supported list](#type-iis-enzymes)) |
| `minimum_melting_temperature` | float | ≥25 | Minimum pad primer Tm (°C) |
| `maximum_melting_temperature` | float | ≤95 | Maximum pad primer Tm (°C) |
| `maximum_repeat_length` | int | 6-20 | Maximum shared repeat length |

**Optional Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `output_file` | str \| None | None | Output CSV path (required for CLI) |
| `random_seed` | int \| None | None | RNG seed for reproducibility |
| `verbose` | bool | True | Print progress output |

**Returns**: `(DataFrame, stats_dict)` - Output contains `5primeSpacer`, `ForwardPrimer`, `<split_column>`, `ReversePrimer`, `3primeSpacer`.

**Important**: The chosen Type IIS recognition site must be absent from all split fragments (in either orientation). `pad` validates this and fails early if internal cut sites are found. To avoid issues, exclude your Type IIS motif (e.g., `GGTCTC`/`GAGACC` for BsaI) from `excluded_motifs` when designing upstream elements (primers, barcodes, motifs, spacers). Note: this only constrains newly designed elements — if your input/core sequences already contain the site, you’ll need to choose a different Type IIS system or redesign those sequences.

**CLI Equivalent**:
```bash
op pad \
    --input-data split_library.csv \
    --oligo-length-limit 200 \
    --split-column Split1 \
    --typeiis-system BsaI \
    --minimum-melting-temperature 52 \
    --maximum-melting-temperature 58 \
    --maximum-repeat-length 10 \
    --output-file padded_split1
```

[↑ Back to TOC](#table-of-contents)

---

### merge

**Purpose**: Concatenate contiguous columns into a single column.

**Signature**:
```python
df, stats = op.merge(
    # Required
    input_data,                    # str | pd.DataFrame
    merge_column,                  # str

    # Optional
    output_file=None,              # str | None
    left_context_column=None,      # str | None
    right_context_column=None,     # str | None
    verbose=True,                  # bool
)
```

**Required Parameters**

| Parameter | Type | Constraints | Description |
|-----------|------|-------------|-------------|
| `input_data` | str \| pd.DataFrame | - | CSV path or DataFrame with `ID` + DNA sequence columns |
| `merge_column` | str | - | Output column name for merged region |

**Optional Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `output_file` | str \| None | None | Output CSV path (required for CLI) |
| `left_context_column` | str \| None | None | First column to merge (defaults to first sequence column) |
| `right_context_column` | str \| None | None | Last column to merge (defaults to last sequence column) |
| `verbose` | bool | True | Print progress output |

**Returns**: `(DataFrame, stats_dict)` - Source columns in range are removed.

**CLI Equivalent**:
```bash
op merge \
    --input-data library.csv \
    --merge-column 5prime_region \
    --left-context-column Primer5 \
    --right-context-column BC1 \
    --output-file merged
```

[↑ Back to TOC](#table-of-contents)

---

### revcomp

**Purpose**: Reverse complement a range of columns and reverse their order.

**Signature**:
```python
df, stats = op.revcomp(
    # Required
    input_data,                    # str | pd.DataFrame

    # Optional
    output_file=None,              # str | None
    left_context_column=None,      # str | None
    right_context_column=None,     # str | None
    verbose=True,                  # bool
)
```

**Required Parameters**

| Parameter | Type | Constraints | Description |
|-----------|------|-------------|-------------|
| `input_data` | str \| pd.DataFrame | - | CSV path or DataFrame with `ID` + DNA sequence columns |

**Optional Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `output_file` | str \| None | None | Output CSV path (required for CLI) |
| `left_context_column` | str \| None | None | First column to revcomp (defaults to first sequence column) |
| `right_context_column` | str \| None | None | Last column to revcomp (defaults to last sequence column) |
| `verbose` | bool | True | Print progress output |

**Returns**: `(DataFrame, stats_dict)`

**CLI Equivalent**:
```bash
op revcomp \
    --input-data library.csv \
    --left-context-column Gene \
    --right-context-column Primer3 \
    --output-file revcomped
```

[↑ Back to TOC](#table-of-contents)

---

### lenstat

**Purpose**: Report length statistics and free space remaining (non-destructive).

**Signature**:
```python
stats = op.lenstat(
    # Required
    input_data,                    # str | pd.DataFrame
    oligo_length_limit,            # int

    # Optional
    verbose=True,                  # bool
)
```

**Required Parameters**

| Parameter | Type | Constraints | Description |
|-----------|------|-------------|-------------|
| `input_data` | str \| pd.DataFrame | - | CSV path or DataFrame with `ID` + DNA sequence columns |
| `oligo_length_limit` | int | ≥4 | Reference length limit for free space calculation |

**Optional Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `verbose` | bool | True | Print progress output |

**Returns**: `stats_dict` (stats-only, no DataFrame, no `output_file`)

**CLI Equivalent**:
```bash
op lenstat \
    --input-data library.csv \
    --oligo-length-limit 200
```

[↑ Back to TOC](#table-of-contents)

---

### verify

**Purpose**: QC check for constraints, architecture, motifs, degeneracy.

**Signature**:
```python
stats = op.verify(
    # Required
    input_data,                    # str | pd.DataFrame

    # Optional
    oligo_length_limit=None,       # int | None
    excluded_motifs=None,          # list | str | pd.DataFrame | None
    verbose=True,                  # bool
)
```

**Required Parameters**

| Parameter | Type | Constraints | Description |
|-----------|------|-------------|-------------|
| `input_data` | str \| pd.DataFrame | - | CSV path or DataFrame with `ID` column (can include metadata/IUPAC columns) |

**Optional Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `oligo_length_limit` | int \| None | None | If provided, checks for length overflow |
| `excluded_motifs` | list \| str \| pd.DataFrame \| None | None | Motifs to scan/report (emergence; junction attribution requires separate sequence columns) |
| `verbose` | bool | True | Print progress output |

**Returns**: `stats_dict` (stats-only, no DataFrame, no `output_file`)

**Column Concatenation**:
- Only **sequence columns** (DNA/IUPAC; may include `'-'`) are concatenated; metadata columns are skipped
- Sequence columns are concatenated **left-to-right in DataFrame column order**
- Gap characters (`'-'`) are **stripped** during concatenation
- If a `CompleteOligo` column exists (from `final()`), it is used directly instead of concatenating
- Junction identification for edge-effect analysis follows this same column order:
  - For columns `[Primer1, BC1, Variant, Primer2]`, junctions are `Primer1|BC1`, `BC1|Variant`, `Variant|Primer2`

**Notes**:
- More permissive than design modules; handles metadata and degenerate/IUPAC columns.
- Reports emergent motifs (occurrences beyond baseline minimum) and attributes them to column junctions.
- Excluded-motif matching is literal substring matching; degenerate/IUPAC bases are not expanded as wildcards.
- Run `verify` **before** `final()` to preserve separate columns for junction attribution.

**CLI Equivalent**:
```bash
op verify \
    --input-data library.csv \
    --oligo-length-limit 200 \
    --excluded-motifs GAATTC,GGATCC
```

[↑ Back to TOC](#table-of-contents)

---

### final

**Purpose**: Concatenate all columns into synthesis-ready oligos.

**Signature**:
```python
df, stats = op.final(
    # Required
    input_data,                    # str | pd.DataFrame

    # Optional
    output_file=None,              # str | None
    verbose=True,                  # bool
)
```

**Required Parameters**

| Parameter | Type | Constraints | Description |
|-----------|------|-------------|-------------|
| `input_data` | str \| pd.DataFrame | - | CSV path or DataFrame with `ID` + DNA sequence columns |

**Optional Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `output_file` | str \| None | None | Output CSV path (required for CLI) |
| `verbose` | bool | True | Print progress output |

**Returns**: `(DataFrame, stats_dict)` - Output contains `ID`, `CompleteOligo`, `OligoLength` only.

**CLI Equivalent**:
```bash
op final \
    --input-data library.csv \
    --output-file synthesis_ready
```

[↑ Back to TOC](#table-of-contents)

---

## Analysis Mode

### index

**Purpose**: Build a searchable index of barcodes (and optional associates) for counting.

**Signature**:
```python
stats = op.index(
    # Required
    barcode_data,                  # str | pd.DataFrame
    barcode_column,                # str
    index_file,                    # str

    # Optional (barcode anchors)
    barcode_prefix_column=None,    # str | None
    barcode_suffix_column=None,    # str | None
    barcode_prefix_gap=0,          # int
    barcode_suffix_gap=0,          # int

    # Optional (associate indexing)
    associate_data=None,           # str | pd.DataFrame | None
    associate_column=None,         # str | None
    associate_prefix_column=None,  # str | None
    associate_suffix_column=None,  # str | None

    # Optional
    verbose=True,                  # bool
)
```

**Required Parameters**

| Parameter | Type | Constraints | Description |
|-----------|------|-------------|-------------|
| `barcode_data` | str \| pd.DataFrame | - | CSV path or DataFrame with `ID` + barcode column |
| `barcode_column` | str | - | Barcode column name to index |
| `index_file` | str | - | Output basename (writes `<name>.oligopool.index`) |

**Optional Parameters (Barcode Anchors)**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `barcode_prefix_column` | str \| None | None | Column for constant prefix anchor (at least one anchor required) |
| `barcode_suffix_column` | str \| None | None | Column for constant suffix anchor (at least one anchor required) |
| `barcode_prefix_gap` | int | 0 | Bases between prefix anchor and barcode in read |
| `barcode_suffix_gap` | int | 0 | Bases between barcode and suffix anchor in read |

**Optional Parameters (Associate Indexing)**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `associate_data` | str \| pd.DataFrame \| None | None | CSV path or DataFrame with associates (can be same as `barcode_data`) |
| `associate_column` | str \| None | None | Column for associate elements |
| `associate_prefix_column` | str \| None | None | Column for constant associate prefix anchor |
| `associate_suffix_column` | str \| None | None | Column for constant associate suffix anchor (at least one anchor required when using associates) |

**Optional Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `verbose` | bool | True | Print progress output |

**Returns**: `stats_dict` (stats-only, no DataFrame)

**Notes**:
- Anchors must be constant (single-unique) sequences, ideally ≥6 bp and adjacent to indexed column.
- At least one of `barcode_prefix_column` or `barcode_suffix_column` is required.

**CLI Equivalent**:
```bash
op index \
    --barcode-data library.csv \
    --barcode-column BC1 \
    --barcode-prefix-column Primer1 \
    --index-file bc1_index
```

[↑ Back to TOC](#table-of-contents)

---

### pack

**Purpose**: Preprocess FastQ files - filter, optionally merge paired-ends, deduplicate.

**Signature**:
```python
stats = op.pack(
    # Required
    r1_fastq_file,                 # str
    r1_read_type,                  # int | str
    pack_type,                     # int | str
    pack_file,                     # str

    # Optional (R1 filters)
    minimum_r1_read_length=1,      # int
    minimum_r1_read_quality=20,    # int

    # Optional (paired-end)
    r2_fastq_file=None,            # str | None
    r2_read_type=None,             # int | str | None
    minimum_r2_read_length=None,   # int | None
    minimum_r2_read_quality=None,  # int | None

    # Optional (performance)
    pack_size=3.0,                 # float
    core_count=0,                  # int
    memory_limit=0.0,              # float

    # Optional
    verbose=True,                  # bool
)
```

**Required Parameters**

| Parameter | Type | Constraints | Description |
|-----------|------|-------------|-------------|
| `r1_fastq_file` | str | - | R1 FastQ path (supports `.gz`) |
| `r1_read_type` | int \| str | 0 or 1 | `0` or `'forward'`=forward, `1` or `'reverse'`=reverse orientation. Also accepts: `'fwd'`, `'f'`, `'rev'`, `'r'` |
| `pack_type` | int \| str | 0 or 1 | `0` or `'concatenate'`=concatenate pairs, `1` or `'merge'`=merge/assemble pairs. Also accepts: `'concatenated'`, `'joined'`, `'join'`, `'merged'`, `'concat'`, `'cat'`, `'assemble'`, `'assembled'`, `'asm'` |
| `pack_file` | str | - | Output basename (writes `<name>.oligopool.pack`) |

**Optional Parameters (R1 Filters)**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `minimum_r1_read_length` | int | 1 | Minimum R1 read length |
| `minimum_r1_read_quality` | int | 20 | Minimum average R1 Phred score |

**Optional Parameters (Paired-End)**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `r2_fastq_file` | str \| None | None | R2 FastQ path |
| `r2_read_type` | int \| str \| None | None | R2 orientation: `0` or `'forward'`=forward, `1` or `'reverse'`=reverse. Also accepts: `'fwd'`, `'f'`, `'rev'`, `'r'` |
| `minimum_r2_read_length` | int \| None | None | Minimum R2 read length |
| `minimum_r2_read_quality` | int \| None | None | Minimum average R2 Phred score |

**Optional Parameters (Performance)**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pack_size` | float | 3.0 | Million unique reads per pack (0.1-5.0) |
| `core_count` | int | 0 | CPU cores (`0`=auto) |
| `memory_limit` | float | 0.0 | GB per core (`0`=auto) |

**Optional Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `verbose` | bool | True | Print progress output |

**Returns**: `stats_dict` (stats-only, no DataFrame)

**CLI Equivalent**:
```bash
op pack \
    --r1-fastq-file reads_R1.fq.gz \
    --r2-fastq-file reads_R2.fq.gz \
    --r1-read-type 0 \
    --r2-read-type 1 \
    --pack-type 1 \
    --minimum-r1-read-quality 30 \
    --minimum-r2-read-quality 30 \
    --pack-file sample
```

[↑ Back to TOC](#table-of-contents)

---

### acount

**Purpose**: Association counting - verify barcode-variant coupling in reads.

**Signature**:
```python
counts_df, stats = op.acount(
    # Required
    index_file,                    # str
    pack_file,                     # str
    count_file,                    # str

    # Optional
    mapping_type=0,                # int | str
    barcode_errors=-1,             # int
    associate_errors=-1,           # int
    callback=None,                 # callable | None
    core_count=0,                  # int
    memory_limit=0.0,              # float
    verbose=True,                  # bool
)
```

**Required Parameters**

| Parameter | Type | Constraints | Description |
|-----------|------|-------------|-------------|
| `index_file` | str | - | Index basename or path (reads `<name>.oligopool.index`) |
| `pack_file` | str | - | Pack basename or path (reads `<name>.oligopool.pack`) |
| `count_file` | str | - | Output basename (writes `<name>.oligopool.acount.csv`) |

**Optional Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mapping_type` | int \| str | 0 | `0` or `'fast'`=fast/near-exact, `1` or `'sensitive'`=slow/sensitive. Also accepts: `'quick'`, `'sens'`, `'accurate'` |
| `barcode_errors` | int | -1 | Max barcode errors (`-1`=auto from index) |
| `associate_errors` | int | -1 | Max associate errors (`-1`=auto from index) |
| `callback` | callable \| None | None | Custom read filter function (Python API only) |
| `core_count` | int | 0 | CPU cores (`0`=auto) |
| `memory_limit` | float | 0.0 | GB per core (`0`=auto) |
| `verbose` | bool | True | Print progress output |

**Returns**: `(counts_DataFrame, stats_dict)` - Output contains `<indexname>.ID`, `BarcodeCounts`, `AssociationCounts`.

**Callback Signature** (Python only):
```python
def callback(read, ID, count, coreid):
    """
    Args:
        read: str - Processed read sequence
        ID: tuple - Identified barcode IDs
        count: int - Read frequency in current pack
        coreid: int - CPU core ID
    Returns:
        bool - True to accept, False to reject
    """
    return True
```

**CLI Equivalent**:
```bash
op acount \
    --index-file bc1_index \
    --pack-file sample \
    --count-file results \
    --mapping-type 1
```

[↑ Back to TOC](#table-of-contents)

---

### xcount

**Purpose**: Barcode-only counting (single or combinatorial) without associate verification.

**Signature**:
```python
counts_df, stats = op.xcount(
    # Required
    index_files,                   # str | list
    pack_file,                     # str
    count_file,                    # str

    # Optional
    mapping_type=0,                # int | str
    barcode_errors=-1,             # int
    callback=None,                 # callable | None
    core_count=0,                  # int
    memory_limit=0.0,              # float
    verbose=True,                  # bool
)
```

**Required Parameters**

| Parameter | Type | Constraints | Description |
|-----------|------|-------------|-------------|
| `index_files` | str \| list | - | Single index (string) or multiple (list) for combinatorial counting |
| `pack_file` | str | - | Pack basename or path (reads `<name>.oligopool.pack`) |
| `count_file` | str | - | Output basename (writes `<name>.oligopool.xcount.csv`) |

**Optional Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mapping_type` | int \| str | 0 | `0` or `'fast'`=fast/near-exact, `1` or `'sensitive'`=slow/sensitive. Also accepts: `'quick'`, `'sens'`, `'accurate'` |
| `barcode_errors` | int | -1 | Max barcode errors (`-1`=auto from index) |
| `callback` | callable \| None | None | Custom read filter function (Python API only) |
| `core_count` | int | 0 | CPU cores (`0`=auto) |
| `memory_limit` | float | 0.0 | GB per core (`0`=auto) |
| `verbose` | bool | True | Print progress output |

**Returns**: `(counts_DataFrame, stats_dict)` - Output contains one `<indexname>.ID` column per index, plus `CombinatorialCounts`. Missing barcodes shown as `'-'`.

**CLI Equivalent**:
```bash
# Single barcode
op xcount \
    --index-files bc1_index \
    --pack-file sample \
    --count-file counts

# Combinatorial (BC1 x BC2)
op xcount \
    --index-files bc1_index,bc2_index \
    --pack-file sample \
    --count-file combo_counts
```

[↑ Back to TOC](#table-of-contents)

---

## Advanced

### vectorDB

**Purpose**: LevelDB-based k-mer storage. Created by `background()`, accessible directly for custom workflows.

**Constructor**:
```python
db = op.vectorDB(
    path,                          # str - Directory path
    maximum_repeat_length,         # int - Max repeat length (k = max_repeat_length + 1)
)
```

**Methods**:
```python
# Add k-mers
db.add(seq, rna=True)              # Add k-mers from single sequence
db.multiadd(seq_list, rna=True)    # Add k-mers from multiple sequences

# Check membership
seq in db                          # True if any k-mer from seq exists
db.multicheck(seq_list, rna=True)  # Returns list[bool]

# Iterate
for kmer in db:
    print(kmer)
len(db)                            # Number of stored k-mers

# Remove k-mers
db.remove(seq, rna=True)
db.multiremove(seq_list, rna=True, clear=False)

# Database management
db.clear(maxreplen=None)           # Wipe DB; optionally set new max repeat length
db.close()                         # Close handle, keep files
db.drop()                          # Delete DB from disk
```

**Parameters**:
- `path`: Directory path (background outputs end with `.oligopool.background`)
- `maximum_repeat_length`: Max repeat length. If reopening existing DB, stored value is used.
- `rna`: If `True`, convert U→T before processing

[↑ Back to TOC](#table-of-contents)

---

### Scry

**Purpose**: 1-nearest-neighbor barcode classifier. Powers `acount`/`xcount` internally.

**Constructor and Methods**:
```python
model = op.Scry()

# Train
model.train(
    X,                             # iterable - Barcode sequences
    Y,                             # iterable - Labels/IDs
    n,                             # int - Barcode length
    k,                             # int - k-mer length for indexing
    t,                             # int - Max errors for model building
    liner=None,                    # callable - Progress printer (optional)
)

# Prime for prediction
model.prime(
    t=0,                           # int - Errors to consider at prediction time
    mode=0,                        # int - 0=fast, 1=sensitive
)

# Predict
label, score = model.predict(x)    # Returns (label, score) or (None, None)
```

**Example**:
```python
model = op.Scry()
model.train(
    X=['ATGCATGC', 'GCTAGCTA', 'TTAATTAA'],
    Y=['bc1', 'bc2', 'bc3'],
    n=8, k=4, t=1
)
model.prime(t=1, mode=0)
label, score = model.predict('ATGCATGC')  # ('bc1', 1.0)
```

[↑ Back to TOC](#table-of-contents)

---

## Appendix

### IUPAC Codes

| Code | Bases | Description |
|------|-------|-------------|
| A | A | Adenine |
| T | T | Thymine |
| G | G | Guanine |
| C | C | Cytosine |
| N | A, T, G, C | Any base |
| W | A, T | Weak (2 H-bonds) |
| S | G, C | Strong (3 H-bonds) |
| M | A, C | Amino |
| K | G, T | Keto |
| R | A, G | Purine |
| Y | C, T | Pyrimidine |
| B | C, G, T | Not A |
| D | A, G, T | Not C |
| H | A, C, T | Not G |
| V | A, C, G | Not T |

[↑ Back to TOC](#table-of-contents)

---

### Common Restriction Sites

| Enzyme | Recognition Sequence |
|--------|---------------------|
| EcoRI | GAATTC |
| BamHI | GGATCC |
| HindIII | AAGCTT |
| XbaI | TCTAGA |
| SalI | GTCGAC |
| PstI | CTGCAG |
| NotI | GCGGCCGC |
| XhoI | CTCGAG |
| NcoI | CCATGG |
| NdeI | CATATG |
| BglII | AGATCT |
| SpeI | ACTAGT |
| AatII | GACGTC |
| KpnI | GGTACC |

[↑ Back to TOC](#table-of-contents)

---

### Type IIS Enzymes

Supported enzymes for `pad()` (34 total):

`AcuI`, `AlwI`, `BbsI`, `BccI`, `BceAI`, `BciVI`, `BcoDI`, `BmrI`, `BpuEI`, `BsaI`, `BseRI`, `BsmAI`, `BsmBI`, `BsmFI`, `BsmI`, `BspCNI`, `BspQI`, `BsrDI`, `BsrI`, `BtgZI`, `BtsCI`, `BtsI`, `BtsIMutI`, `EarI`, `EciI`, `Esp3I`, `FauI`, `HgaI`, `HphI`, `HpyAV`, `MlyI`, `MnlI`, `SapI`, `SfaNI`

For other enzymes, design primers/sites manually with `primer` or `motif`.

[↑ Back to TOC](#table-of-contents)

---

### File Formats

**Input CSV**:
```csv
ID,Variant
v1,ATGCATGCATGCATGC
v2,GCTAGCTAGCTAGCTA
v3,TTAACCGGTTAACCGG
```

**Output File Extensions**:

| Extension | Module | Description |
|-----------|--------|-------------|
| `.oligopool.background` | background | K-mer database directory |
| `.oligopool.barcode.csv` | barcode | Barcode-annotated library |
| `.oligopool.primer.csv` | primer | Primer-annotated library |
| `.oligopool.motif.csv` | motif | Motif/anchor-annotated library |
| `.oligopool.spacer.csv` | spacer | Spacer-annotated library |
| `.oligopool.split.csv` | split | Split-fragment library |
| `.oligopool.pad.csv` | pad | Padded fragment library |
| `.oligopool.merge.csv` | merge | Merged element library |
| `.oligopool.revcomp.csv` | revcomp | Reverse-complemented library |
| `.oligopool.final.csv` | final | Synthesis-ready library |
| `.oligopool.index` | index | Barcode index file |
| `.oligopool.pack` | pack | Packed reads file |
| `.oligopool.acount.csv` | acount | Association counts |
| `.oligopool.xcount.csv` | xcount | Combinatorial counts |

[↑ Back to TOC](#table-of-contents)

---

### CLI Parameter Mapping

Python parameters map to CLI flags by converting `snake_case` to `kebab-case`:

| Python | CLI |
|--------|-----|
| `input_data` | `--input-data` |
| `oligo_length_limit` | `--oligo-length-limit` |
| `barcode_length` | `--barcode-length` |
| `minimum_hamming_distance` | `--minimum-hamming-distance` |
| `maximum_repeat_length` | `--maximum-repeat-length` |
| `barcode_column` | `--barcode-column` |
| `output_file` | `--output-file` |
| `barcode_type` | `--barcode-type` |
| `left_context_column` | `--left-context-column` |
| `right_context_column` | `--right-context-column` |
| `patch_mode` | `--patch-mode` |
| `cross_barcode_columns` | `--cross-barcode-columns` |
| `minimum_cross_distance` | `--minimum-cross-distance` |
| `excluded_motifs` | `--excluded-motifs` |
| `verbose` | `--quiet` (inverted) |
| `random_seed` | `--random-seed` |

**CLI-Specific Flags**:
- `--stats-json`: Print stats dict as JSON to stdout
- `--stats-file PATH`: Write stats JSON to file
- `--quiet`: Suppress verbose output (sets `verbose=False`)

**Sequence Constraint Shorthand** (CLI):
```bash
--primer-sequence-constraint "'N'*20"
--motif-sequence-constraint GCC+N*20+CCG
```

[↑ Back to TOC](#table-of-contents)

---

<p align="center">
<i>Made with molecules and math by the Salis Lab</i>
<br>
<a href="https://github.com/ayaanhossain/oligopool">GitHub</a> ·
<a href="https://pubs.acs.org/doi/10.1021/acssynbio.4c00661">Paper</a> ·
<a href="https://github.com/ayaanhossain/oligopool/issues">Issues</a>
</p>
