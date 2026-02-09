# Oligopool Calculator - API Reference

Complete parameter reference for all `oligopool` modules.

**Related documentation**:
- [README](../README.md) - Overview, installation, and quick start
- [User Guide](docs.md) - Tutorials, examples, and workflows
- [Example Notebook](../examples/OligopoolCalculatorInAction.ipynb) - Interactive design and analysis demo
- [AI Agent Guide](agent-skills.md) - Decision trees, best practices, and gotchas
- [Docker Guide](docker-notes.md) - Run `oligopool` in a container

---

## Table of Contents

**Design Mode**
- [`barcode`](#barcode) - Hamming-distance barcodes
- [`primer`](#primer) - Thermodynamic primers
- [`motif`](#motif) - Sequence motifs & anchors
- [`spacer`](#spacer) - Neutral spacers
- [`background`](#background) - K-mer screening database
- [`merge`](#merge) - Collapse columns
- [`revcomp`](#revcomp) - Reverse complement
- [`final`](#final) - Finalize for synthesis

**Assembly Mode**
- [`split`](#split) - Fragment long oligos
- [`pad`](#pad) - Assembly-ready padding

**Degenerate Mode**
- [`compress`](#compress) - Compress to IUPAC-degenerate oligos
- [`expand`](#expand) - Expand degenerate oligos

**Analysis Mode**
- [`index`](#index) - Build barcode index
- [`pack`](#pack) - Preprocess FastQ
- [`acount`](#acount) - Association counting
- [`xcount`](#xcount) - Combinatorial counting

**QC Mode**
- [`lenstat`](#lenstat) - Length statistics
- [`verify`](#verify) - Conflict detection
- [`inspect`](#inspect) - Inspect artifacts

**Advanced**
- [`vectorDB`](#vectordb) - K-mer database
- [`Scry`](#scry) - Barcode classifier

**Appendix**
- [IUPAC Codes](#iupac-codes)
- [Common Restriction Sites](#common-restriction-sites)
- [Type IIS Enzymes](#type-iis-enzymes)
- [File Formats](#file-formats)
- [CLI Parameter Mapping](#cli-parameter-mapping)

---

## Design Mode

### `barcode`

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
    cross_barcode_columns=None,    # str | list | None
    minimum_cross_distance=None,   # int | None
    excluded_motifs=None,          # list | str | dict | pd.DataFrame | None
    background_directory=None,     # str | list | None
    random_seed=None,              # int | None
    verbose=True,                  # bool
)
```

**Required Parameters**

- `input_data` (str | DataFrame): CSV path or DataFrame with `ID` column + DNA sequence columns
- `oligo_length_limit` (int, ≥4): Maximum allowed oligo length (bp)
- `barcode_length` (int, ≥4): Length of designed barcodes (bp)
- `minimum_hamming_distance` (int, ≥1): Minimum pairwise Hamming distance within newly designed set
- `maximum_repeat_length` (int, ≥4): Maximum shared repeat length between barcode and context/oligos
- `barcode_column` (str): Column name to create/overwrite (or patch-fill when `patch_mode=True`)

**Optional Parameters**

- `output_file` (str | None, default=None): Output CSV path (required for CLI)
- `barcode_type` (int | str, default=0): `0`/`'terminus'`=fast terminus-optimized (targets distinctive 5'/3' ends), `1`/`'spectrum'`=thorough spectrum-optimized (targets k-mer saturation). Also accepts: `'term'`, `'t'`, `'fast'`, `'spec'`, `'s'`, `'slow'`
- `left_context_column` (str | None, default=None): Column for left DNA context
- `right_context_column` (str | None, default=None): Column for right DNA context
- `patch_mode` (bool, default=False): Fill only missing values (`None`/`NaN`/empty/`'-'`); existing values must be valid ATGC of length `barcode_length`
- `cross_barcode_columns` (str | list | None, default=None): Existing barcode column(s) for cross-set separation
- `minimum_cross_distance` (int | None, default=None): Min Hamming distance to cross set (requires `cross_barcode_columns`)
- `excluded_motifs` (list | str | dict | DataFrame | None, default=None): Motifs to exclude. Accepts a list, CSV/FASTA path, DataFrame with `Exmotif` column, comma-string, or multiple sources (list of sources or `{name: source}` dict). Motifs must be strict ATGC only (no IUPAC codes, no dashes).
- `background_directory` (str | list | None, default=None): Background k-mer DB(s) from `background()` (screened against ALL specified DBs; junction-aware when context columns are provided)
- `random_seed` (int | None, default=None): RNG seed for reproducibility
- `verbose` (bool, default=True): Print progress output

**Returns**: `(DataFrame, stats_dict)`

**Notes**:
- At least one of `left_context_column` or `right_context_column` must be specified
- Design order: after primers/motifs, before spacers
- Cross-set mode is global (not per-row): each new barcode must be ≥`minimum_cross_distance` away from every barcode in the union of `cross_barcode_columns`
- If design is challenging: adjust `barcode_length`, `minimum_hamming_distance`, `maximum_repeat_length`, or switch `barcode_type`
- Constant anchors (e.g., index prefix/suffix) are typically designed first via `motif` with `motif_type=1`

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

### `primer`

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
    excluded_motifs=None,          # list | str | dict | pd.DataFrame | None
    background_directory=None,     # str | list | None
    random_seed=None,              # int | None
    verbose=True,                  # bool
)
```

**Required Parameters**

- `input_data` (str | DataFrame): CSV path or DataFrame with `ID` column + DNA sequence columns
- `oligo_length_limit` (int, ≥4): Maximum allowed oligo length (bp)
- `primer_sequence_constraint` (str): IUPAC constraint string (e.g., `'SS' + 'N'*18` for GC clamp)
- `primer_type` (int | str): `0`/`'forward'`=forward, `1`/`'reverse'`=reverse. Also accepts: `'fwd'`, `'f'`, `'rev'`, `'r'`
- `minimum_melting_temperature` (float, ≥25): Minimum primer Tm (°C)
- `maximum_melting_temperature` (float, ≤95): Maximum primer Tm (°C)
- `maximum_repeat_length` (int, ≥6): Maximum shared repeat length between primer and oligos/background
- `primer_column` (str): Output column name

**Optional Parameters**

- `output_file` (str | None, default=None): Output CSV path (required for CLI)
- `left_context_column` (str | None, default=None): Column for left DNA context
- `right_context_column` (str | None, default=None): Column for right DNA context
- `patch_mode` (bool, default=False): Fill only missing values (`None`/`NaN`/empty/`'-'`)
- `oligo_sets` (list | str | DataFrame | None, default=None): Per-row grouping labels for set-specific primers (list, CSV, or DataFrame with `ID` + `OligoSet`)
- `paired_primer_column` (str | None, default=None): Column of paired primer for Tm matching (within 1°C)
- `excluded_motifs` (list | str | dict | DataFrame | None, default=None): Motifs to exclude. Accepts a list, CSV/FASTA path, DataFrame with `Exmotif` column, comma-string, or multiple sources. Strict ATGC only (no IUPAC codes, no dashes).
- `background_directory` (str | list | None, default=None): Background k-mer DB(s) from `background()` (screened against ALL specified DBs)
- `random_seed` (int | None, default=None): RNG seed for reproducibility
- `verbose` (bool, default=True): Print progress output

**Returns**: `(DataFrame, stats_dict)`

**Notes**:
- At least one of `left_context_column` or `right_context_column` must be specified
- Design primers early; for paired primers, design inner primer first, then outer with `paired_primer_column`
- `maximum_repeat_length` screens against `input_data` only; for genome-wide screening, build a background DB with `background()` and pass via `background_directory`
- Chained primer design: design one primer, then design its partner by passing `paired_primer_column`
- With `oligo_sets`: primers are designed per set with cross-set dimer screening; if `paired_primer_column` is provided, Tm matching is applied per set
- Patch mode with `oligo_sets`: existing per-set primers are reused; only missing sets trigger new design

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

### `motif`

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
    excluded_motifs=None,          # list | str | dict | pd.DataFrame | None
    background_directory=None,     # str | list | None
    random_seed=None,              # int | None
    verbose=True,                  # bool
)
```

**Required Parameters**

- `input_data` (str | DataFrame): CSV path or DataFrame with `ID` column + DNA sequence columns
- `oligo_length_limit` (int, ≥4): Maximum allowed oligo length (bp)
- `motif_sequence_constraint` (str): IUPAC constraint string (can be degenerate or constant)
- `maximum_repeat_length` (int, ≥4): Maximum shared repeat length between motif/anchor and oligos
- `motif_column` (str): Output column name

**Optional Parameters**

- `output_file` (str | None, default=None): Output CSV path (required for CLI)
- `motif_type` (int | str, default=0): `0`/`'variable'`=per-variant motifs (unique per row), `1`/`'constant'`=single motif/anchor shared by all rows. Also accepts: `'var'`, `'per-variant'`, `'non-constant'`, `'const'`, `'anchor'`, `'fixed'`
- `left_context_column` (str | None, default=None): Column for left DNA context
- `right_context_column` (str | None, default=None): Column for right DNA context
- `patch_mode` (bool, default=False): Fill only missing values; for `motif_type=1`, existing anchor (must be unique across all rows) is reused
- `excluded_motifs` (list | str | dict | DataFrame | None, default=None): Motifs to exclude. Accepts a list, CSV/FASTA path, DataFrame with `Exmotif` column, comma-string, or multiple sources. Strict ATGC only (no IUPAC codes, no dashes).
- `background_directory` (str | list | None, default=None): Background k-mer DB(s) from `background()` (screened against ALL specified DBs; junction-aware when context columns are provided)
- `random_seed` (int | None, default=None): RNG seed for reproducibility
- `verbose` (bool, default=True): Print progress output

**Returns**: `(DataFrame, stats_dict)`

**Notes**:
- At least one of `left_context_column` or `right_context_column` must be specified
- Design order: before barcodes if designing anchors for indexing
- Use `motif_type=1` to design constant anchors (e.g., barcode prefix/suffix for `index`)
- For anchors, tune `maximum_repeat_length` to control distinctiveness from context
- Constant bases in constraint may conflict with `excluded_motifs` and become impossible to solve
- Patch mode with `motif_type=1`: a compatible existing anchor is reused for new rows

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

### `spacer`

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
    excluded_motifs=None,          # list | str | dict | pd.DataFrame | None
    background_directory=None,     # str | list | None
    random_seed=None,              # int | None
    verbose=True,                  # bool
)
```

**Required Parameters**

- `input_data` (str | DataFrame): CSV path or DataFrame with `ID` column + DNA sequence columns
- `oligo_length_limit` (int, ≥4): Maximum allowed oligo length (bp)
- `maximum_repeat_length` (int, ≥4): Maximum shared repeat length between spacer and oligos
- `spacer_column` (str): Output column name

**Optional Parameters**

- `output_file` (str | None, default=None): Output CSV path (required for CLI)
- `spacer_length` (int | list | str | DataFrame | None, default=None): `None`=auto-fill to `oligo_length_limit`, `int`=fixed length, `list`=per-row lengths, DataFrame/CSV with `ID` + `Length`
- `left_context_column` (str | None, default=None): Column for left DNA context
- `right_context_column` (str | None, default=None): Column for right DNA context
- `patch_mode` (bool, default=False): Fill only missing values (`None`/`NaN`/empty/`'-'`)
- `excluded_motifs` (list | str | dict | DataFrame | None, default=None): Motifs to exclude. Accepts a list, CSV/FASTA path, DataFrame with `Exmotif` column, comma-string, or multiple sources. Strict ATGC only (no IUPAC codes, no dashes).
- `background_directory` (str | list | None, default=None): Background k-mer DB(s) from `background()` (screened against ALL specified DBs; junction-aware when context columns are provided)
- `random_seed` (int | None, default=None): RNG seed for reproducibility
- `verbose` (bool, default=True): Print progress output

**Returns**: `(DataFrame, stats_dict)`

**Notes**:
- At least one of `left_context_column` or `right_context_column` must be specified
- Design order: last, after all other elements (primers, motifs, barcodes)
- If `spacer_length=None`, spacers auto-fill remaining space up to `oligo_length_limit`

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

### `background`

**Purpose**: Build a k-mer database for screening designed sequences against off-target repeats.

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

- `input_data` (list | str | DataFrame): Background sequences: list of DNA strings, CSV/DataFrame with `Sequence` column, or FASTA path
- `maximum_repeat_length` (int, 6-20): Maximum repeat length to screen (k-mer size = `maximum_repeat_length + 1`)
- `output_directory` (str): Output directory basename (writes `<name>.oligopool.background`)

**Optional Parameters**

- `verbose` (bool, default=True): Print progress output

**Returns**: `stats_dict` (stats-only, no DataFrame)

**Notes**:
- Use via `background_directory` in `primer`/`barcode`/`motif`/`spacer` (and as a QC scan in `verify`); supports a single path or a list of paths
- `background()` writes a directory ending in `.oligopool.background`; pass that directory as `background_directory` (single path or list; junction-aware when context columns are provided)
- The DB stores all (k+1)-mers from input sequences where k = `maximum_repeat_length`

**CLI Equivalent**:
```bash
op background \
    --input-data genome.fasta \
    --maximum-repeat-length 15 \
    --output-directory ecoli_bg
```

[↑ Back to TOC](#table-of-contents)

---

### `merge`

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

- `input_data` (str | DataFrame): CSV path or DataFrame with `ID` column + DNA sequence columns
- `merge_column` (str): Output column name for merged region

**Optional Parameters**

- `output_file` (str | None, default=None): Output CSV path (required for CLI)
- `left_context_column` (str | None, default=None): First column to merge (defaults to first sequence column)
- `right_context_column` (str | None, default=None): Last column to merge (defaults to last sequence column)
- `verbose` (bool, default=True): Print progress output

**Returns**: `(DataFrame, stats_dict)` — source columns in range are removed

**Notes**:
- Use to collapse multiple adjacent columns into one before further design steps
- Useful for creating context columns or simplifying architecture mid-pipeline

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

### `revcomp`

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

- `input_data` (str | DataFrame): CSV path or DataFrame with `ID` column + DNA sequence columns

**Optional Parameters**

- `output_file` (str | None, default=None): Output CSV path (required for CLI)
- `left_context_column` (str | None, default=None): First column to revcomp (defaults to first sequence column)
- `right_context_column` (str | None, default=None): Last column to revcomp (defaults to last sequence column)
- `verbose` (bool, default=True): Print progress output

**Returns**: `(DataFrame, stats_dict)`

**Notes**:
- Reverses column order in addition to reverse-complementing sequences
- Use for strand-flipping a region mid-pipeline

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

### `final`

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

- `input_data` (str | DataFrame): CSV path or DataFrame with `ID` column + DNA sequence columns

**Optional Parameters**

- `output_file` (str | None, default=None): Output CSV path (required for CLI)
- `verbose` (bool, default=True): Print progress output

**Returns**: `(DataFrame, stats_dict)` — output contains `ID`, `CompleteOligo`, `OligoLength` only

**Notes**:
- Concatenates all sequence columns left-to-right; gap characters (`'-'`) are stripped
- Run `verify` after design to check for length, motif emergence, and background conflicts

**CLI Equivalent**:
```bash
op final \
    --input-data library.csv \
    --output-file synthesis_ready
```

[↑ Back to TOC](#table-of-contents)

---

## Assembly Mode

Assembly Mode provides tools for fragmenting long oligos that exceed synthesis length limits into overlapping pieces for assembly workflows.

### `split`

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

- `input_data` (str | DataFrame): CSV path or DataFrame with `ID` column + DNA sequence columns
- `split_length_limit` (int): Maximum fragment length (bp)
- `minimum_melting_temperature` (float): Minimum overlap Tm (°C)
- `minimum_hamming_distance` (int): Minimum pairwise Hamming distance between overlaps
- `minimum_overlap_length` (int): Minimum overlap length (bp)
- `maximum_overlap_length` (int): Maximum overlap length (bp)

**Optional Parameters**

- `output_file` (str | None, default=None): Output CSV path (required for CLI)
- `random_seed` (int | None, default=None): RNG seed for reproducibility
- `separate_outputs` (bool, default=False): Return list of per-split DataFrames; if `output_file` set, writes `{base}.SplitN.oligopool.split.csv` files
- `verbose` (bool, default=True): Print progress output

**Returns**:
- `(DataFrame, stats_dict)` when `separate_outputs=False` (default) — output contains `Split1`, `Split2`, ... columns
- `([DataFrame, ...], stats_dict)` when `separate_outputs` is enabled — one DataFrame per `SplitN` column

**Notes**:
- Number of fragments varies per oligo; even-numbered splits (`Split2`, `Split4`, ...) are reverse-complemented
- Original annotation columns are not preserved in output
- Treat each `SplitN` column as its own synthesis pool; run `pad` per `SplitN`, then `final` each padded DataFrame
- Do not run `final()` on raw multi-column `split` output

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

### `pad`

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

- `input_data` (str | DataFrame): Output from `split()` or any DataFrame with fragment column
- `oligo_length_limit` (int, ≥60): Maximum padded fragment length (bp)
- `split_column` (str): Which fragment column to pad (e.g., `Split1`)
- `typeIIS_system` (str): Type IIS enzyme name (see [supported list](#type-iis-enzymes))
- `minimum_melting_temperature` (float, ≥25): Minimum pad primer Tm (°C)
- `maximum_melting_temperature` (float, ≤95): Maximum pad primer Tm (°C)
- `maximum_repeat_length` (int, 6-20): Maximum shared repeat length

**Optional Parameters**

- `output_file` (str | None, default=None): Output CSV path (required for CLI)
- `random_seed` (int | None, default=None): RNG seed for reproducibility
- `verbose` (bool, default=True): Print progress output

**Returns**: `(DataFrame, stats_dict)` — output contains `5primeSpacer`, `ForwardPrimer`, `<split_column>`, `ReversePrimer`, `3primeSpacer`

**Notes**:
- The Type IIS recognition site must be absent from all split fragments (in either orientation); `pad` validates this and fails early if conflicts are found
- To prevent conflicts, add your Type IIS motif and its reverse complement (e.g., `GGTCTC` and `GAGACC` for BsaI) to `excluded_motifs` when designing upstream elements
- This only constrains newly designed elements—if your input sequences already contain the site, choose a different enzyme or redesign those sequences
- Run `pad` separately for each fragment (e.g., `Split1`, `Split2`, ...); you cannot pad all columns in one call
- If a fragment can't fit under `oligo_length_limit`, spacer(s) are set to `'-'`
- Post-synthesis workflow: PCR amplify → Type IIS digest → mung bean nuclease (skip for blunt-cutters like `MlyI`) → assemble via overlaps

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

## Degenerate Mode

Degenerate Mode compresses variant libraries with low mutational diversity into IUPAC-degenerate oligos for cost-efficient synthesis. Best for selection assays where enriched variants are identified by sequencing (no barcode-based readout required). Ideal for selection-based discovery: design → compress → synthesize → select → sequence → map back.

### `compress`

**Purpose**: Compress concrete DNA sequences into IUPAC-degenerate oligos for cost-efficient synthesis.

**Signature**:
```python
mapping_df, synthesis_df, stats = op.compress(
    # Required
    input_data,                    # str | pd.DataFrame

    # Optional
    mapping_file=None,             # str | None
    synthesis_file=None,           # str | None
    rollout_simulations=100,       # int
    rollout_horizon=4,             # int
    random_seed=None,              # int | None
    verbose=True,                  # bool
)
```

**Required Parameters**

- `input_data` (str | DataFrame): CSV path or DataFrame with `ID` column + DNA sequence columns. All non-ID columns are concatenated as the sequence. Must contain only A/T/G/C (no degenerate codes).

**Optional Parameters**

- `mapping_file` (str | None, default=None): Output CSV for variant-to-degenerate mapping
- `synthesis_file` (str | None, default=None): Output CSV for degenerate oligos ready for synthesis
- Output suffixes (if missing): `mapping_file` → `.oligopool.compress.mapping.csv`, `synthesis_file` → `.oligopool.compress.synthesis.csv`
- `rollout_simulations` (int, default=100): Monte Carlo simulations per decision (higher = better compression, slower)
- `rollout_horizon` (int, default=4): Lookahead positions for rollouts
- `random_seed` (int | None, default=None): RNG seed for reproducibility
- `verbose` (bool, default=True): Print progress output

**Returns**: `(mapping_df, synthesis_df, stats_dict)`
- `mapping_df`: Maps original variant IDs to degenerate oligo IDs (`ID`, `Sequence`, `DegenerateID`)
- `synthesis_df`: Degenerate oligos for ordering (`DegenerateID`, `DegenerateSeq`, `Degeneracy`, `OligoLength`)
- `stats_dict`: Includes `compression_ratio`, `min_degeneracy`, `max_degeneracy`, `mean_degeneracy`

**Notes**:
- Core guarantee (lossless): `degeneracy(prefix) <= count(compatible variants)` — no invented sequences
- Sequences of different lengths are compressed independently
- If sequences are too diverse to compress, returns 1:1 mapping
- Monte Carlo rollouts are parallelized across CPU cores automatically via numba

**CLI Equivalent**:
```bash
op compress \
    --input-data variants.csv \
    --mapping-file degenerate_library \
    --synthesis-file degenerate_library \
    --random-seed 42
```

[↑ Back to TOC](#table-of-contents)

---

### `expand`

**Purpose**: Expand IUPAC-degenerate sequences into all concrete A/T/G/C sequences.

**Signature**:
```python
df, stats = op.expand(
    # Required
    input_data,                    # str | pd.DataFrame
    sequence_column,               # str

    # Optional
    mapping_file=None,             # str | pd.DataFrame | None
    output_file=None,              # str | None
    expansion_limit=None,          # int | None
    verbose=True,                  # bool
)
```

**Required Parameters**

- `input_data` (str | DataFrame): CSV path or DataFrame with `ID` column (or `DegenerateID` from `compress` output)
- `sequence_column` (str): Column containing IUPAC-degenerate sequences to expand

**Optional Parameters**

- `mapping_file` (str | DataFrame | None, default=None): CSV path or DataFrame with `ID` and `DegenerateID` columns (the `mapping_df` from `compress`); restores original variant IDs in output
- `mapping_file` is auto-suffixed to `.oligopool.compress.mapping.csv` if missing (so a basename like `degenerate_library` works)
- `output_file` (str | None, default=None): Output CSV path
- `expansion_limit` (int | None, default=None): Safety cap for maximum total expanded sequences
- `verbose` (bool, default=True): Print progress output

**Returns**: `(DataFrame, stats_dict)` — output contains original columns plus `ExpandedSeq` and `OligoLength`; if `mapping_file` provided, includes both `ID` and `DegenerateID` columns

**Notes**:
- Primarily used as a verification tool to confirm `compress` output
- IUPAC codes: N=any, R=A/G, Y=C/T, S=C/G, W=A/T, K=G/T, M=A/C, B=C/G/T, D=A/G/T, H=A/C/T, V=A/C/G
- Expansion can be exponential; use `expansion_limit` for safety with highly degenerate sequences
- Expansion is parallelized across CPU cores when supported (falls back to serial on restricted systems)

**CLI Equivalent**:
```bash
op expand \
    --input-data synthesis.csv \
    --sequence-column DegenerateSeq \
    --mapping-file degenerate_library \
    --output-file expanded
```

[↑ Back to TOC](#table-of-contents)

---

## Analysis Mode

### `index`

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
    associate_prefix_gap=0,        # int
    associate_suffix_gap=0,        # int

    # Optional
    verbose=True,                  # bool
)
```

**Required Parameters**

- `barcode_data` (str | DataFrame): CSV path or DataFrame with `ID` + barcode column
- `barcode_column` (str): Barcode column name to index
- `index_file` (str): Output basename (writes `<name>.oligopool.index`)

**Optional Parameters (Barcode Anchors)**

- `barcode_prefix_column` (str | None, default=None): Column for constant prefix anchor
- `barcode_suffix_column` (str | None, default=None): Column for constant suffix anchor
- `barcode_prefix_gap` (int, default=0): Bases between prefix anchor and barcode in read
- `barcode_suffix_gap` (int, default=0): Bases between barcode and suffix anchor in read

**Optional Parameters (Associate Indexing)**

- `associate_data` (str | DataFrame | None, default=None): CSV path or DataFrame with associates (can be same as `barcode_data`)
- `associate_column` (str | None, default=None): Column for associate elements
- `associate_prefix_column` (str | None, default=None): Column for constant associate prefix anchor
- `associate_suffix_column` (str | None, default=None): Column for constant associate suffix anchor
- `associate_prefix_gap` (int, default=0): Bases between prefix anchor and associate in read
- `associate_suffix_gap` (int, default=0): Bases between associate and suffix anchor in read

**Optional Parameters**

- `verbose` (bool, default=True): Print progress output

**Returns**: `stats_dict` (stats-only, no DataFrame)

**Notes**:
- At least one of `barcode_prefix_column` or `barcode_suffix_column` is required
- Anchors must be constant (single-unique) sequences, ideally ≥6 bp and adjacent to indexed column
- If an anchor appears multiple times in a read, counting keeps the best-scoring placement; ties that yield different barcodes are rejected as `barcode_ambiguous`
- When using associates, at least one of `associate_prefix_column` or `associate_suffix_column` is required

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

### `pack`

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

- `r1_fastq_file` (str): R1 FastQ path (supports `.gz`)
- `r1_read_type` (int | str): `0`/`'forward'`=forward, `1`/`'reverse'`=reverse orientation. Also accepts: `'fwd'`, `'f'`, `'rev'`, `'r'`
- `pack_type` (int | str): `0`/`'concatenate'`=concatenate pairs, `1`/`'merge'`=merge/assemble overlapping pairs. Also accepts: `'concatenated'`, `'concat'`, `'cat'`, `'joined'`, `'join'`, `'merged'`, `'assemble'`, `'assembled'`, `'asm'`
- `pack_file` (str): Output basename (writes `<name>.oligopool.pack`)

**Optional Parameters (R1 Filters)**

- `minimum_r1_read_length` (int, default=1): Minimum R1 read length
- `minimum_r1_read_quality` (int, default=20): Minimum average R1 Phred score

**Optional Parameters (Paired-End)**

- `r2_fastq_file` (str | None, default=None): R2 FastQ path
- `r2_read_type` (int | str | None, default=None): R2 orientation: `0`/`'forward'` or `1`/`'reverse'`
- `minimum_r2_read_length` (int | None, default=None): Minimum R2 read length
- `minimum_r2_read_quality` (int | None, default=None): Minimum average R2 Phred score

**Optional Parameters (Performance)**

- `pack_size` (float, default=3.0): Million unique reads per pack (0.1-5.0)
- `core_count` (int, default=0): CPU cores (`0`=auto)
- `memory_limit` (float, default=0.0): GB per core (`0`=auto)

**Optional Parameters**

- `verbose` (bool, default=True): Print progress output

**Returns**: `stats_dict` (stats-only, no DataFrame)

**Notes**:
- `pack_type='concatenate'`: joins R1+R2 (use when reads don't overlap)
- `pack_type='merge'`: assembles overlapping R1+R2 into consensus (use when reads overlap)
- Pack files store reads as `(r1, r2)` tuples; merged/single-end reads use `r2=None`
- Deduplication is performed automatically

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

### `acount`

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
    failed_reads_file=None,        # str | None
    failed_reads_sample_size=1000, # int
    verbose=True,                  # bool
)
```

**Required Parameters**

- `index_file` (str): Index basename or path (reads `<name>.oligopool.index`)
- `pack_file` (str): Pack basename or path (reads `<name>.oligopool.pack`)
- `count_file` (str): Output basename (writes `<name>.oligopool.acount.csv`)

**Optional Parameters**

- `mapping_type` (int | str, default=0): `0`/`'fast'`=fast/near-exact, `1`/`'sensitive'`=slow/sensitive. Also accepts: `'quick'`, `'near-exact'`, `'sens'`, `'accurate'`, `'slow'`
- `barcode_errors` (int, default=-1): Max barcode errors (`-1`=auto from index)
- `associate_errors` (int, default=-1): Max associate errors (`-1`=auto from index)
- `callback` (callable | None, default=None): Custom read filter function (Python API only)
- `core_count` (int, default=0): CPU cores (`0`=auto)
- `memory_limit` (float, default=0.0): GB per core (`0`=auto)
- `failed_reads_file` (str | None, default=None): Output CSV path for failed read samples (`.oligopool.acount.failed_reads.csv` appended if missing). `None` disables sampling.
- `failed_reads_sample_size` (int, default=1000): Maximum samples per failure category (1-100,000)
- `verbose` (bool, default=True): Print progress output

**Returns**: `(counts_DataFrame, stats_dict)` — output contains `<indexname>.ID`, `BarcodeCounts`, `AssociationCounts`

**Notes**:
- Use `acount` when you need to verify barcode-variant coupling (requires associates in index)
- `mapping_type='sensitive'` is slower but catches more errors
- Multi-anchor reads: if an anchor appears multiple times, the best-scoring one is used; ties with multiple barcodes are rejected as `barcode_ambiguous`
- Failed reads sampling collects representative samples from each failure category for diagnostics:
  - `phix_match`: PhiX contamination detected
  - `low_complexity`: Low-complexity sequence (mono/di/trinucleotide)
  - `anchor_missing`: No anchor found in read
  - `barcode_absent`: Barcode not identified
  - `barcode_ambiguous`: Multiple anchors yield different barcodes with equal best confidence
  - `associate_prefix_missing`: Associate prefix/suffix constants not found
  - `associate_mismatch`: Associate does not match expected variant
  - `callback_false`: Callback function returned False
  - `incalculable`: Other uncategorized failures

**Callback Signature** (Python only):
```python
def callback(r1, r2, ID, count, coreid):
    """
    Args:
        r1: str - Read 1 sequence (always present)
        r2: str | None - Read 2 sequence (None for merged/single-end reads)
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
    --mapping-type 1 \
    --failed-reads-file failures \
    --failed-reads-sample-size 500
```

[↑ Back to TOC](#table-of-contents)

---

### `xcount`

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
    failed_reads_file=None,        # str | None
    failed_reads_sample_size=1000, # int
    verbose=True,                  # bool
)
```

**Required Parameters**

- `index_files` (str | list): Single index (string) or multiple (list) for combinatorial counting
- `pack_file` (str): Pack basename or path (reads `<name>.oligopool.pack`)
- `count_file` (str): Output basename (writes `<name>.oligopool.xcount.csv`)

**Optional Parameters**

- `mapping_type` (int | str, default=0): `0`/`'fast'`=fast/near-exact, `1`/`'sensitive'`=slow/sensitive. Also accepts: `'quick'`, `'near-exact'`, `'sens'`, `'accurate'`, `'slow'`
- `barcode_errors` (int, default=-1): Max barcode errors (`-1`=auto from index)
- `callback` (callable | None, default=None): Custom read filter function (Python API only)
- `core_count` (int, default=0): CPU cores (`0`=auto)
- `memory_limit` (float, default=0.0): GB per core (`0`=auto)
- `failed_reads_file` (str | None, default=None): Output CSV path for failed read samples (`.oligopool.xcount.failed_reads.csv` appended if missing). `None` disables sampling.
- `failed_reads_sample_size` (int, default=1000): Maximum samples per failure category (1-100,000)
- `verbose` (bool, default=True): Print progress output

**Returns**: `(counts_DataFrame, stats_dict)` — output contains one `<indexname>.ID` column per index, plus `CombinatorialCounts`. Missing barcodes shown as `'-'`

**Notes**:
- Use `xcount` for barcode-only counting without associate verification
- For combinatorial counting (BC1 × BC2), pass multiple indexes as a list
- Single index: counts each barcode; multiple indexes: counts barcode combinations
- Multi-anchor reads: if an anchor appears multiple times, the best-scoring one is used; ties with multiple barcodes are rejected as `barcode_ambiguous`
- Failed reads sampling collects representative samples from each failure category for diagnostics:
  - `phix_match`: PhiX contamination detected
  - `low_complexity`: Low-complexity sequence (mono/di/trinucleotide)
  - `anchor_missing`: No anchor found in read
  - `barcode_absent`: Barcode not identified
  - `barcode_ambiguous`: Multiple anchors yield different barcodes with equal best confidence
  - `callback_false`: Callback function returned False
  - `incalculable`: Other uncategorized failures

**CLI Equivalent**:
```bash
# Single barcode
op xcount \
    --index-files bc1_index \
    --pack-file sample \
    --count-file counts \
    --failed-reads-file failures

# Combinatorial (BC1 x BC2)
op xcount \
    --index-files bc1_index,bc2_index \
    --pack-file sample \
    --count-file combo_counts
```

[↑ Back to TOC](#table-of-contents)

---

## QC Mode

### `lenstat`

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

- `input_data` (str | DataFrame): CSV path or DataFrame with `ID` column + DNA sequence columns
- `oligo_length_limit` (int, ≥4): Reference length limit for free space calculation

**Optional Parameters**

- `verbose` (bool, default=True): Print progress output

**Returns**: `stats_dict` (stats-only, no DataFrame, no `output_file`)

**Notes**:
- Non-destructive; use to check remaining space before adding more elements
- Reports min/max/mean lengths and free space distribution

**CLI Equivalent**:
```bash
op lenstat \
    --input-data library.csv \
    --oligo-length-limit 200
```

[↑ Back to TOC](#table-of-contents)

---

### `verify`

**Purpose**: Verify oligo pool for length, motif emergence, and background k-mer conflicts.

**Signature**:
```python
df, stats = op.verify(
    # Required
    input_data,                    # str | pd.DataFrame
    oligo_length_limit,            # int

    # Optional
    output_file=None,              # str | None
    excluded_motifs=None,          # list | str | dict | pd.DataFrame | None
    background_directory=None,     # str | list | None
    verbose=True,                  # bool
)
```

**Required Parameters**

- `input_data` (str | DataFrame): CSV path or DataFrame with `ID` column and at least one DNA column (ATGC only)
- `oligo_length_limit` (int): Maximum allowed oligo length (>= 4)

**Optional Parameters**

- `output_file` (str | None, default=None): Output CSV path (required for CLI)
- `excluded_motifs` (list | str | dict | DataFrame | None, default=None): Motifs to check for emergence. Accepts a list, CSV/FASTA path, DataFrame with `Exmotif` column, comma-string, or multiple sources. Strict ATGC only (no IUPAC codes, no dashes).
- `background_directory` (str | list | None, default=None): Background k-mer DB path(s) from `background()` (checked against all specified DBs)
- `verbose` (bool, default=True): Print progress output

**Returns**: `(DataFrame, stats_dict)` — DataFrame contains per-row conflict flags and details

**Output DataFrame Columns**:
- `ID`: Original row identifier
- `CompleteOligo`: Concatenated oligo sequence
- `OligoLength`: Length of CompleteOligo
- `HasLengthConflict`: True if length exceeds `oligo_length_limit`
- `HasExmotifConflict`: True if motif emergence detected (False if `excluded_motifs=None`)
- `HasBackgroundConflict`: True if k-mer matches background DB (False if `background_directory=None`)
- `HasAnyConflicts`: OR of above three flags
- `LengthConflictDetails`, `ExmotifConflictDetails`, `BackgroundConflictDetails`: Dict or None with violation details

**Notes**:
- Uses `CompleteOligo` column if present; otherwise concatenates all pure ATGC columns left-to-right
- IUPAC/degenerate columns are skipped silently during DNA column detection
- **Emergence**: A motif has emergence when its count exceeds the library-wide minimum (baseline)
- Conflict details are Python dicts in the returned DataFrame; serialized as JSON strings in CSV output
- **Reading CSV back**: Parse JSON for all `*Details` columns (currently the `*ConflictDetails` columns; e.g., loop over columns ending in `Details` and apply `json.loads`).

**CLI Equivalent**:
```bash
# Single source (comma-separated motifs)
op verify \
    --input-data library.csv \
    --oligo-length-limit 200 \
    --excluded-motifs GAATTC,GGATCC \
    --output-file results

# Multiple sources (merged)
op verify \
    --input-data library.csv \
    --oligo-length-limit 200 \
    --excluded-motifs cutsites.csv homopolymers.csv \
    --output-file results
```

[↑ Back to TOC](#table-of-contents)

---

### `inspect`

**Purpose**: Inspect non-CSV artifacts produced by Oligopool Calculator (background DBs, index files, pack files).

**Signature**:
```python
stats = op.inspect(
    # Required
    target,                        # str

    # Optional
    kind='auto',                   # str
    verbose=True,                  # bool
)
```

**Required Parameters**

- `target` (str): Artifact path (background directory, `.oligopool.index`, or `.oligopool.pack`)

**Optional Parameters**

- `kind` (str, default=`auto`): `background`, `index`, `pack`, or `auto`
- `verbose` (bool, default=True): Print progress output

**Returns**: `stats_dict` (stats-only module) — summary is stored under `stats['vars']['meta']` and `stats['vars']['verdict']`

**Notes**:
- Read-only: does not repair artifacts and never unpickles unsafe objects
- Background directories are validated before opening to avoid accidentally creating a new DB
- Index summary loads `meta.map` if present; pack summary loads `packing.stat` if present
- `stats['status']=True` only when `stats['vars']['verdict']=='Valid'`

**CLI Equivalent**:
```bash
# Background directory
op inspect --target demo.oligopool.background --stats-json --quiet

# Index file
op inspect --target BC1.oligopool.index --stats-json --quiet

# Pack file
op inspect --target reads.oligopool.pack --stats-json --quiet
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
db.multicheck(seq_list, rna=True)  # Returns list

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

Why this list: `pad()` uses Type IIS for pad excision (scarless removal after optional blunting). Each supported system is modeled as a recognition motif plus a 3' cut offset into adjacent `N` bases (e.g., BsaI: `GGTCTC` + `N*5`). Enzymes that cut upstream (or require more complex cut models) are not included. The variety is intentional: if one recognition site conflicts with your library, you can often switch enzymes without changing the workflow.

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
| `.oligopool.background` | `background` | K-mer database directory |
| `.oligopool.barcode.csv` | `barcode` | Barcode-annotated library |
| `.oligopool.primer.csv` | `primer` | Primer-annotated library |
| `.oligopool.motif.csv` | `motif` | Motif/anchor-annotated library |
| `.oligopool.spacer.csv` | `spacer` | Spacer-annotated library |
| `.oligopool.split.csv` | `split` | Split-fragment library |
| `.oligopool.pad.csv` | `pad` | Padded fragment library |
| `.oligopool.merge.csv` | `merge` | Merged element library |
| `.oligopool.revcomp.csv` | `revcomp` | Reverse-complemented library |
| `.oligopool.final.csv` | `final` | Synthesis-ready library |
| `.oligopool.index` | `index` | Barcode index file |
| `.oligopool.pack` | `pack` | Packed reads file |
| `.oligopool.acount.csv` | `acount` | Association counts |
| `.oligopool.xcount.csv` | `xcount` | Combinatorial counts |

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

**Special Cases**:
| Python | CLI | Notes |
|--------|-----|-------|
| `typeIIS_system` | `--typeiis-system` | camelCase becomes lowercase |
| `separate_outputs` | `--separate-outputs` / `--no-separate-outputs` | Default inverted: Python=False, CLI=True |
| `index_file` | `--index-file` | Used by `index`, `acount` (singular) |
| `index_files` | `--index-files` | Used by `xcount` (plural, comma/space-separated) |

**CLI-Specific Flags**:
- `--config PATH`: Load YAML config values for a command (CLI flags override; see `docs/docs.md` "Config Files")
- `--stats-json`: Print stats dict as JSON to stdout
- `--stats-file PATH`: Write stats JSON to file
- `--quiet`: Suppress verbose output (sets `verbose=False`)

**CLI Workflow Runner**:
- `op pipeline --config PIPELINE.yaml` (and `--dry-run`): Execute multi-step workflows from a single YAML config (see `docs/docs.md` "Config Files")

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
