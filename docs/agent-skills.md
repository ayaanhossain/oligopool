# Oligopool Calculator - AI Agent Guide

This document is designed for AI assistants to understand and help users with oligopool library design and analysis.

## Operation Policy (Keep This Useful)

- **Source of truth**: `docs/docs.md` + runtime docs (`help(op)` / `help(op.<module>)`) + real CLI behavior.
- **This guide**: Agent-facing contracts, gotchas, and workflow scaffolding (Notebook / Script / CLI).
- **If anything disagrees**: Trust runtime behavior, then update `docs/docs.md` and this guide together.

**Quick facts (the ones that matter in practice):**
- **Entry points**: `op` and `oligopool` are equivalent CLIs.
- **CLI help model**: There is no `--help` flag. Use `op COMMAND` for command options and `op manual COMMAND` for docstrings.
- **Basenames**: Prefer basenames for all outputs (`--output-file`, `--index-file`, `--pack-file`, `--count-file`); suffixes are auto-appended as needed.
- **Return shapes**:
  - Design/transform: `(out_df, stats_dict)`
  - Stats-only: `stats_dict` (`background`, `lenstat`, `verify`, `index`, `pack`)
  - Counting: `(counts_df, stats_dict)` (`acount`, `xcount`)
- **ID handling**: Input requires a unique `ID`. CSV outputs include an explicit `ID` column (no pandas index column).
- **Patch Mode**: `patch_mode=True` / `--patch-mode` fills only missing values (None/NaN/empty/`'-'`) and never overwrites existing designs.
- **Cross-set barcodes**: `cross_barcode_columns` and `minimum_cross_distance` must be provided together.

## Package Overview

**Oligopool Calculator** is a Python library for designing and analyzing oligonucleotide pool (oligopool)
libraries used in massively parallel reporter assays (MPRAs) and similar high-throughput experiments.

**Core workflow:**
1. **Design Mode**: Build a library of oligos with barcodes, primers, motifs, and spacers
2. **Analysis Mode**: Count barcoded reads from NGS data to quantify variant activity

**Installation:**
```bash
pip install oligopool
```

**Import:**
```python
import oligopool as op
```

## Key Concepts

### DataFrame-Centric Design
- Most design/transform modules operate on pandas DataFrames and return `(out_df, stats_dict)`
- Input is a CSV path or DataFrame with an `ID` primary-key column (and DNA sequence columns)
- Some modules return only a `stats` dict (no DataFrame): `background`, `lenstat`, `verify`, `index`, `pack`
- Counting modules return `(counts_df, stats_dict)`: `acount`, `xcount`
- When a DataFrame is passed with `ID` already set as the index, the library preserves that on return;
  otherwise it returns a DataFrame with an explicit `ID` column
- Chain modules to build libraries iteratively (new columns) and use Patch Mode to extend pools (new rows)

### The Stats Dictionary
Every module returns a `stats` dict:
```python
stats = {
    'status': True/False,      # Success or failure
    'basis': 'solved',         # Short reason code (see common values below)
    'step': 1,                 # Pipeline step number
    'step_name': 'designing',  # Step name
    'vars': {...},             # Module-specific variables
    'warns': {...},            # Any warnings
    'module': 'barcode',       # Module name (standardized)
    'input_rows': 100,         # Input row count (standardized)
    'output_rows': 100,        # Output row count (standardized)
    'random_seed': None,       # Present for stochastic modules
}
```

**Common `basis` values:**
- `'solved'` - A valid solution was found (or a deterministic step completed successfully)
- `'complete'` - Early-success (e.g., Patch Mode with no missing values to fill)
- `'unsolved'` - Search finished without reaching the target; check `stats['vars']` for exhaustion flags
- `'infeasible'` - Constraints/inputs make the requested design impossible at a given step

**Important**:
- Invalid argument types/paths/columns typically raise `RuntimeError('Invalid Argument Input(s).')` and do
  not return `stats`. Algorithmic failures return `stats['status'] = False` with a descriptive `step_name`.

### The Oligo Architecture
A typical oligo has this structure:
```
[Primer1] [Barcode1] [Variant] [Barcode2] [Primer2] [Spacer]
```
- **Variants**: The core sequences being tested (user-provided)
- **Barcodes**: Unique identifiers for each variant (Hamming-distance separated)
- **Primers**: Amplification sites (thermodynamically optimized)
- **Motifs/Anchors**: Constant sequences for indexing or restriction sites
- **Spacers**: Neutral filler to reach target oligo length

### Context Columns and Edge Effects
Most modules require `left_context_column` and/or `right_context_column`.

**Edge Effect**: The emergence of an undesired sequence (excluded motif, repeat, etc.) at the fusion boundary
when inserting an element between or next to existing sequence contexts.
For example, inserting barcode `GAATT` next to context ending in `...G` creates `...GGAATT`, which contains
`GAATTC` (EcoRI site) spanning the junction.

Context columns prevent edge effects by:
- Checking that no excluded motifs emerge at element-context junctions
- Ensuring no repeats span the insertion boundaries
- At least one context column is required for proper edge-effect screening

### Patch Mode
Use `patch_mode=True` to extend existing pools:
- Only fills missing values (None/NaN/empty/'-')
- Preserves existing designs
- Useful for iterative pool expansion

---

## Working Modes (Notebook / Script / CLI)

### Notebook mode (Jupyter)
- Prefer the Python API (`import oligopool as op`) so you can pass/inspect DataFrames directly.
- Keep `verbose=True` while iterating; use `lenstat` mid-pipeline and `verify` before ordering.
- Use `random_seed` in stochastic modules (`barcode`, `primer`, `motif`, `spacer`, `split`, `pad`) for reproducible results.
- Read docs inline: `help(op)`, `help(op.barcode)`, etc.

### Scripting mode (Python)
- Treat invalid arguments as exceptions: modules raise `RuntimeError('Invalid Argument Input(s).')` during parsing.
- Treat algorithmic failures as data: check `stats['status']` and `stats['basis']` and branch accordingly.
- For file-based pipelines, pass `output_file=...` to write CSVs; outputs always include an explicit `ID` column
  and omit the pandas index column.
- Use Patch Mode (`patch_mode=True`) to extend pools without overwriting existing designs.

### CLI mode (`op` / `oligopool`)
- The CLI is best for reproducible, file-first pipelines (e.g., in bash scripts or workflows).
- DataFrame-producing commands require `--output-file`; `--quiet` suppresses verbose printing.
- `--stats-json` prints the stats dictionary as JSON to stdout; `--stats-file` writes stats JSON to a file.
- Callbacks for `acount`/`xcount` are Python-only; the CLI always runs with `callback=None`.
- Discover options and docs with `op COMMAND` and `op manual COMMAND`; print citation info with `op cite`;
  enable completion via `op complete --install`.

## Design Mode - Full API Reference

### barcode

**Purpose**: Generate Hamming-distance separated barcodes for unique variant identification.

```python
df, stats = op.barcode(
    # Required
    input_data,                    # str | pd.DataFrame - CSV path or DataFrame with ID + DNA columns
    oligo_length_limit,            # int - Maximum allowed oligo length (≥4)
    barcode_length,                # int - Length of designed barcodes (≥4)
    minimum_hamming_distance,      # int - Minimum pairwise Hamming distance (≥1)
    maximum_repeat_length,         # int - Max shared repeat length with oligos (≥4)
    barcode_column,                # str - Column name for designed barcodes

    # Optional
    output_file=None,              # str - Output CSV path (required for CLI, optional for library)
    barcode_type=0,                # int - 0=terminus-optimized (fast), 1=spectrum-optimized (thorough)
    left_context_column=None,      # str - Column for left DNA context
    right_context_column=None,     # str - Column for right DNA context
    patch_mode=False,              # bool - Fill only missing values in existing column
    cross_barcode_columns=None,    # str | list[str] - Existing barcode column(s) for cross-set separation
    minimum_cross_distance=None,   # int - Min Hamming distance to cross set (use with cross_barcode_columns)
    excluded_motifs=None,          # list | str | pd.DataFrame - Motifs to exclude (list, CSV, DataFrame, or FASTA)
    verbose=True,                  # bool - Print progress
    random_seed=None,              # int - RNG seed for reproducibility
)
```

**Design order**: After primers/motifs, before spacers.

**Tips**:
- Use `barcode_type=0` for large libraries (faster)
- Use `barcode_type=1` for maximum barcode diversity
- For multi-barcode designs, use `cross_barcode_columns` to ensure BC2 is separated from BC1
- Cross-set separation is strict: provide both `cross_barcode_columns` and `minimum_cross_distance`;
  all cross-set barcode values must already be A/T/G/C strings of length `barcode_length`.

---

### primer

**Purpose**: Design thermodynamically optimal primers for amplification.

```python
df, stats = op.primer(
    # Required
    input_data,                    # str | pd.DataFrame - CSV path or DataFrame
    oligo_length_limit,            # int - Maximum allowed oligo length (≥4)
    primer_sequence_constraint,    # str - IUPAC constraint (e.g., 'SS' + 'N'*18 for GC clamp)
    primer_type,                   # int - 0=forward, 1=reverse
    minimum_melting_temperature,   # float - Minimum Tm (≥25°C)
    maximum_melting_temperature,   # float - Maximum Tm (≤95°C)
    maximum_repeat_length,         # int - Max shared repeat length (≥6)
    primer_column,                 # str - Column name for designed primer

    # Optional
    output_file=None,              # str - Output CSV path (required for CLI, optional for library)
    left_context_column=None,      # str - Column for left DNA context
    right_context_column=None,     # str - Column for right DNA context
    patch_mode=False,              # bool - Fill only missing values
    oligo_sets=None,               # list | str | pd.DataFrame - Per-oligo set labels for set-specific primers
    paired_primer_column=None,     # str - Column of paired primer for Tm matching
    excluded_motifs=None,          # list | str | pd.DataFrame - Motifs to exclude (list, CSV, DataFrame, or FASTA)
    background_directory=None,     # str - Background k-mer database for off-target screening
    verbose=True,                  # bool - Print progress
    random_seed=None,              # int - RNG seed
)
```

**Design order**: Design primers early. For paired primers: design inner primer first, then outer with `paired_primer_column`.

**Tips**:
- Use `'SS' + 'N'*18` for 5' GC clamp
- Use `'N'*18 + 'WW'` for 3' AT clamp
- Build background with `op.background()` before primer design for off-target screening
- For multiplexed libraries, use `oligo_sets` to design set-specific primers

---

### motif

**Purpose**: Add sequence motifs or constant anchors.

```python
df, stats = op.motif(
    # Required
    input_data,                    # str | pd.DataFrame - CSV path or DataFrame
    oligo_length_limit,            # int - Maximum allowed oligo length (≥4)
    motif_sequence_constraint,     # str - IUPAC pattern or constant sequence
    maximum_repeat_length,         # int - Max shared repeat length (≥4)
    motif_column,                  # str - Column name for designed motif

    # Optional
    output_file=None,              # str - Output CSV path (required for CLI, optional for library)
    motif_type=0,                  # int - 0=per-variant motifs, 1=constant anchor for all
    left_context_column=None,      # str - Column for left DNA context
    right_context_column=None,     # str - Column for right DNA context
    patch_mode=False,              # bool - Fill only missing values
    excluded_motifs=None,          # list | str | pd.DataFrame - Motifs to exclude (list, CSV, DataFrame, or FASTA)
    verbose=True,                  # bool - Print progress
    random_seed=None,              # int - RNG seed
)
```

**Use cases**:
- Restriction sites: `motif_sequence_constraint='GAATTC'` (EcoRI)
- Barcode anchors: `motif_type=1` with `'N'*10` for constant anchor
- Degenerate regions: `'NNNGGATCCNNN'` (BamHI with flanking Ns)

**Design order**: Before barcodes if designing anchors for indexing.

---

### spacer

**Purpose**: Add neutral filler DNA to reach target oligo length.

```python
df, stats = op.spacer(
    # Required
    input_data,                    # str | pd.DataFrame - CSV path or DataFrame
    oligo_length_limit,            # int - Maximum allowed oligo length (≥4)
    maximum_repeat_length,         # int - Max shared repeat length (≥4)
    spacer_column,                 # str - Column name for designed spacer

    # Optional
    output_file=None,              # str - Output CSV path (required for CLI, optional for library)
    spacer_length=None,            # int | list | str | pd.DataFrame - Fixed, per-variant, or None for auto-fill
    left_context_column=None,      # str - Column for left DNA context
    right_context_column=None,     # str - Column for right DNA context
    patch_mode=False,              # bool - Fill only missing values
    excluded_motifs=None,          # list | str | pd.DataFrame - Motifs to exclude (list, CSV, DataFrame, or FASTA)
    verbose=True,                  # bool - Print progress
    random_seed=None,              # int - RNG seed
)
```

**spacer_length options**:
- `None`: Auto-fill to reach `oligo_length_limit`
- `int`: Fixed length for all oligos
- `list`: Per-variant lengths aligned to input
- `DataFrame`: With 'ID' and 'Length' columns

**Design order**: Last, after all other elements.

---

### background

**Purpose**: Build k-mer database for primer off-target screening.

```python
stats = op.background(
    # Required
    input_data,                    # list | str | pd.DataFrame - Sequences (list, CSV, DataFrame, or FASTA)
    maximum_repeat_length,         # int - k-mer size for screening (6-20)
    output_directory,              # str - Output directory path

    # Optional
    verbose=True,                  # bool - Print progress
)
```

**input_data formats**:
- List of DNA strings: `['ATGC...', 'GCTA...']`
- CSV file: With 'Sequence' column
- DataFrame: With 'Sequence' column
- FASTA file: `.fa`, `.fasta`, `.fna` (optionally gzipped)

**When to use**: Before primer design when screening against host genome/plasmid.

---

### split

**Purpose**: Break long oligos into overlapping fragments for assembly.

```python
df, stats = op.split(
    # Required
    input_data,                    # str | pd.DataFrame - CSV path or DataFrame
    split_length_limit,            # int - Maximum fragment length
    minimum_melting_temperature,   # float - Assembly overlap Tm
    minimum_hamming_distance,      # int - Minimum distance between overlaps
    minimum_overlap_length,        # int - Minimum overlap length
    maximum_overlap_length,        # int - Maximum overlap length

    # Optional
    output_file=None,              # str - Output CSV path (required for CLI, optional for library)
    verbose=True,                  # bool - Print progress
    random_seed=None,              # int - RNG seed for reproducibility
)
```

**When to use**: When oligos exceed synthesis length limits (typically >200 bp).

---

### pad

**Purpose**: Add amplification primers with Type IIS sites for assembly.

```python
df, stats = op.pad(
    # Required
    input_data,                    # str | pd.DataFrame - DataFrame from split()
    oligo_length_limit,            # int - Maximum padded oligo length (≥ 60)
    split_column,                  # str - Which fragment column to pad
    typeIIS_system,                # str - Restriction enzyme (e.g., 'BsaI', 'BbsI', 'BsmBI')
    minimum_melting_temperature,   # float - Pad primer minimum Tm (≥ 25°C)
    maximum_melting_temperature,   # float - Pad primer maximum Tm (≤ 95°C)
    maximum_repeat_length,         # int - Max repeat length (between 6 and 20)

    # Optional
    output_file=None,              # str - Output CSV path (required for CLI, optional for library)
    verbose=True,                  # bool - Print progress
    random_seed=None,              # int - RNG seed for reproducibility
)
```

**When to use**: After split, for each fragment.

**Supported Type IIS enzymes (34 total):**
`AcuI`, `AlwI`, `BbsI`, `BccI`, `BceAI`, `BciVI`, `BcoDI`, `BmrI`, `BpuEI`, `BsaI`, `BseRI`, `BsmAI`,
`BsmBI`, `BsmFI`, `BsmI`, `BspCNI`, `BspQI`, `BsrDI`, `BsrI`, `BtgZI`, `BtsCI`, `BtsI`, `BtsIMutI`,
`EarI`, `EciI`, `Esp3I`, `FauI`, `HgaI`, `HphI`, `HpyAV`, `MlyI`, `MnlI`, `SapI`, `SfaNI`

For unsupported enzymes, design primers/sites manually with `primer` or `motif`.

---

### merge

**Purpose**: Collapse multiple contiguous columns into one.

```python
df, stats = op.merge(
    # Required
    input_data,                    # str | pd.DataFrame - CSV path or DataFrame
    merge_column,                  # str - Name for merged column

    # Optional
    output_file=None,              # str - Output CSV path (required for CLI, optional for library)
    left_context_column=None,      # str - First column to merge (defaults to first column if None)
    right_context_column=None,     # str - Last column to merge (defaults to last column if None)
    verbose=True,                  # bool - Print progress
)
```

**When to use**: Combine multiple elements into single column before `revcomp` or `final`.

---

### revcomp

**Purpose**: Reverse complement a column range.

```python
df, stats = op.revcomp(
    # Required
    input_data,                    # str | pd.DataFrame - CSV path or DataFrame

    # Optional
    output_file=None,              # str - Output CSV path (required for CLI, optional for library)
    left_context_column=None,      # str - First column to revcomp (defaults to first column if None)
    right_context_column=None,     # str - Last column to revcomp (defaults to last column if None)
    verbose=True,                  # bool - Print progress
)
```

**When to use**: When you need reverse complement of designed elements.

---

### lenstat

**Purpose**: Check length statistics and remaining space (non-destructive).

```python
stats = op.lenstat(
    # Required
    input_data,                    # str | pd.DataFrame - CSV path or DataFrame

    oligo_length_limit,            # int - Reference length limit for free space calculation (≥ 4)
    verbose=True,                  # bool - Print progress
)
```

**When to use**: Frequently during design to monitor remaining space.

---

### verify

**Purpose**: QC check before synthesis.

```python
stats = op.verify(
    # Required
    input_data,                    # str | pd.DataFrame - CSV path or DataFrame

    # Optional
    oligo_length_limit=None,       # int - Check for length violations
    excluded_motifs=None,          # list | str | pd.DataFrame - Check for motif violations (list, CSV, DataFrame, or FASTA)
    verbose=True,                  # bool - Print progress
)
```

**When to use**: Before `final()` to catch constraint violations.

---

### final

**Purpose**: Concatenate columns into synthesis-ready oligos.

```python
df, stats = op.final(
    # Required
    input_data,                    # str | pd.DataFrame - CSV path or DataFrame

    # Optional
    output_file=None,              # str - Output CSV path (required for CLI, optional for library)
    verbose=True,                  # bool - Print progress
)
```

**Output**: DataFrame with `CompleteOligo` and `OligoLength` (plus an explicit `ID` column unless the caller
provided `ID` as the DataFrame index).

**When to use**: Final step before synthesis.

---

## Analysis Mode - Full API Reference

### index

**Purpose**: Build barcode index for counting.

```python
stats = op.index(
    # Required
    barcode_data,                  # str | pd.DataFrame - CSV path or DataFrame with barcodes
    barcode_column,                # str - Column containing barcodes
    index_file,                    # str - Output index filename (creates .oligopool.index)

    # Optional
    barcode_prefix_column=None,    # str - Column for constant prefix anchor
    barcode_suffix_column=None,    # str - Column for constant suffix anchor
    barcode_prefix_gap=0,          # int - Distance between prefix and barcode
    barcode_suffix_gap=0,          # int - Distance between suffix and barcode
    associate_data=None,           # str | pd.DataFrame | None - CSV path or DataFrame with associates
    associate_column=None,         # str | None - Column for associate elements
    associate_prefix_column=None,  # str - Column for associate prefix
    associate_suffix_column=None,  # str - Column for associate suffix
    verbose=True,                  # bool - Print progress
)
```

**Tips**:
- Prefix/suffix anchors must be constant (single unique sequence, ≥6 bp) and adjacent to the indexed column.
- For `acount`, specify `associate_data`/`associate_column` and at least one constant adjacent associate prefix/suffix.
- Design anchors with `motif(motif_type=1, ...)` (or use universal primers if they are constant).

---

### pack

**Purpose**: Preprocess and deduplicate FASTQ reads.

```python
stats = op.pack(
    # Required
    r1_fastq_file,                 # str - Path to R1 FASTQ (supports .gz)
    r1_read_type,                  # int - 0=forward, 1=reverse orientation
    pack_type,                     # int - 0=concatenate, 1=merge paired reads
    pack_file,                     # str - Output pack filename (creates .oligopool.pack)

    # Optional
    minimum_r1_read_length=1,      # int - Minimum R1 read length
    minimum_r1_read_quality=20,    # int - Minimum average R1 Phred score
    r2_fastq_file=None,            # str - Path to R2 FASTQ
    r2_read_type=None,             # int - R2 orientation
    minimum_r2_read_length=None,   # int - Minimum R2 read length
    minimum_r2_read_quality=None,  # int - Minimum average R2 Phred score
    pack_size=3.0,                 # float - Million unique reads per pack (0.1-5.0)
    core_count=0,                  # int - CPU cores (0=auto)
    memory_limit=0.0,              # float - GB per core (0=auto)
    verbose=True,                  # bool - Print progress
)
```

**Tips**:
- For single-end reads, use R1 arguments only
- For paired-end merging (`pack_type=1`), provide both `r2_fastq_file` and `r2_read_type`
- For pre-merged reads (e.g., from FLASH), use as single-end

---

### acount

**Purpose**: Association counting - verify barcode-variant coupling.

```python
counts_df, stats = op.acount(
    # Required
    index_file,                    # str - Index filename (from index())
    pack_file,                     # str - Pack filename (from pack())
    count_file,                    # str - Output count filename (creates .oligopool.acount.csv)

    # Optional
    mapping_type=0,                # int - 0=fast, 1=sensitive
    barcode_errors=-1,             # int - Max barcode errors (-1=auto)
    associate_errors=-1,           # int - Max associate errors (-1=auto)
    callback=None,                 # callable - Custom read processing function (Python API only)
    core_count=0,                  # int - CPU cores (0=auto)
    memory_limit=0.0,              # float - GB per core (0=auto)
    verbose=True,                  # bool - Print progress
)
```

**Output columns**: one `<indexname>.ID` column per index (typically one), plus `BarcodeCounts`, `AssociationCounts`

**When to use**:
- Validating synthesis accuracy
- Input library QC
- Post-assay barcode-variant verification

---

### xcount

**Purpose**: Barcode-only counting (single or combinatorial).

```python
counts_df, stats = op.xcount(
    # Required
    index_files,                   # str | list - Single index or list of index filenames
    pack_file,                     # str - Pack filename
    count_file,                    # str - Output count filename (creates .oligopool.xcount.csv)

    # Optional
    mapping_type=0,                # int - 0=fast, 1=sensitive
    barcode_errors=-1,             # int - Max barcode errors (-1=auto)
    callback=None,                 # callable - Custom read processing function (Python API only)
    core_count=0,                  # int - CPU cores (0=auto)
    memory_limit=0.0,              # float - GB per core (0=auto)
    verbose=True,                  # bool - Print progress
)
```

**Output**: one `<indexname>.ID` column per index, plus `CombinatorialCounts`. Missing barcodes shown as `'-'`.

**When to use**:
- Pure barcode counting (single index)
- Multi-barcode combinations (BC1 × BC2)
- Cleaved vs uncleaved quantification (ribozymes, etc.)

---

### Callback Functions (Python API only)

For `acount` and `xcount`, custom callbacks enable additional analysis:

```python
def my_callback(read, ID, count, coreid):
    """
    Args:
        read: str - The processed read sequence
        ID: tuple - Identified barcode IDs (e.g., ('BC1_v1', 'BC2_v1') or ('BC1_v1', '-'))
        count: int - Read/ID frequency in current pack
        coreid: int - CPU core ID processing this read

    Returns:
        bool - True to accept read for counting, False to reject
    """
    # Custom logic here
    return True

counts_df, stats = op.xcount(..., callback=my_callback)
```

**Use cases**:
- Filter reads by additional criteria
- Extract cleavage sites
- Record custom metrics with `multiprocessing.Manager().dict()`

---

## Advanced Modules

### vectorDB

**Purpose**: LevelDB-based scalable on-disk k-mer storage for background databases.

```python
db = op.vectorDB(
    path,                          # str - Path to store/load vectorDB instance
    maximum_repeat_length,         # int - Max shared repeat length (k-mer size = maximum_repeat_length + 1)
)
```

**Methods:**
```python
# Add sequence k-mers to database
db.add(seq, rna=True)              # seq: DNA string, rna: convert U→T if True
db.multiadd(seq_list, rna=True)    # add many sequences

# Check if any k-mer from seq exists in database
seq in db                          # Returns bool (uses __contains__)
db.multicheck(seq_list, rna=True)  # list[bool]

# Iterate over stored k-mers
for kmer in db:
    print(kmer)

# Get number of stored k-mers
len(db)

# Remove k-mers
db.remove(seq, rna=True)
db.multiremove(seq_list, rna=True, clear=False)

# Clear all k-mers (optionally set new maxreplen)
db.clear(maxreplen=None)

# Close vs drop
db.close()                         # close handle, keep files
db.drop()
```

**When to use**: Direct manipulation of background k-mer databases, custom screening pipelines.

**Note**: When reopening an existing vectorDB, `maximum_repeat_length` is ignored and loaded from the instance.

---

### Scry

**Purpose**: 1-nearest-neighbor barcode classifier used internally by `acount`/`xcount`.

```python
model = op.Scry()

# Train the model
model.train(
    X,                             # iterable - Sequences to fit
    Y,                             # iterable - Labels for sequences
    n,                             # int - Length of each sequence
    k,                             # int - k-mer length for sequences
    t,                             # int - Maximum errors in sequence variants
    liner=None,                    # coroutine - For dynamic printing (optional)
)

# Prime for prediction
model.prime(
    t=0,                           # int - Number of errors to consider (default: 0)
    mode=0,                        # int - 0=fast/near-exact, 1=slow/sensitive (default: 0)
)

# Predict label
label, score = model.predict(x)    # returns (label, score) or (None, None)
```

**When to use**: Building custom counting pipelines, debugging barcode classification, custom analysis workflows.

**Example:**
```python
model = op.Scry()
model.train(
    X=['ATGCATGC', 'GCTAGCTA', 'TTAATTAA'],
    Y=['bc1', 'bc2', 'bc3'],
    n=8,
    k=4,
    t=1
)
model.prime(t=1, mode=0)
label, score = model.predict('ATGCATGC')  # Returns ('bc1', 1.0) for exact/near-exact hits
```

---

## Common Workflows

### Basic Library Design
```python
import oligopool as op
import pandas as pd

# 1. Start with variants
df = pd.DataFrame({
    'ID': ['v1', 'v2', 'v3'],
    'Variant': ['ATGCATGC...', 'GCTAGCTA...', 'TTAACCGG...']
})

# 2. Build background (optional, for primer screening)
op.background(
    input_data='plasmid.fasta',
    maximum_repeat_length=12,
    output_directory='plasmid_bg'
)

# 3. Add primers (inner first if paired)
df, _ = op.primer(
    input_data=df,
    oligo_length_limit=200,
    primer_sequence_constraint='SS' + 'N'*18,
    primer_type=0,
    minimum_melting_temperature=53,
    maximum_melting_temperature=55,
    maximum_repeat_length=10,
    primer_column='Primer1',
    left_context_column='Variant',
    background_directory='plasmid_bg',
    excluded_motifs=['GAATTC', 'GGATCC'],
)

df, _ = op.primer(
    input_data=df,
    oligo_length_limit=200,
    primer_sequence_constraint='N'*18 + 'WW',
    primer_type=1,
    minimum_melting_temperature=53,
    maximum_melting_temperature=55,
    maximum_repeat_length=10,
    primer_column='Primer2',
    paired_primer_column='Primer1',
    right_context_column='Variant',
    background_directory='plasmid_bg',
    excluded_motifs=['GAATTC', 'GGATCC'],
)

# 4. Add barcodes
df, _ = op.barcode(
    input_data=df,
    oligo_length_limit=200,
    barcode_length=12,
    minimum_hamming_distance=3,
    maximum_repeat_length=6,
    barcode_column='BC1',
    left_context_column='Primer1',
    right_context_column='Variant',
    excluded_motifs=['GAATTC', 'GGATCC'],
)

# 5. Check length
op.lenstat(input_data=df, oligo_length_limit=200)

# 6. Add spacers to reach target length
df, _ = op.spacer(
    input_data=df,
    oligo_length_limit=200,
    maximum_repeat_length=6,
    spacer_column='Spacer',
    spacer_length=None,  # Auto-fill
    left_context_column='Primer2',
    excluded_motifs=['GAATTC', 'GGATCC'],
)

# 7. QC and finalize
op.verify(input_data=df, oligo_length_limit=200, excluded_motifs=['GAATTC', 'GGATCC'])
final_df, _ = op.final(input_data=df)
final_df.to_csv('library_for_synthesis.csv', index=False)
```

### Multi-Barcode Design (BC1 × BC2)
```python
# Design BC1
df, _ = op.barcode(
    input_data=df,
    barcode_length=11,
    minimum_hamming_distance=3,
    barcode_column='BC1',
    left_context_column='Primer1',
    right_context_column='Variant',
    ...
)

# Design BC2 with cross-set separation from BC1
df, _ = op.barcode(
    input_data=df,
    barcode_length=11,
    minimum_hamming_distance=3,
    barcode_column='BC2',
    left_context_column='Variant',
    right_context_column='Primer2',
    cross_barcode_columns=['BC1'],
    minimum_cross_distance=3,
    ...
)
```

### Analysis Pipeline
```python
# 1. Build index files
#
# IMPORTANT: `index()` requires at least one constant (library-wide) prefix/suffix anchor
# for barcodes (and for associates, if used). These anchors must be:
#   - adjacent to the indexed column in the DataFrame
#   - a single unique sequence across the library (≥ 6 bp)
#
# Example architecture for association counting:
#   Primer1 | BC1 | Variant | Primer2
op.index(
    barcode_data='library.csv',
    barcode_column='BC1',
    barcode_prefix_column='Primer1',
    associate_data='library.csv',
    associate_column='Variant',
    associate_suffix_column='Primer2',
    index_file='bc1_assoc_idx',
)

# Example architecture for combinatorial barcode counting:
#   Primer1 | BC1 | ...payload... | BC2 | Primer2
op.index(
    barcode_data='library.csv',
    barcode_column='BC1',
    barcode_prefix_column='Primer1',
    index_file='bc1_idx',
)

op.index(
    barcode_data='library.csv',
    barcode_column='BC2',
    barcode_suffix_column='Primer2',
    index_file='bc2_idx',
)

# 2. Pack reads
op.pack(
    r1_fastq_file='reads_R1.fq.gz',
    r2_fastq_file='reads_R2.fq.gz',
    r1_read_type=0,
    r2_read_type=1,
    pack_type=1,  # Merge
    minimum_r1_read_quality=30,
    minimum_r2_read_quality=30,
    pack_file='reads',
)

# 3a. Association counting (verify barcode-variant coupling)
ac_df, _ = op.acount(
    index_file='bc1_assoc_idx',
    pack_file='reads',
    count_file='association_counts',
    mapping_type=1,
)

# 3b. Combinatorial counting (BC1 × BC2)
xc_df, _ = op.xcount(
    index_files=['bc1_idx', 'bc2_idx'],
    pack_file='reads',
    count_file='combo_counts',
    mapping_type=1,
)
```

### Extending an Existing Pool (Patch Mode)
```python
# Load existing library
df = pd.read_csv('existing_library.csv')

# Add new variants
new_variants = pd.DataFrame({
    'ID': ['v101', 'v102'],
    'Variant': ['ATGC...', 'GCTA...'],
    'BC1': [None, None],  # Missing - will be designed
    'Primer1': [df['Primer1'].iloc[0], df['Primer1'].iloc[0]],  # Reuse existing
})
df = pd.concat([df, new_variants], ignore_index=True)

# Patch mode fills only missing values
df, _ = op.barcode(
    input_data=df,
    barcode_column='BC1',
    patch_mode=True,  # Only design for v101, v102
    ...
)
```

---

## Troubleshooting Guide

### "Design failed" / status=False
1. Check `stats['basis']` for failure reason
2. Common fixes:
   - Increase `barcode_length` or decrease `minimum_hamming_distance`
   - Relax `excluded_motifs`
   - Increase `maximum_repeat_length`
   - Widen Tm range for primers
   - Try `barcode_type=0` (faster, may find solutions faster)

### Oligo too long
- Use `lenstat` to check remaining space
- Reduce element lengths
- Use `split` + `pad` for assembly

### No barcodes mapping in analysis
- Verify anchor sequences in `index` match designed anchors exactly
- Check `barcode_prefix/suffix_gap` settings
- Use `mapping_type=1` for sensitive mode
- Verify read orientation (`r1_read_type`, `r2_read_type`)

### Cross-contamination in counts
- Increase `minimum_hamming_distance`
- Add `cross_barcode_columns` for multi-barcode designs
- Check for index hopping in sequencing

### Memory issues
- Reduce `pack_size` in `pack()`
- Use `core_count` to limit parallelism
- Set `memory_limit` to constrain per-core memory

---

## Parameter Quick Reference

### Excluded Motifs Formats
All element modules accept `excluded_motifs`:
- **List**: `['GAATTC', 'GGATCC', 'AAAAA']`
- **CSV**: File with 'Exmotif' column
- **DataFrame**: With 'Exmotif' column
- **FASTA**: `.fa`, `.fasta`, `.fna` (optionally gzipped)

### Context Columns
- Always specify at least one of `left_context_column` or `right_context_column`
- Use adjacent columns in the oligo architecture
- Critical for preventing edge effects (undesired motifs/repeats emerging at fusion junctions)

### IUPAC Codes for Constraints
```
A = A          W = A or T
T = T          S = G or C
G = G          M = A or C
C = C          K = G or T
N = any        R = A or G
               Y = C or T
               B = C, G, or T
               D = A, G, or T
               H = A, C, or T
               V = A, C, or G
```

### Common Restriction Sites
```
EcoRI  = GAATTC    BamHI  = GGATCC
HindIII = AAGCTT   XbaI   = TCTAGA
BsaI   = GGTCTC    BbsI   = GAAGAC
BsmBI  = CGTCTC    AatII  = GACGTC
```

---

## File Formats

### Input CSV
```csv
ID,Variant
v1,ATGCATGCATGCATGC
v2,GCTAGCTAGCTAGCTA
v3,TTAACCGGTTAACCGG
```

### Output Files
- `.oligopool.background` - k-mer database directory
- `.oligopool.primer.csv` - primer-annotated library CSV
- `.oligopool.barcode.csv` - barcode-annotated library CSV
- `.oligopool.motif.csv` - motif/anchor-annotated library CSV
- `.oligopool.spacer.csv` - spacer-annotated library CSV
- `.oligopool.split.csv` - split-fragment library CSV
- `.oligopool.pad.csv` - padded fragment library CSV
- `.oligopool.merge.csv` - merged element library CSV
- `.oligopool.revcomp.csv` - reverse-complemented library CSV
- `.oligopool.final.csv` - synthesis-ready concatenated library CSV
- `.oligopool.index` - barcode index file
- `.oligopool.pack` - packed reads file
- `.oligopool.acount.csv` - association counts
- `.oligopool.xcount.csv` - combinatorial counts

---

## CLI Usage

Every module is available via command line:

```bash
# Get help
op                           # List all commands
op barcode                   # Command-specific options (run without args)
op manual barcode            # Detailed documentation
op cite                      # Citation information

# Design
op barcode \
    --input-data variants.csv \
    --oligo-length-limit 200 \
    --barcode-length 12 \
    --minimum-hamming-distance 3 \
    --maximum-repeat-length 6 \
    --barcode-column BC1 \
    --left-context-column Variant \
    --output-file with_barcodes

# Analysis
op index \
    --barcode-data library.csv \
    --barcode-column BC1 \
    --barcode-prefix-column Primer1 \
    --index-file bc1_idx

op pack \
    --r1-fastq-file reads_R1.fq.gz \
    --r2-fastq-file reads_R2.fq.gz \
    --r1-read-type 0 \
    --r2-read-type 1 \
    --pack-type 1 \
    --pack-file reads

op xcount \
    --index-files bc1_idx \
    --pack-file reads \
    --count-file counts

# Tab completion (recommended)
op complete --install
```

### CLI-Specific Notes
- DataFrame-producing commands require `--output-file` (e.g., `barcode`, `primer`, `motif`, `spacer`, `split`,
  `pad`, `merge`, `revcomp`, `final`). Stats-only commands do not (`background`, `lenstat`, `verify`, `index`, `pack`).
- Most output paths auto-append a suffix if missing (e.g., `.oligopool.barcode.csv`), so prefer using a basename like
  `--output-file my_library` to avoid doubled extensions.
- Patch Mode is enabled via `--patch-mode` on element modules (fills missing values without overwriting existing ones).
- Sequence constraints: `--primer-sequence-constraint "'N'*20"` or `--primer-sequence-constraint GCC+N*20+CCG`
- Callbacks are Python-only (not available in CLI)
- Exit codes: `0` success, `1` runtime error, `404` CLI argument parsing error (e.g., missing required args).

---

## Best Practices

1. **Design order**: background → primers → motifs/anchors → barcodes → spacers → verify → final
2. **Use lenstat frequently**: Check remaining space after each element
3. **Save intermediate DataFrames**: Enable rollback if design fails
4. **Set random_seed**: For reproducible designs
5. **Run verify before synthesis**: Catch constraint violations early
6. **Use Patch Mode for extensions**: Don't redesign existing elements
7. **Use cross_barcode_columns**: For multi-barcode designs to ensure separation
8. **Build background first**: If you have host/plasmid sequences to screen against

---

## Decision Guide for Common Questions

**"Which module do I use to...?"**
- Add unique identifiers → `barcode`
- Add amplification sites → `primer`
- Add restriction sites or anchors → `motif`
- Fill to target length → `spacer`
- Screen against genome/plasmid → `background` then `primer` with `background_directory`
- Handle long oligos → `split` + `pad`
- Count reads → `index` + `pack` + `acount` or `xcount`

**"acount or xcount?"**
- Need to verify barcode-variant pairing → `acount`
- Just counting barcodes (no variant verification) → `xcount`
- Combinatorial counting (BC1 × BC2) → `xcount`

**"Design is failing, what do I relax?"**
- Barcodes: increase `barcode_length`, decrease `minimum_hamming_distance`, use `barcode_type=0`
- Primers: widen Tm range, more `N`s in constraint, increase `maximum_repeat_length`
- General: reduce `excluded_motifs`, increase `maximum_repeat_length`

**"What order should I design elements?"**
1. `background` (if screening needed)
2. `primer` (inner primers first if paired)
3. `motif` (anchors for indexing)
4. `barcode`
5. `spacer`
6. `verify` → `final`

**"How do I extend an existing pool?"**
- Use `patch_mode=True` - only fills missing values, preserves existing designs

**"Input requirements?"**
- DataFrame must have unique `ID` column
- For element-design modules (`barcode`, `primer`, `motif`, `spacer`), all non-ID columns must be non-empty DNA strings
  (A/T/G/C) except the target output column in Patch Mode, which may contain missing values (None/NaN/empty/`'-'`).
- `verify` is more permissive and will summarize metadata columns and degenerate/IUPAC columns rather than failing.
- At least one context column required for element design

---

## Links

- Repository: https://github.com/ayaanhossain/oligopool
- Documentation: `docs/docs.md` in repository
- Paper: https://doi.org/10.1021/acssynbio.4c00661
- CLI help: `op manual <command>`
