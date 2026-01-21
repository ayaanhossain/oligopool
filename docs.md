<h1 align="center">
    <a href="https://github.com/ayaanhossain/oligopool/">
        <img src="https://raw.githubusercontent.com/ayaanhossain/repfmt/main/oligopool/img/logo.svg" alt="Oligopool Calculator" width="460"/>
    </a>
</h1>

<p align="center"><i>Your friendly neighborhood DNA library designer</i></p>

---

# Documentation

Welcome to the Oligopool Calculator docs! Whether you're designing your first barcode library or optimizing a million-variant MPRA, you're in the right place.

**TL;DR**: Design oligopool libraries with `barcode`, `primer`, `motif`, `spacer`. Analyze sequencing data with `index`, `pack`, `acount`/`xcount`. Most modules take CSV/DataFrames in and return `(out_df, stats)` (a few return stats only). Chain them together. Ship it.

---

## Table of Contents

- [Quick Start](#quick-start)
- [Core Concepts](#core-concepts)
- [Design Mode](#design-mode)
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
- [Analysis Mode](#analysis-mode)
  - [index](#index) - Build barcode index
  - [pack](#pack) - Preprocess FastQ
  - [acount](#acount) - Association counting
  - [xcount](#xcount) - Combinatorial counting
- [Workflows](#workflows)
- [CLI Reference](#cli-reference)
- [Tips & Tricks](#tips--tricks)
- [Advanced Modules](#advanced-modules)
- [Troubleshooting](#troubleshooting)

---

## Quick Start

[↑ Back to TOC](#table-of-contents)

```python
import pandas as pd
import oligopool as op

# Start with your variants
df = pd.DataFrame({
    'ID': ['V1', 'V2', 'V3'],
    'Promoter': ['ATGCATGCATGCATGC', 'GCTAGCTAGCTAGCTA', 'TTAATTAATTAATTAA']
})

# Add a barcode (yes, it's that easy)
df, stats = op.barcode(
    input_data=df,
    oligo_length_limit=200,
    barcode_length=12,
    minimum_hamming_distance=3,
    maximum_repeat_length=8,
    barcode_column='BC1',
    left_context_column='Promoter'
)

# Check your work
print(df)
print(stats['status'])  # True = success!
```

That's it. You just designed unique barcodes with guaranteed Hamming distance. No PhD required (though it helps with the biology part).

---

## Core Concepts

[↑ Back to TOC](#table-of-contents)

### The DataFrame Flow

Most modules follow one of these patterns:

```python
out_df, stats = op.barcode(input_data=df, ...)
stats = op.verify(input_data=df, ...)
```

- **Input**: CSV path or pandas DataFrame with an `ID` column
- **Output**:
  - **Design/transform modules** return `(out_df, stats)`
  - **Stats-only modules** return `stats` (`background`, `lenstat`, `verify`, `index`, `pack`)
  - **Counting modules** return `(counts_df, stats)` (`acount`, `xcount`)
- **ID handling**: If you pass a DataFrame with `ID` as the index, outputs preserve it; files written to disk always include an explicit `ID` column (no pandas index column).
- **Chainable**: Output of one module feeds into the next

### The Stats Dictionary

Every module returns a `stats` dict with:
- `status`: `True` (success) or `False` (failed)
- `basis`: Why it succeeded/failed
- `step`: Which pipeline step finished
- `vars`: Useful variables from the run
- `warns`: Any warnings you should know about

### Patch Mode: The Secret Weapon

Extending an existing library? Don't redesign everything:

```python
df, stats = op.barcode(..., patch_mode=True)
```

Patch mode fills only missing values (`None`, `NaN`, `'-'`, empty). Your existing designs stay untouched. Your sanity stays intact.

### Context Columns and Edge Effects

Most design modules need to know what's next to the element being designed:

```python
df, stats = op.barcode(
    ...,
    left_context_column='Promoter',   # What's on the left?
    right_context_column='Gene'       # What's on the right?
)
```

**Edge effects** occur when an undesired sequence (excluded motif, repeat) emerges at the fusion boundary when inserting an element. For example, inserting barcode `GAATT` next to context ending in `...G` creates `...GGAATT`, which contains `GAATTC` (EcoRI site) spanning the junction. Context columns let the algorithm check and prevent these.

### Reproducibility

All stochastic design modules support `random_seed` for reproducible results:

```python
df, stats = op.barcode(..., random_seed=42)
```

Same seed + same inputs = same outputs. Great for debugging and publications.

### Verbose Mode

Control output verbosity with `verbose` (default: `True`):

```python
df, stats = op.barcode(..., verbose=False)  # Silent mode
```

---

## Design Mode

[↑ Back to TOC](#table-of-contents)

Design Mode is where you build your library piece by piece. Think of it like molecular LEGO, but with more constraints and fewer stepping-on-brick injuries.

### barcode

[↑ Back to TOC](#table-of-contents)

**What it does**: Designs unique barcodes with guaranteed minimum Hamming distance between all pairs.

**When to use it**: You need to identify individual variants in a pooled sequencing experiment.

```python
df, stats = op.barcode(
    input_data=df,
    oligo_length_limit=200,           # Max oligo length
    barcode_length=12,                # Barcode size (bp)
    minimum_hamming_distance=3,       # Min mismatches between any two barcodes
    maximum_repeat_length=8,          # Avoid repeats longer than this
    barcode_column='BC1',             # Output column name
    left_context_column='Primer5',    # Adjacent left element
    right_context_column='Gene',      # Adjacent right element
    barcode_type=0,                   # 0=fast (terminus), 1=thorough (spectrum)
    excluded_motifs=['GAATTC'],       # Avoid these sequences (e.g., EcoRI)
)
```

**Excluded motifs** (restriction sites, repetitive sequences, etc.):
```python
# Simple list - applied globally to all variants
df, stats = op.barcode(..., excluded_motifs=['GAATTC', 'GGATCC'])

# Via DataFrame - just needs an 'Exmotif' column
exmotifs_df = pd.DataFrame({
    'Exmotif': ['GAATTC', 'GGATCC', 'AAGCTT']  # All excluded from ALL variants
})
df, stats = op.barcode(..., excluded_motifs=exmotifs_df)  # Or 'exmotifs.csv'
```

> **Note**: Excluded motifs are applied globally to all variants (no per-variant exclusion). The DataFrame/CSV only needs an `Exmotif` column, and you can also pass a FASTA path (each record treated as an excluded motif).

**Pro tips**:
- Start with `barcode_type=0` (fast). Switch to `1` only if you're packing barcodes tight.
- Higher Hamming distance = more error tolerance, fewer possible barcodes.
- Use `patch_mode=True` to extend a pool without overwriting existing barcodes.
- Use `cross_barcode_columns` + `minimum_cross_distance` for multi-barcode designs (BC1 vs BC2).
- `excluded_motifs` works the same way in `primer`, `motif`, and `spacer` modules.

**Cross-set separation** (for multiplexed designs):
```python
# Design BC2 that's distinct from BC1
df, stats = op.barcode(
    ...,
    barcode_column='BC2',
    cross_barcode_columns=['BC1'],
    minimum_cross_distance=2,
)
```

<details>
<summary><b>Full parameter reference (Python API)</b></summary>

```python
df, stats = op.barcode(
    # Required
    input_data,                # str | pd.DataFrame
    oligo_length_limit,        # int
    barcode_length,            # int
    minimum_hamming_distance,  # int
    maximum_repeat_length,     # int
    barcode_column,            # str

    # Optional
    output_file=None,          # str | None
    barcode_type=0,            # int (0 terminus, 1 spectrum)
    left_context_column=None,  # str | None
    right_context_column=None, # str | None
    patch_mode=False,          # bool
    cross_barcode_columns=None,# str | list[str] | None
    minimum_cross_distance=None,# int | None
    excluded_motifs=None,      # list | str | pd.DataFrame | None
    verbose=True,              # bool
    random_seed=None,          # int | None
)
```

**Required parameters**
- `input_data`: CSV path or DataFrame with `ID` + DNA sequence columns.
- `oligo_length_limit`: Maximum allowed oligo length (bp).
- `barcode_length`: Barcode length (bp).
- `minimum_hamming_distance`: Minimum pairwise Hamming distance within the newly designed set.
- `maximum_repeat_length`: Maximum shared repeat length between the barcode and its context/oligos.
- `barcode_column`: Column name to create/overwrite (or to patch-fill when `patch_mode=True`).

**Optional parameters**
- `output_file`: Write output CSV (library mode only; required in CLI).
- `barcode_type`: `0` (terminus optimized, fast) or `1` (spectrum optimized, slow).
- `left_context_column` / `right_context_column`: Flanking context for edge-effect screening (at least one required).
- `patch_mode`: If `True`, only fills missing values in an existing `barcode_column` (no overwrites).
- `cross_barcode_columns` + `minimum_cross_distance`: Global cross-set separation against existing barcode columns (must be set together).
- `excluded_motifs`: List, CSV, DataFrame, or FASTA of motifs to exclude.
- `verbose`: Print progress output.
- `random_seed`: Reproducible randomness for stochastic steps.

</details>

---

### primer

[↑ Back to TOC](#table-of-contents)

**What it does**: Designs primers with thermodynamic constraints (Tm), repeat screening, and optional background filtering.

**When to use it**: You need amplification primers, sequencing primers, or assembly overlap regions.

```python
df, stats = op.primer(
    input_data=df,
    oligo_length_limit=200,
    primer_sequence_constraint='NNNNNNNNNNNNNNNNNNNN',  # 20N = any 20-mer
    primer_type=0,                    # 0=forward, 1=reverse
    minimum_melting_temperature=52.0,
    maximum_melting_temperature=58.0,
    maximum_repeat_length=10,
    primer_column='FwdPrimer',
    right_context_column='Gene',
)
```

**Paired primer design** (Tm matching within 1°C):
```python
# First, design the forward primer
df, stats = op.primer(..., primer_column='FwdPrimer', primer_type=0)

# Then design reverse primer matched to forward
df, stats = op.primer(
    ...,
    primer_column='RevPrimer',
    primer_type=1,
    paired_primer_column='FwdPrimer',  # Match Tm to this
)
```

**Per-set primers** (for multiplexed pools):
```python
# Via list (aligned to input rows)
df, stats = op.primer(
    ...,
    oligo_sets=['SetA', 'SetA', 'SetB', 'SetB'],  # Group labels
)

# Via DataFrame - IDs must match input DataFrame!
sets_df = pd.DataFrame({
    'ID': ['V1', 'V2', 'V3', 'V4'],               # Must match input_data IDs
    'OligoSet': ['SetA', 'SetA', 'SetB', 'SetB']
})
df, stats = op.primer(
    ...,
    oligo_sets=sets_df,                          # Or 'sets.csv'
)
```

**Background screening**:
```python
# First build a background database
op.background(
    input_data='genome.csv',
    maximum_repeat_length=15,
    output_directory='my_background'
)

# Then use it
df, stats = op.primer(
    ...,
    background_directory='my_background',
)
```

<details>
<summary><b>Full parameter reference (Python API)</b></summary>

```python
df, stats = op.primer(
    # Required
    input_data,                   # str | pd.DataFrame
    oligo_length_limit,           # int
    primer_sequence_constraint,   # str
    primer_type,                  # int (0 forward, 1 reverse)
    minimum_melting_temperature,  # float
    maximum_melting_temperature,  # float
    maximum_repeat_length,        # int
    primer_column,                # str

    # Optional
    output_file=None,             # str | None
    left_context_column=None,     # str | None
    right_context_column=None,    # str | None
    patch_mode=False,             # bool
    oligo_sets=None,              # list | str | pd.DataFrame | None
    paired_primer_column=None,    # str | None
    excluded_motifs=None,         # list | str | pd.DataFrame | None
    background_directory=None,    # str | None
    verbose=True,                 # bool
    random_seed=None,             # int | None
)
```

**Required parameters**
- `input_data`: CSV path or DataFrame with `ID` + DNA sequence columns.
- `oligo_length_limit`: Maximum allowed oligo length (bp).
- `primer_sequence_constraint`: IUPAC constraint string (also accepts simple Python-style shorthand like `'N'*20` in Python).
- `primer_type`: `0` forward or `1` reverse primer design.
- `minimum_melting_temperature` / `maximum_melting_temperature`: Allowed primer Tm range (°C).
- `maximum_repeat_length`: Maximum shared repeat length between primer and oligos (and background, if provided).
- `primer_column`: Output column name.

**Optional parameters**
- `output_file`: Write output CSV (library mode only; required in CLI).
- `left_context_column` / `right_context_column`: Flanking context for edge-effect screening (at least one required).
- `patch_mode`: If `True`, only fills missing values in an existing `primer_column` (no overwrites).
- `oligo_sets`: Per-row grouping labels (list aligned to `input_data`, or CSV/DataFrame with `ID` + `OligoSet`) to design set-specific primers with cross-set compatibility screening.
- `paired_primer_column`: Column containing the paired primer sequence for Tm matching (supports per-set pairing when `oligo_sets` is used).
- `excluded_motifs`: List, CSV, DataFrame, or FASTA of motifs to exclude.
- `background_directory`: Background k-mer DB created by `background(...)` for off-target screening.
- `verbose`: Print progress output.
- `random_seed`: Reproducible randomness for stochastic steps.

</details>

---

### motif

[↑ Back to TOC](#table-of-contents)

**What it does**: Inserts sequence motifs (per-variant or constant anchors) with constraint satisfaction.

**When to use it**:
- Constant anchors for indexing (barcode prefix/suffix)
- Regulatory elements with sequence constraints
- Any fixed or semi-fixed sequence insertion

```python
# Per-variant motifs (different for each row)
df, stats = op.motif(
    input_data=df,
    oligo_length_limit=200,
    motif_sequence_constraint='NNNNNNNN',  # 8N = any 8-mer
    maximum_repeat_length=6,
    motif_column='Spacer',
    motif_type=0,                          # 0=per-variant
    left_context_column='BC1',
)

# Constant anchor (same for all rows - great for indexing!)
df, stats = op.motif(
    ...,
    motif_type=1,                          # 1=constant anchor
    motif_column='BC1_Prefix',
)
```

**Why anchors matter**: When you index barcodes later, you need constant flanking sequences to locate them in reads. Design anchors with `motif_type=1` *before* designing barcodes.

<details>
<summary><b>Full parameter reference (Python API)</b></summary>

```python
df, stats = op.motif(
    # Required
    input_data,                 # str | pd.DataFrame
    oligo_length_limit,         # int
    motif_sequence_constraint,  # str
    maximum_repeat_length,      # int
    motif_column,               # str

    # Optional
    output_file=None,           # str | None
    motif_type=0,               # int (0 per-variant, 1 constant)
    left_context_column=None,   # str | None
    right_context_column=None,  # str | None
    patch_mode=False,           # bool
    excluded_motifs=None,       # list | str | pd.DataFrame | None
    verbose=True,               # bool
    random_seed=None,           # int | None
)
```

**Required parameters**
- `input_data`: CSV path or DataFrame with `ID` + DNA sequence columns.
- `oligo_length_limit`: Maximum allowed oligo length (bp).
- `motif_sequence_constraint`: IUPAC constraint string (or a constant string).
- `maximum_repeat_length`: Maximum shared repeat length between motif/anchor and surrounding oligos.
- `motif_column`: Output column name.

**Optional parameters**
- `output_file`: Write output CSV (library mode only; required in CLI).
- `motif_type`: `0` per-variant motifs, `1` constant motif/anchor reused across all rows.
- `left_context_column` / `right_context_column`: Flanking context for edge-effect screening (at least one required).
- `patch_mode`: If `True`, only fills missing values in an existing `motif_column` (no overwrites).
- `excluded_motifs`: List, CSV, DataFrame, or FASTA of motifs to exclude.
- `verbose`: Print progress output.
- `random_seed`: Reproducible randomness for stochastic steps.

</details>

---

### spacer

[↑ Back to TOC](#table-of-contents)

**What it does**: Inserts neutral DNA spacers to meet length requirements.

**When to use it**: Your oligos need to hit a specific length (e.g., synthesis array constraints).

```python
# Fixed length spacer
df, stats = op.spacer(
    input_data=df,
    oligo_length_limit=200,
    maximum_repeat_length=8,
    spacer_column='Spacer3',
    spacer_length=15,                      # Fixed 15 bp
    left_context_column='Gene',
)

# Auto-fill to oligo_length_limit
df, stats = op.spacer(
    ...,
    spacer_length=None,                    # Auto-size to fill remaining space
)

# Per-variant lengths (list)
df, stats = op.spacer(
    ...,
    spacer_length=[10, 15, 12, 8],         # Different per row
)

# Per-variant lengths (DataFrame) - IDs must match input DataFrame!
lengths_df = pd.DataFrame({
    'ID': ['V1', 'V2', 'V3', 'V4'],         # Must match input_data IDs
    'Length': [10, 15, 12, 8]
})
df, stats = op.spacer(
    ...,
    spacer_length=lengths_df,              # Or 'lengths.csv'
)
```

<details>
<summary><b>Full parameter reference (Python API)</b></summary>

```python
df, stats = op.spacer(
    # Required
    input_data,              # str | pd.DataFrame
    oligo_length_limit,      # int
    maximum_repeat_length,   # int
    spacer_column,           # str

    # Optional
    output_file=None,        # str | None
    spacer_length=None,      # int | list | str | pd.DataFrame | None
    left_context_column=None,# str | None
    right_context_column=None,# str | None
    patch_mode=False,        # bool
    excluded_motifs=None,    # list | str | pd.DataFrame | None
    verbose=True,            # bool
    random_seed=None,        # int | None
)
```

**Required parameters**
- `input_data`: CSV path or DataFrame with `ID` + DNA sequence columns.
- `oligo_length_limit`: Maximum allowed oligo length (bp).
- `maximum_repeat_length`: Maximum shared repeat length between spacer and surrounding oligos.
- `spacer_column`: Output column name.

**Optional parameters**
- `output_file`: Write output CSV (library mode only; required in CLI).
- `spacer_length`:
  - `None`: auto-size spacer per row to reach `oligo_length_limit`
  - `int`: fixed length for all rows
  - `list`: per-row lengths aligned to `input_data`
  - CSV/DataFrame: must contain `ID` + `Length`
- `left_context_column` / `right_context_column`: Flanking context for edge-effect screening (at least one required).
- `patch_mode`: If `True`, only fills missing values in an existing `spacer_column` (no overwrites).
- `excluded_motifs`: List, CSV, DataFrame, or FASTA of motifs to exclude.
- `verbose`: Print progress output.
- `random_seed`: Reproducible randomness for stochastic steps.

</details>

---

### background

[↑ Back to TOC](#table-of-contents)

**What it does**: Builds a k-mer database for screening primers against off-target sequences.

**When to use it**: You want primers that won't bind to host genome, plasmid backbone, etc.

```python
# From a FASTA file (genome, plasmid, etc.)
stats = op.background(
    input_data='ecoli_genome.fasta',       # Also supports .fa, .fna, .gz
    maximum_repeat_length=15,              # Screen k-mers up to this length
    output_directory='ecoli_bg',           # Creates ecoli_bg.oligopool.background
)

# Or from a CSV/DataFrame with 'Sequence' column
stats = op.background(
    input_data='plasmid_sequences.csv',
    maximum_repeat_length=15,
    output_directory='plasmid_bg',
)

# Now use it in primer design
df, stats = op.primer(
    ...,
    background_directory='ecoli_bg',
)
```

<details>
<summary><b>Full parameter reference (Python API)</b></summary>

```python
stats = op.background(
    input_data,             # list | str | pd.DataFrame
    maximum_repeat_length,  # int
    output_directory,       # str
    verbose=True,           # bool
)
```

**Required parameters**
- `input_data`: Background sequences (list of DNA strings, FASTA path, or CSV/DataFrame with a `Sequence` column).
- `maximum_repeat_length`: Maximum repeat length to screen against the background (bp).
- `output_directory`: Output directory basename (writes `<name>.oligopool.background`).

**Optional parameters**
- `verbose`: Print progress output.

</details>

---

### split

[↑ Back to TOC](#table-of-contents)

**What it does**: Breaks long oligos into overlapping fragments for assembly (Gibson, Golden Gate, etc.).

**When to use it**: Your oligos exceed synthesis length limits (~200 bp).

```python
split_df, stats = op.split(
    input_data=df,
    split_length_limit=150,               # Max fragment length
    minimum_melting_temperature=55.0,     # Overlap Tm
    minimum_hamming_distance=3,           # Overlap uniqueness
    minimum_overlap_length=20,
    maximum_overlap_length=30,
    output_file='split_library',
)
```

Output contains `Split1`, `Split2`, ... columns with fragments ready for assembly.

<details>
<summary><b>Full parameter reference (Python API)</b></summary>

```python
split_df, stats = op.split(
    # Required
    input_data,                 # str | pd.DataFrame
    split_length_limit,         # int
    minimum_melting_temperature,# float
    minimum_hamming_distance,   # int
    minimum_overlap_length,     # int
    maximum_overlap_length,     # int

    # Optional
    output_file=None,           # str | None
    verbose=True,               # bool
    random_seed=None,           # int | None
)
```

**Required parameters**
- `input_data`: CSV path or DataFrame with `ID` + DNA sequence columns.
- `split_length_limit`: Maximum fragment length (bp).
- `minimum_melting_temperature`: Minimum overlap Tm (°C).
- `minimum_hamming_distance`: Minimum pairwise distance between overlap sequences.
- `minimum_overlap_length` / `maximum_overlap_length`: Allowed overlap length range (bp).

**Optional parameters**
- `output_file`: Write output CSV (library mode only; required in CLI).
- `verbose`: Print progress output.
- `random_seed`: Reproducible randomness for stochastic steps.

</details>

---

### pad

[↑ Back to TOC](#table-of-contents)

**What it does**: Adds primers and Type IIS restriction sites to split fragments for scarless assembly.

**When to use it**: After `split`, to make fragments synthesis-ready and assembly-compatible.

```python
pad_df, stats = op.pad(
    input_data=split_df,
    split_column='Split1',                # Which fragment to pad (from split output)
    typeIIS_system='BsaI',                # Enzyme for Golden Gate
    oligo_length_limit=200,
    minimum_melting_temperature=52.0,
    maximum_melting_temperature=58.0,
    maximum_repeat_length=10,
    output_file='padded_split1',
)
```

**Supported Type IIS enzymes (34 total):**
`AcuI`, `AlwI`, `BbsI`, `BccI`, `BceAI`, `BciVI`, `BcoDI`, `BmrI`, `BpuEI`, `BsaI`, `BseRI`, `BsmAI`, `BsmBI`, `BsmFI`, `BsmI`, `BspCNI`, `BspQI`, `BsrDI`, `BsrI`, `BtgZI`, `BtsCI`, `BtsI`, `BtsIMutI`, `EarI`, `EciI`, `Esp3I`, `FauI`, `HgaI`, `HphI`, `HpyAV`, `MlyI`, `MnlI`, `SapI`, `SfaNI`

For other enzymes, design primers manually with the appropriate recognition/cut sites using `primer` or `motif`.

<details>
<summary><b>Full parameter reference (Python API)</b></summary>

```python
pad_df, stats = op.pad(
    # Required
    input_data,                   # str | pd.DataFrame
    oligo_length_limit,           # int
    split_column,                 # str
    typeIIS_system,               # str
    minimum_melting_temperature,  # float
    maximum_melting_temperature,  # float
    maximum_repeat_length,        # int

    # Optional
    output_file=None,             # str | None
    verbose=True,                 # bool
    random_seed=None,             # int | None
)
```

**Required parameters**
- `input_data`: Output from `split` (or any DataFrame with a split fragment column).
- `oligo_length_limit`: Maximum padded fragment length (bp).
- `split_column`: Which fragment column to pad (e.g., `Split1`).
- `typeIIS_system`: Name of the Type IIS enzyme system (see supported list above).
- `minimum_melting_temperature` / `maximum_melting_temperature`: Padding primer Tm range (°C).
- `maximum_repeat_length`: Maximum shared repeat length between primers/spacers and oligos.

**Optional parameters**
- `output_file`: Write output CSV (library mode only; required in CLI).
- `verbose`: Print progress output.
- `random_seed`: Reproducible randomness for stochastic steps.

</details>

---

### merge

[↑ Back to TOC](#table-of-contents)

**What it does**: Concatenates contiguous columns into a single column.

**When to use it**: Mid-pipeline cleanup, combining elements before further processing.

```python
df, stats = op.merge(
    input_data=df,
    merge_column='5prime_region',         # Output column name
    left_context_column='Primer5',
    right_context_column='BC1',           # Merges Primer5 through BC1 (inclusive)
)
```

<details>
<summary><b>Full parameter reference (Python API)</b></summary>

```python
df, stats = op.merge(
    # Required
    input_data,              # str | pd.DataFrame
    merge_column,            # str

    # Optional
    output_file=None,        # str | None
    left_context_column=None,# str | None
    right_context_column=None,# str | None
    verbose=True,            # bool
)
```

**Required parameters**
- `input_data`: CSV path or DataFrame with `ID` + DNA sequence columns.
- `merge_column`: Output column name for the merged region.

**Optional parameters**
- `output_file`: Write output CSV (library mode only; required in CLI).
- `left_context_column` / `right_context_column`: Define the inclusive column range to merge. If omitted, merges from the first to the last sequence column.
- `verbose`: Print progress output.

</details>

---

### revcomp

[↑ Back to TOC](#table-of-contents)

**What it does**: Reverse complements a range of columns and reverses their order.

**When to use it**: Switching strand orientation (e.g., synthesis vs. sequencing direction).

```python
df, stats = op.revcomp(
    input_data=df,
    left_context_column='Gene',
    right_context_column='Primer3',       # RevComp everything from Gene to Primer3 (inclusive)
)
```

<details>
<summary><b>Full parameter reference (Python API)</b></summary>

```python
df, stats = op.revcomp(
    # Required
    input_data,               # str | pd.DataFrame

    # Optional
    output_file=None,         # str | None
    left_context_column=None, # str | None
    right_context_column=None,# str | None
    verbose=True,             # bool
)
```

**Required parameters**
- `input_data`: CSV path or DataFrame with `ID` + DNA sequence columns.

**Optional parameters**
- `output_file`: Write output CSV (library mode only; required in CLI).
- `left_context_column` / `right_context_column`: Define the inclusive column range to reverse-complement and reorder. If omitted, uses the first and last sequence columns.
- `verbose`: Print progress output.

</details>

---

### lenstat

[↑ Back to TOC](#table-of-contents)

**What it does**: Reports length statistics and free space remaining.

**When to use it**: Check progress during design, ensure you're within synthesis limits.

```python
stats = op.lenstat(
    input_data=df,
    oligo_length_limit=200,
)
# Prints a nice table of per-column and total lengths
```

<details>
<summary><b>Full parameter reference (Python API)</b></summary>

```python
stats = op.lenstat(
    input_data,        # str | pd.DataFrame
    oligo_length_limit,# int
    verbose=True,      # bool
)
```

**Required parameters**
- `input_data`: CSV path or DataFrame with `ID` + DNA sequence columns.
- `oligo_length_limit`: Maximum allowed oligo length (bp).

**Optional parameters**
- `verbose`: Print progress output.

</details>

---

### verify

[↑ Back to TOC](#table-of-contents)

**What it does**: QC check for constraints, architecture, motifs, degeneracy.

**When to use it**: Final sanity check before ordering synthesis.

```python
stats = op.verify(
    input_data=df,
    oligo_length_limit=200,
    excluded_motifs=['GAATTC', 'GGATCC'],  # Restriction sites to flag
)
```

Checks:
- Length constraints
- Excluded motif presence
- Degenerate bases (non-ATGC)
- Column architecture

<details>
<summary><b>Full parameter reference (Python API)</b></summary>

```python
stats = op.verify(
    input_data,            # str | pd.DataFrame
    oligo_length_limit=None,# int | None
    excluded_motifs=None,  # list | str | pd.DataFrame | None
    verbose=True,          # bool
)
```

**Required parameters**
- `input_data`: CSV path or DataFrame with an `ID` column (can include metadata and degenerate/IUPAC columns).

**Optional parameters**
- `oligo_length_limit`: If provided, checks for oligo length overflow.
- `excluded_motifs`: Motifs to scan/report (list, CSV/DataFrame with `Exmotif`, or FASTA path).
- `verbose`: Print progress output.

</details>

---

### final

[↑ Back to TOC](#table-of-contents)

**What it does**: Concatenates all columns into synthesis-ready oligos.

**When to use it**: Last step before ordering.

```python
final_df, stats = op.final(
    input_data=df,
    output_file='synthesis_ready',
)
# Output contains 'CompleteOligo' and 'OligoLength' columns
```

> **Tip**: `final` does not preserve annotation columns (it returns only the concatenated oligo + length). Save the annotated design DataFrame separately if you need it for indexing/counting.

<details>
<summary><b>Full parameter reference (Python API)</b></summary>

```python
final_df, stats = op.final(
    input_data,     # str | pd.DataFrame
    output_file=None,# str | None
    verbose=True,   # bool
)
```

**Required parameters**
- `input_data`: CSV path or DataFrame with `ID` + DNA sequence columns.

**Optional parameters**
- `output_file`: Write output CSV (library mode only; required in CLI).
- `verbose`: Print progress output.

</details>

---

## Analysis Mode

[↑ Back to TOC](#table-of-contents)

Analysis Mode is where sequencing data meets your designed library. You've done the experiment, now let's count some barcodes.

### index

[↑ Back to TOC](#table-of-contents)

**What it does**: Builds a searchable index of barcodes (and optional associates) for counting.

**When to use it**: Before counting - you need to index what you're looking for.

```python
stats = op.index(
    barcode_data='library.csv',
    barcode_column='BC1',
    index_file='bc1_index',                # Creates bc1_index.oligopool.index
    barcode_prefix_column='BC1_Prefix',    # Constant anchor before barcode
    barcode_suffix_column='BC1_Suffix',    # Constant anchor after barcode
    barcode_prefix_gap=0,                  # Bases between prefix and barcode
)
```

**With associates** (variant sequences linked to barcodes):
```python
stats = op.index(
    barcode_data='library.csv',
    barcode_column='BC1',
    index_file='bc1_index',
    barcode_prefix_column='Anchor5',
    associate_data='library.csv',          # Can be same or different file
    associate_column='Variant',
    associate_suffix_column='Anchor3',
)
```

<details>
<summary><b>Full parameter reference (Python API)</b></summary>

```python
stats = op.index(
    # Required
    barcode_data,               # str | pd.DataFrame
    barcode_column,             # str
    index_file,                 # str

    # Optional (barcode anchors)
    barcode_prefix_column=None, # str | None
    barcode_suffix_column=None, # str | None
    barcode_prefix_gap=0,       # int
    barcode_suffix_gap=0,       # int

    # Optional (associate indexing)
    associate_data=None,        # str | pd.DataFrame | None
    associate_column=None,      # str | None
    associate_prefix_column=None,# str | None
    associate_suffix_column=None,# str | None

    # Optional
    verbose=True,               # bool
)
```

**Required parameters**
- `barcode_data`: CSV path or DataFrame with `ID` + barcode column.
- `barcode_column`: Barcode column name to index.
- `index_file`: Output basename (writes `<name>.oligopool.index`).

**Optional parameters**
- `barcode_prefix_column` / `barcode_suffix_column`: Constant (single-unique) anchors flanking the barcode (at least one required).
- `barcode_prefix_gap` / `barcode_suffix_gap`: Number of bases between the anchor and the barcode in the read.
- `associate_data` / `associate_column`: Add associate verification for `acount` (barcode↔associate coupling).
- `associate_prefix_column` / `associate_suffix_column`: Constant anchors for the associate (at least one required when associates are used).
- `verbose`: Print progress output.

</details>

---

### pack

[↑ Back to TOC](#table-of-contents)

**What it does**: Preprocesses FastQ files - filters, optionally merges paired-ends, deduplicates.

**When to use it**: After sequencing, before counting.

```python
stats = op.pack(
    r1_fastq_file='sample_R1.fq.gz',
    r2_fastq_file='sample_R2.fq.gz',       # Optional for paired-end
    r1_read_type=0,                        # 0=forward, 1=reverse
    r2_read_type=1,
    minimum_r1_read_length=50,
    minimum_r2_read_length=50,
    minimum_r1_read_quality=30,            # Phred score filter
    minimum_r2_read_quality=30,
    pack_type=1,                           # 1=merge pairs, 0=R1 only
    pack_file='sample',                    # Creates sample.oligopool.pack
)
```

<details>
<summary><b>Full parameter reference (Python API)</b></summary>

```python
stats = op.pack(
    # Required
    r1_fastq_file,            # str
    r1_read_type,             # int (0 forward, 1 reverse)
    pack_type,                # int (0 concatenated, 1 merged)
    pack_file,                # str

    # Optional (R1 filters)
    minimum_r1_read_length=1, # int
    minimum_r1_read_quality=20,# int

    # Optional (paired-end)
    r2_fastq_file=None,       # str | None
    r2_read_type=None,        # int | None
    minimum_r2_read_length=None,# int | None
    minimum_r2_read_quality=None,# int | None

    # Optional (performance)
    pack_size=3.0,            # float (million unique reads per pack)
    core_count=0,             # int (0 auto)
    memory_limit=0.0,         # float (GB per core; 0 auto)

    # Optional
    verbose=True,             # bool
)
```

**Required parameters**
- `r1_fastq_file`: R1 FastQ path (can be gzipped).
- `r1_read_type`: `0` forward, `1` reverse.
- `pack_type`: `0` store concatenated reads, `1` store merged/assembled reads.
- `pack_file`: Output basename (writes `<name>.oligopool.pack`).

**Optional parameters**
- `minimum_r1_read_length` / `minimum_r1_read_quality`: Filters for R1.
- `r2_fastq_file` / `r2_read_type` / `minimum_r2_read_length` / `minimum_r2_read_quality`: Paired-end inputs and filters (set these to enable paired-end/merge workflows).
- `pack_size`: Pack size in millions of unique reads (trade memory vs speed).
- `core_count`: Number of CPU cores (0 = auto).
- `memory_limit`: GB per core (0 = auto).
- `verbose`: Print progress output.

</details>

---

### acount

[↑ Back to TOC](#table-of-contents)

**What it does**: Association counting - verifies that the expected variant is indeed coupled to its target barcode in each read.

**When to use it**: Validating synthesis/cloning accuracy, confirming barcode-variant associations after library construction.

```python
df, stats = op.acount(
    index_file='bc1_index',
    pack_file='sample',
    count_file='counts',                   # Creates counts.oligopool.acount.csv
    barcode_errors=1,                      # Allow 1 mismatch (-1=auto)
)
```

<details>
<summary><b>Full parameter reference (Python API)</b></summary>

```python
counts_df, stats = op.acount(
    # Required
    index_file,         # str
    pack_file,          # str
    count_file,         # str

    # Optional
    mapping_type=0,     # int (0 fast, 1 sensitive)
    barcode_errors=-1,  # int (-1 auto)
    associate_errors=-1,# int (-1 auto)
    callback=None,      # callable | None (Python only)
    core_count=0,       # int (0 auto)
    memory_limit=0.0,   # float (GB per core; 0 auto)
    verbose=True,       # bool
)
```

**Required parameters**
- `index_file`: Index basename or path (reads `<name>.oligopool.index`).
- `pack_file`: Pack basename or path (reads `<name>.oligopool.pack`).
- `count_file`: Output basename (writes `<name>.oligopool.acount.csv`).

**Optional parameters**
- `mapping_type`: `0` fast/near-exact, `1` slow/sensitive.
- `barcode_errors` / `associate_errors`: Max errors allowed (use `-1` to auto-infer from index).
- `callback`: Python-only hook for custom filtering (CLI always uses `None`).
- `core_count` / `memory_limit`: Parallelization controls (0 = auto).
- `verbose`: Print progress output.

</details>

---

### xcount

[↑ Back to TOC](#table-of-contents)

**What it does**: Barcode-focused counting - maps reads to one or more barcode indices without associate verification.

**When to use it**: Pure barcode counting (single or multiple indices), combinatorial barcode designs (BC1 × BC2).

```python
# Single barcode index
df, stats = op.xcount(
    index_files=['bc1_index'],
    pack_file='sample',
    count_file='barcode_counts',
)

# Multiple barcode indices (combinatorial)
df, stats = op.xcount(
    index_files=['bc1_index', 'bc2_index'],
    pack_file='sample',
    count_file='combo_counts',
    mapping_type=1,                          # 0=fast, 1=sensitive
    barcode_errors=-1,                       # Auto-infer
)
```

Output includes all observed combinations, with `'-'` for missing barcodes in partial reads.

**acount vs xcount**: Use `acount` when you need barcode+variant association verification; use `xcount` for barcode-only counting (single or combinatorial).

**Custom callbacks** (Python API only):
```python
def my_filter(read, ID, count, coreid):
    # Custom logic here
    return True  # Accept read

df, stats = op.xcount(..., callback=my_filter)
```

<details>
<summary><b>Full parameter reference (Python API)</b></summary>

```python
counts_df, stats = op.xcount(
    # Required
    index_files,        # str | list
    pack_file,          # str
    count_file,         # str

    # Optional
    mapping_type=0,     # int (0 fast, 1 sensitive)
    barcode_errors=-1,  # int (-1 auto)
    callback=None,      # callable | None (Python only)
    core_count=0,       # int (0 auto)
    memory_limit=0.0,   # float (GB per core; 0 auto)
    verbose=True,       # bool
)
```

**Required parameters**
- `index_files`: One index (string) or many (list) for combinatorial counting.
- `pack_file`: Pack basename or path (reads `<name>.oligopool.pack`).
- `count_file`: Output basename (writes `<name>.oligopool.xcount.csv`).

**Optional parameters**
- `mapping_type`: `0` fast/near-exact, `1` slow/sensitive.
- `barcode_errors`: Max barcode errors (use `-1` to auto-infer from index).
- `callback`: Python-only hook for custom filtering (CLI always uses `None`).
- `core_count` / `memory_limit`: Parallelization controls (0 = auto).
- `verbose`: Print progress output.

</details>

---

## Workflows

[↑ Back to TOC](#table-of-contents)

### Basic Library Design

```python
import pandas as pd
import oligopool as op

# 1. Start with variants
df = pd.read_csv('variants.csv')  # Must have 'ID' column

# 2. Add forward primer
df, _ = op.primer(
    input_data=df,
    oligo_length_limit=200,
    primer_sequence_constraint='N' * 20,
    primer_type=0,
    minimum_melting_temperature=55,
    maximum_melting_temperature=60,
    maximum_repeat_length=10,
    primer_column='FwdPrimer',
    right_context_column='Variant',
)

# 3. Add barcode
df, _ = op.barcode(
    input_data=df,
    oligo_length_limit=200,
    barcode_length=12,
    minimum_hamming_distance=3,
    maximum_repeat_length=8,
    barcode_column='BC1',
    left_context_column='FwdPrimer',
    right_context_column='Variant',
)

# 4. Add reverse primer (Tm-matched)
df, _ = op.primer(
    input_data=df,
    oligo_length_limit=200,
    primer_sequence_constraint='N' * 20,
    primer_type=1,
    minimum_melting_temperature=55,
    maximum_melting_temperature=60,
    maximum_repeat_length=10,
    primer_column='RevPrimer',
    left_context_column='Variant',
    paired_primer_column='FwdPrimer',
)

# 5. Check length
op.lenstat(input_data=df, oligo_length_limit=200)

# 6. Verify and finalize
op.verify(input_data=df, oligo_length_limit=200)
# Save the annotated library (for indexing/counting) before `final`, which drops annotation columns.
df.to_csv('library_design.csv', index=False)
final_df, _ = op.final(input_data=df, output_file='library_for_synthesis')
```

### Analysis Pipeline

```python
import oligopool as op

# 1. Index your barcodes
op.index(
    barcode_data='library_design.csv',
    barcode_column='BC1',
    barcode_prefix_column='FwdPrimer',
    associate_data='library_design.csv',
    associate_column='Variant',
    associate_suffix_column='RevPrimer',
    index_file='bc1_assoc',
)

# 2. Pack your reads
op.pack(
    r1_fastq_file='experiment_R1.fq.gz',
    r2_fastq_file='experiment_R2.fq.gz',
    r1_read_type=0,
    r2_read_type=1,
    minimum_r1_read_quality=30,
    minimum_r2_read_quality=30,
    pack_type=1,
    pack_file='experiment',
)

# 3. Count!
counts_df, stats = op.acount(
    index_file='bc1_assoc',
    pack_file='experiment',
    count_file='results',
)

print(counts_df)
```

---

## CLI Reference

[↑ Back to TOC](#table-of-contents)

Every module is available via command line:

```bash
# See all commands
op

# Get help on a specific command
op barcode

# Read the manual
op manual barcode
op cite

# Run a command
op barcode \
    --input-data variants.csv \
    --oligo-length-limit 200 \
    --barcode-length 12 \
    --minimum-hamming-distance 3 \
    --maximum-repeat-length 8 \
    --barcode-column BC1 \
    --left-context-column Variant \
    --output-file with_barcodes
```

**Install tab completion** (highly recommended):
```bash
op complete --install
```

### CLI-Specific Notes

**`--output-file` is required for DataFrame-producing commands**: Unlike the Python API where `output_file=None` returns results in-memory, CLI commands that produce a DataFrame require `--output-file`.

**Output filenames are auto-suffixed**: Most commands append a suffix if missing (e.g., `.oligopool.barcode.csv`), so prefer basenames like `--output-file my_library` to avoid doubled extensions.

**Stats output**: Use `--stats-json` to print the stats dict as JSON to stdout, or `--stats-file path.json` to write it to disk.

**Quiet mode**: Use `--quiet` to suppress verbose module output (sets `verbose=False`).

**Patch mode**: For `barcode`/`primer`/`motif`/`spacer`, use `--patch-mode` to fill only missing values in an existing column (no overwrites).

**Cross-set barcodes**: For multi-barcode designs (BC1 vs BC2), use `--cross-barcode-columns BC1` with `--minimum-cross-distance N` on the subsequent barcode designs.

**Sequence constraint shorthand**: For `--primer-sequence-constraint` and `--motif-sequence-constraint`, you can use:
```bash
# Quoted Python-style expression
--primer-sequence-constraint "'N'*20"

# Shorthand concatenation (no quotes needed)
--motif-sequence-constraint GCC+N*20+CCG
```

**Callbacks are Python-only**: The `callback` parameter in `acount`/`xcount` is not available via CLI. For custom read processing logic, use the Python API.

---

## Tips & Tricks

[↑ Back to TOC](#table-of-contents)

### Design Tips

1. **Design anchors before barcodes**: Use `motif(motif_type=1)` for constant anchors, then design barcodes adjacent to them.

2. **Check length early and often**: Run `lenstat` after each element to avoid surprises.

3. **Start loose, tighten constraints**: Begin with relaxed parameters, then increase stringency if design succeeds easily.

4. **Use Patch Mode for iterations**: When adding variants to an existing pool, use `patch_mode=True` to preserve existing designs.

5. **Mind the edge effects**: Always specify context columns - repeat collisions at boundaries are sneaky.

### Analysis Tips

1. **Index anchors matter**: The quality of your anchors determines counting accuracy. Design them well.

2. **Quality filter aggressively**: Low-quality reads just add noise. Use `minimum_*_read_quality=30` or higher.

3. **Check your stats**: The returned stats dict tells you exactly what happened. Read it.

4. **Start with sensitive mapping**: Use `mapping_type=1` initially to catch more reads, then compare with `mapping_type=0`.

### Performance Tips

1. **Parallelize**: Most modules auto-detect cores. For analysis, ensure you have RAM (memory_limit parameter).

2. **Use fast barcode mode**: `barcode_type=0` is usually sufficient and much faster.

3. **Pack once, count many**: The pack file can be reused across multiple counting runs.

---

## Advanced Modules

[↑ Back to TOC](#table-of-contents)

For power users who want to peek under the hood:

### vectorDB

LevelDB-based k-mer storage. Created by `background()`, but you can access it directly:

```python
db = op.vectorDB(
    path='ecoli_bg.oligopool.background',
    maximum_repeat_length=15,  # ignored if reopening an existing DB
)
```

Useful for inspecting or manipulating background databases.

<details>
<summary><b>vectorDB reference</b></summary>

```python
db = op.vectorDB(path, maximum_repeat_length)

db.add(seq, rna=True)
db.multiadd(seq_list, rna=True)

seq in db                 # True if any k-mer from seq exists in the DB
db.multicheck(seq_list)   # list[bool] for a list of sequences

db.remove(seq, rna=True)
db.multiremove(seq_list, rna=True, clear=False)

db.clear(maxreplen=None)  # wipe DB; optionally change max repeat length
db.close()                # close handle (keeps files)
db.drop()                 # delete DB from disk
```

**Constructor parameters**
- `path` (`str`): Directory path where the DB lives (background outputs end with `.oligopool.background`).
- `maximum_repeat_length` (`int`): Max repeat length (`K = maximum_repeat_length + 1`). If reopening an existing DB, the stored `K` is used.

</details>

### Scry

1-nearest-neighbor barcode classifier. Powers `acount`/`xcount` internally:

```python
model = op.Scry()
model.train(X=['ATGC...','GCTA...'], Y=[0,1], n=8, k=6, t=1)
model.prime(t=1, mode=0)   # 0 fast, 1 sensitive
label, score = model.predict('ATGC....')
```

Useful for building custom counting pipelines or debugging classification issues.

<details>
<summary><b>Scry reference</b></summary>

```python
model = op.Scry()

model.train(X, Y, n, k, t, liner=None)
model.prime(t=0, mode=0)
model.predict(x)
```

**Key parameters**
- `train(X, Y, n, k, t, liner=None)`:
  - `X`: iterable of fixed-length barcode strings
  - `Y`: iterable of labels/IDs (same length as `X`)
  - `n`: barcode length
  - `k`: k-mer length used for indexing/screening
  - `t`: max errors used to build the model
  - `liner`: optional progress printer (used internally by oligopool)
- `prime(t=0, mode=0)`:
  - `t`: errors assumed at prediction time
  - `mode`: `0` fast/near-exact, `1` slow/sensitive
- `predict(x)`: returns `(label, score)` or `(None, None)` if ambiguous/unresolvable.

</details>

---

## Troubleshooting

[↑ Back to TOC](#table-of-contents)

### Design Failures

**"Barcode Design Infeasible"**
- Reduce `minimum_hamming_distance`
- Increase `barcode_length`
- Switch to `barcode_type=0` (terminus optimized)
- Relax `maximum_repeat_length`
- Remove some `excluded_motifs`

**"Primer Design Infeasible"**
- Widen the Tm range (`minimum/maximum_melting_temperature`)
- Use more degenerate constraint (more `N`s)
- Relax `maximum_repeat_length`
- Check if background is too restrictive

**"Not enough free space"**
- Your oligos are too long for `oligo_length_limit`
- Run `lenstat` to see where space is being used
- Consider splitting with `split` + `pad`

### Analysis Issues

**Low barcode mapping rate**
- Check your anchors are correct and present in reads
- Try `mapping_type=1` (sensitive)
- Increase `barcode_errors` tolerance
- Verify read orientation (`r1_read_type`, `r2_read_type`)

**Missing combinations in xcount**
- Partial reads get `'-'` for missing barcodes (this is expected)
- Check if anchors for all indexes are present in reads
- Verify read length covers all barcode positions

### General

**"ID column missing/invalid"**
- Input DataFrame must have a unique `ID` column
- IDs must be unique (no duplicates)
- IDs should be strings

**Memory issues**
- Use `memory_limit` parameter in analysis functions
- Process in batches for very large datasets
- Pack files are chunked for efficient processing

---

<p align="center">
<i>Made with molecules and math by the Salis Lab</i>
<br>
<a href="https://github.com/ayaanhossain/oligopool">GitHub</a> ·
<a href="https://pubs.acs.org/doi/10.1021/acssynbio.4c00661">Paper</a> ·
<a href="https://github.com/ayaanhossain/oligopool/issues">Issues</a>
</p>
