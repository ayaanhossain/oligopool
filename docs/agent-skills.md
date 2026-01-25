# Oligopool Calculator - AI Agent Guide

This document is designed for AI assistants to understand and help users with oligopool library design and analysis.

## Using This Guide

Share this document with your AI assistant for better oligopool design help:
- **Claude Code / Cursor / Copilot:** The AI can read `docs/agent-skills.md` directly from the repo
- **ChatGPT / Claude.ai:** Upload this file, or share the raw URL (AIs can fetch raw URLs directly)
- **Custom agents:** Include this document in your system prompt or fetch from the raw URL

```
https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/agent-skills.md
```

**Related documentation**:
- [README](https://github.com/ayaanhossain/oligopool/blob/master/README.md) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/README.md)] - Overview, installation, and quick start
- [User Guide](https://github.com/ayaanhossain/oligopool/blob/master/docs/docs.md) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/docs.md)] - Tutorials, examples, and workflows for end users
- [API Reference](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] - Complete parameter documentation for all modules
- [Docker Guide](https://github.com/ayaanhossain/oligopool/blob/master/docs/docker-notes.md) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/docker-notes.md)] - Run `oligopool` in a container

## Operation Policy (Keep This Useful)

- **Source of truth**: `docs/docs.md` + `docs/api.md` + runtime docs (`help(op)` / `help(op.<module>)`) + real CLI behavior.
- **This guide**: Agent-facing contracts, gotchas, and workflow scaffolding (Notebook / Script / CLI).
- **If anything disagrees**: Trust runtime behavior, then update `docs/docs.md`, `docs/api.md`, and this guide together.

**Quick facts (the ones that matter in practice):**
- **Entry points**: `op` and `oligopool` are equivalent CLIs.
- **CLI help model**: There is no `--help` flag. Use `op COMMAND` for command options and `op manual COMMAND` for docstrings.
- **Basenames**: Prefer basenames for all outputs (`--output-file`, `--index-file`, `--pack-file`, `--count-file`); suffixes are auto-appended as needed.
- **String type parameters**: Type parameters accept integers OR descriptive strings (case-insensitive, fuzzy-matched):
  - `barcode_type`: `0`/`'terminus'`/`'fast'` or `1`/`'spectrum'`/`'slow'`
  - `primer_type`: `0`/`'forward'` or `1`/`'reverse'`
  - `motif_type`: `0`/`'variable'`/`'per-variant'` or `1`/`'constant'`/`'anchor'`
  - `pack_type`: `0`/`'concatenate'` or `1`/`'merge'`
  - `r1/r2_read_type`: `0`/`'forward'` or `1`/`'reverse'`
  - `mapping_type`: `0`/`'fast'` or `1`/`'sensitive'`
- **Return shapes**:
  - Design/transform: `(out_df, stats_dict)`
  - Stats-only: `stats_dict` (`background`, `lenstat`, `verify`, `index`, `pack`)
  - Counting: `(counts_df, stats_dict)` (`acount`, `xcount`)
  - Split with `separate_outputs=True`: `([df_Split1, df_Split2, ...], stats_dict)`
- **ID handling**: Input requires a unique `ID`. CSV outputs include an explicit `ID` column (no pandas index column).
- **Patch Mode**: `patch_mode=True` / `--patch-mode` fills only missing values (None/NaN/empty/`'-'`) and never overwrites existing designs.
- **Cross-set barcodes**: `cross_barcode_columns` and `minimum_cross_distance` must be provided together.

## Package Overview

`Oligopool Calculator` is a Python library for designing and analyzing oligonucleotide pool (oligopool)
libraries used in massively parallel reporter assays (MPRAs) and similar high-throughput experiments.

**Common applications**:
- **MPRAs**: Promoter/enhancer activity measurement
- **Saturation mutagenesis**: Systematic variant testing
- **CRISPR screens**: Pooled guide libraries
- **Ribozyme/mRNA studies**: Cleavage and stability quantification

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

### Docker mode
`Oligopool Calculator` can be run in a Docker container, useful for:
- Cross-platform consistency (works identically on macOS, Windows, Linux)
- Isolated environments without affecting your system Python
- Reproducible analysis pipelines in CI/CD or HPC environments

**Quick start:**
```bash
# Clone and build the image
git clone https://github.com/ayaanhossain/oligopool.git
cd oligopool
docker build -t oligopool-docker .

# Run from your project directory
cd /path/to/your/project
docker run -it -v $(pwd):/op-workspace --name op-container oligopool-docker

# Inside the container, use Python API or CLI
python -c "import oligopool as op; print(op.__version__)"
op barcode  # Show barcode command options
```

**With Jupyter:**
```bash
# Map port for Jupyter access
docker run -p 8888:8080 -it -v $(pwd):/op-workspace --name op-container oligopool-docker

# Inside container, start Jupyter
jupyter notebook --ip=0.0.0.0 --port=8080 --no-browser --allow-root
# Access at http://localhost:8888
```

For detailed Docker instructions, see [docker-notes.md](https://github.com/ayaanhossain/oligopool/blob/master/docs/docker-notes.md) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/docker-notes.md)].

---

## Design Mode - Module Reference

For complete parameter documentation, see [api.md](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)].

### barcode

**Purpose**: Generate Hamming-distance separated barcodes for unique variant identification.

**Design order**: After primers/motifs, before spacers.

**API**: See [`barcode`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#barcode) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for parameters.

**Tips**:
- Use `barcode_type='terminus'` (or `0`) for large libraries (faster)
- Use `barcode_type='spectrum'` (or `1`) for maximum barcode diversity
- For multi-barcode designs, use `cross_barcode_columns` to ensure BC2 is separated from BC1
- Cross-set separation is strict: provide both `cross_barcode_columns` and `minimum_cross_distance`;
  all cross-set barcode values must already be A/T/G/C strings of length `barcode_length`.

---

### primer

**Purpose**: Design thermodynamically optimal primers for amplification.

**Design order**: Design primers early. For paired primers: design inner primer first, then outer with `paired_primer_column`.

**API**: See [`primer`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#primer) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for parameters.

**Tips**:
- Use `'SS' + 'N'*18` for 5' GC clamp
- Use `'N'*18 + 'WW'` for 3' AT clamp
- Build background with `op.background()` before primer design for off-target screening
- For multiplexed libraries, use `oligo_sets` to design set-specific primers

**`oligo_sets` for multiplexed libraries**:
When your library has distinct groups (e.g., different promoters, different genes), use `oligo_sets` to design group-specific primers while ensuring cross-compatibility:
```python
# Each set gets its own primers, but all primers are screened against each other
df, _ = op.primer(
    input_data=df,
    oligo_sets=['setA', 'setA', 'setB', 'setB', 'setC'],  # Per-row labels
    # Or: oligo_sets='sets.csv' with ID + OligoSet columns
    # Or: oligo_sets=pd.DataFrame({'ID': [...], 'OligoSet': [...]})
    ...
)
```
- All sets share the same Tm/constraint requirements
- Primers are cross-screened to prevent amplification of wrong sets
- In patch mode, existing per-set primers are reused for their sets

---

### motif

**Purpose**: Add sequence motifs or constant anchors.

**Design order**: Before barcodes if designing anchors for indexing.

**API**: See [`motif`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#motif) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for parameters.

**Use cases**:
- Restriction sites: `motif_sequence_constraint='GAATTC'` (EcoRI)
- Barcode anchors: `motif_type='constant'` (or `1`) with `'N'*10` for constant anchor
- Degenerate regions: `'NNNGGATCCNNN'` (BamHI with flanking Ns)

---

### spacer

**Purpose**: Add neutral filler DNA to reach target oligo length.

**Design order**: Last, after all other elements.

**API**: See [`spacer`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#spacer) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for parameters.

**spacer_length options**:
- `None`: Auto-fill to reach `oligo_length_limit`
- `int`: Fixed length for all oligos
- `list`: Per-variant lengths aligned to input
- `DataFrame`: With 'ID' and 'Length' columns

---

### background

**Purpose**: Build k-mer database for primer off-target screening.

**When to use**: Before primer design when screening against host genome/plasmid.

**API**: See [`background`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#background) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for parameters.

**input_data formats**:
- List of DNA strings: `['ATGC...', 'GCTA...']`
- CSV file: With 'Sequence' column
- DataFrame: With 'Sequence' column
- FASTA file: `.fa`, `.fasta`, `.fna` (optionally gzipped)

---

### split

**Purpose**: Break long oligos into overlapping fragments for assembly.

**When to use**: When oligos exceed synthesis length limits (typically >200 bp).

**API**: See [`split`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#split) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for parameters.

**Key concept**: Each `SplitN` column = a separate oligo pool to synthesize. If `split` returns `Split1`, `Split2`, `Split3`, you order **three separate pools** from your vendor.

**Return formats** (different defaults for library vs CLI):
- **Library default** (`separate_outputs=False`): Single DataFrame with `Split1`, `Split2`, ... columns
- **Library with** `separate_outputs=True`: List of DataFrames `[df_Split1, df_Split2, ...]`
- **CLI default**: Separate files (`out.Split1.oligopool.split.csv`, etc.)
- **CLI with** `--no-separate-outputs`: Single combined file

**Tips**:
- Raw split output is NOT synthesis-ready - use `pad` to add primers/Type IIS sites
- Even-numbered splits (`Split2`, `Split4`, ...) are reverse-complemented for PCR assembly
- Fragment count is auto-determined and can vary per oligo
- Use `separate_outputs=True` in library mode for cleaner pad workflows (see Long Oligo Assembly)

---

### pad

**Purpose**: Add amplification primers with Type IIS sites for assembly.

**When to use**: After split, for each fragment.

**API**: See [`pad`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#pad) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for parameters.

**Critical workflow**: Run `pad` **once per split column**, then `final` on each:

```python
# Using separate_outputs=True (recommended)
split_dfs, _ = op.split(..., separate_outputs=True)  # Returns [df_Split1, df_Split2, ...]

for i, split_df in enumerate(split_dfs, start=1):
    pad_df, _ = op.pad(split_df, split_column=f'Split{i}', typeIIS_system='BsaI', ...)
    final_df, _ = op.final(pad_df, output_file=f'synthesis_Split{i}')
```

**Tips**:
- You cannot pad all columns in one call - iterate
- Output: `5primeSpacer`, `ForwardPrimer`, `<split_column>`, `ReversePrimer`, `3primeSpacer`
- **Exclude your Type IIS motif from upstream elements** (e.g., `GGTCTC`/`GAGACC` for BsaI in `excluded_motifs`) - `pad` validates fragments and fails early if internal sites found

**Post-synthesis**: PCR amplify → Type IIS digest (removes pads, leaves enzyme-specific overhangs) → mung bean nuclease (blunts overhangs; skip for blunt-cutters like `MlyI`) → assemble via **split-designed overlaps** (Gibson, overlap-extension PCR). Type IIS removes pads; the 15–30 bp overlaps from `split` drive assembly.

**Supported Type IIS enzymes (34 total):**
`AcuI`, `AlwI`, `BbsI`, `BccI`, `BceAI`, `BciVI`, `BcoDI`, `BmrI`, `BpuEI`, `BsaI`, `BseRI`, `BsmAI`,
`BsmBI`, `BsmFI`, `BsmI`, `BspCNI`, `BspQI`, `BsrDI`, `BsrI`, `BtgZI`, `BtsCI`, `BtsI`, `BtsIMutI`,
`EarI`, `EciI`, `Esp3I`, `FauI`, `HgaI`, `HphI`, `HpyAV`, `MlyI`, `MnlI`, `SapI`, `SfaNI`

For unsupported enzymes, design primers/sites manually with `primer` or `motif`.

---

### merge

**Purpose**: Collapse multiple contiguous columns into one.

**Design order**: Mid-pipeline maneuver. Use after designing elements you want to combine, before further processing.

**API**: See [`merge`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#merge) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for parameters.

**When to use**:
- Simplify DataFrame structure by combining adjacent elements (e.g., `Primer1 + BC1` → `5prime_region`)
- Prepare for `revcomp` when you want to reverse-complement a logical unit
- Reduce column count before `final` for cleaner output

**Note**: Source columns in the merge range are removed from the output DataFrame.

---

### revcomp

**Purpose**: Reverse complement a column range and reverse column order.

**Design order**: Mid-pipeline maneuver. Use when switching strand orientation after designing elements.

**API**: See [`revcomp`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#revcomp) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for parameters.

**When to use**:
- Design in "readout orientation" but synthesize in opposite orientation
- Verify `split` fragment orientation (even-numbered splits are reverse-complemented by design)
- Switch between sense/antisense strand representation mid-pipeline

**Note**: Both the sequence content AND the column order are reversed within the specified range.

---

### lenstat

**Purpose**: Check length statistics and remaining space (non-destructive).

**When to use**: Frequently during design to monitor remaining space.

**API**: See [`lenstat`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#lenstat) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for parameters.

---

### verify

**Purpose**: QC check before synthesis.

**When to use**: Before `final()` to catch constraint violations.

**API**: See [`verify`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#verify) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for parameters.

**Column concatenation**:
- Only **sequence columns** (DNA/IUPAC) concatenated; metadata columns skipped
- Sequence columns joined **left-to-right in DataFrame column order**
- Gap characters (`'-'`) stripped during concatenation
- If `CompleteOligo` exists, used directly (junction attribution skipped)
- Junctions follow column order: `[A, B, C]` → `A|B`, `B|C`

---

### final

**Purpose**: Concatenate columns into synthesis-ready oligos.

**When to use**: Final step before synthesis.

**Output**: DataFrame with `CompleteOligo` and `OligoLength` (plus an explicit `ID` column unless the caller
provided `ID` as the DataFrame index).

**API**: See [`final`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#final) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for parameters.

---

## Analysis Mode - Module Reference

For complete parameter documentation, see [api.md](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)].

### index

**Purpose**: Build barcode index for counting.

**API**: See [`index`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#index) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for parameters.

**Tips**:
- Prefix/suffix anchors must be constant (single unique sequence, ≥6 bp) and adjacent to the indexed column.
- For `acount`, specify `associate_data`/`associate_column` and at least one constant adjacent associate prefix/suffix.
- Design anchors with `motif(motif_type='constant', ...)` (or use universal primers if they are constant).

---

### pack

**Purpose**: Preprocess and deduplicate FASTQ reads.

**API**: See [`pack`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#pack) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for parameters.

**Tips**:
- For single-end reads, use R1 arguments only
- For paired-end merging (`pack_type='merge'`), provide both `r2_fastq_file` and `r2_read_type`
- For pre-merged reads (e.g., from FLASH), use as single-end

**`pack_size` parameter**:
Controls memory usage by splitting reads into chunks of N million unique reads (default: 3.0, range: 0.1-5.0).
- Lower values (0.5-1.0): Less memory, more pack files, useful for memory-constrained systems
- Higher values (3.0-5.0): More memory, fewer pack files, faster counting
- Each pack file is processed independently during counting

---

### acount

**Purpose**: Association counting - verify barcode-variant coupling.

**API**: See [`acount`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#acount) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for parameters.

**Output columns**: one `<indexname>.ID` column per index (typically one), plus `BarcodeCounts`, `AssociationCounts`

**When to use**:
- Validating synthesis accuracy
- Input library QC
- Post-assay barcode-variant verification

**Error tolerance parameters**:
- `barcode_errors` (`int`, default `-1`): Max mismatches allowed in barcode matching. `-1` auto-infers from index (typically `minimum_hamming_distance // 2`).
- `associate_errors` (`int`, default `-1`): Max mismatches allowed in associate (variant) matching. `-1` auto-infers from index.
- Set explicit values (e.g., `0` or `1`) for stricter/looser matching

---

### xcount

**Purpose**: Barcode-only counting (single or combinatorial).

**API**: See [`xcount`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#xcount) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for parameters.

**Output**: one `<indexname>.ID` column per index, plus `CombinatorialCounts`. Missing barcodes shown as `'-'`.

**When to use**:
- Pure barcode counting (single index)
- Multi-barcode combinations (BC1 × BC2)
- Cleaved vs uncleaved quantification (ribozymes, etc.)

**Error tolerance**:
- `barcode_errors` (`int`, default `-1`): Max mismatches allowed in barcode matching. `-1` auto-infers from index.

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

For complete method documentation, see [api.md](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#advanced) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)].

### vectorDB

**Purpose**: LevelDB-based scalable on-disk k-mer storage for background databases.

**When to use**: Direct manipulation of background k-mer databases, custom screening pipelines.

**API**: See [`vectorDB`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#vectordb) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for methods.

**Note**: When reopening an existing vectorDB, `maximum_repeat_length` is ignored and loaded from the instance.

---

### Scry

**Purpose**: 1-nearest-neighbor barcode classifier used internally by `acount`/`xcount`.

**When to use**: Building custom counting pipelines, debugging barcode classification, custom analysis workflows.

**API**: See [`Scry`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#scry) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for methods.

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
    primer_type='forward',
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
    primer_type='reverse',
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

---

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

---

### Long Oligo Assembly (split → pad → final)

When oligos exceed synthesis length limits (~200 bp), use `split` + `pad` to break them into assembly-ready fragments:

```python
import oligopool as op

# 1. Split long oligos into separate DataFrames (one per fragment)
split_dfs, stats = op.split(
    input_data='long_oligos.csv',
    split_length_limit=150,
    minimum_melting_temperature=55.0,
    minimum_hamming_distance=3,
    minimum_overlap_length=20,
    maximum_overlap_length=30,
    separate_outputs=True,  # Enable to return [df_Split1, df_Split2, ...]
)

# 2. Pad each fragment, then finalize
# Each df in split_dfs has ID + one SplitN column
for i, split_df in enumerate(split_dfs, start=1):
    split_col = f'Split{i}'

    # Add primers + Type IIS sites
    pad_df, _ = op.pad(
        input_data=split_df,
        split_column=split_col,
        typeIIS_system='BsaI',
        oligo_length_limit=200,
        minimum_melting_temperature=52.0,
        maximum_melting_temperature=58.0,
        maximum_repeat_length=10,
    )

    # Generate synthesis-ready oligos
    final_df, _ = op.final(
        input_data=pad_df,
        output_file=f'synthesis_{split_col}',
    )

# You now have N synthesis files to order from your vendor
# After synthesis, combine fragments via overlap-based assembly (Gibson, overlap-extension PCR, etc.)
```

**Key points**:
- Enable `separate_outputs` to get a list of DataFrames directly
- Each `SplitN` column = one separate oligo pool to order
- Run `pad` once per fragment (cannot batch)
- Even-numbered splits are reverse-complemented (for PCR assembly orientation)
- Raw split output is NOT synthesis-ready - always use `pad` + `final`
- Exclude your Type IIS motif from upstream elements to prevent internal cut sites

---

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
    r1_read_type='forward',
    r2_read_type='reverse',
    pack_type='merge',
    minimum_r1_read_quality=30,
    minimum_r2_read_quality=30,
    pack_file='reads',
)

# 3a. Association counting (verify barcode-variant coupling)
ac_df, _ = op.acount(
    index_file='bc1_assoc_idx',
    pack_file='reads',
    count_file='association_counts',
    mapping_type='sensitive',
)

# 3b. Combinatorial counting (BC1 × BC2)
xc_df, _ = op.xcount(
    index_files=['bc1_idx', 'bc2_idx'],
    pack_file='reads',
    count_file='combo_counts',
    mapping_type='sensitive',
)
```

---

### Saturation Mutagenesis Library
Saturation mutagenesis systematically tests every possible mutation at each position in a sequence. `Oligopool Calculator` handles the oligo design; you provide the variant sequences.

```python
import oligopool as op
import pandas as pd

# 1. Generate saturation mutagenesis variants (user-provided)
# Example: all single-nucleotide mutations of a 20bp region
wild_type = 'ATGCATGCATGCATGCATGC'
variants = []
for pos in range(len(wild_type)):
    for nt in 'ATGC':
        if nt != wild_type[pos]:
            mutant = wild_type[:pos] + nt + wild_type[pos+1:]
            variants.append({
                'ID': f'pos{pos+1}_{wild_type[pos]}to{nt}',
                'Variant': mutant
            })

# Add wild-type control
variants.insert(0, {'ID': 'WT', 'Variant': wild_type})
df = pd.DataFrame(variants)

# 2. Design oligo architecture (standard workflow)
df, _ = op.primer(
    input_data=df,
    oligo_length_limit=200,
    primer_sequence_constraint='SS' + 'N'*18,
    primer_type='forward',
    minimum_melting_temperature=53,
    maximum_melting_temperature=55,
    primer_column='Primer1',
    left_context_column='Variant',
)

df, _ = op.barcode(
    input_data=df,
    oligo_length_limit=200,
    barcode_length=12,
    minimum_hamming_distance=3,
    barcode_column='BC1',
    left_context_column='Primer1',
    right_context_column='Variant',
)

# Continue with primer2, spacers, verify, final...
```

**Key points**:
- Variant generation is your responsibility (`oligopool` designs the oligo architecture)
- For protein saturation mutagenesis, generate codon variants and provide as DNA sequences
- Use `patch_mode=True` to add new variants to an existing saturation library
- Standard analysis pipeline (index → pack → xcount) quantifies variant activity

---

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

### CRISPR Guide Library
Design a barcoded CRISPR guide library for pooled screens.

```python
import oligopool as op
import pandas as pd

# 1. Start with guide sequences (user-provided)
df = pd.DataFrame({
    'ID': ['gene1_g1', 'gene1_g2', 'gene2_g1'],
    'Guide': ['ATGCATGCATGCATGCATGC', 'GCTAGCTAGCTAGCTAGCTA', 'TTAATTAATTAATTAATTAA']
})

# 2. Add scaffold and flanking sequences
# Architecture: U6_promoter | Guide | Scaffold | BC | Primer
df['Scaffold'] = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCG'  # tracrRNA

# 3. Add barcode for guide identification
df, _ = op.barcode(
    input_data=df,
    oligo_length_limit=200,
    barcode_length=12,
    minimum_hamming_distance=3,
    barcode_column='BC',
    left_context_column='Scaffold',
    excluded_motifs=['TTTT'],  # Avoid Pol III terminators
)

# 4. Add primer for amplification/sequencing
df, _ = op.primer(
    input_data=df,
    oligo_length_limit=200,
    primer_sequence_constraint='N'*20,
    primer_type='reverse',
    primer_column='Primer',
    right_context_column='BC',
)

# 5. Finalize
op.verify(input_data=df, oligo_length_limit=200)
final_df, _ = op.final(input_data=df)
```

**Key points**:
- Exclude Pol III terminators (`TTTT`) from barcodes/spacers
- Guide sequences are typically 20bp for SpCas9
- Consider excluding PAM sequences from barcodes if needed

---

### Ribozyme Cleavage Quantification
Quantify self-cleaving ribozyme activity using dual barcodes (cleaved vs uncleaved).

```python
import oligopool as op
import pandas as pd

# Architecture: Primer1 | BC_5prime | Ribozyme | BC_3prime | Primer2
# After cleavage: [Primer1 | BC_5prime | ...] and [... | BC_3prime | Primer2]

# 1. Design library with dual barcodes
df = pd.DataFrame({
    'ID': ['rz1', 'rz2', 'rz3'],
    'Ribozyme': ['ATGC...', 'GCTA...', 'TTAA...']
})

df, _ = op.barcode(input_data=df, barcode_column='BC_5prime', ...)
df, _ = op.barcode(
    input_data=df,
    barcode_column='BC_3prime',
    cross_barcode_columns=['BC_5prime'],  # Ensure separation
    minimum_cross_distance=3,
    ...
)

# 2. After sequencing, count barcode combinations
op.index(barcode_data=df, barcode_column='BC_5prime', index_file='bc5_idx', ...)
op.index(barcode_data=df, barcode_column='BC_3prime', index_file='bc3_idx', ...)

# 3. Combinatorial counting reveals cleavage
xc_df, _ = op.xcount(
    index_files=['bc5_idx', 'bc3_idx'],
    pack_file='reads',
    count_file='cleavage_counts',
)

# Interpretation:
# - Matched pairs (BC_5prime_v1, BC_3prime_v1): Uncleaved ribozyme
# - Unmatched pairs (BC_5prime_v1, '-') or ('-', BC_3prime_v1): Cleaved products
# - Cleavage rate = cleaved / (cleaved + uncleaved)
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
- Use `mapping_type='sensitive'` for sensitive mode
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

## Quick Reference

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

See [api.md#iupac-codes](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#iupac-codes) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for the complete table.

### Common Restriction Sites

See [api.md#common-restriction-sites](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#common-restriction-sites) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for the complete table.

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
    --r1-read-type forward \
    --r2-read-type reverse \
    --pack-type merge \
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

> **CLI Parameter Mapping**: See [api.md#cli-parameter-mapping](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#cli-parameter-mapping) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for the complete mapping.

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
- Barcodes: increase `barcode_length`, decrease `minimum_hamming_distance`, use `barcode_type='terminus'`
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
- README: [README.md](https://github.com/ayaanhossain/oligopool/blob/master/README.md) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/README.md)]
- User Guide: [docs.md](https://github.com/ayaanhossain/oligopool/blob/master/docs/docs.md) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/docs.md)]
- API Reference: [api.md](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)]
- Docker Guide: [docker-notes.md](https://github.com/ayaanhossain/oligopool/blob/master/docs/docker-notes.md) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/docker-notes.md)]
- Paper: https://doi.org/10.1021/acssynbio.4c00661
- CLI help: `op manual <command>`
