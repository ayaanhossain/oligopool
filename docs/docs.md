<h1 align="center">
    <a href="https://github.com/ayaanhossain/oligopool/">
        <img src="https://raw.githubusercontent.com/ayaanhossain/repfmt/main/oligopool/img/logo.svg" alt="Oligopool Calculator" width="460"/>
    </a>
</h1>

---

# Documentation

Welcome to the `Oligopool Calculator` docs. Whether you're designing your first barcode library or optimizing a million-variant MPRA, you're in the right place.

**TL;DR**: Design oligo pool libraries with `barcode`, `primer`, `motif`, `spacer`. Analyze sequencing data with `index`, `pack`, `acount`/`xcount`. Most modules take CSV/DataFrames in and return `(out_df, stats)`, and a few return stats only. Chain them together. Ship it.

**Other resources**:
- [README](../README.md) - Overview, installation, and quick start
- [API Reference](api.md) - Complete parameter documentation for all modules
- [Example Notebook](../examples/OligopoolCalculatorInAction.ipynb) - Interactive design and analysis demo
- [AI Agent Guide](agent-skills.md) - Decision trees, best practices, and gotchas for AI assistants
- [Docker Guide](docker-notes.md) - Run `oligopool` in a container for cross-platform consistency

---

## Table of Contents

- [Quick Start](#quick-start)
- [Core Concepts](#core-concepts)
  - [The DataFrame Flow](#the-dataframe-flow)
  - [The Stats Dictionary](#the-stats-dictionary)
  - [Patch Mode: The Secret Weapon](#patch-mode-the-secret-weapon)
  - [String-Friendly Type Selectors](#string-friendly-type-selectors)
  - [Context Columns and Edge Effects](#context-columns-and-edge-effects)
  - [Reproducibility](#reproducibility)
  - [Verbose Mode](#verbose-mode)
- [Design Mode](#design-mode)
  - [`barcode`](#barcode) - Hamming-distance barcodes
  - [`primer`](#primer) - Thermodynamic primers
  - [`motif`](#motif) - Sequence motifs & anchors
  - [`spacer`](#spacer) - Neutral spacers
  - [`background`](#background) - K-mer screening database
  - [`merge`](#merge) - Collapse columns
  - [`revcomp`](#revcomp) - Reverse complement
  - [`join`](#join) - Join two tables
  - [`final`](#final) - Finalize for synthesis
- [Assembly Mode](#assembly-mode)
  - [`split`](#split) - Fragment long oligos
  - [`pad`](#pad) - Assembly-ready padding
- [Degenerate Mode](#degenerate-mode)
  - [`compress`](#compress) - Compress to IUPAC-degenerate
  - [`expand`](#expand) - Expand degenerate oligos
- [Analysis Mode](#analysis-mode)
  - [`index`](#index) - Build barcode index
  - [`pack`](#pack) - Preprocess FastQ
  - [`acount`](#acount) - Association counting
  - [`xcount`](#xcount) - Combinatorial counting
- [QC Mode](#qc-mode)
  - [`lenstat`](#lenstat) - Length statistics
  - [`verify`](#verify) - Conflict detection
  - [`inspect`](#inspect) - Inspect artifacts
- [Workflows](#workflows)
  - [Basic Library Design](#basic-library-design)
  - [Analysis Pipeline](#analysis-pipeline)
- [CLI Reference](#cli-reference)
  - [CLI-Specific Notes](#cli-specific-notes)
- [Config Files](#config-files)
  - [Why YAML?](#why-yaml)
  - [Single Command Config](#single-command-config)
  - [Pipeline Execution](#pipeline-execution)
  - [Parallel Pipeline Execution](#parallel-pipeline-execution)
  - [Dry Run Validation](#dry-run-validation)
  - [Config Precedence](#config-precedence)
  - [Parameter Sharing](#parameter-sharing)
  - [Config Tips](#config-tips)
- [Tips & Tricks](#tips-tricks)
  - [Design Tips](#design-tips)
  - [Analysis Tips](#analysis-tips)
  - [Performance Tips](#performance-tips)
- [Composability](#composability)
  - [The Decomposition Principle](#the-decomposition-principle)
  - [Design Recipes](#design-recipes)
  - [Constraint Composition](#constraint-composition)
  - [Analysis Recipes](#analysis-recipes)
  - [Degenerate Recipes](#degenerate-recipes)
  - [Application Templates](#application-templates)
  - [Composability Cheat Sheet](#composability-cheat-sheet)
  - [Stuff to Note](#stuff-to-note)
- [Advanced Modules](#advanced-modules)
  - [vectorDB](#vectordb)
  - [Scry](#scry)
- [Troubleshooting](#troubleshooting)
  - [Design Failures](#design-failures)
  - [Assembly Failures](#assembly-failures)
  - [Analysis Issues](#analysis-issues)
  - [General](#general)

---

## Quick Start

[^ Back to TOC](#table-of-contents)

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

[^ Back to TOC](#table-of-contents)

### The DataFrame Flow

[^ Back to TOC](#table-of-contents)

Most modules return both a DataFrame and a stats dictionary:

```python
out_df, stats = op.barcode(input_data=df, ...)
```

A few modules only return stats (no DataFrame to speak of):

```python
stats = op.background(input_data=[...], output_directory='ref_bg')
```

- **Input**: CSV path or pandas DataFrame with an `ID` column
- **Output**:
  - **Design/transform modules** return `(out_df, stats)`
  - **Stats-only modules** return `stats` (`background`, `lenstat`, `inspect`, `index`, `pack`)
  - **Counting modules** return `(counts_df, stats)` (`acount`, `xcount`)
- **Chainable**: Output of one module feeds into the next
- **Placeholders**: `'-'` is a conventional gap value. Some modules treat it as missing in Patch Mode, and `final` ignores it when concatenating.

### The Stats Dictionary

[^ Back to TOC](#table-of-contents)

Every module returns a `stats` dict with:
- `status`: `True` (success) or `False` (failed)
- `basis`: Why it succeeded/failed
- `step`: Numeric step index for the last completed stage
- `step_name`: Human-readable name of the last completed stage
- `vars`: Useful variables from the run
- `warns`: Any warnings you should know about
- `module`: Module name that emitted the stats
- `input_rows`: Number of input rows processed
- `output_rows`: Number of output rows produced

When debugging, start with `status`, `basis`, `step`, and `step_name`, then inspect `vars` for the concrete values that drove the decision.

### Patch Mode: The Secret Weapon

[^ Back to TOC](#table-of-contents)

Extending an existing library? Don't redesign everything:

```python
df, stats = op.barcode(..., patch_mode=True)
```

Patch mode fills only missing values (`None`, `NaN`, `'-'`, empty). Your existing designs stay untouched. Your sanity stays intact.

Concrete example (append rows, fill only the new ones):

```python
import pandas as pd
import oligopool as op

df = pd.DataFrame({
    'ID': ['V1', 'V2'],
    'Variant': ['ATGC' * 10, 'GCTA' * 10],
})

# Initial design
df, _ = op.barcode(
    input_data=df,
    oligo_length_limit=200,
    barcode_length=12,
    minimum_hamming_distance=3,
    maximum_repeat_length=8,
    barcode_column='BC',
    left_context_column='Variant',
)

# If your pipeline keeps ID as the index, normalize back to an ID column:
if 'ID' not in df.columns:
    df = df.reset_index()

# Later: append new variants; mark missing barcodes with '-'
df = pd.concat([df, pd.DataFrame({
    'ID': ['V3'],
    'Variant': ['TTAA' * 10],
    'BC': ['-'],
})], ignore_index=True)

# Patch mode fills only the '-' row; existing BC values are preserved
df, _ = op.barcode(
    input_data=df,
    oligo_length_limit=200,
    barcode_length=12,
    minimum_hamming_distance=3,
    maximum_repeat_length=8,
    barcode_column='BC',
    left_context_column='Variant',
    patch_mode=True,
)
```

Two common orchestration patterns:

- **Pre-seed placeholders**: In scripted pipelines, it's fine to initialize "future" element columns to `'-'` as placeholders, then fill them later with Patch Mode.
- **Drop + redesign**: If you intentionally want to redesign a whole element column, drop it first (otherwise most modules treat it as an output column that already exists).

Note: `'-'` is not meaningful biological context - avoid using placeholder-only columns as `left_context_column`/`right_context_column`.

### String-Friendly Type Selectors

[^ Back to TOC](#table-of-contents)

Many categorical selector parameters accept either integer codes or descriptive strings (case-insensitive).
This includes most `*_type` and `*_policy` parameters, e.g. `barcode_type='terminus'`,
`primer_type='forward'`, `motif_type='anchor'`, `pack_type='merge'`, `mapping_type='sensitive'`,
and `join_policy='left'`. See the [API Reference](api.md) for the full alias list.

### Context Columns and Edge Effects

[^ Back to TOC](#table-of-contents)

Most design modules need to know what's next to the element being designed:

```python
df, stats = op.barcode(
    ...,
    left_context_column='Promoter',   # What's on the left?
    right_context_column='Gene'       # What's on the right?
)
```

For the main element-design modules (`barcode`, `primer`, `motif`, `spacer`), you need at least one of
`left_context_column` or `right_context_column`.

When both are provided, newly designed columns are inserted between those context columns in the DataFrame.
This means dependency/context specification also resolves final linear column ordering for downstream `merge`.

**Edge effects** happen when an undesired sequence (excluded motif, repeat) appears at an insertion boundary. Example: inserting barcode `AATTC` next to context ending in `...G` creates `...GAATTC`, which contains the EcoRI site `GAATTC` across the junction. Context columns let the algorithm check for this and avoid it.

**Background screening** is the same idea, but against a whole reference (genome, transcriptome, plasmid backbone, ...). Build one or more k-mer DBs with `background()`. Then pass `background_directory` (single path or list of paths) to `barcode`/`primer`/`motif`/`spacer`, or `verify` for QC. Screening is junction-aware when context columns are provided.

### Reproducibility

[^ Back to TOC](#table-of-contents)

All stochastic design modules support `random_seed` for reproducible results:

```python
df, stats = op.barcode(..., random_seed=42)
```

Same seed + same inputs = same outputs. Great for debugging and publications.

### Verbose Mode

[^ Back to TOC](#table-of-contents)

Control output verbosity with `verbose` (default: `True`):

```python
df, stats = op.barcode(..., verbose=False)  # Silent mode
```

---

## Design Mode

[^ Back to TOC](#table-of-contents)

Design Mode is where you build your library piece by piece. Think of it like molecular LEGO, but with more constraints and fewer stepping-on-brick injuries.

**When to use Design Mode**: You're building a new oligo pool library from scratch, or extending an existing one. Design Mode handles molecular constraints (Hamming distance for barcodes, primer thermodynamics, excluded-motif screening, and `oligo_length_limit`) so you can focus on the biology.

The typical workflow:
1. Start with your core sequences (variants, promoters, genes, etc.)
2. Add functional elements: `primer` → `motif` → `barcode` → `spacer`
3. Validate in QC Mode with `lenstat` and `verify`
4. Finalize with `final` to get synthesis-ready oligos

### `barcode`

[^ Back to TOC](#table-of-contents)

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

**Column placement behavior sketch:**
```text
Input:  Left, Core, Right
barcode(..., barcode_column='BC', left_context_column='Left', right_context_column='Core')
Output: Left, BC, Core, Right
```

**Excluded motifs** let you ban restriction sites, repetitive sequences, or anything else you don't want appearing in your barcodes.

The simplest approach is to pass a list:

```python
df, stats = op.barcode(..., excluded_motifs=['GAATTC', 'GGATCC'])
```

You can also pass a CSV path, a DataFrame with an `Exmotif` column, or a FASTA file. To use
multiple motif sets, pass a list of sources or a `{name: source}` dict:

```python
df, stats = op.barcode(..., excluded_motifs=['cutsites.csv', 'homopolymers.csv'])
df, stats = op.barcode(..., excluded_motifs={'cutsites': 'cutsites.csv', 'homo': ['AAAA', 'TTTT']})
```

All sources are merged; motifs must be strict ATGC (no IUPAC, no dashes). See the
[API Reference](api.md#barcode) for details.

> **Note**: Excluded motifs are applied globally to all variants (no per-variant exclusion).

**Stuff to Note:**
- `terminus` vs `spectrum`: `terminus` optimized barcodes enforce distinctive 5'/3' ends; `spectrum` optimized barcodes saturate k-mers (slower but tighter). Higher `minimum_hamming_distance` buys you error tolerance but shrinks the design space.
- If you plan to index barcodes from reads, design constant anchors first with `motif(motif_type=1, ...)`.
- `patch_mode` preserves existing barcodes; any pre-existing values must be valid ATGC strings of length `barcode_length`.
- Cross-set separation is strict: set `cross_barcode_columns` and `minimum_cross_distance` together; the cross-set barcodes must be strict ATGC strings of length `barcode_length`.
- `background_directory` screens designed barcodes against one or more background k-mer DBs (junction-aware when context columns are provided).

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

> **API Reference**: See [`barcode`](api.md#barcode) for complete parameter documentation.

---

### `primer`

[^ Back to TOC](#table-of-contents)

**What it does**: Designs primers with thermodynamic and structural constraints. That includes Tm targeting, repeat screening, anti-structure checks (hairpin, homodimer, heterodimer, cross-primer dimer), and optional background filtering.

**When to use it**: You need ultra-high-quality amplification primers, sequencing primers, or assembly overlap regions.

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

**Column placement behavior sketch:**
```text
Input:  Left, Core
primer(..., primer_column='FwdPrimer', right_context_column='Core')
Output: Left, FwdPrimer, Core
```

**Paired primer design** ensures your forward and reverse primers have matched Tm (within 1 °C). Design the inner primer first:

```python
df, stats = op.primer(..., primer_column='FwdPrimer', primer_type=0)
```

Then design the outer primer, telling it to match Tm with the first while maintaining favorable thermodynamic coupling:

```python
df, stats = op.primer(
    ...,
    primer_column='RevPrimer',
    primer_type=1,
    paired_primer_column='FwdPrimer',
)
```

**Per-set primers** let you design different primers for different groups in your pool. Pass group labels as a list (aligned to input rows):

```python
df, stats = op.primer(
    ...,
    oligo_set=['SetA', 'SetA', 'SetB', 'SetB'],
)
```

You can also pass a CSV path or a DataFrame with `ID` and `OligoSet` columns. See the
[API Reference](api.md#primer) for details.

**Background screening** keeps your primers from binding to off-target sequences (host genome, plasmid backbone, etc.). First, build a background database:

```python
op.background(
    input_data='genome.fasta',
    maximum_repeat_length=15,
    output_directory='my_background'
)
```

Then point to it during primer design:

```python
df, stats = op.primer(
    ...,
    background_directory=['my_background', 'my_other_background'],
)
```

**Stuff to Note:**
- `maximum_repeat_length` controls non-repetitiveness against `input_data` only; screening against a background requires `background_directory`.
- If `paired_primer_column` is provided, the paired primer type is inferred and Tm matching is applied within 1 °C.
- When `oligo_set` is provided, primers are designed per set and screened for cross-set compatibility; if `paired_primer_column` is also provided, it must be constant within each set.
- `primer_sequence_constraint` is flexible by design; mix fixed and degenerate bases to enforce patterns (for example, GC clamp-like starts such as `'SS' + 'N'*18`).
- Structural screening is first-class: candidates with strong hairpin/homodimer/heterodimer behavior (including cross-primer dimer risks when relevant) are rejected during optimization.
- `patch_mode` preserves existing primers; in `oligo_set` mode, existing per-set primers are reused and only missing sets trigger new primer design.

> **API Reference**: See [`primer`](api.md#primer) for complete parameter documentation.

---

### `motif`

[^ Back to TOC](#table-of-contents)

**What it does**: Inserts sequence motifs (per-variant or constant anchors) with constraint satisfaction.

**When to use it**:
- Constant anchors for indexing (barcode prefix/suffix)
- Regulatory elements with sequence constraints
- Any fixed or semi-fixed sequence insertion

For per-variant motifs (different sequence for each row), use `motif_type=0`:

```python
df, stats = op.motif(
    input_data=df,
    oligo_length_limit=200,
    motif_sequence_constraint='NNNNNNNN',
    maximum_repeat_length=6,
    motif_column='Spacer',
    motif_type=0,
    left_context_column='BC1',
)
```

For constant anchors (same sequence for all rows - great for indexing!), use `motif_type=1`:

```python
df, stats = op.motif(
    ...,
    motif_type=1,
    motif_column='BC1_Prefix',
)
```

**Column placement behavior sketch:**
```text
Input:  Left, Core
motif(..., motif_column='Anchor', left_context_column='Left')
Output: Left, Anchor, Core
```

**Why anchors matter**: When you index barcodes later, you need constant flanking sequences to locate them in reads. Design anchors with `motif_type=1` *before* designing barcodes.

**Stuff to Note:**
- Constant bases in `motif_sequence_constraint` can accidentally force an excluded motif (or repeat) and make the design infeasible.
- For anchors, tune `maximum_repeat_length` to control how distinct the anchor is from surrounding sequence.
- `patch_mode` fills only missing values; for `motif_type=1`, an existing compatible constant anchor is reused for new rows.
- `background_directory` screens designed motifs/anchors against one or more background k-mer DBs (junction-aware when context columns are provided).

> **API Reference**: See [`motif`](api.md#motif) for complete parameter documentation.

---

### `spacer`

[^ Back to TOC](#table-of-contents)

**What it does**: Inserts neutral DNA spacers to meet length requirements.

**When to use it**: Your oligos need to hit a specific length (e.g., synthesis array constraints).

Need a fixed-length spacer? Just say how long:

```python
df, stats = op.spacer(
    input_data=df,
    oligo_length_limit=200,
    maximum_repeat_length=8,
    spacer_column='Spacer3',
    spacer_length=15,
    left_context_column='Gene',
)
```

**Column placement behavior sketch:**
```text
Input:  Left, Core
spacer(..., spacer_column='Spacer', left_context_column='Left')
Output: Left, Spacer, Core
```

Want spacers to auto-fill whatever space remains? Leave `spacer_length` as `None`:

```python
df, stats = op.spacer(..., spacer_length=None)
```

In auto-fill mode, the spacer length is computed per row as the remaining free space under your
`oligo_length_limit` after concatenating all existing sequence columns (gaps `'-'` are ignored). Rows
with no free space get a spacer of `'-'`.

Need different lengths for different variants? Pass a list (aligned to input rows):

```python
df, stats = op.spacer(..., spacer_length=[10, 15, 12, 8])
```

You can also pass a CSV path or a DataFrame with `ID` and `SpacerLength` columns. See the
[API Reference](api.md#spacer) for details.

**Stuff to Note:**
- If a row is already at (or over) `oligo_length_limit`, spacer becomes `'-'` (an intentional "no-space" sentinel).
- If `spacer_length` comes from CSV/DataFrame, include `ID` and `SpacerLength` columns aligned to your input IDs.
- `patch_mode` fills only missing values and never overwrites existing spacers (some rows may still end up with `'-'` if no spacer can fit).
- `background_directory` screens designed spacers against one or more background k-mer DBs (useful when spacers must avoid matching a reference).

> **API Reference**: See [`spacer`](api.md#spacer) for complete parameter documentation.

---

### `background`

[^ Back to TOC](#table-of-contents)

**What it does**: Builds a k-mer database for screening designed sequences (primers, barcodes, motifs, spacers) against off-target k-mers.

**When to use it**: You want primers that won't bind to host genome, plasmid backbone, etc.

Build a background from a FASTA file (genome, plasmid, whatever you want to avoid):

```python
stats = op.background(
    input_data='ecoli_genome.fasta',
    maximum_repeat_length=15,
    output_directory='ecoli_bg',
)
```

Works with `.fa`, `.fna`, and gzipped files too. Prefer CSV? Just include a `Sequence` column:

```python
stats = op.background(
    input_data='plasmid_sequences.csv',
    maximum_repeat_length=15,
    output_directory='plasmid_bg',
)
```

Then use your background database during design (or QC):

```python
df, stats = op.primer(..., background_directory=['ecoli_bg', 'plasmid_bg'])
```

**Stuff to Note:**
- The background output directory ends with `.oligopool.background` and is what you pass as `background_directory` (single path or list of paths).
- `background(maximum_repeat_length=...)` screens against the background only; `maximum_repeat_length` in design modules screens against your input oligos.
- If you need to inspect or modify the DB, use `vectorDB` (see the Advanced Modules section).

> **API Reference**: See [`background`](api.md#background) for complete parameter documentation.

---

### `merge`

[^ Back to TOC](#table-of-contents)

**What it does**: Concatenates contiguous columns into a single column.

**When to use it**: Mid-pipeline cleanup, combining elements before further processing.

```python
df, stats = op.merge(
    input_data=df,
    merge_column='5prime_region',         # Output column name
    left_context_column='Primer5',
    right_context_column='BC1',           # Merges Primer5 through BC1
)
```

**Column order behavior sketch:**
```text
Input:  A, B, C, D
merge(B..C -> X)
Output: A, X, D
```

**Stuff to Note:**
- `left_context_column` and `right_context_column` do not need to be adjacent; `merge` collapses the full inclusive range.
- If you omit the bounds, `merge` collapses the first → last sequence columns.
- The merged source columns are removed from the output DataFrame.

> **API Reference**: See [`merge`](api.md#merge) for complete parameter documentation.

---

### `revcomp`

[^ Back to TOC](#table-of-contents)

**What it does**: Reverse complements a range of columns and reverses their order.

**When to use it**: Switching strand orientation (e.g., synthesis vs. sequencing direction).

```python
df, stats = op.revcomp(
    input_data=df,
    left_context_column='Gene',
    right_context_column='Primer3',       # RevComp everything from Gene to Primer3 (inclusive)
)
```

**Column order behavior sketch:**
```text
Input:  A, B, C, D
revcomp(B..C)
Output: A, C', B', D   (sequences are reverse-complemented; order is reversed)
```

**Stuff to Note:**
- `left_context_column` and `right_context_column` do not need to be adjacent; `revcomp` acts on the full inclusive range.
- If you omit the bounds, `revcomp` reverse-complements the first → last sequence columns and reverses their order.
- Useful mid-pipeline when you design in "readout orientation" but must synthesize in the opposite orientation (and for sanity-checking `split` fragment orientation).

> **API Reference**: See [`revcomp`](api.md#revcomp) for complete parameter documentation.

---

### `join`

[^ Back to TOC](#table-of-contents)

**What it does**: Joins two DataFrames/CSVs on `ID`, keeps `input_data` column order as backbone, and inserts only new columns from `other_data`.

**When to use it**: Recombining parallel design branches (for example, from a YAML CLI DAG) back into one design table, or reconciling two independently designed DataFrames.

```python
df, stats = op.join(
    input_data=df_backbone,
    other_data=df_branch,
    oligo_length_limit=200,
    join_policy='left',                # or 1 / 'right'
)
```

**Column order behavior sketch:**
```text
Input: A, B, C, D, E
Other: A, B, X, D, E

join_policy='left'  -> A, B, X, C, D, E
join_policy='right' -> A, B, C, X, D, E
```

Another common pattern: if `Input` has `Left, Core, Right` and `Other` has `Left, Spacer, Core, Right`,
then `Spacer` is inserted between `Left` and `Core` (this is unambiguous, so `join_policy` does not matter).

**Stuff to Note:**
- `join` is an inter-table operation (contrast with `merge`/`revcomp`, which operate within one DataFrame).
- `input_data` and `other_data` must contain the same `ID` set; `join` never creates or drops rows (ID mismatches are an error).
- Overlapping non-`ID` column names from `other_data` must match `input_data` exactly; they are verified then ignored (backbone is preserved).
- New columns from `other_data` are inserted near their nearest shared anchors; `join_policy` resolves ambiguous placements: `0`/`left` is left-biased, `1`/`right` is right-biased.
- Post-join QC: if any joined oligo exceeds `oligo_length_limit`, `join` is infeasible and returns no output table.
- For a CLI/YAML example, see [Parallel Pipeline Execution](#parallel-pipeline-execution) in the Config Files section.

> **API Reference**: See [`join`](api.md#join) for complete parameter documentation.

---

### `final`

[^ Back to TOC](#table-of-contents)

**What it does**: Concatenates all columns into synthesis-ready oligos.

**When to use it**: Last step before ordering.

```python
final_df, stats = op.final(
    input_data=df,
    output_file='synthesis_ready',
)
# Output contains 'CompleteOligo' and 'OligoLength' columns
```

**Output columns behavior sketch:**
```text
Input:  A, B, C
final()
Output: CompleteOligo, OligoLength   (ID is preserved)
```

**Stuff to Note:**
- `final` gives you only `CompleteOligo` and `OligoLength` (annotation columns are dropped).
- Save your annotated design DataFrame/CSV before `final` if you plan to `index`, `pack`, and count later.
- In most workflows, `final` is your last design step before synthesis.

> **API Reference**: See [`final`](api.md#final) for complete parameter documentation.

---

## Assembly Mode

[^ Back to TOC](#table-of-contents)

Assembly Mode provides tools for fragmenting long oligos that exceed synthesis length limits into overlapping pieces for assembly workflows.

**When to use Assembly Mode**: Your designed oligos exceed synthesis length limits (typically ~200-300 bp, vendor-dependent). Assembly Mode splits them into overlapping fragments that can be synthesized separately, then reassembled post-synthesis via Gibson assembly, overlap-extension PCR, or similar methods.

The typical workflow:
1. Design your full-length library with Design Mode
2. Use `split` to fragment oligos into overlapping pieces
3. Use `pad` to add Type IIS primer sites for each fragment
4. Synthesize fragments separately, then assemble via overlaps

### `split`

[^ Back to TOC](#table-of-contents)

**What it does**: Breaks long oligos into overlapping fragments for overlap-based assembly (Gibson, overlap-extension PCR, etc.).

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

**Output columns behavior sketch:**
```text
Input:  (annotated design table)
split()
Output: ID, Split1, Split2, ...   (annotation columns are dropped)
```

**Key concept: Each column = one oligo pool to synthesize**

When `split` returns `Split1`, `Split2`, `Split3` columns, you will **order three separate oligo pools** from your synthesis vendor (one per column). After synthesis, you combine those fragments via overlap-based assembly (Gibson, overlap-extension PCR, etc.) to reconstruct the full-length oligo.

**Stuff to Note:**
- The number of fragments is automatically chosen and can vary per oligo (some rows may not have all `Split*` columns populated).
- Fragments are returned in PCR/assembly order. Even-numbered split columns (`Split2`, `Split4`, ...) are reverse-complemented by design. This alternating orientation supports efficient PCR-based assembly where adjacent fragments anneal through overlaps.
- `split` returns only fragment columns (your original annotation columns are dropped). Keep a saved annotated library if you will need it later.
- As a rule of thumb, keep `minimum_overlap_length` comfortably larger than `minimum_hamming_distance`.
- **Raw split output is not assembly-ready yet**: run `pad` per `SplitN` to add primers/Type IIS sites, then run `final` on each padded DataFrame.
- For separate per-fragment outputs: in Python, enable `separate_outputs`. In CLI, this is already the default (writes `out.Split1.oligopool.split.csv`, etc.). Use `--no-separate-outputs` in CLI for a single combined file.

> **API Reference**: See [`split`](api.md#split) for complete parameter documentation.

---

### `pad`

[^ Back to TOC](#table-of-contents)

**What it does**: Adds primers and Type IIS restriction sites to split fragments for scarless assembly.

**When to use it**: After `split`, to make fragments synthesis-ready and assembly-compatible.

**Run `pad` once per split column** - each call produces a synthesis-ready DataFrame for that fragment:

```python
# After split returns Split1, Split2, Split3 columns:
split_df, _ = op.split(input_data=df, ...)

# Pad each fragment separately
pad1_df, _ = op.pad(input_data=split_df, split_column='Split1', typeIIS_system='BsaI', ...)
pad2_df, _ = op.pad(input_data=split_df, split_column='Split2', typeIIS_system='BsaI', ...)
pad3_df, _ = op.pad(input_data=split_df, split_column='Split3', typeIIS_system='BsaI', ...)

# Finalize each for synthesis
final1_df, _ = op.final(input_data=pad1_df, output_file='synthesis_split1')
final2_df, _ = op.final(input_data=pad2_df, output_file='synthesis_split2')
final3_df, _ = op.final(input_data=pad3_df, output_file='synthesis_split3')

# You now have 3 CSV files to send to your synthesis vendor
```

**Output columns behavior sketch:**
```text
Input:  ID, Split1, Split2, ... (from split)
pad(split_column='Split1')
Output: ID, 5primeSpacer, ForwardPrimer, Split1, ReversePrimer, 3primeSpacer
```

**Stuff to Note:**
- **The chosen Type IIS recognition site must be absent from your split fragments** (in either orientation). If fragments contain internal cut sites, they will be cleaved during digest. `pad` checks this and fails early if conflicts are found.
- **Exclude your Type IIS motif from upstream design elements.** If you plan to use `BsaI` for padding, add its recognition site (`GGTCTC`) and reverse complement (`GAGACC`) to `excluded_motifs` when designing primers, barcodes, motifs, and spacers. This prevents new internal cut sites in designed elements. It does not clean a core/variant sequence that already contains the site; in that case, pick another enzyme or redesign those sequences.
- Output columns are `5primeSpacer`, `ForwardPrimer`, `<split_column>`, `ReversePrimer`, `3primeSpacer` (other `split_df` columns are not preserved).
- If a fragment can't fit under `oligo_length_limit`, the spacer(s) for that row are set to `'-'` (a deliberate "no-space" sentinel).
- The supported Type IIS list is the "batteries included" set for `pad`, not an arbitrary list. These 34 enzymes are included because `pad` has a clear recognition motif plus 3' cut-offset model for each one (for example, BsaI: `GGTCTC` + `N*5`). That gives predictable excision and scarless primer-pad removal (optionally followed by blunting). Upstream/complex cutters are not included.
- Choose an enzyme whose recognition site is absent from your fragments and that matches your downstream handling preferences (e.g., sticky-end cutters vs blunt cutters like `MlyI`).

**Post-synthesis workflow** (after you get your oligos back):
1. **PCR amplify** the padded fragments
2. **Type IIS digest** → excises primers and spacers, leaving enzyme-specific sticky overhangs (the overhangs are primer-derived, not from your fragment)
3. **Mung bean nuclease** (optional) → blunts the overhangs; skip this step if using a blunt-cutter like `MlyI`
4. **Assemble** the cleaned fragments via their **split-designed overlaps** (Gibson, overlap-extension PCR, etc.)

The Type IIS enzyme is for **pad removal**, not fragment-to-fragment ligation - the 15-30 bp overlaps from `split` drive the actual assembly.

**Supported Type IIS enzymes (34 total):**
`AcuI`, `AlwI`, `BbsI`, `BccI`, `BceAI`, `BciVI`, `BcoDI`, `BmrI`, `BpuEI`, `BsaI`, `BseRI`, `BsmAI`, `BsmBI`, `BsmFI`, `BsmI`, `BspCNI`, `BspQI`, `BsrDI`, `BsrI`, `BtgZI`, `BtsCI`, `BtsI`, `BtsIMutI`, `EarI`, `EciI`, `Esp3I`, `FauI`, `HgaI`, `HphI`, `HpyAV`, `MlyI`, `MnlI`, `SapI`, `SfaNI`

For other enzymes, build custom pads from primitives with `primer` + `motif` + `spacer` (i.e., reproduce the same primitive composition that `pad` automates for built-in Type IIS systems).

> **API Reference**: See [`pad`](api.md#pad) for complete parameter documentation.

---

## Degenerate Mode

[^ Back to TOC](#table-of-contents)

Degenerate Mode enables cost-efficient synthesis of variant libraries with low mutational diversity. Instead of synthesizing each variant individually, similar sequences are compressed into IUPAC-degenerate oligos that expand to cover multiple variants.

**When to use Degenerate Mode**: You have a large library of similar sequences and want to reduce synthesis costs by grouping compatible variants. It works best for selection assays where you identify enriched variants by sequencing (no barcode-based readout required). ML-generated libraries and saturation mutagenesis libraries often compress well.

This enables **selection-based discovery workflows**:
1. Design variants computationally (e.g., promoters, UTRs)
2. Compress similar variants into IUPAC-degenerate oligos (`compress`)
3. Synthesize, clone, and run selection (growth, FACS, survival screens)
4. Sequence survivors to identify enriched variants
5. Map sequences back to original variant IDs using `mapping_df`

### `compress`

[^ Back to TOC](#table-of-contents)

**What it does**: Compresses concrete DNA sequences into IUPAC-degenerate representation.

**When to use it**: You have many similar sequences and want fewer oligos to synthesize.

```python
mapping_df, synthesis_df, stats = op.compress(
    input_data=df,                        # DataFrame with ID + DNA columns
    # Any basename/path is fine; `compress` appends distinct suffixes for each output.
    mapping_file='degenerate_library',    # Variant-to-degenerate map
    synthesis_file='degenerate_library',  # Degenerate oligos for ordering
    random_seed=42,
)

# Check compression results
print(f"Compressed {stats['vars']['input_variants']} variants "
      f"into {stats['vars']['degenerate_oligos']} oligos "
      f"({stats['vars']['compression_ratio']}x compression)")
```

**Output columns behavior sketch:**
```text
compress()
mapping_df:   ID, Sequence, DegenerateID
synthesis_df: DegenerateID, DegenerateSeq, Degeneracy, OligoLength
```

**Stuff to Note:**
- Input must be concrete DNA (A/T/G/C only, no degenerate codes)
- All non-ID columns are concatenated (left-to-right) to form the full sequence
- Include only columns you intend to make degenerate (e.g., diverse barcodes can prevent good compression)
- Similar sequences compress well; diverse sequences may not compress at all (1:1 mapping)
- The algorithm guarantees no invented sequences (lossless compression)

> **API Reference**: See [`compress`](api.md#compress) for complete parameter documentation.

---

### `expand`

[^ Back to TOC](#table-of-contents)

**What it does**: Expands IUPAC-degenerate sequences into all concrete A/T/G/C sequences.

**When to use it**: Verifying that `compress` output covers all your original variants.

```python
# Expand the synthesis output to verify coverage
expanded_df, stats = op.expand(
    input_data=synthesis_df,
    sequence_column='DegenerateSeq',
)

# Check all original sequences are recovered
original_seqs = set(df['Sequence'])
expanded_seqs = set(expanded_df['ExpandedSeq'])
assert original_seqs == expanded_seqs, "Lossless guarantee violated!"
```

**Output columns behavior sketch:**
```text
expand()
Output: ExpandedSeq, OligoLength   (+ ID and optionally DegenerateID if mapping_file is provided)
```

**Stuff to Note:**
- Expansion can be exponential (N positions = 4^N for all-N sequence)
- Use `expansion_limit` as a safety cap for highly degenerate sequences
- Does NOT recover original variant IDs; use `mapping_df` from `compress` for that

> **API Reference**: See [`expand`](api.md#expand) for complete parameter documentation.

---

## Analysis Mode

[^ Back to TOC](#table-of-contents)

Analysis Mode is where sequencing data meets your designed library. You've done the experiment, now let's count some barcodes.

**When to use Analysis Mode**: You've run your experiment (MPRA, CRISPR screen, etc.) and have FastQ files ready for barcode quantification. Analysis Mode handles the pipeline from raw reads to count matrices.

The typical workflow:
1. Build one index per barcode column with `index` (multi-barcode designs use multiple index files)
2. Preprocess FastQ files with `pack` (length/quality filtering, optional paired-end merging, deduplication)
3. Count associations with `acount` (barcode <-> associate) or `xcount` (combinatorial)
4. Export count matrices for downstream analysis

### `index`

[^ Back to TOC](#table-of-contents)

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
    associate_suffix_gap=0,                # Bases between variant and suffix anchor in read
)
```

**Stuff to Note:**
- You need at least one of `barcode_prefix_column` or `barcode_suffix_column`; anchors should be constant (single-unique) sequences and ideally adjacent to the barcode/associate.
- `barcode_prefix_gap`/`barcode_suffix_gap` specify how many bases separate the anchor and barcode in the read (real sequencing is messy; this makes it configurable).
- If your associate anchors are not directly adjacent, use `associate_prefix_gap`/`associate_suffix_gap` to specify the exact gaps in the read.
- For association counting, partial presence of the associate sequence can be sufficient, but the associate anchors must be adjacent and fully present.
- You can use multiple index files with `xcount` for combinatorial counting (associate info is ignored in that mode).

> **API Reference**: See [`index`](api.md#index) for complete parameter documentation.

---

### `pack`

[^ Back to TOC](#table-of-contents)

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
    pack_type=1,                           # 0=concatenate pairs, 1=merge pairs
    pack_file='sample',                    # Creates sample.oligopool.pack
)
```

**Stuff to Note:**
- For paired-end workflows, both reads must pass the length/quality filters for the pair to be retained.
- If reads are merged externally, pass the merged reads as single-end (R1 only) and leave all R2 args as `None`.
- `pack_type=0` (concatenated) tends to be IO-bound; `pack_type=1` (merged) is more compute-heavy. Pack files are reusable across counting runs.
- Internally, pack files store reads as `(r1, r2)` tuples; merged/single-end reads use `r2=None` (relevant for Python callbacks).

> **API Reference**: See [`pack`](api.md#pack) for complete parameter documentation.

---

### `acount`

[^ Back to TOC](#table-of-contents)

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

**Output columns behavior sketch:**
```text
acount()
Output columns: <index_name>.ID, BarcodeCount, AssociationCount
```

**Stuff to Note:**
- Reads with unresolved associates are excluded from the counts (that's the point of `acount`).
- `callback` is Python-only; the CLI always runs with `callback=None`.
- If you want to inspect *examples* of discarded reads, set `failed_reads_file` to write a small per-category diagnostic CSV.
- `acount` operates on a single index + pack pair (for combinatorial counting, use `xcount`).

> **API Reference**: See [`acount`](api.md#acount) for complete parameter documentation.

---

### `xcount`

[^ Back to TOC](#table-of-contents)

**What it does**: Barcode-focused counting - maps reads to one or more barcode indices without associate verification.

**When to use it**: Pure barcode counting (single or multiple indices), combinatorial barcode designs (BC1 x BC2).

For single barcode counting, point to one index:

```python
df, stats = op.xcount(
    index_files=['bc1_index'],
    pack_file='sample',
    count_file='barcode_counts',
)
```

For combinatorial counting (BC1 x BC2), list multiple indices:

```python
df, stats = op.xcount(
    index_files=['bc1_index', 'bc2_index'],
    pack_file='sample',
    count_file='combo_counts',
    mapping_type=1,
)
```

**Output columns behavior sketch:**
```text
xcount()
Output columns: <index1_name>.ID, <index2_name>.ID, ..., CombinatorialCount
```

Output includes all observed combinations. Reads missing a barcode show `'-'` for that position.

**`acount` vs `xcount`**: Use `acount` when you need barcode+variant association verification; use `xcount` for barcode-only counting (single or combinatorial).

**Stuff to Note:**
- Reads are retained if at least one barcode maps; missing barcodes are represented as `'-'` in the output combination.
- `callback` is Python-only; the CLI always runs with `callback=None`.
- If you want to inspect *examples* of discarded reads, set `failed_reads_file` to write a small per-category diagnostic CSV.
- Associate information in the index is ignored in `xcount` mode.

**Custom callbacks** (Python API only):
```python
def my_filter(r1, r2, ID, count, coreid):
    # r2 is None for merged reads; otherwise you get (r1, r2)
    return True  # Accept read(s)

df, stats = op.xcount(..., callback=my_filter)
```

> **API Reference**: See [`xcount`](api.md#xcount) for complete parameter documentation.

---

## QC Mode

[^ Back to TOC](#table-of-contents)

QC Mode provides utilities for validating designed libraries and inspecting non-CSV artifacts before synthesis or downstream analysis.

The typical workflow:
1. `lenstat` to check length budget and free space
2. `verify` to detect length, excluded-motif emergence, and background conflicts
3. `inspect` to sanity-check background/index/pack artifacts

### `lenstat`

[^ Back to TOC](#table-of-contents)

**What it does**: Reports length statistics and free space remaining.

**When to use it**: Check progress during design, ensure you're within synthesis limits.

```python
stats = op.lenstat(
    input_data=df,
    oligo_length_limit=200,
)
# Prints a nice table of per-column and total lengths
```

**Stuff to Note:**
- `lenstat` is stats-only: it won't modify your DataFrame and it doesn't write an `output_file`.
- It assumes all non-ID columns are DNA strings; if you have metadata columns or degenerate/IUPAC bases, use `verify` instead.
- Run it early and often-especially before `spacer`, `split`, and `pad`.

> **API Reference**: See [`lenstat`](api.md#lenstat) for complete parameter documentation.

---

### `verify`

[^ Back to TOC](#table-of-contents)

**What it does**: Detects integrity, length, motif emergence, and background k-mer conflicts in your oligo pool.

**When to use it**: QC check after design, before ordering synthesis.

```python
df, stats = op.verify(
    input_data=df,
    oligo_length_limit=200,
    excluded_motifs=['GAATTC', 'GGATCC'],  # Restriction sites to flag
    output_file='verify_results',
)
```

**Output columns behavior sketch:**
```text
verify()
Output columns:
  - CompleteOligo, OligoLength
  - HasIntegrityConflict, HasLengthConflict, HasExmotifConflict, HasBackgroundConflict, HasAnyConflicts
  - IntegrityConflictDetails, LengthConflictDetails, ExmotifConflictDetails, BackgroundConflictDetails

(When written to CSV, the *Details columns are JSON-serialized values; parse them back with json.loads().)
```

**Checks:**
- **Integrity conflicts**: `CompleteOligo` vs constituent column mismatch or `OligoLength` vs actual length mismatch
- **Length conflicts**: Oligos exceeding `oligo_length_limit`
- **Excluded motif conflicts**: Motif emergence (count exceeds library-wide baseline)
- **Background conflicts**: K-mer matches in background DB(s)

**Output DataFrame columns:**
- `HasIntegrityConflict`, `HasLengthConflict`, `HasExmotifConflict`, `HasBackgroundConflict`: Boolean flags
- `HasAnyConflicts`: Combined OR of all four flags
- `IntegrityConflictDetails`, `LengthConflictDetails`, `ExmotifConflictDetails`, `BackgroundConflictDetails`: Dict with violation details (or None)

**How columns are concatenated:**
- Uses `CompleteOligo` if present; otherwise concatenates all **pure ATGC columns** left-to-right
- IUPAC/degenerate columns are skipped silently
- Gap characters (`'-'`) are stripped

**Reading conflict details from CSV:**
```python
import json
import pandas as pd

df = pd.read_csv('verify_results.oligopool.verify.csv')
# Parse JSON-serialized dicts in the *Details columns
for col in ['IntegrityConflictDetails', 'LengthConflictDetails', 'ExmotifConflictDetails', 'BackgroundConflictDetails']:
    df[col] = df[col].apply(lambda x: json.loads(x) if pd.notna(x) else None)
```

**Stuff to Note:**
- `verify` returns `(DataFrame, stats)` like other design/transform modules.
- `oligo_length_limit` is mandatory for `verify`.
- Motif **emergence** = count exceeds library-wide minimum; flagged even if baseline >= 1.
- Motif matching is literal substring matching; IUPAC bases are not expanded as wildcards.
- Rows with a `CompleteOligo`-vs-constituent mismatch attribute exmotif/background hits to `CompleteOligo` instead of constituent columns.

> **API Reference**: See [`verify`](api.md#verify) for complete parameter documentation.

---

### `inspect`

[^ Back to TOC](#table-of-contents)

`inspect` is a read-only, stats-only utility that pulls key facts from complex non-CSV artifacts produced by `Oligopool Calculator`.
Use `inspect` when generation-step stats are unavailable, or when you want to quickly confirm what a background/index/pack artifact contains before reusing it.

- Accepted targets: `.oligopool.background` directories, `.oligopool.index` files, `.oligopool.pack` files.
- `kind='auto'` (default) infers artifact type from the target path.
- No output file is written; the module returns `stats` only.
- `stats['vars']['verdict']` is `Valid`, `Corrupted`, or `Invalid`.
- `stats['status']=True` only when verdict is `Valid`.
- Artifact summaries are reported under `stats['vars']['meta']` (for example k-mer counts for backgrounds or packing stats for packs).

```python
stats = op.inspect(target='demo.oligopool.background')
```

```bash
op inspect --target demo.oligopool.background --stats-json --quiet
```

> **API Reference**: See [`inspect`](api.md#inspect) for complete parameter documentation.

---

## Workflows

[^ Back to TOC](#table-of-contents)

### Basic Library Design

[^ Back to TOC](#table-of-contents)

Start with your variants in a CSV (must have an `ID` column):

```python
import pandas as pd
import oligopool as op

df = pd.read_csv('variants.csv')
```

Add a forward primer. It'll sit to the left of your variants:

```python
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
```

Add barcodes between the primer and variants:

```python
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
```

Add a reverse primer on the other side, Tm-matched to the forward:

```python
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
```

Check that everything fits within your length budget:

```python
op.lenstat(input_data=df, oligo_length_limit=200)
```

Run QC, save your annotated design, and finalize for synthesis:

```python
verify_df, verify_stats = op.verify(input_data=df, oligo_length_limit=200)
print(f"Conflicts: {verify_stats['vars']['any_conflict']}")

df.to_csv('library_design.csv', index=False)  # Keep this for indexing later!
final_df, _ = op.final(input_data=df, output_file='library_for_synthesis')
```

### Analysis Pipeline

[^ Back to TOC](#table-of-contents)

First, index your barcodes. The anchors (constant flanking sequences) help locate barcodes in reads:

```python
import oligopool as op

op.index(
    barcode_data='library_design.csv',
    barcode_column='BC1',
    barcode_prefix_column='FwdPrimer',
    associate_data='library_design.csv',
    associate_column='Variant',
    associate_suffix_column='RevPrimer',
    index_file='bc1_assoc',
)
```

Pack your sequencing reads (filters, deduplicates, and optionally merges paired-ends):

```python
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
```

Count! This verifies barcode-variant associations:

```python
counts_df, stats = op.acount(
    index_file='bc1_assoc',
    pack_file='experiment',
    count_file='results',
)

print(counts_df)
```

---

## CLI Reference

[^ Back to TOC](#table-of-contents)

Every module is available via command line. See all commands:

```bash
op
```

Get help on a specific command (shows options and defaults):

```bash
op barcode
```

Read the full manual for a command:

```bash
op manual barcode
```

Run a command:

```bash
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

Install tab completion (you'll thank yourself later):

```bash
op complete --install
```

### CLI-Specific Notes

[^ Back to TOC](#table-of-contents)

**CLI outputs require an output basename**: unlike the Python API (where many modules can return results in memory), CLI commands that write files require an output basename (for example, `--output-file`, `--index-file`, `--pack-file`, `--count-file`, `--mapping-file`, `--synthesis-file`).

**Output filenames are auto-suffixed**: commands append a suffix if missing (for example, `.oligopool.barcode.csv`). Prefer basenames like `--output-file output_basename` to avoid doubled extensions.
For example, `output_file: "my_library.csv"` becomes `my_library.csv.oligopool.barcode.csv`, while `output_file: "my_library.oligopool.barcode.csv"` is left unchanged.

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

**Callbacks are Python-only**: The `callback` parameter in `acount`/`xcount` is not available via CLI (including YAML pipelines, which execute via CLI). For custom read processing logic, use the Python API.

> **API Reference**: See [api.md](api.md#cli-parameter-mapping) for complete CLI parameter mapping.

---

## Config Files

[^ Back to TOC](#table-of-contents)

The CLI supports YAML config files for repeatable, documented workflows. Config files eliminate long command lines and enable multi-step pipeline execution.

### Why YAML?

[^ Back to TOC](#table-of-contents)

YAML is useful here because it gives you one readable, versioned run spec:

1. **It keeps long command lines manageable**: You can review parameters as structured data instead of parsing huge CLI strings.
2. **It improves reproducibility**: The exact settings live in a file you can commit, diff, and reuse.
3. **It supports multi-step workflows cleanly**: Pipeline step ordering/dependencies stay explicit in one place.
4. **It reduces repetition**: Lists, comments, and anchors/aliases make repetitive configs easier to maintain.
5. **It still stays flexible**: You can override any value at runtime with CLI flags.

### Single Command Config

[^ Back to TOC](#table-of-contents)

Use `--config` with any command to load parameters from a YAML file:

```bash
op barcode --config barcode_design.yaml
```

Example config file (`barcode_design.yaml`):

```yaml
# Barcode design config
barcode:
  input_data: "library.csv"
  oligo_length_limit: 250
  barcode_length: 16
  minimum_hamming_distance: 3
  maximum_repeat_length: 8
  barcode_column: "BC"
  left_context_column: "Variant"
  output_file: "library_barcoded"
  barcode_type: "spectrum"
  random_seed: 42  # Strongly recommended for reproducible, documented barcode design runs.
  excluded_motifs:
    - "GGATCC"                 # Inline motif
    - "cloning_blacklist.csv"  # CSV with Exmotif column
    - "TCTAGA"
  # Multi-value args like excluded_motifs and background_directory are YAML lists.
  # CLI equivalent: --background-directory host_bg plasmid_bg
  background_directory:
    - "host_bg"
    - "plasmid_bg"
```

Config keys use snake_case (same as Python API). CLI arguments override config values:

```bash
# Override barcode_length from config
op barcode --config barcode_design.yaml --barcode-length 20
```

### Pipeline Execution

[^ Back to TOC](#table-of-contents)

Run multi-step workflows from a single config with `op pipeline`:

```bash
op pipeline --config mpra_design_serial.yaml
```

Example pipeline config (condensed):

```yaml
# MPRA Design Pipeline
pipeline:
  name: "MPRA Library Design"
  steps:
    - primer
    - barcode
    - spacer
    - final

primer:
  input_data: "variants.csv"
  output_file: "step1_primer"
  oligo_length_limit: 250
  primer_sequence_constraint: "N*22"
  primer_type: "forward"
  primer_column: "FwdPrimer"
  minimum_melting_temperature: 52
  maximum_melting_temperature: 58
  maximum_repeat_length: 10
  # You can stack multiple backgrounds by listing them here.
  # CLI equivalent: --background-directory host_bg plasmid_bg
  background_directory: ["host_bg", "plasmid_bg"]

barcode:
  input_data: "step1_primer"
  output_file: "step2_barcode"
  oligo_length_limit: 250
  barcode_length: 16
  minimum_hamming_distance: 3
  maximum_repeat_length: 8
  barcode_column: "BC"
  left_context_column: "FwdPrimer"
  # Lists map to multi-value CLI args (nargs='+'): pass multiple values after the flag.
  # CLI equivalent: --excluded-motifs GGATCC TCTAGA
  excluded_motifs: ["GGATCC", "TCTAGA"]

spacer:
  input_data: "step2_barcode"
  output_file: "step3_spacer"
  oligo_length_limit: 250
  maximum_repeat_length: 8
  spacer_column: "Spacer"
  left_context_column: "BC"

final:
  input_data: "step3_spacer"
  output_file: "final_library"
```

For a runnable file version, see `examples/cli-yaml-pipeline/mpra_design_serial.yaml`.

Pipeline input resolution supports both explicit filenames and basename chaining. Example:
- `output_file: "step1_primer"` writes `step1_primer.oligopool.primer.csv`
- a downstream `input_data: "step1_primer"` resolves to that produced file
- if `input_data` is an explicit existing path, it is used as-is
- if multiple steps produce the same basename alias, pipeline validation fails with a config error

### Parallel Pipeline Execution

[^ Back to TOC](#table-of-contents)

Parallel execution helps most when steps are naturally independent. In practice, that is usually analysis workflows (`index`, `pack`, `xcount`, `acount`), while design workflows are often dependency-dense and mostly sequential.

The fundamental YAML change is in `pipeline.steps`:
- **Serial format**: a list of command names (`- primer`, `- barcode`, ...), which runs in listed order.
- **Parallel/DAG format**: a list of step objects with `name`, `command`, and optional `after`, which defines dependencies explicitly.
- **Mixed format (supported)**: you can mix strings and step objects in one list; if any dict-style step appears, DAG parsing is used and string steps are treated as `name=command` with no dependencies.

```yaml
# Serial
pipeline:
  steps:
    - primer
    - barcode
    - spacer

# Parallel/DAG
pipeline:
  steps:
    - name: index_bc1
      command: index
    - name: count
      command: xcount
      after: [index_bc1]

# Mixed (also valid)
pipeline:
  steps:
    - primer
    - name: barcode_step
      command: barcode
      after: [primer]
```

**Step fields:**
- `name`: Step identifier (required)
- `command`: Oligopool command to run (defaults to `name`)
- `after`: List of step names this step depends on (optional)
- `config`: Config section name to use (defaults to `name`)

**Execution model:**
- Steps with no `after` dependencies form level 1 (eligible to run concurrently).
- Each subsequent level waits for dependencies, and independent steps are grouped automatically for concurrent execution.
- `--dry-run` shows execution levels and parallelism.

```yaml
# analysis_multi.yaml
pipeline:
  name: "Multi-Sample Analysis DAG Demo"
  steps:
    - name: index_bc1
      command: index
    - name: index_bc2
      command: index
    - name: pack_sample_a
      command: pack
    - name: xcount_sample_a
      command: xcount
      after: [index_bc1, index_bc2, pack_sample_a]
    - name: acount_sample_a
      command: acount
      after: [index_bc2, pack_sample_a]
    # Repeat pack/xcount/acount step groups for sample_b, sample_c, ...
```

`--dry-run` should show:
- Level 1 (parallel): `index_bc1`, `index_bc2`, `pack_sample_a`, ...
- Level 2 (parallel): `xcount_sample_a`, `acount_sample_a`, ...

On real execution, level 1 builds shared index artifacts once and packs each sample concurrently. Level 2 then runs per-sample counting branches in parallel (`xcount` and/or `acount`) using those artifacts.

For runnable DAG examples, see:
- `examples/cli-yaml-pipeline/analysis_multi.yaml` (index once, then per-sample `pack → (xcount and/or acount)` branches)
- `examples/cli-yaml-pipeline/analysis_single.yaml` (same pattern reduced to one sample)

**Rare parallel design branches with `join`:**

Most Design Mode workflows are sequential, because element constraints tend to couple steps. If two design steps are truly independent (they depend only on already-fixed context and do not depend on each other), you can run them as parallel branches and then recombine their outputs with `join`.

```yaml
pipeline:
  name: "Parallel Branch Join (Design)"
  steps:
    - name: branch_a
      command: barcode
    - name: branch_b
      command: barcode
    - name: rejoin
      command: join
      after: [branch_a, branch_b]
    - name: spacer
      command: spacer
      after: [rejoin]
    - name: final
      command: final
      after: [spacer]

common_design: &common_design
  input_data: "variants.csv"
  oligo_length_limit: 200
  maximum_repeat_length: 8
  random_seed: 42

branch_a:
  <<: *common_design
  output_file: "bc_a"
  barcode_column: "BC_A"
  # ... barcode params (length, distance, context, motifs, backgrounds)

branch_b:
  <<: *common_design
  output_file: "bc_b"
  barcode_column: "BC_B"
  # ... barcode params

rejoin:
  input_data: "bc_a"          # resolves to bc_a.oligopool.barcode.csv
  other_data: "bc_b"          # resolves to bc_b.oligopool.barcode.csv
  oligo_length_limit: 200
  join_policy: left           # 0/left or 1/right
  output_file: "bc_ab"

spacer:
  input_data: "bc_ab"
  output_file: "bc_ab_spacer"
  oligo_length_limit: 200
  maximum_repeat_length: 8
  spacer_column: "Spacer"
  left_context_column: "BC_A"

final:
  input_data: "bc_ab_spacer"
  output_file: "bc_ab_final"
```

For a runnable branch-then-join example, see `examples/cli-yaml-pipeline/mpra_design_parallel.yaml`.

### Dry Run Validation

[^ Back to TOC](#table-of-contents)

Validate a pipeline config without executing:

```bash
op pipeline --config mpra_design_serial.yaml --dry-run
op pipeline --config analysis_multi.yaml --dry-run
```

This checks that all steps are valid commands and displays the parameters for each step.

Common validation errors you might see:

| Message (or symptom) | Typical cause | Fix |
|---|---|---|
| ambiguous output alias / duplicate basename | Multiple steps write the same logical basename | Give each producing step a unique basename |
| missing config section for step | A step in `pipeline.steps` has no matching config block | Add the missing section or point `config:` to the correct section |
| unknown command in pipeline step | Typo or non-existent command name | Use a valid `op` command name (`primer`, `barcode`, `xcount`, ...) |

### Config Precedence

[^ Back to TOC](#table-of-contents)

Parameter values are resolved in this order (highest priority first):

1. **CLI arguments** - Always win
2. **Config file values** - Fill in unset args
3. **Command defaults** - Built-in defaults

### Parameter Sharing

[^ Back to TOC](#table-of-contents)

When multiple steps share the same arguments, define a common block once and merge it into each step.
This keeps large pipelines readable and reduces copy/paste errors.

```yaml
common_count: &common_count
  mapping_type: sensitive
  failed_reads_sample_size: 1000

xcount_sample_a:
  <<: *common_count
  index_files: ["ms_bc1", "ms_bc2"]
  pack_file: "ms_s1_reads"
  count_file: "ms_s1_xcount"

xcount_sample_b:
  <<: *common_count
  index_files: ["ms_bc1", "ms_bc2"]
  pack_file: "ms_s2_reads"
  count_file: "ms_s2_xcount"

xcount_sample_c:
  <<: *common_count
  index_files: ["ms_bc1", "ms_bc2"]
  pack_file: "ms_s3_reads"
  count_file: "ms_s3_xcount"
```

The `<<` merge key is standard YAML and works with this CLI because configs are loaded with `PyYAML`.
If you preprocess YAML with stricter tooling, verify merge-key support there too.

### Config Tips

[^ Back to TOC](#table-of-contents)

1. **Use comments liberally**: YAML supports `#` comments - document your design choices.
2. **Keep configs in version control**: Reproducibility for future you.
3. **Choose explicit or shorthand**: Use explicit paths for full control, or basename chaining (e.g., `output_file: step1`, next step `input_data: step1`) for concise configs.
4. **Lists in YAML**: For multi-value parameters like `excluded_motifs` and `background_directory`, use YAML lists. This corresponds to the CLI style of passing multiple values after the same flag.
5. **Type aliases work**: Use `"spectrum"` or `1` for `barcode_type` - same as CLI.
6. **Avoid alias collisions**: Reusing the same output basename in multiple steps creates ambiguous aliases, which pipeline validation rejects as a config error.

---

## Tips & Tricks

[^ Back to TOC](#table-of-contents)

### Design Tips

[^ Back to TOC](#table-of-contents)

1. **Design in dependency order**: Start with constrained context-defining regions (fewest degenerate bases), then design dependent regions (anchors before barcodes; inner/primary primers before paired outer/reverse primers).

2. **Check length early and often**: Run `lenstat` after each element to avoid surprises.

3. **Read `stats` after every step**: All modules return `stats`; use them to catch infeasibility, conflicts, and weak settings early.

4. **Layer constraints**: Set core constraints first, add optional layers incrementally, and if infeasible remove optional layers one at a time while keeping core layers fixed.

5. **Start loose, then tighten**: Begin with permissive settings (algorithm mode, distance/Tm/repeat thresholds, optional constraints), then tighten only after a feasible baseline is found.

6. **Use Patch Mode for iterations**: When adding variants to an existing pool, use `patch_mode=True` to preserve existing designs.

7. **Mind edge effects**: Always specify context columns; boundary motif/repeat emergence is common.

8. **Use metadata-driven design inputs**: Keep a per-oligo `metadata_df` keyed by `ID` (for example `OligoSet`, `SpacerLength`, condition labels), then pass module-specific views as needed:

```python
# Reuse one metadata table across modules
df, _ = op.primer(..., oligo_set=metadata_df[['ID', 'OligoSet']])
df, _ = op.spacer(..., spacer_length=metadata_df[['ID', 'SpacerLength']])
```

### Analysis Tips

[^ Back to TOC](#table-of-contents)

1. **Index anchors matter**: The quality of your anchors determines counting accuracy. Design them well.

2. **Quality filter aggressively**: Low-quality reads just add noise. Use `minimum_*_read_quality=30` or higher.

3. **Check run summaries**: For counting runs, inspect returned stats and sample `failed_reads_file` (when enabled) to diagnose filtering, mapping, and ambiguity.

4. **Start with sensitive mapping**: Use `mapping_type=1` initially to catch more reads, then compare with `mapping_type=0`.

5. **Index once, pack once, iterate counts**: Index files and pack files are reusable; iterate `acount`/`xcount` settings without rebuilding artifacts.

6. **Scale knobs**: For analysis, tune `core_count` and `memory_limit` together when runs are slow or memory constrained.

### Performance Tips

[^ Back to TOC](#table-of-contents)

1. **Parallelize**: Most modules auto-detect cores. For analysis, ensure you have RAM (memory_limit parameter).

2. **Use fast barcode mode**: `barcode_type=0` is usually sufficient and much faster.

---

## Composability

[^ Back to TOC](#table-of-contents)

Most complex oligo library designs don't require new features - they require composing existing primitives. This section teaches you how.

### The Decomposition Principle

[^ Back to TOC](#table-of-contents)

You can build any compound oligo element from the primitive modules (`motif`, `barcode`, `primer`, `spacer`) using `merge`, `revcomp`, and context columns. You can build complex analyses from independent `index` files, reusable `pack` files, `acount`/`xcount`, and `callback` functions. Constraint composition (`excluded_motifs` + `background`) threads across the whole pipeline. Because modules accept and return DataFrames, small `pandas` transforms fit naturally between steps. If an ask sounds like a new feature, try decomposing it into existing primitives first.

### Design Recipes

[^ Back to TOC](#table-of-contents)

#### Barcode with Embedded Restriction Site

**The ask**: "I need a barcode that contains an EcoRI site in the middle."

**The decomposition**: Conceptually split the barcode into two halves: place a constant restriction-site motif as the center core, design the left and right extensions with `barcode`, then merge all three pieces into one compound element column.

```python
import pandas as pd
import oligopool as op

df = pd.DataFrame({
    'ID': ['V1', 'V2', 'V3'],
    'Variant': ['ATGCATGCATGCATGC', 'GCTAGCTAGCTAGCTA', 'TTAATTAATTAATTAA']
})

# Step 1: Insert the constant EcoRI site (center core)
df, _ = op.motif(
    input_data=df,
    oligo_length_limit=200,
    motif_sequence_constraint='GAATTC',
    maximum_repeat_length=8,
    motif_column='EcoRI_Site',
    motif_type=1,                     # Constant - same for all rows
    left_context_column='Variant',
)

# Step 2: Design the first barcode half (left extension)
df, _ = op.barcode(
    input_data=df,
    oligo_length_limit=200,
    barcode_length=6,
    minimum_hamming_distance=2,
    maximum_repeat_length=8,
    barcode_column='BC1_Left',
    left_context_column='Variant',
    right_context_column='EcoRI_Site',
)

# Step 3: Design the second barcode half (right extension)
df, _ = op.barcode(
    input_data=df,
    oligo_length_limit=200,
    barcode_length=6,
    minimum_hamming_distance=2,
    maximum_repeat_length=8,
    barcode_column='BC1_Right',
    left_context_column='EcoRI_Site',
)

# Step 4: Merge into a single compound barcode column
df, _ = op.merge(
    input_data=df,
    merge_column='BC1',
    left_context_column='BC1_Left',
    right_context_column='BC1_Right',
)

# Result: df now has a 'BC1' column containing 6bp + GAATTC + 6bp = 18bp compound barcodes
# Each compound barcode has a guaranteed EcoRI site in the middle
```

**Why it works**: `merge` concatenates the column range from `BC1_Left` through `BC1_Right` into a single `BC1` column. The source columns are removed. The two barcode halves provide Hamming separation on either side of the constant site.

#### Barcode with Flanking Anchors (Bridge to Analysis)

**The ask**: "I need constant sequences flanking my barcode so I can find it in reads later."

**The decomposition**: Design both constant anchor motifs first, then design the barcode between them using both context columns. These anchors serve double duty: they prevent edge effects during design and become the prefix/suffix anchors for `index`.

```python
# Design constant prefix anchor (>=6 bp recommended for indexing).
# Note: 'NNNNNNNN' is an IUPAC constraint - the anchor is *designed*, then held constant across rows.
df, _ = op.motif(
    input_data=df,
    oligo_length_limit=200,
    motif_sequence_constraint='NNNNNNNN',   # 8N - enough for robust anchoring
    maximum_repeat_length=8,
    motif_column='BC_Prefix',
    motif_type=1,                           # Constant anchor
    left_context_column='Variant',
)

# Design constant suffix anchor
df, _ = op.motif(
    input_data=df,
    oligo_length_limit=200,
    motif_sequence_constraint='NNNNNNNN',
    maximum_repeat_length=8,
    motif_column='BC_Suffix',
    motif_type=1,
    left_context_column='Variant',
)

# Design barcode between anchors (use both contexts)
df, _ = op.barcode(
    input_data=df,
    oligo_length_limit=200,
    barcode_length=12,
    minimum_hamming_distance=3,
    maximum_repeat_length=8,
    barcode_column='BC1',
    left_context_column='BC_Prefix',
    right_context_column='BC_Suffix',
)

# Later, at analysis time - these anchors locate the barcode in reads:
# op.index(
#     barcode_data=df,
#     barcode_column='BC1',
#     barcode_prefix_column='BC_Prefix',
#     barcode_suffix_column='BC_Suffix',
#     index_file='bc1_index',
# )
```

#### Cross-Talk-Free Barcode Sets

**The ask**: "I need BC1 and BC2 in the same pool but they must never be confused with each other."

**The decomposition**: Design BC1 first, then design BC2 with `cross_barcode_columns` pointing at BC1.

```python
# Design first barcode set
df, _ = op.barcode(
    input_data=df,
    oligo_length_limit=200,
    barcode_length=12,
    minimum_hamming_distance=3,
    maximum_repeat_length=8,
    barcode_column='BC1',
    left_context_column='Variant',
)

# Design second barcode set - globally separated from BC1
df, _ = op.barcode(
    input_data=df,
    oligo_length_limit=200,
    barcode_length=12,
    minimum_hamming_distance=3,
    maximum_repeat_length=8,
    barcode_column='BC2',
    left_context_column='Variant',
    cross_barcode_columns=['BC1'],     # Stay far from every BC1 barcode
    minimum_cross_distance=2,          # At least 2 mismatches to any BC1
)
```

Cross-set separation is global: each new BC2 barcode must be at least `minimum_cross_distance` away from **every** barcode in BC1 (not just the barcode in the same row).

#### General Compound Element Pattern

Any compound element follows the same recipe: design sub-regions in dependency/context order, starting with the most constrained pieces (the ones with the fewest degenerate bases), then `merge`.

```python
# Example: constant tag + variable region + constant tag
# Design order (dependency/context-driven, non-linear):
df, _ = op.motif(..., motif_column='Tag5', motif_type=1)          # Constant 5' tag (more constrained)
df, _ = op.motif(..., motif_column='Tag3', motif_type=1,
                 left_context_column='Tag5')                      # Constant 3' tag (equally constrained)
df, _ = op.motif(..., motif_column='VarRegion', motif_type=0,
                 left_context_column='Tag5',
                 right_context_column='Tag3')                     # Variable middle region (less constrained)

df, _ = op.merge(..., merge_column='CompoundElement',
                 left_context_column='Tag5',
                 right_context_column='Tag3')
```

#### Long Constructs via `split` → `pad` → `final`

If an oligo is too long to synthesize directly, split it into overlapping fragments (`split`), pad each fragment into an assembly-ready layout (`pad`), then finalize each fragment for synthesis (`final`).

```python
import oligopool as op

# 1) Split into overlapping fragments
frags, _ = op.split(
    input_data='full_length_library.csv',
    split_length_limit=180,
    minimum_melting_temperature=55.0,
    minimum_hamming_distance=2,
    minimum_overlap_length=20,
    maximum_overlap_length=30,
    separate_outputs=True,  # list of per-fragment DataFrames
)

# 2) Pad + finalize each fragment as its own synthesis pool
for i, frag_df in enumerate(frags, 1):
    split_col = next(col for col in frag_df.columns if col != 'ID')
    padded_df, _ = op.pad(
        input_data=frag_df,
        oligo_length_limit=200,
        split_column=split_col,
        typeIIS_system='BsaI',
        minimum_melting_temperature=55.0,
        maximum_melting_temperature=75.0,
        maximum_repeat_length=12,
        output_file=f'frag_{i}.padded',
    )
    final_df, _ = op.final(
        input_data=padded_df,
        output_file=f'frag_{i}.final',
    )
```

#### Flip Orientation Mid-Pipeline (`revcomp`)

Sometimes a region is designed in the "readout" orientation but must be installed/synthesized in the opposite orientation (or you want to sanity-check split fragment orientation). Use `revcomp` to reverse-complement sequences and reverse the column order across a column range.

```python
import pandas as pd

df = pd.DataFrame({
    'ID': ['V1'],
    'Left': ['ATGC'],
    'Middle': ['GATTACA'],
    'Right': ['CCGG'],
})

# Reverse-complement Left..Right, and reverse their column order:
df, _ = op.revcomp(
    input_data=df,
    left_context_column='Left',
    right_context_column='Right',
)
```

### Constraint Composition

[^ Back to TOC](#table-of-contents)

#### Cut-Site-Free Library with Host Genome Screening

**The ask**: "I need a library with no EcoRI or BamHI sites anywhere, and no homology to the E. coli genome."

**The decomposition**: Build a background database for the host genome. Then pass the **same** `excluded_motifs` list and `background_directory` to **every** design module and to `verify` for final QC.

```python
import oligopool as op

# Restriction sites to exclude - for non-palindromic motifs, include the reverse complement too.
cut_sites = [
    'GAATTC',  # EcoRI (palindromic - reverse complement is itself)
    'GGATCC',  # BamHI (palindromic - reverse complement is itself)
    'CGTCTC', 'GAGACG',  # BsmBI (non-palindromic - MUST include both)
]

# Build host genome background DB
op.background(
    input_data='ecoli_genome.fasta',
    maximum_repeat_length=15,
    output_directory='ecoli_bg',
)
bg_dir = 'ecoli_bg'

# Now pass BOTH constraints to every design step:
df, _ = op.primer(
    ...,
    excluded_motifs=cut_sites,
    background_directory=bg_dir,
)

df, _ = op.motif(
    ...,
    excluded_motifs=cut_sites,
    background_directory=bg_dir,
)

df, _ = op.barcode(
    ...,
    excluded_motifs=cut_sites,
    background_directory=bg_dir,
)

df, _ = op.spacer(
    ...,
    excluded_motifs=cut_sites,
    background_directory=bg_dir,
)

# Final QC - verify catches anything that slipped through at junctions
verify_df, verify_stats = op.verify(
    input_data=df,
    oligo_length_limit=200,
    excluded_motifs=cut_sites,
    background_directory=bg_dir,
)
print(f"Any conflicts? {verify_stats['vars']['any_conflict']}")
```

**Why reverse complements matter**: `excluded_motifs` does literal substring matching on the designed oligo sequence (strict ATGC only). If your downstream construct is double-stranded (e.g., cloning) and you want to avoid a restriction site regardless of which strand/orientation it appears on, include the motif **and its reverse complement** for each non-palindromic site (e.g., BsmBI: `CGTCTC` + `GAGACG`). For palindromic sites (EcoRI `GAATTC`, BamHI `GGATCC`), the reverse complement is identical, so listing it once is enough.

**Why pipeline-wide matters**: Each design module avoids creating the motif within its own element. But a motif could emerge at the **junction** between two independently designed elements. Context columns help within a single module call, but `verify` checks the fully concatenated oligo for any remaining conflicts.

#### Constraint Layering Strategy (Add/Remove Layers)

Use multiple `excluded_motifs` sources and multiple background DBs as explicit constraint layers that you can toggle while troubleshooting feasibility.

```python
constraint_layers = {
    'cloning_sites': 'cutsites.csv',
    'platform_rules': 'platform_forbidden.csv',
    'homopolymers': ['AAAAA', 'CCCCC', 'GGGGG', 'TTTTT'],
}
background_layers = ['ecoli_bg', 'vector_bg', 'phix_bg']

df, stats = op.barcode(
    ...,
    excluded_motifs=constraint_layers,
    background_directory=background_layers,
)
```

Recommended workflow:
1. Start with core mandatory layers only (for example cloning sites + primary host background).
2. Add optional hardening layers incrementally (extra host DBs, platform-specific motif sets).
3. If infeasible, remove one optional layer at a time and rerun.
4. Lock your final layer set and thread it unchanged through all design modules and `verify`.

This gives you better documentation (`stats` can report per-set motif attribution) and a controlled path to recover feasibility without losing track of which constraints are active.

### Analysis Recipes

[^ Back to TOC](#table-of-contents)

#### Multi-Index Combinatorial Counting

**The ask**: "I have BC1 and BC2 in each oligo. I need a count matrix of all BC1 x BC2 combinations."

**The decomposition**: Build a separate `index` for each barcode (each with its own anchors), then pass both to `xcount`.

```python
# Index BC1 using its flanking anchors
op.index(
    barcode_data='library.csv',
    barcode_column='BC1',
    barcode_prefix_column='BC1_Prefix',
    barcode_suffix_column='BC1_Suffix',
    index_file='bc1_idx',
)

# Index BC2 using its flanking anchors
op.index(
    barcode_data='library.csv',
    barcode_column='BC2',
    barcode_prefix_column='BC2_Prefix',
    barcode_suffix_column='BC2_Suffix',
    index_file='bc2_idx',
)

# Pack reads (once - reusable)
op.pack(
    r1_fastq_file='reads_R1.fq.gz',
    r2_fastq_file='reads_R2.fq.gz',
    r1_read_type=0,
    r2_read_type=1,
    pack_type=1,
    pack_file='experiment',
)

# Combinatorial counting - xcount handles multiple indexes natively
counts_df, stats = op.xcount(
    index_files=['bc1_idx', 'bc2_idx'],
    pack_file='experiment',
    count_file='combo_counts',
)
# Output: columns bc1_idx.ID, bc2_idx.ID, CombinatorialCount
# Reads where one barcode is missing show '-' for that position
```

This scales to N-way matrices: just add more index files to the list.

#### Custom Callback Filtering

**The ask**: "I only want to count reads that meet custom criteria (e.g., minimum length, specific subsequence)."

**The decomposition**: Write a callback function and pass it to `acount` or `xcount` (Python API only - not available via CLI).

```python
def my_filter(r1, r2, ID, count, coreid):
    """
    r1:     str       - Read 1 sequence (always present)
    r2:     str|None  - Read 2 (None for merged/single-end)
    ID:     tuple     - Identified barcode ID(s)
    count:  int       - Read frequency in current pack
    coreid: int       - CPU core ID
    Returns: bool     - True to accept, False to reject
    """
    # Example: reject short reads
    if len(r1) < 100:
        return False
    # Example: require a specific subsequence
    if 'GCTAGC' not in r1:
        return False
    return True

counts_df, stats = op.xcount(
    index_files=['bc1_idx'],
    pack_file='experiment',
    count_file='filtered_counts',
    callback=my_filter,
)
```

### Degenerate Recipes

[^ Back to TOC](#table-of-contents)

Degenerate mode (`compress`/`expand`) supports multiple workflows: (1) cost-efficient synthesis of large variant sets via IUPAC-degenerate oligos, (2) selection-based screens where you sequence survivors and map them back using `mapping_df`, and (3) "compression for analysis" where a single degenerate library (plus its mapping) becomes a reusable bridge between sequencing results and the original variant IDs.

#### Saturation Mutagenesis Compression

**The ask**: "I have 1000 single-amino-acid substitution variants. Can I synthesize fewer oligos?"

**The decomposition**: Generate all substitution variants, then `compress` into IUPAC-degenerate oligos. Similar variants collapse into shared degenerate representations.

```python
import pandas as pd
import oligopool as op

# Your substitution variants (e.g., alanine scanning, full saturation)
df = pd.DataFrame({
    'ID': [f'mut_{i}' for i in range(1000)],
    'Sequence': [...]  # 1000 single-substitution variants, strict ATGC
})

# Compress - similar sequences share degenerate representations
mapping_df, synthesis_df, stats = op.compress(
    input_data=df,
    mapping_file='satmut_library',
    synthesis_file='satmut_library',
    random_seed=42,
)

print(f"Compressed {stats['vars']['input_variants']} variants "
      f"into {stats['vars']['degenerate_oligos']} oligos "
      f"({stats['vars']['compression_ratio']}x compression)")

# Verify losslessness
expanded_df, _ = op.expand(
    input_data=synthesis_df,
    sequence_column='DegenerateSeq',
    mapping_file='satmut_library',
)

# Order synthesis_df for synthesis - these are your degenerate oligos
```

#### Selection-Based Discovery Without Barcodes

**The ask**: "I want to screen variants by selection (growth, FACS) without barcodes."

**The decomposition**: `compress` → synthesize → select → sequence survivors → map back to original variant IDs using `mapping_df`.

```python
# 1. Compress variants
mapping_df, synthesis_df, stats = op.compress(
    input_data='variants.csv',
    mapping_file='discovery_lib',
    synthesis_file='discovery_lib',
)

# 2. Order synthesis_df for synthesis (6-20x fewer oligos than individual variants)
# 3. Clone, select, sequence survivors

# 4. Map surviving sequences back to variant IDs
#    mapping_df has columns: ID, Sequence, DegenerateID
#    Match sequenced survivors to 'Sequence' column to recover original 'ID'
survivors = pd.read_csv('sequenced_survivors.csv')
hits = mapping_df[mapping_df['Sequence'].isin(survivors['Sequence'])]
```

**Key property**: Compression is lossless by construction. `expand` recovers exactly the input set, with no invented sequences. Sequences of different lengths are compressed into separate groups automatically.

**Important caveat**: treat `compress` output (`synthesis_df`) as a terminal branch in design mode. You usually don't chain it into `barcode`, `primer`, or `spacer`, because barcodes are intentionally maximally different and that works against IUPAC compression. Use compression for selection-style workflows, and use the standard design pipeline for barcode-based readouts.

### Application Templates

[^ Back to TOC](#table-of-contents)

These compact templates show how the primitives compose for common biological applications. Adapt parameters to your specific constraints.

#### Promoter MPRA Library

```
pipeline/
|-- background(host.fasta)
|-- primer(fwd)
|-- motif(constant, anchor)            # BC prefix for indexing
|-- motif(constant, anchor)            # BC suffix for indexing
|-- barcode                            # Designed with both motif flanks as context
|-- primer(rev, paired)
|-- spacer(auto)
|-- lenstat
|-- verify
`-- final
```

Thread `excluded_motifs` (cloning sites + reverse complements) and `background_directory` through every design module.
This same single-barcode MPRA pattern also applies to other regulatory/stability-element libraries (e.g., UTR variants, degrons, structured elements).

#### CRISPR Guide Library

```
pipeline/
|-- background(host.fasta)
|-- primer(fwd)
|-- barcode
|-- spacer
|-- lenstat
|-- verify
`-- final
```

Guides are the variant column. Add `excluded_motifs` for cloning-system recognition sites (e.g., BsmBI: `CGTCTC`/`GAGACG`, BbsI: `GAAGAC`/`GTCTTC`) to prevent internal cleavage.

#### Ribozyme Library (Pre/Post Cleavage Readout)

```
pipeline/
|-- background(host.fasta)
|-- motif(constant, conserved_core)                    # Catalytic/structural region
|-- motif(constant, BCmolecule_anchor)                 # Anchor for indexing
|-- barcode(BC_molecule)                               # Tracks variant/molecule identity (present in all outcomes)
|-- motif(constant, BCstate_anchor)                    # Anchor for indexing
|-- barcode(BC_state, cross_barcode_columns=['BC_molecule'])  # "State" barcode designed to avoid cross-talk
|-- primer(fwd)
|-- primer(rev, paired)
|-- spacer
|-- lenstat
|-- verify
`-- final
```

Variants flank the conserved core. Use a dual-barcode strategy when cleavage changes what part of the construct is sequenced in your assay:

- `BC_molecule`: on a region present in both cleaved/uncleaved products (per-molecule ID).
- `BC_state`: placed so it is only observed in one state (assay-dependent capture/amplification design).

At analysis time, build two indexes (`BC_molecule`, `BC_state`) and use `xcount(index_files=[...])` to count barcode combinations. Because `xcount` represents missing barcodes as `'-'`, you can quantify state by comparing counts for `(BC_molecule, BC_state)` pairs vs `BC_state='-'` for the same `BC_molecule`.

### Composability Cheat Sheet

[^ Back to TOC](#table-of-contents)

| # | Ask | Compose |
|---|-----|---------|
| 1 | Embed restriction site in barcode | `motif(const, site)` → `barcode(left_ext, right_context=site)` → `barcode(right_ext, left_context=site)` → `merge` |
| 2 | Flanking anchors for barcode | `motif(const, prefix)` → `motif(const, suffix)` → `barcode` (use both as contexts) |
| 3 | Cross-talk-free barcode sets | `barcode(BC1)` → `barcode(BC2, cross_barcode_columns)` |
| 4 | Tm-matched primer pair | `primer(fwd)` → `primer(rev, paired_primer_column)` |
| 5 | Extend existing pool | append rows → `patch_mode=True` on each element |
| 6 | Multi-genome screening | `background` per genome → `background_directory=[bg1, bg2]` |
| 7 | Auto-fill to length | `spacer(spacer_length=None)` |
| 8 | Compound element | chain `motif`/`barcode`/`primer`/`spacer` → `merge`/`revcomp` → repeat |
| 9 | Flip orientation mid-pipeline | `revcomp(left_context_column=..., right_context_column=...)` |
| 10 | Long constructs / assembly | `split` → `pad` (per fragment) → `final` (per fragment) |
| 11 | Mid-pipeline length check | `lenstat` → adjust element lengths → rerun |
| 12 | Inspect artifacts quickly | `inspect(background/index/pack)` before reuse |
| 13 | Cut-site-free library | `excluded_motifs` (motif + reverse complement) on all modules + `verify` |
| 14 | No host homology | `background_directory` on all modules + `verify` |
| 15 | Multi-barcode counting | separate `index` per barcode → `xcount(index_files=[...])` |
| 16 | Custom read filtering | `callback(r1, r2, ID, count, coreid)` (Python only) |
| 17 | Pack once, count many | single `pack` → reuse with different indexes/modes |
| 18 | Cost-efficient mutagenesis | `compress` → order `synthesis_df` |
| 19 | Selection without barcodes | `compress` → select → map back via `mapping_df` |
| 20 | Recombine parallel branch outputs | parallel branches → `join` → `verify` |

### Stuff to Note

[^ Back to TOC](#table-of-contents)

Composability works because each module enforces local constraints against explicit context, then passes structured outputs forward. Most surprises come from where constraints actually apply (junctions, reverse complements, terminal outputs) and from features that are Python-only.

- **Always set context columns during design**: Use `left_context_column`/`right_context_column` at each design step to prevent motif/repeat artifacts across junction boundaries.
- **Treat `merge` as an inclusive range operation**: It consumes `left_context_column` through `right_context_column` (inclusive) and removes the source columns.
- **Treat `compress` as a terminal branch**: `synthesis_df` is not chainable back into design modules like `barcode`, `primer`, or `spacer`.
- **Know the interface boundary**: `callback` is Python-only for `acount`/`xcount`; YAML pipelines run via CLI too, so use post-run filtering on count outputs there as well.

---

## Advanced Modules

[^ Back to TOC](#table-of-contents)

For power users who want to peek under the hood:

### vectorDB

[^ Back to TOC](#table-of-contents)

ShareDB-based k-mer storage. Created by `background()`, but you can access it directly.

Open or create a database:

```python
db = op.vectorDB(
    path='ecoli_bg.oligopool.background',
    maximum_repeat_length=15,
)
```

Check if a sequence has k-mers in the database:

```python
if 'ATGCATGCATGC' in db:
    print("Found matching k-mers!")
```

Add or remove sequences:

```python
db.add('ATGCATGCATGC')
db.remove('ATGCATGCATGC')
```

Useful for inspecting or manipulating background databases.

> **API Reference**: See [`vectorDB`](api.md#vectordb) for complete method documentation.

### Scry

[^ Back to TOC](#table-of-contents)

1-nearest-neighbor barcode classifier. Powers `acount`/`xcount` internally.

Create and train a model with your barcodes:

```python
model = op.Scry()
model.train(X=['ATGCATGC', 'GCTAGCTA'], Y=['bc1', 'bc2'], n=8, k=6, t=1)
```

Prime for prediction (set error tolerance and mode):

```python
model.prime(t=1, mode=0)  # 0=fast, 1=sensitive
```

Predict labels:

```python
label, score = model.predict('ATGCATGC')  # Returns ('bc1', 1.0)
```

Useful for building custom counting pipelines or debugging classification issues.

> **API Reference**: See [`Scry`](api.md#scry) for complete method documentation.

---

## Troubleshooting

[^ Back to TOC](#table-of-contents)

### Design Failures

[^ Back to TOC](#table-of-contents)

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

### Assembly Failures

[^ Back to TOC](#table-of-contents)

**"Split Design Infeasible"**
- Check `stats['step_name']` and `stats['vars']` first (especially contig counts and infeasibility flags).
- If `step_name='parsing-split-limit'`: increase `split_length_limit`; widen overlap bounds (`minimum_overlap_length`, `maximum_overlap_length`); lower `minimum_melting_temperature`; lower `minimum_hamming_distance`.
- If `step_name='parsing-variable-contig'`: your pool may be too conserved; loosen overlap/distance constraints or increase `split_length_limit`.
- If `step_name='computing-split'`: tune the same core constraints; for borderline cases, try a different `random_seed`.

### Analysis Issues

[^ Back to TOC](#table-of-contents)

**Debugging failed reads**
- Use `failed_reads_file` to sample discarded reads by failure category (anchor missing, barcode absent, barcode ambiguous, callback rejected, etc.)
- Inspect the CSV to understand why reads are being filtered out
- Common issues: wrong anchors, incorrect read orientation, barcode errors too strict

**Barcode ambiguous failures**
- Occurs when an anchor appears multiple times and different valid barcodes are found with equal confidence
- To avoid: design unique anchors using `motif(motif_type=1, ...)` with appropriate `maximum_repeat_length`

**Low barcode mapping rate**
- Check your anchors are correct and present in reads
- Try `mapping_type=1` (sensitive)
- Increase `barcode_errors` tolerance
- Verify read orientation (`r1_read_type`, `r2_read_type`)
- Use `failed_reads_file` to see examples of why reads are failing

**Missing combinations in `xcount`**
- Partial reads get `'-'` for missing barcodes (this is expected)
- Check if anchors for all indexes are present in reads
- Verify read length covers all barcode positions

### General

[^ Back to TOC](#table-of-contents)

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
<a href="https://github.com/ayaanhossain/oligopool">GitHub</a> -
<a href="https://pubs.acs.org/doi/10.1021/acssynbio.4c00661">Paper</a> -
<a href="https://github.com/ayaanhossain/oligopool/issues">Issues</a>
</p>
