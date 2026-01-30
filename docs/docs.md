<h1 align="center">
    <a href="https://github.com/ayaanhossain/oligopool/">
        <img src="https://raw.githubusercontent.com/ayaanhossain/repfmt/main/oligopool/img/logo.svg" alt="Oligopool Calculator" width="460"/>
    </a>
</h1>

<p align="center"><i>Your friendly neighborhood DNA library designer</i></p>

---

# Documentation

Welcome to the `Oligopool Calculator` docs! Whether you're designing your first barcode library or optimizing a million-variant MPRA, you're in the right place.

**TL;DR**: Design oligopool libraries with `barcode`, `primer`, `motif`, `spacer`. Analyze sequencing data with `index`, `pack`, `acount`/`xcount`. Most modules take CSV/DataFrames in and return `(out_df, stats)` (a few return stats only). Chain them together. Ship it.

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
- [Design Mode](#design-mode)
  - [barcode](#barcode) - Hamming-distance barcodes
  - [primer](#primer) - Thermodynamic primers
  - [motif](#motif) - Sequence motifs & anchors
  - [spacer](#spacer) - Neutral spacers
  - [background](#background) - K-mer screening database
  - [merge](#merge) - Collapse columns
  - [revcomp](#revcomp) - Reverse complement
  - [lenstat](#lenstat) - Length statistics
  - [verify](#verify) - QC check
  - [final](#final) - Finalize for synthesis
- [Assembly Mode](#assembly-mode)
  - [split](#split) - Fragment long oligos
  - [pad](#pad) - Assembly-ready padding
- [Degenerate Mode](#degenerate-mode)
  - [compress](#compress) - Compress to IUPAC-degenerate
  - [expand](#expand) - Expand degenerate oligos
- [Analysis Mode](#analysis-mode)
  - [index](#index) - Build barcode index
  - [pack](#pack) - Preprocess FastQ
  - [acount](#acount) - Association counting
  - [xcount](#xcount) - Combinatorial counting
- [Workflows](#workflows)
- [CLI Reference](#cli-reference)
- [Config Files](#config-files)
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

Most modules return both a DataFrame and a stats dictionary:

```python
out_df, stats = op.barcode(input_data=df, ...)
```

A few modules only return stats (no DataFrame to speak of):

```python
stats = op.verify(input_data=df, ...)
```

- **Input**: CSV path or pandas DataFrame with an `ID` column
- **Output**:
  - **Design/transform modules** return `(out_df, stats)`
  - **Stats-only modules** return `stats` (`background`, `lenstat`, `verify`, `index`, `pack`)
  - **Counting modules** return `(counts_df, stats)` (`acount`, `xcount`)
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

### String-Friendly `*_type` Parameters

Many `*_type` parameters accept either integer codes or descriptive strings (case-insensitive), e.g.
`barcode_type='terminus'`, `primer_type='forward'`, `motif_type='anchor'`, `pack_type='merge'`,
`mapping_type='sensitive'`. See the [API Reference](api.md) for the full alias list.

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

**Excluded motifs** let you ban restriction sites, repetitive sequences, or anything else you don't want appearing in your barcodes.

The simplest approach is to pass a list:

```python
df, stats = op.barcode(..., excluded_motifs=['GAATTC', 'GGATCC'])
```

You can also pass a CSV path, a DataFrame with an `Exmotif` column, or a FASTA file. See the
[API Reference](api.md#barcode) for details.

> **Note**: Excluded motifs are applied globally to all variants (no per-variant exclusion).

**Notes (the stuff that bites people):**
- You must provide at least one context column (`left_context_column` or `right_context_column`) so edge effects can be screened.
- Terminus vs spectrum: terminus optimized barcodes enforce distinctive 5'/3' ends; spectrum optimized barcodes saturate k-mers (slower but tighter). Higher `minimum_hamming_distance` buys you error tolerance but shrinks the design space.
- If you plan to index barcodes from reads, design constant anchors first with `motif(motif_type=1, ...)`.
- Patch mode (`patch_mode=True`) fills only missing values (`None`, `NaN`, `''`, `'-'`) and never overwrites existing barcodes (existing values must already be valid ATGC of length `barcode_length`).
- Cross-set separation is strict: set `cross_barcode_columns` and `minimum_cross_distance` together; the cross-set barcodes must be strict ATGC strings of length `barcode_length`.
- `excluded_motifs` works the same way in `primer`, `motif`, and `spacer`.

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

**Paired primer design** ensures your forward and reverse primers have matched Tm (within 1°C). Design the inner primer first:

```python
df, stats = op.primer(..., primer_column='FwdPrimer', primer_type=0)
```

Then design the outer primer, telling it to match Tm with the first:

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
    oligo_sets=['SetA', 'SetA', 'SetB', 'SetB'],
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
    background_directory='my_background',
)
```

**Notes (the stuff that bites people):**
- You must provide at least one context column (`left_context_column` or `right_context_column`) so edge effects can be screened.
- `maximum_repeat_length` controls non-repetitiveness against `input_data` only; screening against a background requires `background_directory`.
- If `paired_primer_column` is provided, the paired primer type is inferred and Tm matching is applied within 1°C.
- When `oligo_sets` is provided, primers are designed per set and screened for cross-set compatibility; if `paired_primer_column` is also provided, it must be constant within each set.
- Patch mode (`patch_mode=True`) preserves existing primer sequences and fills only missing values; in `oligo_sets` mode, existing per-set primers are reused and only missing-only sets trigger new primer design.

> **API Reference**: See [`primer`](api.md#primer) for complete parameter documentation.

---

### motif

[↑ Back to TOC](#table-of-contents)

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

**Why anchors matter**: When you index barcodes later, you need constant flanking sequences to locate them in reads. Design anchors with `motif_type=1` *before* designing barcodes.

**Notes (the stuff that bites people):**
- You must provide at least one context column (`left_context_column` or `right_context_column`) so edge effects can be screened.
- Constant bases in `motif_sequence_constraint` can force an excluded motif (or repeat) and make the design infeasible.
- For anchors, tune `maximum_repeat_length` to control how distinct the anchor is from surrounding sequence.
- Patch mode (`patch_mode=True`) fills only missing values; for `motif_type=1`, an existing compatible constant anchor (must be unique across existing rows) is reused for new rows.

> **API Reference**: See [`motif`](api.md#motif) for complete parameter documentation.

---

### spacer

[↑ Back to TOC](#table-of-contents)

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

You can also pass a CSV path or a DataFrame with `ID` and `Length` columns. See the
[API Reference](api.md#spacer) for details.

**Notes (the stuff that bites people):**
- You must provide at least one context column (`left_context_column` or `right_context_column`) so edge effects can be screened.
- If a row is already at (or over) `oligo_length_limit`, its spacer is set to `'-'` (a deliberate "no-space" sentinel).
- If `spacer_length` is a CSV/DataFrame, it must contain `ID` and `Length` columns aligned to your input IDs.
- Patch mode (`patch_mode=True`) fills only missing values and never overwrites existing spacers (some rows may still end up with `'-'` if no spacer can fit).

> **API Reference**: See [`spacer`](api.md#spacer) for complete parameter documentation.

---

### background

[↑ Back to TOC](#table-of-contents)

**What it does**: Builds a k-mer database for screening primers against off-target sequences.

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

Then use your background database during primer design:

```python
df, stats = op.primer(..., background_directory='ecoli_bg')
```

**Notes (the stuff that bites people):**
- `background(maximum_repeat_length=...)` screens against the background only; `primer(maximum_repeat_length=...)` screens against your input oligos.
- The background output directory ends with `.oligopool.background` and is what you pass as `background_directory`.
- If you need to inspect or modify the DB, use `vectorDB` (see the Advanced Modules section).

> **API Reference**: See [`background`](api.md#background) for complete parameter documentation.

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

**Notes (the stuff that bites people):**
- `left_context_column` and `right_context_column` do not need to be adjacent; `merge` collapses the full inclusive range.
- If you omit the bounds, `merge` collapses the first → last sequence columns.
- The merged source columns are removed from the output DataFrame.

> **API Reference**: See [`merge`](api.md#merge) for complete parameter documentation.

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

**Notes (the stuff that bites people):**
- `left_context_column` and `right_context_column` do not need to be adjacent; `revcomp` acts on the full inclusive range.
- If you omit the bounds, `revcomp` reverse-complements the first → last sequence columns and reverses their order.
- Useful mid-pipeline when you design in "readout orientation" but must synthesize in the opposite orientation (and for sanity-checking `split` fragment orientation).

> **API Reference**: See [`revcomp`](api.md#revcomp) for complete parameter documentation.

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

**Notes (the stuff that bites people):**
- `lenstat` is stats-only and does not modify or write your DataFrame (no `output_file`).
- It assumes all non-ID columns are DNA strings; if you have metadata columns or degenerate/IUPAC bases, use `verify` instead.
- Run it early and often—especially before `spacer`, `split`, and `pad`.

> **API Reference**: See [`lenstat`](api.md#lenstat) for complete parameter documentation.

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
- Degenerate/IUPAC bases
- Column architecture

**How columns are concatenated:**
- Only **sequence columns** (DNA/IUPAC; may include `'-'`) are concatenated; metadata columns are skipped
- Sequence columns are joined **left-to-right in DataFrame column order**
- Gap characters (`'-'`) are stripped (so `ABC-DEF` becomes `ABCDEF`)
- If `CompleteOligo` exists (from `final()`), it's used directly instead of concatenating

**Notes (the stuff that bites people):**
- `verify` is stats-only and never modifies or writes your DataFrame.
- Metadata columns are tracked and excluded from sequence-only checks; degenerate/IUPAC columns are flagged (not treated as hard errors).
- Excluded-motif checks report motif "emergence" in assembled oligos (a motif occurring more times than its minimum occurrence across the library), which is often what you care about in practice.
- Excluded-motif matching is literal substring matching; degenerate/IUPAC bases are not treated as wildcards (so degenerate columns can hide potential motifs).
- When emergent motifs are detected, `verify` also reports which column junctions contribute to motif emergence (helpful for debugging edge effects), attributing only occurrences beyond the baseline minimum. This requires separate sequence columns (run `verify` before `final`).
- Junction attribution follows column order: for `[Primer1, BC1, Variant, Primer2]`, junctions are `Primer1|BC1`, `BC1|Variant`, `Variant|Primer2`.

> **API Reference**: See [`verify`](api.md#verify) for complete parameter documentation.

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

**Notes (the stuff that bites people):**
- `final` returns only `CompleteOligo` and `OligoLength` (annotation columns are not preserved).
- Save your annotated design DataFrame/CSV before `final` if you plan to `index`, `pack`, and count later.
- `final` is typically the terminal design step before synthesis.

> **API Reference**: See [`final`](api.md#final) for complete parameter documentation.

---

## Assembly Mode

[↑ Back to TOC](#table-of-contents)

Assembly Mode provides tools for fragmenting long oligos that exceed synthesis length limits into overlapping pieces for assembly workflows.

### split

[↑ Back to TOC](#table-of-contents)

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

**Key concept: Each column = one oligo pool to synthesize**

When `split` returns `Split1`, `Split2`, `Split3` columns, you will **order three separate oligo pools** from your synthesis vendor (one per column). After synthesis, these fragments are combined via overlap-based assembly (Gibson, overlap-extension PCR, etc.) to reconstruct the full-length oligo.

**Notes (the stuff that bites people):**
- The number of fragments is automatically chosen and can vary per oligo (some rows may not have all `Split*` columns populated).
- Fragments are returned in PCR/assembly order; even-numbered split columns (`Split2`, `Split4`, ...) are reverse-complemented by design. This alternating orientation enables efficient PCR-based assembly where adjacent fragments anneal via their overlapping regions.
- `split` returns only the split fragments (your original annotation columns are not preserved). Save the annotated library separately if you need it later.
- As a rule of thumb, keep `minimum_overlap_length` comfortably larger than `minimum_hamming_distance`.
- **Raw split output is not synthesis-ready for assembly** — run `pad` per `SplitN` to add primers/Type IIS sites, then `final` each padded DataFrame.
- For separate per-fragment outputs: in Python, enable `separate_outputs`; in CLI, this is already the default (writes `out.Split1.oligopool.split.csv`, etc.). Use `--no-separate-outputs` in CLI for a single combined file.

> **API Reference**: See [`split`](api.md#split) for complete parameter documentation.

---

### pad

[↑ Back to TOC](#table-of-contents)

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

**Notes (the stuff that bites people):**
- **Run `pad` separately for each fragment** (e.g., `Split1`, `Split2`, ...). You cannot pad all columns in one call.
- **The chosen Type IIS recognition site must be absent from your split fragments** (in either orientation). If fragments contain internal cut sites, they will be cleaved during digest. `pad` checks this and fails early if conflicts are found.
- **Exclude your Type IIS motif from upstream design elements.** If you plan to use `BsaI` for padding, add its recognition site (`GGTCTC`) and reverse complement (`GAGACC`) to `excluded_motifs` when designing primers, barcodes, motifs, and spacers. This helps prevent internal cut sites from showing up in newly designed elements (it doesn't "clean" a core/variant sequence that already contains the site — in that case, choose a different enzyme or redesign the offending sequences).
- Output columns are `5primeSpacer`, `ForwardPrimer`, `<split_column>`, `ReversePrimer`, `3primeSpacer` (other `split_df` columns are not preserved).
- If a fragment can't fit under `oligo_length_limit`, the spacer(s) for that row are set to `'-'` (a deliberate "no-space" sentinel).
- The supported Type IIS enzymes list is the "batteries included" set for pad removal.

**Post-synthesis workflow** (after you get your oligos back):
1. **PCR amplify** the padded fragments
2. **Type IIS digest** → excises primers and spacers, leaving enzyme-specific sticky overhangs (the overhangs are primer-derived, not from your fragment)
3. **Mung bean nuclease** (optional) → blunts the overhangs; skip this step if using a blunt-cutter like `MlyI`
4. **Assemble** the cleaned fragments via their **split-designed overlaps** (Gibson, overlap-extension PCR, etc.)

The Type IIS enzyme is for **pad removal**, not fragment-to-fragment ligation — the 15–30 bp overlaps from `split` drive the actual assembly.

**Supported Type IIS enzymes (34 total):**
`AcuI`, `AlwI`, `BbsI`, `BccI`, `BceAI`, `BciVI`, `BcoDI`, `BmrI`, `BpuEI`, `BsaI`, `BseRI`, `BsmAI`, `BsmBI`, `BsmFI`, `BsmI`, `BspCNI`, `BspQI`, `BsrDI`, `BsrI`, `BtgZI`, `BtsCI`, `BtsI`, `BtsIMutI`, `EarI`, `EciI`, `Esp3I`, `FauI`, `HgaI`, `HphI`, `HpyAV`, `MlyI`, `MnlI`, `SapI`, `SfaNI`

For other enzymes, design primers manually with the appropriate recognition/cut sites using `primer` or `motif`.

> **API Reference**: See [`pad`](api.md#pad) for complete parameter documentation.

---

## Degenerate Mode

[↑ Back to TOC](#table-of-contents)

Degenerate Mode enables cost-efficient synthesis of variant libraries with low mutational diversity. Instead of synthesizing each variant individually, similar sequences are compressed into IUPAC-degenerate oligos that expand to cover multiple variants.

**When to use Degenerate Mode**: You have a large library of similar sequences and want to reduce synthesis costs by grouping compatible variants. Best for selection assays where you identify winners by sequencing (no barcode counting required). ML-generated libraries and saturation mutagenesis libraries often compress well.

### compress

[↑ Back to TOC](#table-of-contents)

**What it does**: Compresses concrete DNA sequences into IUPAC-degenerate representation.

**When to use it**: You have many similar sequences and want fewer oligos to synthesize.

```python
mapping_df, synthesis_df, stats = op.compress(
    input_data=df,                    # DataFrame with ID + DNA columns
    mapping_file='mapping',           # Links variants to degenerate oligos
    synthesis_file='synthesis',       # Degenerate oligos for ordering
    random_seed=42,
)

# Check compression results
print(f"Compressed {stats['vars']['input_variants']} variants "
      f"into {stats['vars']['degenerate_oligos']} oligos "
      f"({stats['vars']['compression_ratio']}x compression)")
```

**Notes (the stuff that matters):**
- Input must be concrete DNA (A/T/G/C only, no degenerate codes)
- All non-ID columns are concatenated (left-to-right) to form the full sequence
- Include only columns you intend to make degenerate (e.g., diverse barcodes can prevent good compression)
- Similar sequences compress well; diverse sequences may not compress at all (1:1 mapping)
- The algorithm guarantees no invented sequences (lossless compression)

> **API Reference**: See [`compress`](api.md#compress) for complete parameter documentation.

---

### expand

[↑ Back to TOC](#table-of-contents)

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

**Notes (the stuff that matters):**
- Expansion can be exponential (N positions = 4^N for all-N sequence)
- Use `expansion_limit` as a safety cap for highly degenerate sequences
- Does NOT recover original variant IDs; use `mapping_df` from `compress` for that

> **API Reference**: See [`expand`](api.md#expand) for complete parameter documentation.

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
    associate_suffix_gap=0,                # Bases between variant and suffix anchor in read
)
```

**Notes (the stuff that bites people):**
- You must specify at least one of `barcode_prefix_column` or `barcode_suffix_column`; anchors should be constant (single-unique) sequences and ideally adjacent to the barcode/associate.
- `barcode_prefix_gap`/`barcode_suffix_gap` specify how many bases separate the anchor and barcode in the read (real sequencing is messy; this makes it configurable).
- If your associate anchors are not directly adjacent, use `associate_prefix_gap`/`associate_suffix_gap` to specify the exact gaps in the read.
- For association counting, partial presence of the associate sequence can be sufficient, but the associate anchors must be adjacent and fully present.
- You can use multiple index files with `xcount` for combinatorial counting (associate info is ignored in that mode).

> **API Reference**: See [`index`](api.md#index) for complete parameter documentation.

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
    pack_type=1,                           # 0=concatenate pairs, 1=merge pairs
    pack_file='sample',                    # Creates sample.oligopool.pack
)
```

**Notes (the stuff that bites people):**
- For paired-end workflows, both reads must pass the length/quality filters for the pair to be retained.
- If reads are merged externally, pass the merged reads as single-end (R1 only) and leave all R2 args as `None`.
- `pack_type=0` (concatenated) tends to be IO-bound; `pack_type=1` (merged) is more compute-heavy. Pack files are reusable across counting runs.
- Internally, pack files store reads as `(r1, r2)` tuples; merged/single-end reads use `r2=None` (relevant for Python callbacks).

> **API Reference**: See [`pack`](api.md#pack) for complete parameter documentation.

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

**Notes (the stuff that bites people):**
- Reads with unresolved associates are excluded from the counts (that's the point of `acount`).
- `callback` is Python-only; the CLI always runs with `callback=None`.
- If you want to inspect *examples* of discarded reads, set `failed_reads_file` to write a small per-category diagnostic CSV.
- `acount` operates on a single index + pack pair (for combinatorial counting, use `xcount`).

> **API Reference**: See [`acount`](api.md#acount) for complete parameter documentation.

---

### xcount

[↑ Back to TOC](#table-of-contents)

**What it does**: Barcode-focused counting - maps reads to one or more barcode indices without associate verification.

**When to use it**: Pure barcode counting (single or multiple indices), combinatorial barcode designs (BC1 × BC2).

For single barcode counting, point to one index:

```python
df, stats = op.xcount(
    index_files=['bc1_index'],
    pack_file='sample',
    count_file='barcode_counts',
)
```

For combinatorial counting (BC1 × BC2), list multiple indices:

```python
df, stats = op.xcount(
    index_files=['bc1_index', 'bc2_index'],
    pack_file='sample',
    count_file='combo_counts',
    mapping_type=1,
)
```

Output includes all observed combinations. Reads missing a barcode show `'-'` for that position.

**acount vs xcount**: Use `acount` when you need barcode+variant association verification; use `xcount` for barcode-only counting (single or combinatorial).

**Notes (the stuff that bites people):**
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

## Workflows

[↑ Back to TOC](#table-of-contents)

### Basic Library Design

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
op.verify(input_data=df, oligo_length_limit=200)

df.to_csv('library_design.csv', index=False)  # Keep this for indexing later!
final_df, _ = op.final(input_data=df, output_file='library_for_synthesis')
```

### Analysis Pipeline

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

[↑ Back to TOC](#table-of-contents)

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

**CLI outputs require an output basename**: Unlike the Python API (where many modules can return results in-memory), CLI commands that write files require an output basename (e.g., `--output-file`, `--index-file`, `--pack-file`, `--count-file`, `--mapping-file`, `--synthesis-file`).

**Output filenames are auto-suffixed**: Commands append a suffix if missing (e.g., `.oligopool.barcode.csv`), so prefer basenames like `--output-file my_library` to avoid doubled extensions.

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

> **API Reference**: See [api.md](api.md#cli-parameter-mapping) for complete CLI parameter mapping.

---

## Config Files

[↑ Back to TOC](#table-of-contents)

The CLI supports YAML config files for repeatable, documented workflows. Config files eliminate long command lines and enable multi-step pipeline execution.

### Single Command Config

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
  output_file: "library_barcoded.csv"
  barcode_type: "spectrum"
  excluded_motifs:
    - "GGATCC"
    - "TCTAGA"
```

Config keys use snake_case (same as Python API). CLI arguments override config values:

```bash
# Override barcode_length from config
op barcode --config barcode_design.yaml --barcode-length 20
```

### Pipeline Execution

Run multi-step workflows from a single config with `op pipeline`:

```bash
op pipeline --config mpra_pipeline.yaml
```

Example pipeline config (`mpra_pipeline.yaml`):

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
  output_file: "step1_primer.csv"
  oligo_length_limit: 250
  primer_sequence_constraint: "N*22"
  primer_type: "forward"
  primer_column: "FwdPrimer"
  minimum_melting_temperature: 52
  maximum_melting_temperature: 58
  maximum_repeat_length: 10

barcode:
  input_data: "step1_primer.csv"
  output_file: "step2_barcode.csv"
  oligo_length_limit: 250
  barcode_length: 16
  minimum_hamming_distance: 3
  maximum_repeat_length: 8
  barcode_column: "BC"
  left_context_column: "FwdPrimer"

spacer:
  input_data: "step2_barcode.csv"
  output_file: "step3_spacer.csv"
  oligo_length_limit: 250
  maximum_repeat_length: 8
  spacer_column: "Spacer"
  left_context_column: "BC"

final:
  input_data: "step3_spacer.csv"
  output_file: "final_library.csv"
```

Each step's config section specifies explicit `input_data` and `output_file` paths, giving you full control over the data flow.

### Parallel Pipeline Execution

For workflows with independent branches, use the parallel (DAG) format to run steps concurrently:

```yaml
# parallel_pipeline.yaml
pipeline:
  name: "Parallel Design"
  steps:
    - name: fwd_primer
      command: primer
    - name: rev_primer
      command: primer
      # No 'after' - runs in parallel with fwd_primer
    - name: add_barcode
      command: barcode
      after: [fwd_primer]  # Waits for fwd_primer only
    - name: finalize
      command: final
      after: [add_barcode, rev_primer]  # Waits for both

fwd_primer:
  input_data: "variants.csv"
  output_file: "fwd_primer.csv"
  primer_type: "forward"
  primer_sequence_constraint: "N*20"
  # ...

rev_primer:
  input_data: "variants.csv"
  output_file: "rev_primer.csv"
  primer_type: "reverse"
  primer_sequence_constraint: "N*20"
  # ...

add_barcode:
  input_data: "fwd_primer.csv"
  output_file: "with_barcode.csv"
  # ...

finalize:
  input_data: "with_barcode.csv"
  output_file: "final.csv"
```

**Step fields:**
- `name`: Step identifier (required)
- `command`: Oligopool command to run (defaults to `name`)
- `after`: List of step names this step depends on (optional)
- `config`: Config section name to use (defaults to `name`)

**Execution model:**
- Steps with no `after` dependencies form level 1 (run in parallel)
- Each subsequent level waits for its dependencies
- `--dry-run` shows execution levels and parallelism

```bash
op pipeline --config parallel_pipeline.yaml --dry-run
# Output shows:
#   Level 1: fwd_primer, rev_primer (parallel)
#   Level 2: add_barcode
#   Level 3: finalize
```

### Dry Run Validation

Validate a pipeline config without executing:

```bash
op pipeline --config mpra_pipeline.yaml --dry-run
```

This checks that all steps are valid commands and displays the parameters for each step.

### Config Precedence

Parameter values are resolved in this order (highest priority first):

1. **CLI arguments** - Always win
2. **Config file values** - Fill in unset args
3. **Command defaults** - Built-in defaults

### Config Tips

1. **Use comments liberally**: YAML supports `#` comments - document your design choices.
2. **Keep configs in version control**: Reproducibility for future you.
3. **Explicit paths over magic**: Pipeline steps use explicit `input_data`/`output_file` paths - no implicit chaining.
4. **Lists in YAML**: Use YAML list syntax for multi-value parameters like `excluded_motifs`.
5. **Type aliases work**: Use `"spectrum"` or `1` for `barcode_type` - same as CLI.

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

LevelDB-based k-mer storage. Created by `background()`, but you can access it directly.

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
