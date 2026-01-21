<h1 align="center">
    <a href="https://github.com/ayaanhossain/oligopool/">
        <img src="https://raw.githubusercontent.com/ayaanhossain/repfmt/main/oligopool/img/logo.svg" alt="Oligopool Calculator" width="460"/>
    </a>
</h1>

<p align="center"><i>Your friendly neighborhood DNA library designer</i></p>

---

# Documentation

Welcome to the Oligopool Calculator docs! Whether you're designing your first barcode library or optimizing a million-variant MPRA, you're in the right place.

**TL;DR**: Design oligopool libraries with `barcode`, `primer`, `motif`, `spacer`. Analyze sequencing data with `index`, `pack`, `acount`/`xcount`. All modules take DataFrames in, spit DataFrames out. Chain them together. Ship it.

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

---

## Quick Start

[↑ Back to TOC](#table-of-contents)

```python
import pandas as pd
import oligopool as op

# Start with your variants
df = pd.DataFrame({
    'ID': ['V1', 'V2', 'V3'],
    'Promoter': ['ATGCATGC...', 'GCTAGCTA...', 'TTAATTAA...']
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

Every module follows the same pattern:

```python
output_df, stats = op.module_name(input_data=input_df, ...)
```

- **Input**: CSV path or pandas DataFrame with an `ID` column
- **Output**: Updated DataFrame + stats dictionary
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

### Context Columns

Most design modules need to know what's next to the element being designed:

```python
df, stats = op.barcode(
    ...,
    left_context_column='Promoter',   # What's on the left?
    right_context_column='Gene'       # What's on the right?
)
```

This prevents repeat collisions at element boundaries. Edge effects are real, folks.

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

**Pro tips**:
- Start with `barcode_type=0` (fast). Switch to `1` only if you're packing barcodes tight.
- Higher Hamming distance = more error tolerance, fewer possible barcodes.
- Use `cross_barcode_columns` + `minimum_cross_distance` for multi-barcode designs (BC1 vs BC2).

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
df, stats = op.primer(
    ...,
    oligo_sets=['SetA', 'SetA', 'SetB', 'SetB'],  # Group labels
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

# Per-variant lengths
df, stats = op.spacer(
    ...,
    spacer_length=[10, 15, 12, 8],         # Different per row
)
```

---

### background

[↑ Back to TOC](#table-of-contents)

**What it does**: Builds a k-mer database for screening primers against off-target sequences.

**When to use it**: You want primers that won't bind to host genome, plasmid backbone, etc.

```python
stats = op.background(
    input_data='ecoli_genome.csv',         # CSV with ID, Sequence columns
    maximum_repeat_length=15,              # Screen k-mers up to this length
    output_directory='ecoli_bg',           # Creates ecoli_bg.oligopool.background
)

# Now use it in primer design
df, stats = op.primer(
    ...,
    background_directory='ecoli_bg',
)
```

---

### split

[↑ Back to TOC](#table-of-contents)

**What it does**: Breaks long oligos into overlapping fragments for assembly (Gibson, Golden Gate, etc.).

**When to use it**: Your oligos exceed synthesis length limits (~200 bp).

```python
df, stats = op.split(
    input_data=df,
    split_length_limit=150,               # Max fragment length
    minimum_melting_temperature=55.0,     # Overlap Tm
    minimum_hamming_distance=3,           # Overlap uniqueness
    minimum_overlap_length=20,
    maximum_overlap_length=30,
    output_file='split_library.csv',
)
```

Output contains `Split-1`, `Split-2`, etc. columns with fragments ready for assembly.

---

### pad

[↑ Back to TOC](#table-of-contents)

**What it does**: Adds primers and Type IIS restriction sites to split fragments for scarless assembly.

**When to use it**: After `split`, to make fragments synthesis-ready and assembly-compatible.

```python
df, stats = op.pad(
    input_data=split_df,
    split_column='Split-1',               # Which fragment to pad
    typeIIS_system='BsaI',                # Enzyme for Golden Gate
    oligo_length_limit=200,
    minimum_melting_temperature=52.0,
    maximum_melting_temperature=58.0,
    maximum_repeat_length=10,
    output_file='padded_split1.csv',
)
```

Supports 34 Type IIS enzymes: BsaI, BbsI, BsmBI, SapI, and many more.

---

### merge

[↑ Back to TOC](#table-of-contents)

**What it does**: Concatenates contiguous columns into a single column.

**When to use it**: Mid-pipeline cleanup, combining elements before further processing.

```python
df, stats = op.merge(
    input_data=df,
    left_column='Primer5',
    right_column='BC1',                   # Merges Primer5 through BC1
    merge_column='5prime_region',         # Output column name
)
```

---

### revcomp

[↑ Back to TOC](#table-of-contents)

**What it does**: Reverse complements a range of columns and reverses their order.

**When to use it**: Switching strand orientation (e.g., synthesis vs. sequencing direction).

```python
df, stats = op.revcomp(
    input_data=df,
    left_column='Gene',
    right_column='Primer3',               # RevComp everything from Gene to Primer3
)
```

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

---

### final

[↑ Back to TOC](#table-of-contents)

**What it does**: Concatenates all columns into synthesis-ready oligos.

**When to use it**: Last step before ordering.

```python
df, stats = op.final(
    input_data=df,
    output_file='synthesis_ready.csv',
)
# Adds 'CompleteOligo' and 'OligoLength' columns
```

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

---

### acount

[↑ Back to TOC](#table-of-contents)

**What it does**: Association counting - maps reads to barcodes AND their associated variants.

**When to use it**: When you need barcode-to-variant mapping (e.g., variant activity assays).

```python
df, stats = op.acount(
    index_file='bc1_index',
    pack_file='sample.oligopool.pack',
    count_file='counts',                   # Creates counts.oligopool.acount.csv
    barcode_errors=1,                      # Allow 1 mismatch (-1=auto)
)
```

---

### xcount

[↑ Back to TOC](#table-of-contents)

**What it does**: Combinatorial counting - maps reads to multiple barcode combinations.

**When to use it**: Multi-barcode designs (BC1 × BC2), combinatorial libraries.

```python
df, stats = op.xcount(
    index_files=['bc1_index', 'bc2_index'],  # Multiple indices
    pack_file='sample.oligopool.pack',
    count_file='combo_counts',
    mapping_type=1,                          # 0=fast, 1=sensitive
    barcode_errors=-1,                       # Auto-infer
)
```

Output includes all observed combinations, with `'-'` for missing barcodes in partial reads.

**Custom callbacks** (Python API only):
```python
def my_filter(read, ID, count, coreid):
    # Custom logic here
    return True  # Accept read

df, stats = op.xcount(..., callback=my_filter)
```

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
    left_context_column='Variant',
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
    left_context_column='BC1',
    paired_primer_column='FwdPrimer',
)

# 5. Check length
op.lenstat(input_data=df, oligo_length_limit=200)

# 6. Verify and finalize
op.verify(input_data=df, oligo_length_limit=200)
df, _ = op.final(input_data=df, output_file='library.csv')
```

### Analysis Pipeline

```python
import oligopool as op

# 1. Index your barcodes
op.index(
    barcode_data='library.csv',
    barcode_column='BC1',
    barcode_prefix_column='FwdPrimer',
    index_file='bc1',
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
    index_file='bc1',
    pack_file='experiment.oligopool.pack',
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
op barcode --help

# Read the manual
op manual barcode

# Run a command
op barcode \
    --input-data variants.csv \
    --oligo-length-limit 200 \
    --barcode-length 12 \
    --minimum-hamming-distance 3 \
    --maximum-repeat-length 8 \
    --barcode-column BC1 \
    --left-context-column Variant \
    --output-file with_barcodes.csv
```

**Install tab completion** (highly recommended):
```bash
op complete --install
```

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

<p align="center">
<i>Made with molecules and math by the Salis Lab</i>
<br>
<a href="https://github.com/ayaanhossain/oligopool">GitHub</a> ·
<a href="https://pubs.acs.org/doi/10.1021/acssynbio.4c00661">Paper</a> ·
<a href="https://github.com/ayaanhossain/oligopool/issues">Issues</a>
</p>
