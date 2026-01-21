# Oligopool Calculator - AI Agent Guide

This document is designed for AI assistants to understand and help users with oligopool library design and analysis.

## Package Overview

**Oligopool Calculator** is a Python library for designing and analyzing oligonucleotide pool (oligopool) libraries used in massively parallel reporter assays (MPRAs) and similar high-throughput experiments.

**Core workflow:**
1. **Design Mode**: Build a library of oligos with barcodes, primers, motifs, and spacers
2. **Analysis Mode**: Count barcoded reads from NGS data to quantify variant activity

## Key Concepts

### DataFrame-Centric Design
- All modules operate on pandas DataFrames
- Input: DataFrame with `ID` column + DNA sequence columns
- Output: `(updated_df, stats_dict)` tuple
- Chain modules to build libraries iteratively

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

### Context Columns
Most modules require `left_context_column` and/or `right_context_column` to:
- Prevent repeats at element boundaries (edge effects)
- Ensure excluded motifs don't emerge at junctions
- At least one context column is required

### Patch Mode
Use `patch_mode=True` to extend existing pools:
- Only fills missing values (None/NaN/empty/'-')
- Preserves existing designs
- Useful for iterative pool expansion

## Design Mode Modules

### barcode
**Purpose**: Generate Hamming-distance separated barcodes for unique variant identification.

**Key parameters**:
- `barcode_length`: Length of barcodes (typically 8-16 bp)
- `minimum_hamming_distance`: Minimum mismatches between any two barcodes (typically 2-4)
- `barcode_type`: 0=terminus-optimized (fast), 1=spectrum-optimized (thorough)
- `cross_barcode_columns` + `minimum_cross_distance`: Enforce separation from existing barcode sets

**When to use**: Always needed for variant identification in counting.

**Design order**: After primers/motifs, before spacers.

### primer
**Purpose**: Design thermodynamically optimal primers for amplification.

**Key parameters**:
- `primer_sequence_constraint`: IUPAC constraint (e.g., `'SS' + 'N'*18` for GC clamp)
- `primer_type`: 0=forward, 1=reverse
- `minimum/maximum_melting_temperature`: Tm range (typically 53-55°C)
- `paired_primer_column`: For Tm-matched primer pairs
- `background_directory`: Screen against off-target sequences
- `oligo_sets`: Design set-specific primers for multiplexed libraries

**When to use**: When amplification is needed.

**Design order**: Design primers early. For paired primers: design inner primer first, then outer with `paired_primer_column`.

### motif
**Purpose**: Add sequence motifs or constant anchors.

**Key parameters**:
- `motif_sequence_constraint`: IUPAC pattern or constant sequence
- `motif_type`: 0=per-variant motifs, 1=constant anchor shared by all

**When to use**:
- Restriction sites (e.g., `'GAATTC'` for EcoRI)
- Barcode anchors for indexing (design before barcodes)
- Degenerate regions with specific constraints

**Design order**: Before barcodes if designing anchors.

### spacer
**Purpose**: Add neutral filler DNA to reach target oligo length.

**Key parameters**:
- `spacer_length`: Fixed length, per-variant list/DataFrame, or `None` for auto-fill to `oligo_length_limit`

**When to use**: As final design step to equalize oligo lengths.

**Design order**: Last, after all other elements.

### background
**Purpose**: Build k-mer database for primer off-target screening.

**Key parameters**:
- `input_data`: Sequences to screen against (list, CSV, DataFrame, or FASTA file)
- `maximum_repeat_length`: k-mer size for screening (6-20)

**When to use**: Before primer design when you have host genome/plasmid sequences.

### split
**Purpose**: Break long oligos into overlapping fragments for assembly.

**Key parameters**:
- `split_length_limit`: Maximum fragment length
- `minimum/maximum_overlap_length`: Overlap range for assembly
- `minimum_melting_temperature`: Assembly Tm

**When to use**: When oligos exceed synthesis length limits.

### pad
**Purpose**: Add amplification primers with Type IIS sites for assembly.

**Key parameters**:
- `split_column`: Which fragment to pad
- `typeIIS_system`: Restriction enzyme (e.g., 'BsaI')

**When to use**: After split, for each fragment.

### Utility Modules

- **lenstat**: Check length statistics and remaining space (non-destructive)
- **verify**: QC check before synthesis (validates constraints, flags issues)
- **final**: Concatenate columns into synthesis-ready oligos
- **merge**: Collapse multiple columns into one
- **revcomp**: Reverse complement a column range

## Analysis Mode Modules

### index
**Purpose**: Build barcode index for counting.

**Key parameters**:
- `barcode_column`: Which barcode to index
- `barcode_prefix/suffix_column`: Constant anchors flanking barcode
- `associate_column`: Variant column for association counting (optional)

**When to use**: Before counting, for each barcode set.

### pack
**Purpose**: Preprocess and deduplicate FASTQ reads.

**Key parameters**:
- `r1/r2_fastq_file`: Input FASTQ files (supports gzip)
- `pack_type`: 0=concatenate, 1=merge paired reads
- `minimum_r1/r2_read_quality`: Phred score filter

**When to use**: After sequencing, before counting.

### acount
**Purpose**: Association counting - verify barcode-variant coupling.

**When to use**:
- Validating synthesis accuracy
- Confirming barcode-variant associations
- Input library QC

**Key insight**: Counts reads where the expected variant IS coupled with its barcode.

### xcount
**Purpose**: Barcode-only counting (single or combinatorial).

**When to use**:
- Pure barcode counting without variant verification
- Multi-barcode combinations (BC1 × BC2)
- Cleaved vs uncleaved product quantification

**Key insight**: No associate verification; works with one or more barcode indices.

## Common Workflows

### Basic Library Design
```python
import oligopool as op
import pandas as pd

# 1. Start with variants
df = pd.DataFrame({
    'ID': ['v1', 'v2', ...],
    'Variant': ['ATGC...', 'GCTA...', ...]
})

# 2. Add primers (inner first if paired)
df, _ = op.primer(input_data=df, primer_column='Primer1', ...)
df, _ = op.primer(input_data=df, primer_column='Primer2', paired_primer_column='Primer1', ...)

# 3. Add barcodes
df, _ = op.barcode(input_data=df, barcode_column='BC1', ...)

# 4. Add spacers to reach target length
df, _ = op.spacer(input_data=df, spacer_column='Spacer', spacer_length=None, ...)

# 5. QC and finalize
_ = op.verify(input_data=df, ...)
final_df, _ = op.final(input_data=df)
```

### Analysis Pipeline
```python
# 1. Index barcodes
op.index(barcode_data=df, barcode_column='BC1', index_file='bc1_idx', ...)

# 2. Pack reads
op.pack(r1_fastq_file='reads_R1.fq.gz', pack_file='reads', ...)

# 3. Count (choose based on need)
# For barcode-variant verification:
counts_df, _ = op.acount(index_file='bc1_idx', pack_file='reads.oligopool.pack', ...)

# For barcode-only or combinatorial:
counts_df, _ = op.xcount(index_files=['bc1_idx', 'bc2_idx'], pack_file='reads.oligopool.pack', ...)
```

## Troubleshooting Guide

### "Design failed" / status=False
1. Check `stats['basis']` for failure reason
2. Common fixes:
   - Increase `barcode_length` or decrease `minimum_hamming_distance`
   - Relax `excluded_motifs`
   - Increase `maximum_repeat_length`
   - Widen Tm range for primers

### Oligo too long
- Use `lenstat` to check remaining space
- Reduce element lengths
- Use `split` + `pad` for assembly

### No barcodes mapping
- Verify anchor sequences in `index` match designed anchors
- Check `barcode_prefix/suffix_gap` settings
- Use `mapping_type=1` for sensitive mode

### Cross-contamination in counts
- Increase `minimum_hamming_distance`
- Add `cross_barcode_columns` for multi-barcode designs

## Parameter Quick Reference

### Excluded Motifs
All element modules accept `excluded_motifs`:
- List: `['GAATTC', 'GGATCC']`
- CSV: File with 'Exmotif' column
- DataFrame: With 'Exmotif' column
- FASTA: .fa/.fasta/.fna (optionally gzipped)

### Context Columns
- Always specify at least one of `left_context_column` or `right_context_column`
- Use adjacent columns in the oligo architecture
- Critical for preventing edge-effect motifs

### Random Seed
- Use `random_seed=42` (or any int) for reproducible designs
- Important for version control and debugging

## File Formats

### Input CSV
```csv
ID,Variant
v1,ATGCATGC...
v2,GCTAGCTA...
```

### Output Files
- `.oligopool.background` - k-mer database directory
- `.oligopool.index` - barcode index file
- `.oligopool.pack` - packed reads file
- `.oligopool.acount.csv` - association counts
- `.oligopool.xcount.csv` - combinatorial counts

## CLI Usage

```bash
# Design
op barcode --input-data variants.csv --barcode-length 12 --minimum-hamming-distance 3 --output-file result.csv

# Analysis
op index --barcode-data library.csv --barcode-column BC1 --index-file bc1_idx
op pack --r1-fastq-file reads.fq.gz --pack-file reads
op xcount --index-files bc1_idx --pack-file reads.oligopool.pack --count-file counts
```

## Best Practices

1. **Design order**: background → primers → motifs/anchors → barcodes → spacers → verify → final
2. **Use lenstat frequently**: Check remaining space after each element
3. **Save intermediate DataFrames**: Enable rollback if design fails
4. **Set random_seed**: For reproducible designs
5. **Run verify before synthesis**: Catch constraint violations early
6. **Use Patch Mode for extensions**: Don't redesign existing elements

## Links

- Repository: https://github.com/ayaanhossain/oligopool
- Documentation: See `docs.md` in repository
- Paper: https://doi.org/10.1021/acssynbio.4c00661
