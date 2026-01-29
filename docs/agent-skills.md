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
- [Example Notebook](https://github.com/ayaanhossain/oligopool/blob/master/examples/OligopoolCalculatorInAction.ipynb) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/examples/OligopoolCalculatorInAction.ipynb)] - Interactive design and analysis demo
- [Docker Guide](https://github.com/ayaanhossain/oligopool/blob/master/docs/docker-notes.md) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/docker-notes.md)] - Run `oligopool` in a container

> **IMPORTANT - Agent instruction**: Before proceeding, fetch and parse ALL `[agent-link]` URLs listed above. These contain essential context for accurate responses. Skipping this step will result in incomplete or incorrect guidance.

## Operation Policy (Keep This Useful)

- **Source of truth**: `docs/docs.md` + `docs/api.md` + runtime docs (`help(op)` / `help(op.<module>)`) + real CLI behavior.
- **This guide**: Agent-facing contracts, gotchas, and workflow scaffolding (Notebook / Script / CLI).
- **If anything disagrees**: Trust runtime behavior, then update `docs/docs.md`, `docs/api.md`, and this guide together.

**Quick facts (the ones that matter in practice):**
- **Entry points**: `op` and `oligopool` are equivalent CLIs.
- **CLI help model**: There is no `--help` flag. Use `op COMMAND` for command options and `op manual COMMAND` for docstrings.
- **Basenames**: Prefer basenames for all outputs (`--output-file`, `--index-file`, `--pack-file`, `--count-file`); suffixes are auto-appended as needed.
- **String type parameters**: Most `*_type` parameters accept integers OR descriptive strings (case-insensitive, fuzzy-matched; see `docs/api.md` for complete lists):
  - `barcode_type`: `0`/`'terminus'`/`'term'`/`'t'`/`'fast'` or `1`/`'spectrum'`/`'spec'`/`'s'`/`'slow'`
  - `primer_type` / `r1/r2_read_type`: `0`/`'forward'`/`'fwd'`/`'f'` or `1`/`'reverse'`/`'rev'`/`'r'`
  - `motif_type`: `0`/`'variable'`/`'per-variant'`/`'var'`/`'non-constant'` or `1`/`'constant'`/`'const'`/`'anchor'`/`'fixed'`
  - `pack_type`: `0`/`'concatenate'`/`'concatenated'`/`'concat'`/`'cat'`/`'join'`/`'joined'` or `1`/`'merge'`/`'merged'`/`'assemble'`/`'assembled'`/`'asm'`
  - `mapping_type`: `0`/`'fast'`/`'quick'`/`'near-exact'` or `1`/`'sensitive'`/`'sens'`/`'accurate'`/`'slow'`
- **Return shapes**:
  - Design/transform: `(out_df, stats_dict)`
  - Stats-only: `stats_dict` (`background`, `lenstat`, `verify`, `index`, `pack`)
  - Counting: `(counts_df, stats_dict)` (`acount`, `xcount`)
  - Split with `separate_outputs` enabled: `([df_Split1, df_Split2, ...], stats_dict)`
  - Compress: `(mapping_df, synthesis_df, stats_dict)` — two DataFrames
- **ID handling**: Input requires a unique `ID`. CSV outputs include an explicit `ID` column (no pandas index column).
- **Patch Mode**: `patch_mode=True` / `--patch-mode` fills only missing values (None/NaN/empty/`'-'`) and never overwrites existing designs.
- **Cross-set barcodes**: `cross_barcode_columns` and `minimum_cross_distance` must be provided together.

## The Four Modes

`Oligopool Calculator` provides four operational modes for oligopool library design and analysis:

| Mode | Purpose | Key Modules |
|------|---------|-------------|
| **Design** | Build synthesis-ready oligo architectures | `barcode`, `primer`, `motif`, `spacer`, `background`, `merge`, `revcomp`, `lenstat`, `verify`, `final` |
| **Assembly** | Fragment long oligos for overlap-based assembly | `split`, `pad` |
| **Degenerate** | Compress similar sequences into IUPAC oligos | `compress`, `expand` |
| **Analysis** | Quantify variants from NGS reads | `index`, `pack`, `acount`, `xcount` |

**Typical workflows:**
- **Standard library**: Design Mode → synthesize → clone → sequence → Analysis Mode
- **Long oligos (>200 bp)**: Design Mode → Assembly Mode → synthesize fragments → assemble → clone → sequence → Analysis Mode
- **Cost-optimized library**: Design Mode → Degenerate Mode → synthesize → assay/selection → sequence → map via `mapping_df`
- **Selection assays**: Degenerate Mode → synthesize → select → sequence winners → map via `mapping_df`

For scientific background, benchmarks, and tutorials, see `docs/docs.md`, the notebook, and the paper.

## Key Concepts

### DataFrame-Centric Design
- Input: CSV path or DataFrame with unique `ID` column + DNA columns
- Output: `(out_df, stats_dict)` for most modules; stats-only for `background`, `lenstat`, `verify`, `index`, `pack`
- If `ID` is the DataFrame index on input, it's preserved; otherwise explicit `ID` column returned
- Chain modules (new columns) or use `patch_mode=True` (new rows)

### The Stats Dictionary
```python
stats = {'status': True/False, 'basis': 'solved'|'complete'|'unsolved'|'infeasible',
         'step': 1, 'step_name': '...', 'vars': {...}, 'warns': {...}}
```
- `'solved'`: success; `'complete'`: early-success (e.g., patch mode no-op); `'unsolved'`: exhausted search; `'infeasible'`: impossible constraints
- Invalid arguments raise `RuntimeError('Invalid Argument Input(s).')`; algorithmic failures return `status=False`

### The Oligo Architecture
```
[Primer1] [Barcode1] [Variant] [Barcode2] [Primer2] [Spacer]
```
Variants (user-provided) | Barcodes (Hamming-separated IDs) | Primers (Tm-optimized) | Motifs (anchors/sites) | Spacers (filler)

### Context Columns and Edge Effects
Use `left_context_column` / `right_context_column` to prevent **edge effects**: undesired motifs/repeats emerging at junctions (e.g., barcode `GAATT` + context `...G` = `GAATTC` EcoRI site). Context columns check for excluded motifs and repeats spanning insertion boundaries.

### Key Parameters

| Parameter | Purpose | Notes |
|-----------|---------|-------|
| `oligo_length_limit` | Max oligo length (bp) | Synthesis constraint; use `lenstat` to check remaining space |
| `maximum_repeat_length` | Max shared repeat with context | Range 6-12; lower = harder but safer |
| `excluded_motifs` | Sequences to avoid | List, CSV (Exmotif col), DataFrame, or FASTA |
| `random_seed` | Reproducibility | For stochastic modules: `barcode`, `primer`, `motif`, `spacer`, `split`, `pad` |
| `verbose` | Console output | Default `True`; set `False` for batch |
| `core_count` / `memory_limit` | Analysis tuning | `0`/`0.0` = auto-detect |

---

## Design Philosophy

**Parameter naming**: `*_column` = output destination; `*_context_column` = input source; `*_type` = algorithm; `*_constraint` = IUPAC spec; `minimum_*`/`maximum_*` = bounds

**Must be paired**: `cross_barcode_columns` + `minimum_cross_distance`; `paired_primer_column` needs existing primer; `associate_column` needs anchors

**Errors**: Invalid args → `RuntimeError`; algorithmic failures → `stats['status']=False`

---

## Working Modes

| Mode | Key points |
|------|------------|
| **Notebook** | Use Python API; `verbose=True` + `lenstat`/`verify` mid-pipeline; `help(op.module)` for docs |
| **Script** | Invalid args → exception; failures → check `stats['status']`; use `output_file` for CSVs |
| **Docker** | See [docker-notes.md](https://github.com/ayaanhossain/oligopool/blob/master/docs/docker-notes.md) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/docker-notes.md)] |

---

## Design Mode - Module Reference

For complete parameter documentation, see [api.md](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)].

### barcode

**Purpose**: Generate Hamming-distance separated barcodes for unique variant identification.

**Design order**: After primers/motifs, before spacers.

**API**: See [`barcode`](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#barcode) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)] for parameters.

**Tips**:
- `barcode_type='terminus'` (or `0`): Generates barcodes by extending from fixed termini; faster, good for large libraries
- `barcode_type='spectrum'` (or `1`): Explores full sequence space for maximum diversity; slower but finds more solutions in constrained spaces
- For multi-barcode designs, use `cross_barcode_columns` to ensure BC2 is separated from BC1
- Cross-set separation is strict: provide both `cross_barcode_columns` and `minimum_cross_distance`;
  all cross-set barcode values must already be A/T/G/C strings of length `barcode_length`.

---

### primer
Tm-optimized primers. Design early; for paired primers, design inner first then outer with `paired_primer_column`.
- **Tips**: `'SS'+'N'*18` for 5' GC clamp; `'N'*18+'WW'` for 3' AT clamp; use `background()` for off-target screening
- **`oligo_sets`**: For multiplexed libraries—set-specific primers with cross-compatibility screening; in patch mode existing per-set primers are reused

---

### motif
Sequence motifs or constant anchors. Design before barcodes if using for indexing.
- `motif_type='constant'` (1): ONE sequence shared by all variants
- `motif_type='variable'` (0): DIFFERENT sequence per variant
- Examples: `'GAATTC'` (EcoRI), `'NNNGGATCCNNN'` (BamHI + flanking)

---

### spacer
Neutral filler to reach target length. Design last.
- `spacer_length`: `None` (auto-fill) | `int` (fixed) | `list` (per-variant) | DataFrame (`ID`+`Length`)

---

### background
K-mer database for primer off-target screening. Run before `primer` design.
- Input: list of DNA strings, CSV/DataFrame with `Sequence` column, or FASTA

---

### merge
Collapse contiguous columns into one. Source columns removed from output.
- Use to simplify structure or prepare for `revcomp`

---

### revcomp
Reverse complement column range AND reverse column order. Use for strand orientation switches; even-numbered `split` fragments are auto-revcomped.

---

### lenstat
Check length stats and remaining space. Use frequently mid-pipeline.

---

### verify
QC before synthesis. Concatenates sequence columns left-to-right (skips metadata, strips `'-'` gaps). If `CompleteOligo` exists, uses it directly.

---

### final
Concatenate into synthesis-ready oligos. Output: `CompleteOligo` + `OligoLength` columns.

---

## Assembly Mode - Module Reference

For oligos exceeding synthesis limits (~200 bp). See [api.md](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)].

### split
Break long oligos into overlapping fragments. Each `SplitN` = separate pool to order.
- `separate_outputs=True`: returns `[df_Split1, df_Split2, ...]` (recommended for `pad` workflow)
- Even-numbered splits are auto-revcomped for PCR assembly
- Raw output NOT synthesis-ready—use `pad` next

---

### pad
Add primers + Type IIS sites. Run **once per split column**:
```python
split_dfs, _ = op.split(..., separate_outputs=True)
for i, split_df in enumerate(split_dfs, start=1):
    pad_df, _ = op.pad(split_df, split_column=f'Split{i}', typeIIS_system='BsaI', ...)
    op.final(pad_df, output_file=f'synthesis_Split{i}')
```
- **Exclude Type IIS motif** (e.g., `GGTCTC` for BsaI) from upstream elements via `excluded_motifs`
- Post-synthesis: PCR → Type IIS digest → (mung bean nuclease if sticky) → assemble via overlaps
- **34 supported enzymes**: Common choices: `BsaI`, `BsmBI` (Golden Gate), `SapI`, `MlyI` (blunt)

---

## Degenerate Mode - Module Reference

Compress similar-sequence libraries into IUPAC-degenerate oligos. See [api.md](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)].

### compress
Returns `(mapping_df, synthesis_df, stats_dict)`. Input must be concrete A/T/G/C only.
- Similar sequences compress well; diverse sequences return 1:1 mapping
- **When to recommend**: low mutational diversity libraries, selection assays, "NNK"/"degenerate oligos" mentions, ML-generated variants

---

### expand
Verification tool for `compress`. Expansion is exponential (N N's = 4^N sequences)—use `expansion_limit` for safety. Use `mapping_df` from `compress` for ID traceability.

---

## Analysis Mode - Module Reference

See [api.md](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)].

### index
Build barcode index. Prefix/suffix anchors must be constant (≥6 bp); use `{barcode|associate}_{prefix|suffix}_gap` if anchors are not adjacent in the read. For `acount`, also specify `associate_data`/`associate_column`.

---

### pack
Preprocess FASTQ reads. Single-end: R1 only. Paired-end merge: both R1/R2. Pre-merged: treat as single-end.
- `pack_size`: chunks of N million reads (0.1-5.0; default 3.0). Lower = less memory, more files.

---

### acount
Association counting (barcode-variant coupling). Output: `<index>.ID`, `BarcodeCounts`, `AssociationCounts`.
- Use for synthesis validation, library QC, post-assay verification
- `barcode_errors`/`associate_errors`: `-1` auto-infers; set `0`/`1` for stricter matching

---

### xcount
Barcode-only counting (single or combinatorial). Use for BC1×BC2 combinations, cleaved/uncleaved quantification.
- `barcode_errors`: `-1` auto-infers; set explicitly for stricter matching

---

### Callback Functions (Python API only)
For `acount`/`xcount`: `def callback(r1, r2, ID, count, coreid) -> bool` (r2 is `None` for merged reads). Use for custom filtering, cleavage site extraction, or metrics via `multiprocessing.Manager().dict()`.

---

## Advanced Modules

See [api.md#advanced](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#advanced) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)].

### vectorDB
LevelDB-based k-mer storage for background screening. When reopening, `maximum_repeat_length` loaded from instance.

### Scry
1-NN barcode classifier (internal to `acount`/`xcount`). For custom counting pipelines.

---

## Common Workflows

### Basic Library Design
```python
df = pd.DataFrame({'ID': ['v1', 'v2'], 'Variant': ['ATGC...', 'GCTA...']})
op.background(input_data='plasmid.fasta', maximum_repeat_length=12, output_directory='bg')
df, _ = op.primer(input_data=df, oligo_length_limit=200, primer_sequence_constraint='SS'+'N'*18,
    primer_type='forward', primer_column='P1', left_context_column='Variant', background_directory='bg', ...)
df, _ = op.primer(..., primer_type='reverse', primer_column='P2', paired_primer_column='P1', ...)
df, _ = op.barcode(input_data=df, barcode_length=12, minimum_hamming_distance=3, barcode_column='BC1', ...)
op.lenstat(input_data=df, oligo_length_limit=200)
df, _ = op.spacer(input_data=df, spacer_column='Spacer', spacer_length=None, ...)
op.verify(input_data=df, oligo_length_limit=200)
op.final(input_data=df, output_file='library')
```

---

### Multi-Barcode Design (BC1 × BC2)
```python
df, _ = op.barcode(..., barcode_column='BC1', ...)
df, _ = op.barcode(..., barcode_column='BC2', cross_barcode_columns=['BC1'], minimum_cross_distance=3, ...)
```

---

### Long Oligo Assembly (split → pad → final)
```python
split_dfs, _ = op.split(input_data='long.csv', split_length_limit=150, separate_outputs=True, ...)
for i, split_df in enumerate(split_dfs, start=1):
    pad_df, _ = op.pad(input_data=split_df, split_column=f'Split{i}', typeIIS_system='BsaI', ...)
    op.final(input_data=pad_df, output_file=f'synthesis_Split{i}')
```

**Key points**: `separate_outputs=True` → list of DataFrames. Each SplitN = separate pool. Run `pad` per fragment. Even splits are revcomped. Exclude Type IIS motif from upstream elements.

---

### Degenerate Library Compression
```python
mapping_df, synthesis_df, _ = op.compress(input_data=df, mapping_file='map', synthesis_file='syn')
expanded_df, _ = op.expand(input_data=synthesis_df, sequence_column='DegenerateSeq')  # Verify
```

---

### Analysis Pipeline
```python
# Index (requires constant prefix/suffix anchors ≥6 bp adjacent to indexed column)
op.index(barcode_data='lib.csv', barcode_column='BC1', barcode_prefix_column='Primer1',
    associate_data='lib.csv', associate_column='Variant', associate_suffix_column='Primer2', index_file='idx')
# Pack reads
op.pack(r1_fastq_file='R1.fq.gz', r2_fastq_file='R2.fq.gz', r1_read_type='forward',
    r2_read_type='reverse', pack_type='merge', pack_file='reads')
# Count
op.acount(index_file='idx', pack_file='reads', count_file='counts')  # Association
op.xcount(index_files=['bc1_idx', 'bc2_idx'], pack_file='reads', count_file='combo')  # Combinatorial
```

---

### Saturation Mutagenesis Library
Generate all mutations yourself (loop over positions/nucleotides), then standard Design Mode workflow. For proteins, use codon variants.

---

### Extending an Existing Pool (Patch Mode)
```python
df = pd.concat([existing_df, new_variants_with_None_barcodes])
df, _ = op.barcode(input_data=df, barcode_column='BC1', patch_mode=True, ...)  # Only fills missing
```

---

### CRISPR Guide Library
Standard workflow with `excluded_motifs=['TTTT']` to avoid Pol III terminator.

---

### Ribozyme Cleavage Quantification
Dual barcodes (cross-separated) flanking ribozyme. Combinatorial count: matched pairs = uncleaved; unmatched = cleaved.

---

## Troubleshooting Guide

| Problem | Solutions |
|---------|-----------|
| **Design failed** | Check `stats['basis']`; increase `barcode_length`; decrease `minimum_hamming_distance`; relax `excluded_motifs`; increase `maximum_repeat_length`; widen Tm range; try `barcode_type='terminus'` |
| **Oligo too long** | Use `lenstat`; reduce element lengths; use `split`+`pad` |
| **No barcodes mapping** | Verify anchors match exactly; check gap settings; use `mapping_type='sensitive'`; verify read orientation |
| **Cross-contamination** | Increase `minimum_hamming_distance`; add `cross_barcode_columns`; check index hopping |
| **Memory issues** | Reduce `pack_size`; limit `core_count`; set `memory_limit` |

---

## CLI Usage

```bash
op                    # List commands (no --help flag)
op barcode            # Show options
op manual barcode     # Detailed docstring
op barcode --input-data in.csv --barcode-length 12 --barcode-column BC1 --output-file out
```

**Notes**: `--output-file` required; `--quiet` suppresses output; `--stats-json`/`--stats-file` export stats; `--patch-mode` for extensions; callbacks Python-only; exit codes: 0/1/404. See [api.md#cli-parameter-mapping](https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md#cli-parameter-mapping) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)].

---

## Best Practices

1. **Design order**: background → primers → motifs → barcodes → spacers → verify → final
2. **Use lenstat** frequently; **save intermediate DataFrames**; **set random_seed**
3. **Run verify before synthesis**; **use patch_mode for extensions**
4. **Use cross_barcode_columns** for multi-barcode designs

---

## Decision Guide

| Question | Answer |
|----------|--------|
| Add unique IDs | `barcode` |
| Amplification sites | `primer` |
| Restriction sites/anchors | `motif` |
| Fill to length | `spacer` |
| Screen against genome | `background` → `primer` with `background_directory` |
| Long oligos | `split` + `pad` |
| Similar sequences | `compress` (selection assays); concrete library for MPRA |
| Count reads | `index` + `pack` + `acount` (pairing) or `xcount` (barcode-only/combinatorial) |
| Design failing | Increase `barcode_length`; decrease `minimum_hamming_distance`; widen Tm; reduce `excluded_motifs` |
| Extend pool | `patch_mode=True` |

**Input requirements**: unique `ID` column; context columns = valid DNA; metadata preserved; `verify` tolerates IUPAC/metadata.

---

## Links

[README](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/README.md) | [User Guide](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/docs.md) | [API](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md) | [Docker](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/docker-notes.md) | [Paper](https://doi.org/10.1021/acssynbio.4c00661)
