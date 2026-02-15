# Oligopool Calculator - AI Agent Guide

This is the machine-optimized guide for helping users design, assemble, compress,
and analyze oligo pool libraries with `oligopool` (Python + CLI).

It intentionally does **not** reprint the full API. Use `docs/api.md` for
parameter-by-parameter reference and `docs/docs.md` for tutorials/examples.

## Agent Links

- AI agent guide [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/agent-skills.md)]
- README [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/README.md)]
- User guide (`docs/docs.md`) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/docs.md)]
- API reference (`docs/api.md`) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)]
- Example notebook [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/examples/OligopoolCalculatorInAction.ipynb)]
- CLI YAML pipeline example (repo) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/examples/cli-yaml-pipeline/README.md)]
- Docker notes [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/docker-notes.md)]

Agent policy: start with this guide, then explore only the minimum additional docs needed for the task:
`docs/api.md` for exact signatures/allowed values,
`docs/docs.md` for workflows,
and then module docstrings (`help(op.<module>)` in Python or `op manual <COMMAND>` in CLI) for runtime truth.

## Surface Area

Five modes (library + CLI):
- Design Mode: `background`, `primer`, `motif`, `barcode`, `spacer`, `merge`, `revcomp`, `join`, `final`.
- Assembly Mode: `split`, `pad` (for synthesis limits / post-synthesis assembly).
- Degenerate Mode: `compress`, `expand` (cost optimization via IUPAC-degenerate oligos, selection-based discovery workflows).
- Analysis Mode: `index`, `pack`, `acount`, `xcount` (counting from sequencing reads).
- QC Mode: `lenstat`, `verify`, `inspect` (validate outputs + inspect non-CSV artifacts).

## Interface Map

- Python entry point: `import oligopool as op`
- Python module docs: `help(op.<module>)` for full runtime function details
- CLI entry points: `op` and `oligopool` are equivalent.
- CLI help model:
  - `op --help` / `op COMMAND --help` (quick reference)
  - `op manual TOPIC` and `op manual topics` (full docs)
  - `op cite` (citation)
  - `op complete` (argcomplete)
- CLI file outputs: prefer basenames (no suffix). The CLI appends the appropriate
  `.oligopool.<module>...` suffix(es). For commands that produce a CSV, an
  output basename (e.g., `--output-file out`) is required in CLI mode.

## Core Data Contract

- `input_data` is either a CSV path or a `pd.DataFrame`.
- A unique `ID` column is required across the library (primary key).
- Sequence columns are DNA strings by default (A/T/G/C). Some workflows may
  introduce IUPAC-degenerate columns (notably Degenerate Mode); check
  `docs/api.md` for each module's accepted alphabet.
- Library/CLI CSV outputs include an explicit `ID` column (not a pandas index
  column).

## Cross-Cutting Parameters

- Context columns: `left_context_column` / `right_context_column` prevent *edge
  effects* (undesired motifs/repeats spanning an insertion boundary).
  When both are set, new design columns are inserted between those context
  columns, which also resolves linear column order for later `merge`.
- `excluded_motifs`: global "do not create these" list (restriction sites, Type
  IIS sites, primers, etc.; many modules accept it). Accepts single source or
  multiple sources (list of sources or `{name: source}` dict). Strict ATGC only.
- `background_directory`: screen designed elements against one or more background
  k-mer DB(s) created by `background()` (single path or list of paths; supported
  by multiple Design/QC modules).
- `*_type` and `*_policy` parameters usually accept either ints or descriptive strings
  (case-insensitive; often fuzzy-matched). See `docs/api.md` for each module's
  accepted values.
- `random_seed`: reproducibility for stochastic modules (see `docs/api.md` for
  which modules accept it).
- `verbose`: enables the step-by-step console report. CLI also supports `--quiet`.

## Return Shapes

- Most design/transform modules: `(out_df, stats_dict)`
- Stats-only modules: `background`, `lenstat`, `inspect`, `index`, `pack`
- Counting modules: `(counts_df, stats_dict)` for `acount` and `xcount`
- `compress`: `(mapping_df, synthesis_df, stats_dict)` (two DataFrames)
- `split`: returns either `(df, stats)` or `([df_Split1, df_Split2, ...], stats)`
  depending on `separate_outputs` (CLI defaults to separate outputs enabled).
- `join`: `(out_df, stats_dict)`; joins two tables on `ID` and inserts new columns
   from `other_data` in resolved order, for recombining parallel design branches.

## Special Contracts

### Patch Mode (`patch_mode`)

Goal: extend an existing library without overwriting already-designed elements.

- `patch_mode=True` fills *only* missing values in the output column.
- Missing is any of: `None`, `NaN`, empty string, or `'-'`.
- Existing non-missing values are kept as-is; they must still be valid for the
  module (at minimum: correct type/length/alphabet) or validation will fail.

### Cross-set barcodes (`barcode`)

Goal: design BC2/BC3 so they stay far from previously designed `barcode` columns.

- `cross_barcode_columns` and `minimum_cross_distance` must be provided together.
- Cross-set separation is global (not per-row): each new barcode must be at
  least `minimum_cross_distance` away from every barcode in the union of the
  cross columns.

### Per-set primers (`primer` + `oligo_set`)

Goal: multiplexing / selective amplification with set-specific primers.

- `oligo_set` (CLI: `--oligo-set`) groups rows by label; one primer is designed
  per set and assigned to all rows in that set.
- Primers are screened for cross-set compatibility (to reduce cross-dimer risk).
- With `paired_primer_column`, Tm matching is applied within each set and the
  paired primer must be constant within each set.
- With `patch_mode=True`, existing per-set primers are reused; only missing
  rows/sets trigger new design.

### Rejoining branch outputs (`join`)

Goal: recombine two independent design branches back into one table (most commonly
after a YAML CLI DAG fan-out).

- `join` requires the exact same `ID` set in `input_data` and `other_data` (row order
  may differ); it never creates or drops rows.
- Column order is preserved from `input_data` (backbone). Only new, non-overlapping
  columns from `other_data` are inserted into that order.
- If insertion placement is ambiguous, `join_policy` resolves it (`left` vs `right`).

### Split outputs (`split`)

Goal: avoid "one giant CSV with Split1..SplitN columns" when you actually need
separate pools per fragment.

- With `separate_outputs` enabled, `split` returns a list of DataFrames and, if
  writing to disk, emits one CSV per split fragment (suffixing Split1, Split2,
  ...).
- Even-numbered split fragments are auto-`revcomp`ed for assembly workflows.
- `pad` is typically run **once per split fragment**.

### Type IIS compatibility (`pad`)

`pad` adds Type IIS sites for downstream assembly. It also checks that the
chosen Type IIS system is compatible with the pool (i.e., you didn't accidentally
embed the recognition site where it breaks the workflow). Upstream design can
also exclude the recognition site via `excluded_motifs`. The sites are meant to
facilitate scarless removal of the installed amplification pads for overlap-based
assembly strategies identified via `split`.

The built-in Type IIS set is intentionally limited to systems with a clearly
defined motif + 3' cut-offset model in `pad`; for non-built-ins, manually compose
padding with `primer` + `motif` + `spacer`.
Think of this set as the "batteries included" Type IIS options for predictable
3' excision and scarless primer-pad removal.

### Counting callbacks (Python only)

`acount`/`xcount` support a `callback` function in Python API mode only.
CLI does not run callbacks.

### Index anchors (`index` / counting)

`index` relies on constant prefix/suffix anchors (typically >=6 bp). If anchors
are not adjacent in the read, use the gap parameters; mismatched anchors/gaps
silently reduce usable reads downstream.

Counting detail: if an anchor appears multiple times in a read, the engine keeps
the best-scoring placement(s); if multiple placements tie and yield different
barcodes, the read is rejected as ambiguous. Design anchors to be unique (e.g.,
`motif(motif_type=1, maximum_repeat_length=...)`) to avoid this.

### Failed reads sampling (`acount` / `xcount`)

`failed_reads_file` writes a small per-category diagnostic CSV (sampling is off
by default). It is not returned as a DataFrame.

### Sequence constraint shorthand (CLI)

`--primer-sequence-constraint` / `--motif-sequence-constraint` accept either an
IUPAC string (`NNNN...`) or shorthand expressions like `"'N'*20"` /
`GCC+N*20+CCG` (quote as one argument when needed).

## CLI Config Files + Pipelines

The CLI supports YAML config files:
- Any command: `op COMMAND --config config.yaml`
- Override config at the CLI: `op barcode --config x.yaml --barcode-length 20`

Pipeline execution:
- `op pipeline --config pipeline.yaml`
- `op pipeline --config pipeline.yaml --dry-run`

Key rules:
- Config keys use Python parameter names (snake_case).
- CLI flags override config values.
- Pipeline configs can express sequential steps and parallel groups; see
  `docs/docs.md` examples for the supported schema.
- Recommended pattern: keep design workflows mostly serial; use parallel DAGs
  primarily for analysis (`index`/`pack`/`acount`/`xcount` branches).
- Pipeline shorthand is supported: downstream `input_data` can reference
  a prior step `output_file` basename; explicit existing paths are preserved.
- If the same basename is produced by multiple steps, that alias is treated as
  ambiguous and pipeline validation raises a config error.
- See `docs/docs.md` Config Files for the exact schema + examples.

## Module Cheat Sheet

Design Mode:
- `background`: build a background k-mer DB for screening designed elements.
- `primer`: design primers (supports paired primers + `oligo_set` multiplexing).
- `motif`: insert motifs or design constant anchors (`motif_type='constant'`).
- `barcode`: design Hamming-separated barcodes (supports cross-set constraints).
- `spacer`: fill length (supports per-ID lengths; can auto-fill).
- `merge`/`revcomp`/`join`: mid-pipeline architecture maneuvers.
- `final`: concatenate into synthesis-ready `CompleteOligo` (+ length).

Assembly Mode:
- `split`: fragment long oligos into overlapping pieces.
- `pad`: add primers + Type IIS sites for post-synthesis assembly (per split).

Degenerate Mode:
- `compress`: reduce similar sequences into IUPAC-degenerate oligos.
- `expand`: sanity-check expansions (be mindful of exponential blow-up).

Analysis Mode:
- `index`: build barcode(+associate) index; requires constant anchors.
- `pack`: preprocess FastQ into packs for counting.
- `acount`: association counting (barcode <-> associate).
- `xcount`: barcode-only counting (single or combinatorial).
  - both `acount` and `xcount` support `failed_reads_file` sampling
    for assessing analysis failures

QC Mode:
- `lenstat`: ruler/checkpoint for oligo lengths mid-pipeline.
- `verify`: detect length, motif emergence, and background conflicts (CLI: "motif").
- `inspect`: inspect non-CSV artifacts (background/index/pack) and summarize their metadata.

## Workflow Templates

Design:
`background → primer → motif → barcode → spacer → lenstat → verify → final`

QC:
`lenstat → verify` (and `inspect` for background/index/pack artifacts as needed)

Assembly (long constructs):
`split (separate outputs) → pad (per fragment) → final (per fragment)`

Degenerate (cost optimization):
`compress → (optional: expand sanity-check) → order synthesis_df`

Analysis:
`index → pack → (acount or xcount)`

Notes:
- `index` files and `pack` files are reusable artifacts; iterate `acount`/`xcount` settings without rebuilding.
- For scale, tune `core_count` and `memory_limit` in analysis modules.

Extension (add new rows later):
use `patch_mode=True` on the element you need to fill (e.g., barcodes/primers).

## Composability Recipes

Dense lookup for common "compose it from primitives" patterns. For worked examples and caveats, see `docs/docs.md#composability`.

Note: `'-'` is a conventional placeholder value; Patch Mode treats `'-'` as missing and fills it.

### Recipe Table

| Ask | Compose |
|-----|---------|
| Embed restriction site inside a barcode | `motif(constant, 'GAATTC')` → `barcode(N1, right_context='site')` → `barcode(N2, left_context='site')` → `merge` |
| Two barcodes that don't cross-talk | `barcode(col='BC1')` → `barcode(col='BC2', cross_barcode_columns=['BC1'])` |
| Tm-matched primer pair | `primer(forward)` → `primer(reverse, paired_primer_column='Fwd')` |
| Extend pool with new variants | append rows → `barcode(patch_mode=True)` → `spacer(patch_mode=True)` |
| Auto-fill to target length | `spacer(spacer_length=None)` - auto-computes per-row fill |
| Compound element (mixed fixed/variable) | chain `motif`/`barcode`/`primer`/`spacer` sub-regions → `merge`/`revcomp` → repeat |
| Flip a region's orientation | `revcomp(left_context_column=..., right_context_column=...)` |
| Long constructs / assembly | `split(separate_outputs=True)` → `pad` (per `SplitN`) → `final` (per fragment) |
| Mid-pipeline length check | `lenstat` → adjust element lengths → rerun |
| Inspect artifacts quickly | `inspect(background/index/pack)` before reuse |
| Cut-site-free library | same `excluded_motifs` (motifs + reverse complements) on every design module + `verify` |
| No host-genome homology | `background(host.fasta)` → `background_directory` on every design module + `verify` |
| Count multi-barcode combos (BC1 x BC2) | separate `index` per barcode → `xcount(index_files=[idx1, idx2])` |
| N-way combinatorial matrix | N independent `index` files → `xcount(index_files=[idx1,...,idxN])` |
| Verify barcode-variant coupling | `index` with associate columns → `acount` |
| Filter reads by custom criteria | `acount`/`xcount` with `callback(r1, r2, ID, count, coreid) → bool` (Python only) |
| Pack once, count many ways | single `pack` → reuse with different indexes/callbacks/counting modes |
| Debug discarded reads | `failed_reads_file` → sample reads per failure category |
| Design-to-analysis anchor bridge | Design: `motif(const, prefix)` → `motif(const, suffix)` → `barcode(left_context=prefix, right_context=suffix)`; Analysis: `index(barcode_prefix_column=..., barcode_suffix_column=...)` |
| Cost-efficient saturation mutagenesis | generate all substitutions → `compress` → order `synthesis_df` (6-20x fewer oligos) |
| Compression for analysis mapping | `compress` once → reuse `mapping_df` to map sequenced survivors/readouts back to variant IDs |
| Selection-based discovery (no barcodes) | `compress` → synthesize → select → sequence → map back via `mapping_df` |
| Verify compression before synthesis | `compress` → `expand` → confirm expanded set == original set |
| Promoter MPRA library | `background` → `primer(fwd)` → `motif(const, BC_prefix_anchor)` → `motif(const, BC_suffix_anchor)` → `barcode` → `primer(rev, paired)` → `spacer(auto)` → `verify` → `final` |
| CRISPR guide library | guides as variant column → `background` → `primer` → `barcode` → `spacer` → `verify` → `final` - with `excluded_motifs` for cloning sites (BsmBI/BbsI) |
| Saturation mutagenesis (barcoded) | generate all substitutions → standard design pipeline for individual tracking |
| Saturation mutagenesis (degenerate) | generate all substitutions → `compress` → selection-based discovery at fraction of cost |

### Context Threading Rules

1. `excluded_motifs` is strict ATGC only - for **non-palindromic** motifs, list both orientations (motif + reverse complement), e.g. BsmBI: `CGTCTC` + `GAGACG`.
2. Pass `left_context_column`/`right_context_column` at every design step - prevents motif/repeat emergence at element boundaries.
3. If both context columns are specified, insertion order is resolved automatically (new column lands between contexts), so merge ranges usually line up without manual column reordering.
4. Thread `excluded_motifs` and `background_directory` identically through every design module AND `verify` - constraint composition is pipeline-wide, not per-module.
5. Treat motif/background inputs as named layers (core vs optional); add optional layers incrementally, and if infeasible, remove optional layers one at a time while keeping core layers fixed.
6. `verify` checks the full concatenated oligo + junctions; run it as final QC to catch anything that slipped through.

### Pointers

- Examples + explanations: `docs/docs.md#composability`
- Exact parameter semantics: `docs/api.md` (module sections for `barcode`, `primer`, `motif`, `spacer`, `index`, `pack`, `acount`, `xcount`, `compress`, `expand`, `verify`)

## Troubleshooting

- Design fails: read `stats['basis']`; relax constraints (shorter repeats, wider
  Tm), increase design space (longer barcode/primer), or change algorithm type.
- Too long: use `lenstat`, then reduce element lengths or use `split` + `pad`.
- Counting fails: verify anchors/gaps/orientation, then increase mapping sensitivity;
  use `failed_reads_file` to inspect why reads are being discarded.

## CLI Exit Codes

- `0`: success
- `1`: runtime error / design failure / invalid args
- `404`: unknown command
