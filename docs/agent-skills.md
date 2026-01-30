# Oligopool Calculator - AI Agent Guide (token-efficient)

This is the machine-optimized guide for helping users design, assemble, compress,
and analyze oligopool libraries with `oligopool` (Python + CLI).

It intentionally does **not** reprint the full API. Use `docs/api.md` for
parameter-by-parameter reference and `docs/docs.md` for tutorials/examples.

## Agent Links (raw)

- README [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/README.md)]
- User guide (`docs/docs.md`) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/docs.md)]
- API reference (`docs/api.md`) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/api.md)]
- Example notebook [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/examples/OligopoolCalculatorInAction.ipynb)]
- CLI pipeline example (repo) [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/examples/cli-pipeline/README.md)]
- Docker notes [[agent-link](https://raw.githubusercontent.com/ayaanhossain/oligopool/refs/heads/master/docs/docker-notes.md)]

Agent policy: don't fetch everything by default. Fetch the minimum needed:
`docs/api.md` for exact signatures/allowed values, `docs/docs.md` for workflows,
and then module docstrings (`help(op.<module>)`) for runtime truth.

## Surface Area (what users can do)

Four modes (library + CLI):
- Design Mode: `background`, `primer`, `motif`, `barcode`, `spacer`, plus `lenstat`,
  `verify`, `final`, `merge`, `revcomp`.
- Assembly Mode: `split`, `pad` (for synthesis limits / post-synthesis assembly).
- Degenerate Mode: `compress`, `expand` (cost optimization via IUPAC-degenerate synthesis).
- Analysis Mode: `index`, `pack`, `acount`, `xcount` (counting from sequencing reads).

## Interface Map (Python + CLI)

- Python entry point: `import oligopool as op`
- CLI entry points: `op` and `oligopool` are equivalent.
- CLI help model: there is no `--help` flag; use:
  - `op` (main menu)
  - `op COMMAND` (command options)
  - `op manual TOPIC` and `op manual topics`
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

## Cross-Cutting Parameters (know these once)

- Context columns: `left_context_column` / `right_context_column` prevent *edge
  effects* (undesired motifs/repeats spanning an insertion boundary).
- `excluded_motifs`: global "do not create these" list (restriction sites, Type
  IIS sites, primers, etc.; many modules accept it).
- `*_type` parameters usually accept either ints or descriptive strings
  (case-insensitive; often fuzzy-matched). See `docs/api.md` for each module's
  accepted values.
- `random_seed`: reproducibility for stochastic modules (see `docs/api.md` for
  which modules accept it).
- `verbose`: enables the step-by-step console report. CLI also supports `--quiet`.

## Return Shapes (high-signal)

- Most design/transform modules: `(out_df, stats_dict)`
- Stats-only modules: `background`, `lenstat`, `verify`, `index`, `pack`
- Counting modules: `(counts_df, stats_dict)` for `acount` and `xcount`
- `compress`: `(mapping_df, synthesis_df, stats_dict)` (two DataFrames)
- `split`: returns either `(df, stats)` or `([df_Split1, df_Split2, ...], stats)`
  depending on `separate_outputs` (CLI defaults to separate outputs enabled).

## Special Contracts (the gotchas)

### Patch Mode (`patch_mode`)

Goal: extend an existing library without overwriting already-designed elements.

- `patch_mode=True` fills *only* missing values in the output column.
- Missing is any of: `None`, `NaN`, empty string, or `'-'`.
- Existing non-missing values are kept as-is; they must still be valid for the
  module (at minimum: correct type/length/alphabet) or validation will fail.

### Cross-set barcodes (`barcode`)

Goal: design BC2/BC3 so they stay far from previously designed barcode columns.

- `cross_barcode_columns` and `minimum_cross_distance` must be provided together.
- Cross-set separation is global (not per-row): each new barcode must be at
  least `minimum_cross_distance` away from every barcode in the union of the
  cross columns.

### Per-set primers (`primer` + `oligo_sets`)

Goal: multiplexing / selective amplification with set-specific primers.

- `oligo_sets` (CLI: `--oligo-sets`) groups rows by label; one primer is designed
  per set and assigned to all rows in that set.
- Primers are screened for cross-set compatibility (to reduce cross-dimer risk).
- With `paired_primer_column`, Tm matching is applied within each set and the
  paired primer must be constant within each set.
- With `patch_mode=True`, existing per-set primers are reused; only missing
  rows/sets trigger new design.

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
also exclude the recognition site via `excluded_motifs`.

### Counting callbacks (Python only)

`acount`/`xcount` support a `callback` function in Python API mode only.
CLI does not run callbacks.

### Index anchors (`index` / counting)

`index` relies on constant prefix/suffix anchors (typically >=6 bp). If anchors
are not adjacent in the read, use the gap parameters; mismatched anchors/gaps
silently reduce usable reads downstream.

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
- See `docs/docs.md` Config Files for the exact schema + examples.

## Module Cheat Sheet (one-liners)

Design Mode:
- `background`: build a background k-mer DB for primer screening.
- `primer`: design primers (supports paired primers + `oligo_sets` multiplexing).
- `motif`: insert motifs or design constant anchors (`motif_type='constant'`).
- `barcode`: design Hamming-separated barcodes (supports cross-set constraints).
- `spacer`: fill length (supports per-ID lengths; can auto-fill).
- `lenstat`: ruler/telemetry for oligo lengths mid-pipeline.
- `verify`: QC: columns, lengths, excluded motifs, and edge effects.
- `final`: concatenate into synthesis-ready `CompleteOligo` (+ length).
- `merge`/`revcomp`: mid-pipeline architecture maneuvers.

Assembly Mode:
- `split`: fragment long oligos into overlapping pieces.
- `pad`: add primers + Type IIS sites for post-synthesis assembly (per split).

Degenerate Mode:
- `compress`: reduce similar sequences into IUPAC-degenerate oligos.
- `expand`: sanity-check expansions (be mindful of exponential blow-up).

Analysis Mode:
- `index`: build barcode(+associate) index; requires constant anchors.
- `pack`: preprocess FastQ into packs for counting.
- `acount`: association counting (barcode <-> associate); supports `failed_reads_file` sampling.
- `xcount`: barcode-only counting (single or combinatorial); supports `failed_reads_file` sampling.

## Workflow Templates (minimal)

Design:
`background -> primer -> motif -> barcode -> spacer -> verify -> final`

Assembly (long constructs):
`split (separate outputs) -> pad (per fragment) -> final (per fragment)`

Degenerate (cost optimization):
`compress -> (optional: expand sanity-check) -> order synthesis_df`

Analysis:
`index -> pack -> (acount or xcount)`

Extension (add new rows later):
use `patch_mode=True` on the element you need to fill (e.g., barcodes/primers).

## Troubleshooting (fast heuristics)

- Design fails: read `stats['basis']`; relax constraints (shorter repeats, wider
  Tm), increase design space (longer barcode/primer), or change algorithm type.
- Too long: use `lenstat`, then reduce element lengths or use `split`+`pad`.
- Counting fails: verify anchors/gaps/orientation, then increase mapping sensitivity.

## CLI Exit Codes (useful for automation)

- `0`: success
- `1`: runtime error / design failure / invalid args
- `404`: unknown command
