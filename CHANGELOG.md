2026.02.09
1. Docs: added Composability section to `docs/docs.md` with worked Python examples for design composition, constraint composition, analysis composition, degenerate recipes, application templates (MPRA, CRISPR, ribozyme dual-barcode readout), and cheat sheet.
2. Docs: added Composability Recipes section to `docs/agent-skills.md` with dense recipe lookup table, context threading rules, and analysis workflow notes.
3. Docs: added concrete Patch Mode example to Core Concepts in `docs/docs.md`.
4. Docs: added back-to-TOC navigation links to all subsections in `docs/docs.md` and `docs/api.md`.
5. Docs: added stats-dict debugging hint, analysis tips (artifact reuse, scale knobs), and clarified ribozyme dual-barcode narrative in notebook.
6. Bugfix: `inspect` now correctly converts `pack_size` from raw reads-per-pack to millions for display.

2026.02.08
1. Feature: `excluded_motifs` now accepts multiple motif sets (list of sources or `{name: source}` dict) in `barcode`, `primer`, `motif`, `spacer`, and `verify`; all sources are merged for screening, with per-set attribution in stats.
2. CLI: `--excluded-motifs` now accepts multiple values (space-separated) for barcode/primer/motif/spacer/verify.
3. Validation/Core: added `validation_parsing.get_parsed_exmotifs()` polymorphic wrapper; added `utils.is_strict_DNA()` and `utils.stamp_exmotif_stats()`.
4. Strict ATGC enforcement: excluded motifs must be concrete ATGC strings (no IUPAC codes, no dashes).
5. Utility: added `inspect` (stats-only) to summarize non-CSV artifacts (background DBs, `.oligopool.index`, `.oligopool.pack`) and added `op inspect` CLI command.
6. Inspect: reports a `verdict` (`Valid`/`Corrupted`/`Invalid`) plus artifact metadata (background K/LEN, index anchors/gaps, pack type/size/read counts).

2026.02.04
1. Feature: `background_directory` now accepts multiple background DBs (list of `.oligopool.background` paths or `vectorDB` instances) in `barcode`, `primer`, `motif`, `spacer`, and `verify`; screening is applied across ALL specified DBs.
2. CLI: `--background-directory` now accepts multiple values (space- or comma-separated) for barcode/primer/motif/spacer/verify.
3. Validation/Core: added `validation_parsing.get_parsed_backgrounds()` and updated background feasibility logic in `core_background.py` and `core_primer.py` to accept lists.
4. Tests: `smoketest/smoke_multiple_backgrounds.py` covers single + multi-background workflows.
5. Breaking: redesigned `verify` — now requires `oligo_length_limit`, returns `(DataFrame, stats)`, and can write `.oligopool.verify.csv` with per-row conflict flags (`HasLengthConflict`, `HasExmotifConflict`, `HasBackgroundConflict`, `HasAnyConflicts`) plus `*ConflictDetails` dict columns (serialized as JSON strings in CSV).
6. Validation: added `validation_parsing.get_parsed_verify_indata_info()` for verify-specific input validation (requires ID column + at least one DNA column); error shown in Required Arguments section if no DNA columns detected.
7. CLI: stats JSON output now works for tuple-return modules like `compress` (stats are extracted from the dict element rather than assuming `(df, stats)`), including in pipeline mode.
8. Style: `verify` now follows design module conventions — type hints without spaces (`str|pd.DataFrame`), `id_from_index` handling for caller intent preservation.
9. Docs: clarified that `json.loads` applies to all `*Details` columns in `verify` CSV output (currently the `*ConflictDetails` columns).
10. CLI: `op verify` now requires `--output-file` (Python API still allows `output_file=None`).

2026.02.01
1. Background integration: added `background_directory` parameter to `barcode`, `motif`, `spacer`, and `verify` modules for off-target k-mer screening.
2. Background integration: added `is_background_feasible()` to `core_background.py` with context-aware junction checking (checks k-mers spanning element boundaries).
3. Background integration: added `background_fail` stat counter to barcode, motif, spacer, and primer modules (separated from `repeat_fail`).
4. Background integration: added `--background-directory` CLI flag to barcode, motif, spacer, and verify commands.
5. Background integration: `verify` module now scans assembled oligos for background k-mer violations and reports failing IDs.
6. Multi-anchor handling: added `get_anchored_read_candidates()` and `get_best_barcode_from_candidates()` to `core_count.py` for reads with multiple anchor occurrences.
7. Multi-anchor handling: added `barcode_ambiguous` failure category to `acount` and `xcount` for reads where different barcodes tie on Scry confidence.
8. Multi-anchor handling: added per-anchor segment splitting via `_segments_from_spans()` to evaluate each anchor occurrence independently.
9. Multi-anchor handling: added fast-path optimization using `read.count()` for exact anchor matches before edlib alignment.
10. Bugfix: barcode background screening now checks each candidate junction against background DB (was only checking first context).
11. Bugfix: motif/spacer modules now initialize `prefixdict`/`suffixdict` when using background-only paths (fixes crash when no excluded_motifs specified).
12. Docs: added background integration documentation to `api.md`, `docs.md`, and `agent-skills.md`.
13. Docs: added multi-anchor handling and `barcode_ambiguous` documentation to `api.md` and `agent-skills.md`.
14. Tests: added `smoketest/smoke_background_integration.py` for background k-mer screening across all design modules.
15. Tests: added `smoketest/smoke_multi_anchor.py` for multi-anchor read handling in acount/xcount.
16. Docs: standardized terminology—"oligo pool" (two words) for the concept, "oligopool" for the tool/package name—across README, docstrings, docs, and notebook.
17. Cleanup: removed unused functions (`get_hdist`, `get_edist`, `get_categorical_validity`, `merge_config_with_args`, `_nb_popcount4`, `get_validated_concrete_sequence`).
18. Docs: expanded mode introductions in `docs.md`—added "When to use" guidance and workflow steps for Design Mode, Assembly Mode, and Analysis Mode (matching Degenerate Mode style).
19. Breaking: `compress` output suffixes are now distinct—`mapping_file` writes `.oligopool.compress.mapping.csv`, and `synthesis_file` writes `.oligopool.compress.synthesis.csv` (previously both used `.oligopool.compress.csv`).
20. Docs: updated `compress` CLI help text and API docs to reflect the new mapping/synthesis output suffixes.
21. Docs: tightened module docstrings and mode workflow text for accuracy (e.g., `pack` does not perform adapter trimming; `index` builds one index per barcode column).
22. Notebook: updated `OligopoolCalculatorInAction.ipynb` to reflect Analysis Mode terminology (replaced “Counting Mode”), documented new `compress` output suffixes, and standardized string quoting in code cells.
23. Expand: `mapping_file` now auto-suffixes to `.oligopool.compress.mapping.csv` (matching index/pack input behavior); mapping parsing moved into `validation_parsing.py`.
24. Consistency: standardized `random_seed` validation + display across all RNG-using design modules via `validation_parsing.get_parsed_random_seed_info()`.

2026.01.31
1. CLI: added `--config` flag for YAML config file support on all commands.
2. CLI: added `op pipeline` subcommand for multi-step workflow execution from YAML configs.
3. CLI: pipeline supports both sequential steps and parallel DAG execution with `after:` dependencies.
4. CLI: added `--dry-run` flag for pipeline config validation without execution.
5. Core: added `oligopool/base/config_loader.py` with YAML loading, validation, and topological sort for DAG pipelines.
6. Docs: added "Config Files" section to `docs/docs.md` with single command, pipeline, and parallel examples.
7. Docs: rewrote `agent-skills.md` as compact machine-first "contracts + patterns" format.
8. Docs: added `--config` and `op pipeline` to CLI appendix in `api.md`.
9. Examples: added `examples/cli-pipeline/` with sequential and parallel YAML pipeline examples.
10. Examples: added README.md files to all example directories.
11. Pack format: migrated from delimiter-based concatenation (`r1 + '-----' + r2`) to tuple-based storage `(r1, r2)` for cleaner semantics and type safety.
12. Core: `get_concatenated_reads()` now returns `(r1, r2)` tuple instead of delimiter-joined string.
13. Core: `get_merged_reads()` now returns `(merged, None)` tuple on success for consistent format.
14. Analysis: `acount_engine()` and `xcount_engine()` extract R1 from tuple (`read_tuple[0]`) instead of stripping delimiter.
15. Analysis: `stream_packed_reads()` updated to handle tuple-based pack format for callback validation.
16. Cleanup: removed `CONCAT_DELIMITER`, `CONCAT_DELIMITER_BYTES` constants and `strip_concat_delimiter()` function from utils.py.
17. Tests: updated `smoke_concat_delimiter.py` to verify tuple-based storage format.
18. Merging: added `MAX_MISMATCH_DENSITY = 0.10` threshold for adaptive alignment k-values.
19. Merging: added numba-compiled `get_exact_innie_overlap()` and `get_exact_outie_overlap()` for fast-path exact overlap detection.
20. Merging: added numba-compiled `get_innie_consensus_fast()` and `get_outie_consensus_fast()` using byte arrays for ~10x faster consensus building.
21. Merging: `get_innie_merged()` and `get_outie_merged()` now use fast-path exact overlap before falling back to alignment.
22. Merging: adaptive seed length (`max(10, min(20, len(r) // 10))`) and density-based k-values for better accuracy on variable-length reads.
23. Merging: full bytes support throughout merge pipeline for reduced memory allocations.
24. Utils: added `complement_bytes` translation table and `get_revcomp_bytes()` for fast byte-level reverse complement.
25. Utils: added `N_BYTE` constant for efficient byte-level N filtering.
26. Utils: added `readbytes` parameter to `stream_fastq_engine()` for returning reads as bytes, enabling faster downstream processing.
27. Counting: dual-read processing - anchor search now tries R1 first, then R2 if anchor not found.
28. Counting: associate matching now searches both R1 and R2 for associate sequences.
29. Counting: **BREAKING** callback signature changed from `callback(read=..., ID=..., count=..., coreid=...)` to `callback(r1=..., r2=..., ID=..., count=..., coreid=...)` where r2 is None for merged reads.
30. Counting: callback validation updated to test with r1/r2 signature and mixed merged/concat read tuples.
31. Pack storage: `get_concatenated_reads()` and `get_merged_reads()` now decode bytes to strings, ensuring consistent string storage in pack files.
32. Degenerate: `expand` now accepts optional `mapping_file` parameter to restore original variant IDs from `compress` output.
33. Index: added `associate_prefix_gap` / `associate_suffix_gap` to support non-adjacent associate anchors in `index()` for `acount`.
34. Counting: added `failed_reads_file` and `failed_reads_sample_size` parameters to `acount()` and `xcount()` for diagnostic sampling of failed reads by category (phix, low_complexity, anchor_missing, barcode_absent, etc.).

2026.01.27
1. Docs: unified Degenerate Mode messaging across all documentation - reframed as selection/sequence-identity workflow for low mutational diversity libraries.
2. Docs: added Four Modes overview (Design, Assembly, Degenerate, Analysis) to README and agent-skills.md with workflow diagrams and integration points.
3. Docs: optimized agent-skills.md for token efficiency - reduced from ~1000 to ~430 lines while preserving all essential information.
4. Docstrings: aligned `compress` and `expand` docstrings with `barcode`/`primer` style (condensed parameters, added Notes sections).
5. CLI: differentiated `split` and `compress` descriptions from their mode headers in help output.

2026.01.26
1. Degenerate Mode: added `compress` module to compress ML-generated variant libraries into IUPAC-degenerate oligos for cost-efficient synthesis.
2. Degenerate Mode: added `expand` module to expand IUPAC-degenerate oligos into concrete sequences for verification.
3. Core: added `base/core_degenerate.py` with Monte Carlo rollout compression algorithm using bitmask-encoded IUPAC codes.
4. CLI: added `compress` and `expand` subcommands with full parameter support.
5. Docs: added Degenerate Mode sections to `api.md`, `docs.md`, and `agent-skills.md` with usage examples.
6. Notebook: enhanced Degenerate Mode section in `OligopoolCalculatorInAction.ipynb` with scientific motivation, IUPAC code table, realistic 500-variant example, and post-selection workflow.
7. Notebook: added consistent problem-solution framing, comparison tables, and key takeaways across all sections (Design Mode, Assembly Mode, Degenerate Mode, Analysis Mode).
8. README: added Degenerate Mode to Features section and updated `help(op)` and CLI command list examples.

2026.01.25
1. Docs: added Example Notebook links to all documentation files (agent-skills.md, docs.md, api.md, docker-notes.md).
2. Docs: added agent instruction to fetch related documentation in agent-skills.md.
3. Docs: refactored `api.md` parameter tables to bulleted lists for readability.
4. Docs: added Notes sections to all modules in `api.md` with key information from docstrings.

2026.01.24
1. Type parameters: all type parameters now accept both integers and intelligent string specifications with fuzzy matching.
2. Design mode: `barcode_type`, `primer_type`, and `motif_type` accept string aliases (e.g., `'terminus'`, `'forward'`, `'variable'`).
3. Analysis mode: `r1_read_type`, `r2_read_type`, `pack_type`, and `mapping_type` accept string aliases (e.g., `'forward'`, `'concatenate'`, `'fast'`).
4. Validation: added `resolve_type_parameter()`, `get_typed_categorical_validity()`, and `get_optional_typed_categorical_validity()` to `validation_parsing.py` with fuzzy matching via `difflib.get_close_matches(cutoff=0.6)`.
5. CLI: type parameters now accept both integer and string values via `_parse_type_param()` helper.
6. Docs: updated `docs/api.md` to document string aliases for all type parameters.
7. Notebook: updated `OligopoolCalculatorInAction.ipynb` to use string type aliases (e.g., `primer_type='forward'`, `barcode_type='spectrum'`, `mapping_type='sensitive'`).
8. Repo hygiene: notebooks are stored output-free (no execution outputs/counts) to keep diffs clean.

2026.01.23
1. Pad module: added Type IIS compatibility check as separate pipeline step (Step 2) that validates the chosen enzyme's recognition site is absent from split fragments.
2. Pad module: refactored compatibility check into `core_primer.get_typeIIS_compatibility()` alongside existing `get_parsed_typeIIS_constraint()`.
3. Docs: documented Type IIS site exclusion requirement in `docs/docs.md`, `docs/api.md`, and `docs/agent-skills.md` (exclude enzyme motif from upstream design elements).
4. Notebook: updated `OligopoolCalculatorInAction.ipynb` to exclude BsaI motif (`GGTCTC`/`GAGACC`) from all design elements when using BsaI for padding.

2026.01.22
1. Docs: extracted API documentation into `docs/api.md` as single source of truth for parameter references.
2. Docs: refactored `docs/docs.md` to link to `api.md` for parameter details (removed inline `<details>` blocks).
3. Docs: refactored `docs/agent-skills.md` to link to `api.md` (removed inline function signatures).
4. Docs: added strategic cross-references between all documentation files (docs.md, api.md, agent-skills.md, docker-notes.md).
5. Docs: clarified `merge` and `revcomp` as mid-pipeline maneuvers in agent-skills.md.
6. Docs: added Docker mode section to agent-skills.md.
7. Docs: documented `verify` column concatenation behavior (sequence columns only, left-to-right order, gap stripping).
8. README: added API Reference link to header navigation and Documentation section.

2026.01.18
1. Development in progress.
2. Packaging: migrated to `pyproject.toml` (PEP 621 / PEP 517).
3. Barcode module: added cross-set Hamming constraints (`cross_barcode_columns`, `minimum_cross_distance`) and CLI flags.
4. Barcode module: refactored cross-set validation + precomputation into shared helpers (`validation_parsing.py`, `utils.py`).
5. CLI: allow multi-value `--cross-barcode-columns` and `--index-files` (space- or comma-separated).
6. Docs: tightened module docstring descriptions and added high-signal Notes bullets (e.g., primer chaining, motif anchors).
7. Notebook: updated `OligopoolCalculatorInAction.ipynb` to explain cross-barcode sets, introduce `lenstat`/`verify` usage, mention `merge`, and fix minor text issues.
8. CLI: `op manual topic` is now an alias of `op manual topics`.
9. CLI: `op complete --install` now installs a lazy hook that activates after `conda activate` / env activation.
10. CLI/Docs: minor help text refinements (main menu footer, motif help includes anchors).
11. Added `verify` module (stats-only QC) and CLI subcommand for constraint/architecture checks.
12. Verify/Lenstat: refactored shared length-stat formatting and improved `verify` reporting (module-style verdicts, aligned fields, and actionable column/motif summaries).
13. Primer module: added per-set primer design via `oligo_sets` (cross-set dimer checks + per-set Tm pairing) and CLI flag `--oligo-sets`.
14. Element modules: added Patch Mode (`patch_mode` / `--patch-mode`) for barcode/primer/motif/spacer, including reuse of existing per-set primers.
15. Validation/CLI: small robustness fixes for Patch Mode fills (context concatenation tolerates missing values; Notes epilog handles nested bullets).
16. Bugfix: Patch Mode context extraction now slices to missing rows before edge parsing (avoids list/boolean indexing crash in barcode/motif/spacer).
17. CLI help: hide argparse multi-value `[...]` token for cleaner alignment (e.g., `--cross-barcode-columns`).
18. DataFrame UX: normalize `ID` to string, return DataFrames with an explicit `ID` column by default (preserve `ID` as index only when the input DataFrame used it), and write CSVs without an extra numeric index column.
19. CLI/Docs: reordered modules to list `barcode` first; synced module ordering and descriptions across CLI, `__init__.py`, and README.
20. Docs: streamlined `help(op)` docstring for scannable Jupyter/REPL usage (compact module listing with one-line descriptions).
21. Docs: added comprehensive `docs/docs.md` with TOC, examples, workflows, and tips.
22. API: `excluded_motifs` DataFrame/CSV now only requires an `Exmotif` column (removed unnecessary `ID` column requirement).
23. API: `background` and `excluded_motifs` now accept FASTA files (.fa/.fasta/.fna, optionally gzipped).
24. Docstrings: clarified `acount` (barcode+variant verification) vs `xcount` (barcode-only counting, single or combinatorial).
25. Docs: added `docs/agent-skills.md` for AI assistant integration.

2026.01.17
1. Development in progress.
2. Docstring/comment spelling fixes.
3. Internal variable spelling fixes (no behavior change).
4. CLI main menu: added short description and refined command summaries.
5. CLI subcommand help: added one-line descriptions.
6. CLI subcommand help: added footer notes (from module docstrings).
7. CLI subcommand help: fixed notes bullet formatting.
8. CLI subcommand help: show full notes and normalize whitespace.
9. CLI manual: clarified usage and examples.
10. CLI manual: expanded topics listing to include meta topics.
11. CLI manual: improved topics formatting.
12. CLI subcommand help: translate docstring parameter names to CLI flags.
13. `help(oligopool)`: show package manual only.
14. CLI: allow simple Python-style sequence-constraint expressions (e.g. "'N'*20") for primer/motif constraints.
15. Docstrings: corrected minor typos and clarified a few parameter descriptions.
16. CLI: allow shorthand sequence-constraint concatenation (e.g. "GCC+N*20+CCG") for primer/motif constraints.
17. CLI/argcomplete: reduced completion latency by deferring heavy imports until after argument parsing.

2026.01.16
1. Version bump.
2. Standardized RNG seeding for stochastic design/assembly modules.
3. Added baseline stats keys across modules (module, input_rows, output_rows).
4. Improved validation messages with row-context examples.
5. Clarified public API exports via `__all__`.
6. Added CLI entry points (`oligopool`, `op`) including `manual` command and subcommands for major modules.
7. CLI UX: help formatting, "did you mean" suggestions, and `op complete` autocompletion installer (argcomplete required).
8. Documentation updates: README citation/version, CLI usage, completion instructions, and Docker notes.
9. Docker updates: added `Dockerfile` and refreshed container setup notes.
10. CLI help command listing: `manual` first, `complete` last.
11. CLI `xcount` options: removed unsupported associate error flag.
12. Documented that `acount`/`xcount` callbacks are Python-only (not supported in CLI mode).
13. Added `cite` CLI command to print citation information and paper link.
14. `op complete` now prints an inline hint about `--install` when run without flags.
15. Improved `op manual` autocompletion for manual topics.
16. `op complete --install [bash|zsh|fish]` now works (shell arg is optional).
17. Updated library docs examples to match current API parameter names.
18. `op complete` output is banner-free for `eval`/sourcing.

2026.01.15
1. Version bump.

2024.12.02
1. Fixed a small `core_count` print gap.
2. Version bump.

2024.12.01
1. Update version.
2. Show Manager dict in notebook.

2024.11.24
1. Update `__init__` description.
2. Auto-clamped primers in `primer` module.
3. Fix verbosity issue in `lenstat`.
4. Fixed degenerate constraint issue in `motif`.
5. Automated robust `homology` initialization.
6. Adjusted edge-effect length calculation for repeats.
7. Adding `maximum_repeat_length` to `motif` and `spacer`.
8. Refactor non-repetitiveness of primer, spacer and motif,

2024.11.03
1. Changed automatic padding constraints.
2. Introduced design of constants as a motif type.
3. Added Orphan Oligo indicators to `motif` and `spacer`.
4. Updated the oligo length limit parsing.
5. Updated edge-effect calculation in `pad`.
6. Removed `bounter` since deprecated.
7. Switched to `ShareDB`.
8. Added robust parallelization for POSIX systems.

2024.10.29
1. Fix edge extraction and edge-effect length calc.
2. Update `barcode` to include returning orphan oligos on failure.

2024.10.27
1. Drop two-bit encoding idea.
2. Use set-based non-repetitiveness per context.

2024.10.24
1. Updated stats variables in `final` and `lenstat`.
2. Added `revcomp` and `merge` auxiliary modules.
3. Bugfix in `utils`.

2024.10.20
1. Updated `OligopoolCalculatorInAction.ipynb`.
2. Updated `acount` and `xcount`.

2024.10.08
1. Updated `README.md`.
2. Updated `OligopoolCalculatorInAction.ipynb`.
3. Removed shell scripts -- don't need 'em.
4. Updated `index`, `background`, `lenstat`, `pack`, `acount`, and `xcount`.
5. Updated `design_parser.py`.
6. Updated `analysis_pipeline.py`.

2024.10.07
1. Updated `README.md`.
2. Finalize and test `design_parser.py`.

2024.10.05
1. Updated docstrings in `lenstat`, `vectorDB`.
2. Updated `spacer`, `motif`, `final`.
3. Reverting to stats-only for `lenstat`.
4. Removed `deprecated` directory (will use commit history).
5. Updated docstrings for `split` and `pad`.
6. `Design Mode` is completed.

2024.10.02
1. Updated `__init__.py`
2. Updated `OligopoolCalculatorInAction.ipynb`.
3. Updated `README.md`.
4. Mitigated `background` open issue on error in `primer`.
5. Updated docstrings in `barcode`, `primer`, and `background`.

2024.09.30
1. Updated `README.md`.
2. Started `OligopoolCalculatorInAction.ipynb`.
3. Expose `vectorDB` and fixed wrapping.
4. Updated `background` and finished docstring and stats dict.
5. Updated `primer` and finished docstring and stats dict.

2024.09.29
1. Switched to repeat elimination via edge-effects in `barcode` design.
2. Update `__init__.py`.
3. Update `README.md`.
4. Update `__version__`.
5. Update `setup.py`.

2024.09.27
1. Updated `setup.py`.
2. Switching to `plyvel` because `LevelDB` no longer a thing.
3. Updated `vectorDB`.
4. Removed TODOs.
5. Restructured packaging.
6. Updated interfaces.
