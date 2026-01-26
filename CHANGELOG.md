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
