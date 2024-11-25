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