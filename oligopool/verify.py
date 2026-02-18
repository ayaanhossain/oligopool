import json
import time as tt

import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp
from .base import vectordb as db

from typing import Tuple


def verify(
    input_data:str|pd.DataFrame,
    oligo_length_limit:int,
    output_file:str|None=None,
    excluded_motifs:list|str|dict|pd.DataFrame|None=None,
    background_directory:str|list|None=None,
    verbose:bool=True) -> Tuple[pd.DataFrame, dict]:
    '''
    Verify an oligo pool for integrity, length, excluded motif emergence, and background k-mer conflicts.
    Returns a DataFrame with per-row conflict flags and details.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with an 'ID' column
            and at least one DNA column (ATGC only).
        - `oligo_length_limit` (`int`): Maximum allowed oligo length (>= 4).

    Optional Parameters:
        - `output_file` (`str` / `None`): Output CSV filename; required in CLI usage, optional in
            library usage (default: `None`). A `.oligopool.verify.csv` suffix is added if missing.
        - `excluded_motifs` (`list` / `str` / `dict` / `pd.DataFrame` / `None`): Excluded motif source(s) (default: `None`).
        - `background_directory` (`str` / `list` / `None`): Background k-mer DB(s) from `background()`
            (default: `None`).
        - `verbose` (`bool`): If `True`, logs progress to stdout (default: `True`).

    Returns:
        - A pandas DataFrame with conflict flags and details for each oligo.
        - A dictionary of stats from the verification.

    Notes:
        - If `CompleteOligo` column exists, it is used directly; otherwise, all concrete DNA columns
            (ATGC only) are concatenated left-to-right.
        - IUPAC/degenerate columns are skipped silently during DNA column detection.
        - Integrity conflict: `CompleteOligo` vs constituent mismatch or `OligoLength` vs actual mismatch.
        - Length conflict: oligo length exceeds `oligo_length_limit`.
        - Exmotif conflict: motif count exceeds library-wide minimum (baseline).
            `excluded_motifs`: one or more sources (list/dict), merged; strict ATGC only; CSV/DataFrame requires 'Exmotif';
            FASTA sources are supported.
        - Background conflict: any k-mer in an oligo matches a background DB; with multiple DBs, any match flags a conflict.
        - Conflict details are stored as dicts in the DataFrame; when written to CSV, they are
            serialized as JSON strings (parse back with `json.loads()` for any `*Details` column; currently
            the `*ConflictDetails` columns).
        - When constituent DNA columns are present alongside `CompleteOligo`, conflict positions are
            attributed to their originating column. The `columns` key in `ExmotifConflictDetails` and
            `BackgroundConflictDetails`, and the `column_lengths` key in `LengthConflictDetails`, provide
            this attribution. Rows with a `CompleteOligo`-vs-constituent mismatch attribute exmotif/background
            hits to `CompleteOligo` instead of constituent columns. Stats include `excluded_motif_column_conflicts`
            and `background_column_conflicts` dicts keyed by column name.
    '''

    # Preserve return style when the caller intentionally used ID as index.
    id_from_index = ut.get_id_index_intent(input_data)

    # Argument Aliasing
    indata     = input_data
    oligolimit = oligo_length_limit
    outfile    = output_file
    exmotifs   = excluded_motifs
    background = background_directory
    verbose    = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Verify Verbiage Print
    liner.send('\n[Oligopool Calculator: QC Mode - Verify]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # Full indata Parsing and Validation (requires ID + at least one DNA column)
    (indf,
    data_name,
    dna_columns,
    constituent_columns,
    indata_valid) = vp.get_parsed_verify_indata_info(
        data=indata,
        data_field='      Input Data    ',
        liner=liner)
    input_rows = len(indf.index) if isinstance(indf, pd.DataFrame) else 0

    # Show update on successful parse
    if indata_valid:
        liner.send(
            '      Input Data    : {} w/ {:,} Record(s)\n'.format(
                data_name,
                len(indf.index)))

    # Full oligolimit Validation
    oligolimit_valid = vp.get_numeric_validity(
        numeric=oligolimit,
        numeric_field='      Oligo Limit   ',
        numeric_pre_desc=' At most ',
        numeric_post_desc=' Base Pair(s)',
        minval=4,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full outfile Validation
    outfile_valid = vp.get_outdf_validity(
        outdf=outfile,
        outdf_suffix='.oligopool.verify.csv',
        outdf_field='     Output File    ',
        liner=liner)

    # Adjust outfile Suffix
    if outfile is not None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.verify.csv')

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full exmotifs Parsing and Validation (optional)
    (exmotifs,
    exmotifs_valid,
    exmotif_inputs) = vp.get_parsed_exmotifs(
        exmotifs=exmotifs,
        exmotifs_field='   Excluded Motifs  ',
        liner=liner)

    # Full background Parsing and Validation
    (background_valid,
    backgrounds_info,
    _) = vp.get_parsed_backgrounds(
        backgrounds=background,
        backgrounds_field=' Background Database',
        liner=liner)

    # First Pass Validation
    if not all([
        indata_valid,
        oligolimit_valid,
        outfile_valid,
        exmotifs_valid,
        background_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Open Backgrounds
    backgrounds = []
    opened_backgrounds = []
    if backgrounds_info:
        for bg_ref, bg_type, _ in backgrounds_info:
            if bg_type == 'path':
                bg_path = ut.get_adjusted_path(
                    path=bg_ref,
                    suffix='.oligopool.background')
                bg = db.vectorDB(
                    path=bg_path,
                    maximum_repeat_length=None)
                opened_backgrounds.append(bg)
            else:
                bg = bg_ref
            backgrounds.append(bg)

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Parameters
    oligolimit = round(oligolimit)

    # Book-keeping
    outdf = None
    stats = None
    warns = {}

    # Step 1: Scan Oligo Pool
    liner.send('\n[Step 1: Scanning Oligo Pool]\n')

    # Check if CompleteOligo was detected
    used_complete_oligo = (dna_columns == ['CompleteOligo'])

    # Column list for position attribution:
    # - Scenario A (CompleteOligo + constituents): use constituent columns
    # - Scenario B (CompleteOligo only): use ['CompleteOligo'] as single column
    # - Scenario C (no CompleteOligo): use constituent columns (same as dna_columns)
    attribution_columns = constituent_columns if constituent_columns else ['CompleteOligo']

    # Report all DNA columns
    all_dna_columns = list(constituent_columns)
    if used_complete_oligo:
        all_dna_columns.append('CompleteOligo')
    liner.send(' DNA Column(s): {} Column(s)\n'.format(len(all_dna_columns)))
    for col in all_dna_columns:
        liner.send('   - {}\n'.format(col))

    # Padded row count width for progress lines
    num_rows = len(indf)
    plen = ut.get_printlen(value=num_rows)

    # Step 2: Build concatenated oligos
    oligo_sequences = []
    oligo_lengths = []

    for idx in range(num_rows):
        if (idx + 1) % 1000 == 0 or idx + 1 == num_rows:
            liner.send('     Oligos Constructed: Scanned {:>{},} Variant(s)'.format(idx + 1, plen))
        if used_complete_oligo:
            oligo = str(indf['CompleteOligo'].iloc[idx]).replace('-', '').upper()
        else:
            parts = [str(indf[col].iloc[idx]).replace('-', '').upper() for col in dna_columns]
            oligo = ''.join(parts)
        oligo_sequences.append(oligo)
        oligo_lengths.append(len(oligo))
    liner.send('     Oligos Constructed: Scanned {:>{},} Variant(s)\n'.format(num_rows, plen))

    # Consistency check: CompleteOligo vs constituents (Scenario A only)
    complete_oligo_conflicts = [False] * num_rows
    if used_complete_oligo and constituent_columns:
        for idx in range(num_rows):
            if (idx + 1) % 1000 == 0 or idx + 1 == num_rows:
                liner.send(' CompleteOligo Verified: Scanned {:>{},} Variant(s)'.format(idx + 1, plen))
            concat = ''.join(
                str(indf[col].iloc[idx]).replace('-', '').upper()
                for col in constituent_columns)
            if oligo_sequences[idx] != concat:
                complete_oligo_conflicts[idx] = True
        liner.send(' CompleteOligo Verified: Scanned {:>{},} Variant(s)\n'.format(num_rows, plen))

    # OligoLength check (all scenarios, when column present)
    oligo_length_conflicts = [False] * num_rows
    oligo_length_series = None
    if 'OligoLength' in indf.columns:
        oligo_length_series = pd.to_numeric(indf['OligoLength'], errors='coerce')
        if oligo_length_series.notna().all():
            for idx in range(num_rows):
                if (idx + 1) % 1000 == 0 or idx + 1 == num_rows:
                    liner.send('   OligoLength Verified: Scanned {:>{},} Variant(s)'.format(idx + 1, plen))
                expected = int(oligo_length_series.iloc[idx])
                actual = oligo_lengths[idx]
                if expected != actual:
                    oligo_length_conflicts[idx] = True
            liner.send('   OligoLength Verified: Scanned {:>{},} Variant(s)\n'.format(num_rows, plen))

    # Build per-row integrity conflict data
    integrity_conflicts = []
    integrity_details = []
    for idx in range(num_rows):
        co = complete_oligo_conflicts[idx]
        ol = oligo_length_conflicts[idx]
        if co or ol:
            integrity_conflicts.append(True)
            integrity_details.append({
                'complete_oligo_conflict': co,
                'oligo_length_conflict': ol,
                'expected_oligo_length': int(oligo_length_series.iloc[idx]) if ol else None,
                'actual_oligo_length': oligo_lengths[idx] if ol else None,
            })
        else:
            integrity_conflicts.append(False)
            integrity_details.append(None)

    integrity_conflict_count = sum(integrity_conflicts)
    complete_oligo_conflict_count = sum(complete_oligo_conflicts)
    oligo_length_conflict_count = sum(oligo_length_conflicts)

    min_oligo_len = min(oligo_lengths)
    max_oligo_len = max(oligo_lengths)

    # Precompute per-row column boundaries for attribution
    row_boundaries = []
    for idx in range(num_rows):
        if (idx + 1) % 1000 == 0 or idx + 1 == num_rows:
            liner.send('     Bounds Constructed: Scanned {:>{},} Variant(s)'.format(idx + 1, plen))
        row_boundaries.append(
            ut.get_column_boundaries(idx, attribution_columns, indf))
    liner.send('     Bounds Constructed: Scanned {:>{},} Variant(s)\n'.format(num_rows, plen))

    # Override boundaries for rows with CompleteOligo integrity conflicts
    # (constituent boundaries are unreliable when CompleteOligo != concat)
    for idx in range(num_rows):
        if complete_oligo_conflicts[idx]:
            row_boundaries[idx] = [('CompleteOligo', 0, oligo_lengths[idx])]

    # Step 3: Conflict Detection
    # 3a. Length Conflict Detection
    length_conflicts = []
    length_details = []

    for idx, length in enumerate(oligo_lengths):
        if (idx + 1) % 1000 == 0 or idx + 1 == num_rows:
            liner.send('    Lengths Evaluated  : Scanned {:>{},} Variant(s)'.format(idx + 1, plen))
        if length > oligolimit:
            length_conflicts.append(True)
            # Build column_lengths from boundaries
            column_lengths = {}
            for col_name, start, end in row_boundaries[idx]:
                column_lengths[col_name] = end - start
            length_details.append({
                'oligo_length': length,
                'length_limit': oligolimit,
                'excess_bp': length - oligolimit,
                'column_lengths': column_lengths
            })
        else:
            length_conflicts.append(False)
            length_details.append(None)
    liner.send('    Lengths Evaluated  : Scanned {:>{},} Variant(s)\n'.format(num_rows, plen))

    # 3b. Exmotif Emergence Detection
    exmotif_conflicts = []
    exmotif_details = []
    baselines = {}
    emergence_summary = {}
    excluded_motif_column_conflicts = None

    if exmotifs:
        # Pre-compute baselines
        for motif in exmotifs:
            counts = [oligo.count(motif) for oligo in oligo_sequences]
            baselines[motif] = min(counts) if counts else 0

        # Per-row detection
        for row_idx, oligo in enumerate(oligo_sequences):
            if (row_idx + 1) % 1000 == 0 or row_idx + 1 == num_rows:
                liner.send('   Exmotifs Evaluated  : Scanned {:>{},} Variant(s)'.format(row_idx + 1, plen))
            row_details = []
            boundaries = row_boundaries[row_idx]
            for motif in exmotifs:
                count = oligo.count(motif)
                if count > baselines[motif]:
                    # Find all positions
                    positions = []
                    columns = []
                    start = 0
                    while True:
                        pos = oligo.find(motif, start)
                        if pos < 0:
                            break
                        positions.append(pos)
                        columns.append(ut.get_boundary_column(pos, boundaries))
                        # Non-overlapping semantics (match str.count()).
                        start = pos + len(motif)
                    row_details.append({
                        'motif': motif,
                        'occurrences': count,
                        'library_baseline': baselines[motif],
                        'excess_occurrences': count - baselines[motif],
                        'positions': positions,
                        'columns': columns
                    })
            if row_details:
                exmotif_conflicts.append(True)
                exmotif_details.append(row_details)
            else:
                exmotif_conflicts.append(False)
                exmotif_details.append(None)

        # Build emergence summary for stats
        for motif in exmotifs:
            emergent_rows = sum(
                1 for d in exmotif_details
                if d and any(e['motif'] == motif for e in d))
            max_emergence = max(
                (e['excess_occurrences'] for d in exmotif_details
                 if d for e in d if e['motif'] == motif),
                default=0)
            emergence_summary[motif] = {
                'emergent_row_count': emergent_rows,
                'max_excess_occurrences': max_emergence
            }

        # Compute excluded_motif_column_conflicts: per-column count of oligos with exmotif hits
        excluded_motif_column_conflicts = {col: 0 for col in attribution_columns}
        for detail in exmotif_details:
            if detail is None:
                continue
            # Collect all columns hit across all motifs in this row
            cols_hit = set()
            for entry in detail:
                cols_hit.update(entry['columns'])
            for col in cols_hit:
                if col in excluded_motif_column_conflicts:
                    excluded_motif_column_conflicts[col] += 1
        liner.send('   Exmotifs Evaluated  : Scanned {:>{},} Variant(s)\n'.format(num_rows, plen))
    else:
        exmotif_conflicts = [False] * len(oligo_sequences)
        exmotif_details = [None] * len(oligo_sequences)

    # 3c. Background K-mer Conflict Detection
    background_conflicts = []
    background_details = []
    background_column_conflicts = None

    if backgrounds:
        for row_idx, oligo in enumerate(oligo_sequences):
            if (row_idx + 1) % 1000 == 0 or row_idx + 1 == num_rows:
                liner.send(' Background Evaluated  : Scanned {:>{},} Variant(s)'.format(row_idx + 1, plen))
            row_details = {
                'kmers': [],
                'positions': [],
                'columns': [],
                'match_count': 0
            }
            has_conflict = False
            boundaries = row_boundaries[row_idx]

            for bg in backgrounds:
                K = bg.K
                if len(oligo) < K:
                    continue
                for i in range(len(oligo) - K + 1):
                    kmer = oligo[i:i+K]
                    if kmer in bg:
                        has_conflict = True
                        row_details['kmers'].append(kmer)
                        row_details['positions'].append(i)
                        row_details['columns'].append(
                            ut.get_boundary_column(i, boundaries))
                        row_details['match_count'] += 1

            if has_conflict:
                background_conflicts.append(True)
                background_details.append(row_details)
            else:
                background_conflicts.append(False)
                background_details.append(None)

        # Compute background_column_conflicts
        background_column_conflicts = {col: 0 for col in attribution_columns}
        for detail in background_details:
            if detail is None:
                continue
            cols_hit = set(detail['columns'])
            for col in cols_hit:
                if col in background_column_conflicts:
                    background_column_conflicts[col] += 1
        liner.send(' Background Evaluated  : Scanned {:>{},} Variant(s)\n'.format(num_rows, plen))
    else:
        background_conflicts = [False] * len(oligo_sequences)
        background_details = [None] * len(oligo_sequences)

    # Build HasAnyConflicts
    has_any_conflicts = [
        l or e or b or i
        for l, e, b, i in zip(length_conflicts, exmotif_conflicts,
                              background_conflicts, integrity_conflicts)
    ]

    # Step 4: Build Output DataFrame (index from indf, ID column added at return)
    outdf = pd.DataFrame(index=indf.index)
    outdf['CompleteOligo'] = oligo_sequences
    outdf['OligoLength'] = oligo_lengths
    outdf['HasIntegrityConflict'] = integrity_conflicts
    outdf['HasLengthConflict'] = length_conflicts
    outdf['HasExmotifConflict'] = exmotif_conflicts
    outdf['HasBackgroundConflict'] = background_conflicts
    outdf['HasAnyConflicts'] = has_any_conflicts
    outdf['IntegrityConflictDetails'] = integrity_details
    outdf['LengthConflictDetails'] = length_details
    outdf['ExmotifConflictDetails'] = exmotif_details
    outdf['BackgroundConflictDetails'] = background_details

    liner.send(' Time Elapsed: {:.2f} sec\n'.format(tt.time() - t0))

    # Step 5: Build Stats & Write Output
    length_conflict_count = sum(length_conflicts)
    excluded_motif_conflict_count = sum(exmotif_conflicts)
    background_conflict_count = sum(background_conflicts)
    any_conflict_count = sum(has_any_conflicts)
    clean_count = len(outdf) - any_conflict_count
    total_count = len(outdf)

    # Compute any_column_conflicts: per-column count of rows with any
    # position-level conflict (exmotif or background) in that column
    any_column_conflicts = None
    if excluded_motif_column_conflicts is not None or background_column_conflicts is not None:
        any_column_conflicts = {col: 0 for col in attribution_columns}
        for row_idx in range(num_rows):
            cols_hit = set()
            ex_detail = exmotif_details[row_idx]
            if ex_detail:
                for entry in ex_detail:
                    cols_hit.update(entry['columns'])
            bg_detail = background_details[row_idx]
            if bg_detail:
                cols_hit.update(bg_detail['columns'])
            for col in cols_hit:
                if col in any_column_conflicts:
                    any_column_conflicts[col] += 1

    stats = {
        'status': any_conflict_count == 0,
        'basis': 'verified' if any_conflict_count == 0 else 'conflicts',
        'step': 1,
        'step_name': 'scanning-oligo-pool',
        'vars': {
            'total_count':                   total_count,
            'dna_columns':                   dna_columns,
            'has_complete_oligo':            used_complete_oligo,
            'constituent_columns':           constituent_columns if constituent_columns else None,
            'length_limit':                  oligolimit,
            'min_oligo_length':              min_oligo_len,
            'max_oligo_length':              max_oligo_len,
            'integrity_conflict_count':      integrity_conflict_count,
            'complete_oligo_conflict_count': complete_oligo_conflict_count,
            'oligo_length_conflict_count':   oligo_length_conflict_count,
            'length_conflict_count':         length_conflict_count,
            'excluded_motif_conflict_count':        excluded_motif_conflict_count,
            'background_conflict_count':     background_conflict_count,
            'any_conflict_count':            any_conflict_count,
            'clean_count':                   clean_count,
            'excluded_motifs':               list(exmotifs) if exmotifs else None,
            'library_baselines':             baselines if baselines else None,
            'motif_emergence':               emergence_summary if emergence_summary else None,
            'background_kmer_lengths':       [bg.K for bg in backgrounds] if backgrounds else None,
            'excluded_motif_column_conflicts':      excluded_motif_column_conflicts,
            'background_column_conflicts':   background_column_conflicts,
            'any_column_conflicts':          any_column_conflicts,
        },
        'warns': warns
    }

    # Excluded-motif attribution for stats (uniform across modules)
    ut.stamp_exmotif_stats(
        stats=stats,
        exmotif_inputs=exmotif_inputs)

    # Write output CSV if specified
    if outfile is not None:
        # Serialize dict columns as JSON strings for CSV
        csv_df = outdf.copy()
        for col in ['IntegrityConflictDetails', 'LengthConflictDetails', 'ExmotifConflictDetails', 'BackgroundConflictDetails']:
            csv_df[col] = csv_df[col].apply(lambda x: json.dumps(x) if x is not None else None)
        # Add ID column for CSV output
        csv_df = ut.get_df_with_id_column(csv_df)
        ut.write_df_csv(df=csv_df, outfile=outfile, sep=',')

    # Verification Statistics
    liner.send('\n[Verification Statistics]\n')

    verify_status = 'Successful' if any_conflict_count == 0 else 'Conflicts Detected'
    liner.send('      Verify Status   : {}\n'.format(verify_status))
    liner.send('       Total Count    : {:{},d} Variant(s)\n'.format(total_count, plen))
    liner.send('   Integrity Conflicts: {:{},d} Variant(s) ({:6.2f} %)\n'.format(
        integrity_conflict_count, plen,
        ut.safediv(A=integrity_conflict_count * 100., B=total_count)))
    liner.send('      Length Conflicts: {:{},d} Variant(s) ({:6.2f} %)\n'.format(
        length_conflict_count, plen,
        ut.safediv(A=length_conflict_count * 100., B=total_count)))
    liner.send('     Exmotif Conflicts: {:{},d} Variant(s) ({:6.2f} %)\n'.format(
        excluded_motif_conflict_count, plen,
        ut.safediv(A=excluded_motif_conflict_count * 100., B=total_count)))
    liner.send('  Background Conflicts: {:{},d} Variant(s) ({:6.2f} %)\n'.format(
        background_conflict_count, plen,
        ut.safediv(A=background_conflict_count * 100., B=total_count)))
    liner.send('         Any Conflicts: {:{},d} Variant(s) ({:6.2f} %)\n'.format(
        any_conflict_count, plen,
        ut.safediv(A=any_conflict_count * 100., B=total_count)))
    liner.send('       Clean Count    : {:{},d} Variant(s) ({:6.2f} %)\n'.format(
        clean_count, plen,
        ut.safediv(A=clean_count * 100., B=total_count)))

    # Column-wise conflict breakdown (any = exmotif or background)
    if any_column_conflicts is not None and any_conflict_count > 0:
        liner.send(' Column-wise Conflicts:\n')
        ac_len = max(len(c) for c in attribution_columns)
        for col in attribution_columns:
            count = any_column_conflicts[col]
            liner.send('   - {:>{}} {:{},d} Variant(s) ({:6.2f} %)\n'.format(
                col, ac_len, count, plen,
                ut.safediv(A=count * 100., B=total_count)))

    # Show Time Elapsed
    liner.send(' Time Elapsed: {:.2f} sec\n'.format(tt.time() - t0))

    # Close Backgrounds
    for bg in opened_backgrounds:
        bg.close()

    # Close Liner
    liner.close()

    # Return Results
    stats = ut.stamp_stats(
        stats=stats,
        module='verify',
        input_rows=input_rows,
        output_rows=len(outdf))
    outdf_return = outdf
    if not id_from_index:
        outdf_return = ut.get_df_with_id_column(outdf)
    return (outdf_return, stats)
