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
    Verify an oligo pool for length, excluded motif emergence, and background k-mer conflicts.
    Returns a DataFrame with per-row conflict flags and details.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with an 'ID' column
            and at least one DNA column (ATGC only).
        - `oligo_length_limit` (`int`): Maximum allowed oligo length (≥ 4).

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
        - Length conflict: oligo length exceeds `oligo_length_limit`.
        - Exmotif conflict: motif count exceeds library-wide minimum (baseline).
            `excluded_motifs`: one or more sources (list/dict), merged; strict ATGC only; CSV/DataFrame requires 'Exmotif';
            FASTA sources are supported.
        - Background conflict: any k-mer in an oligo matches a background DB; with multiple DBs, any match flags a conflict.
        - Conflict details are stored as dicts in the DataFrame; when written to CSV, they are
            serialized as JSON strings (parse back with `json.loads()` for any `*Details` column; currently
            the `*ConflictDetails` columns).
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

    # Report DNA columns
    liner.send(' DNA Column(s): {} Column(s)\n'.format(len(dna_columns)))
    for col in dna_columns:
        liner.send('    - {}\n'.format(col))

    if used_complete_oligo:
        liner.send(' Complete Oligo : Present (using directly)\n')
    else:
        liner.send(' Complete Oligo : Not Present (concatenating above)\n')

    # Step 2: Build concatenated oligos
    oligo_sequences = []
    oligo_lengths = []

    for idx in range(len(indf)):
        if used_complete_oligo:
            oligo = str(indf['CompleteOligo'].iloc[idx]).replace('-', '').upper()
        else:
            parts = [str(indf[col].iloc[idx]).replace('-', '').upper() for col in dna_columns]
            oligo = ''.join(parts)
        oligo_sequences.append(oligo)
        oligo_lengths.append(len(oligo))

    min_oligo_len = min(oligo_lengths)
    max_oligo_len = max(oligo_lengths)

    if min_oligo_len == max_oligo_len:
        liner.send('    Oligo Length: {:,} Base Pair(s)\n'.format(min_oligo_len))
    else:
        liner.send('    Oligo Length: {:,} to {:,} Base Pair(s)\n'.format(
            min_oligo_len, max_oligo_len))

    # Step 3: Conflict Detection
    # 3a. Length Conflict Detection
    length_conflicts = []
    length_details = []

    for length in oligo_lengths:
        if length > oligolimit:
            length_conflicts.append(True)
            length_details.append({
                'length': length,
                'limit': oligolimit,
                'overflow': length - oligolimit
            })
        else:
            length_conflicts.append(False)
            length_details.append(None)

    # 3b. Exmotif Emergence Detection
    exmotif_conflicts = []
    exmotif_details = []
    baselines = {}
    emergence_summary = {}

    if exmotifs:
        # Pre-compute baselines
        for motif in exmotifs:
            counts = [oligo.count(motif) for oligo in oligo_sequences]
            baselines[motif] = min(counts) if counts else 0

        # Per-row detection
        for oligo in oligo_sequences:
            row_details = {}
            for motif in exmotifs:
                count = oligo.count(motif)
                if count > baselines[motif]:
                    # Find all positions
                    positions = []
                    start = 0
                    while True:
                        pos = oligo.find(motif, start)
                        if pos < 0:
                            break
                        positions.append(pos)
                        start = pos + 1
                    row_details[motif] = {
                        'count': count,
                        'baseline': baselines[motif],
                        'emergence': count - baselines[motif],
                        'positions': positions
                    }
            if row_details:
                exmotif_conflicts.append(True)
                exmotif_details.append(row_details)
            else:
                exmotif_conflicts.append(False)
                exmotif_details.append(None)

        # Build emergence summary for stats
        for motif in exmotifs:
            emergent_rows = sum(1 for d in exmotif_details if d and motif in d)
            max_emergence = max(
                (d[motif]['emergence'] for d in exmotif_details if d and motif in d),
                default=0)
            emergence_summary[motif] = {
                'rows': emergent_rows,
                'max': max_emergence
            }
    else:
        exmotif_conflicts = [False] * len(oligo_sequences)
        exmotif_details = [None] * len(oligo_sequences)

    # 3c. Background K-mer Conflict Detection
    background_conflicts = []
    background_details = []

    if backgrounds:
        for oligo in oligo_sequences:
            row_details = {
                'matching_kmers': [],
                'kmer_count': 0,
                'first_match_position': None
            }
            has_conflict = False

            for bg in backgrounds:
                K = bg.K
                if len(oligo) < K:
                    continue
                for i in range(len(oligo) - K + 1):
                    kmer = oligo[i:i+K]
                    if kmer in bg:
                        has_conflict = True
                        row_details['matching_kmers'].append(kmer)
                        if row_details['first_match_position'] is None:
                            row_details['first_match_position'] = i
                        row_details['kmer_count'] += 1

            if has_conflict:
                background_conflicts.append(True)
                background_details.append(row_details)
            else:
                background_conflicts.append(False)
                background_details.append(None)
    else:
        background_conflicts = [False] * len(oligo_sequences)
        background_details = [None] * len(oligo_sequences)

    # Build HasAnyConflicts
    has_any_conflicts = [
        l or e or b
        for l, e, b in zip(length_conflicts, exmotif_conflicts, background_conflicts)
    ]

    # Step 4: Build Output DataFrame (index from indf, ID column added at return)
    outdf = pd.DataFrame(index=indf.index)
    outdf['CompleteOligo'] = oligo_sequences
    outdf['OligoLength'] = oligo_lengths
    outdf['HasLengthConflict'] = length_conflicts
    outdf['HasExmotifConflict'] = exmotif_conflicts
    outdf['HasBackgroundConflict'] = background_conflicts
    outdf['HasAnyConflicts'] = has_any_conflicts
    outdf['LengthConflictDetails'] = length_details
    outdf['ExmotifConflictDetails'] = exmotif_details
    outdf['BackgroundConflictDetails'] = background_details

    liner.send(' Scanning Completed\n')
    liner.send(' Time Elapsed: {:.2f} sec\n'.format(tt.time() - t0))

    # Step 5: Build Stats & Write Output
    length_conflict_count = sum(length_conflicts)
    exmotif_conflict_count = sum(exmotif_conflicts)
    background_conflict_count = sum(background_conflicts)
    any_conflict_count = sum(has_any_conflicts)
    clean_count = len(outdf) - any_conflict_count

    stats = {
        'status': any_conflict_count == 0,
        'basis': 'verified' if any_conflict_count == 0 else 'conflicts',
        'step': 1,
        'step_name': 'scanning-oligo-pool',
        'vars': {
            'total_count': len(outdf),
            'dna_columns': dna_columns,
            'complete_oligo': used_complete_oligo,
            'oligo_limit': oligolimit,
            'min_oligo_len': min_oligo_len,
            'max_oligo_len': max_oligo_len,
            'length_conflict': length_conflict_count,
            'exmotif_conflict': exmotif_conflict_count,
            'background_conflict': background_conflict_count,
            'any_conflict': any_conflict_count,
            'clean_count': clean_count,
            'excluded_motifs': list(exmotifs) if exmotifs else None,
            'motif_baselines': baselines if baselines else None,
            'motif_emergence': emergence_summary if emergence_summary else None,
            'background_k': [bg.K for bg in backgrounds] if backgrounds else None,
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
        for col in ['LengthConflictDetails', 'ExmotifConflictDetails', 'BackgroundConflictDetails']:
            csv_df[col] = csv_df[col].apply(lambda x: json.dumps(x) if x is not None else None)
        # Add ID column for CSV output
        csv_df = ut.get_df_with_id_column(csv_df)
        ut.write_df_csv(df=csv_df, outfile=outfile, sep=',')

    # Verification Statistics
    liner.send('\n[Verification Statistics]\n')

    total_count = len(outdf)
    plen = ut.get_printlen(value=total_count)

    verify_status = 'Successful' if any_conflict_count == 0 else 'Conflicts Detected'
    liner.send('     Verify Status   : {}\n'.format(verify_status))
    liner.send('      Total Count    : {:{},d} Row(s)\n'.format(total_count, plen))
    liner.send('     Length Conflicts: {:{},d} Row(s) ({:6.2f} %)\n'.format(
        length_conflict_count, plen,
        ut.safediv(A=length_conflict_count * 100., B=total_count)))
    liner.send('    Exmotif Conflicts: {:{},d} Row(s) ({:6.2f} %)\n'.format(
        exmotif_conflict_count, plen,
        ut.safediv(A=exmotif_conflict_count * 100., B=total_count)))
    liner.send(' Background Conflicts: {:{},d} Row(s) ({:6.2f} %)\n'.format(
        background_conflict_count, plen,
        ut.safediv(A=background_conflict_count * 100., B=total_count)))
    liner.send('        Any Conflicts: {:{},d} Row(s) ({:6.2f} %)\n'.format(
        any_conflict_count, plen,
        ut.safediv(A=any_conflict_count * 100., B=total_count)))
    liner.send('      Clean Count    : {:{},d} Row(s) ({:6.2f} %)\n'.format(
        clean_count, plen,
        ut.safediv(A=clean_count * 100., B=total_count)))

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
