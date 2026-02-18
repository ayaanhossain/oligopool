import json
import time as tt

import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp
from .base import vectordb as db
from .base import core_verify as cv

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

    # Run Engine
    result = cv.verify_engine(
        indf=indf,
        dna_columns=dna_columns,
        constituent_columns=constituent_columns,
        oligolimit=oligolimit,
        exmotifs=exmotifs,
        backgrounds=backgrounds,
        liner=liner)

    # Unpack engine results
    used_complete_oligo             = result['used_complete_oligo']
    attribution_columns             = result['attribution_columns']
    excluded_motif_column_conflicts = result['excluded_motif_column_conflicts']
    background_column_conflicts     = result['background_column_conflicts']
    any_column_conflicts            = result['any_column_conflicts']

    # DataFrame Assembly
    liner.send('  DataFrame Assembly   ...')
    outdf = pd.DataFrame(index=indf.index)
    outdf['CompleteOligo']              = result['oligo_sequences']
    outdf['OligoLength']                = result['oligo_lengths']
    outdf['HasIntegrityConflict']       = result['integrity_conflicts']
    outdf['HasLengthConflict']          = result['length_conflicts']
    outdf['HasExmotifConflict']         = result['exmotif_conflicts']
    outdf['HasBackgroundConflict']      = result['background_conflicts']
    outdf['HasAnyConflicts']            = result['has_any_conflicts']
    outdf['IntegrityConflictDetails']   = result['integrity_details']
    outdf['LengthConflictDetails']      = result['length_details']
    outdf['ExmotifConflictDetails']     = result['exmotif_details']
    outdf['BackgroundConflictDetails']  = result['background_details']
    liner.send('  DataFrame Assembly   : Completed\n')

    liner.send(' Time Elapsed: {:.2f} sec\n'.format(tt.time() - t0))

    # Step 5: Build Stats & Write Output
    integrity_conflict_count       = result['integrity_conflict_count']
    complete_oligo_conflict_count  = result['complete_oligo_conflict_count']
    oligo_length_conflict_count    = result['oligo_length_conflict_count']
    length_conflict_count          = sum(result['length_conflicts'])
    excluded_motif_conflict_count  = sum(result['exmotif_conflicts'])
    background_conflict_count      = sum(result['background_conflicts'])
    any_conflict_count             = sum(result['has_any_conflicts'])
    total_count = len(outdf)
    clean_count = total_count - any_conflict_count
    baselines         = result['baselines']
    emergence_summary = result['emergence_summary']

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
            'min_oligo_length':              result['min_oligo_length'],
            'max_oligo_length':              result['max_oligo_length'],
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

    plen = ut.get_printlen(value=total_count)
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
        ac_len = max(len(c) for c in any_column_conflicts)
        for col in any_column_conflicts:
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
