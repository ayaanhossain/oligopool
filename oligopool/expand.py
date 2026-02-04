import time as tt

import atexit as ae

import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp
from .base import core_degenerate as cd

from typing import Tuple


def expand(
    input_data:str|pd.DataFrame,
    sequence_column:str,
    mapping_file:str|pd.DataFrame|None=None,
    output_file:str|None=None,
    expansion_limit:int|None=None,
    verbose:bool=True) -> Tuple[pd.DataFrame, dict]:
    '''
    Expand IUPAC-degenerate sequences into all concrete A/T/G/C sequences.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with degenerate sequences.
        - `sequence_column` (`str`): Column name containing IUPAC-degenerate sequences to expand.

    Optional Parameters:
        - `mapping_file` (`str` / `pd.DataFrame` / `None`): Mapping output from `compress` (default: `None`). See Notes.
        - `output_file` (`str` / `None`): Filename for output DataFrame (default: `None`).
            A `.oligopool.expand.csv` suffix is added if missing.
        - `expansion_limit` (`int` / `None`): Safety cap for maximum total expanded sequences;
            if estimated expansion exceeds this limit, expansion fails (default: `None` for no limit).
        - `verbose` (`bool`): If `True`, logs progress to stdout (default: `True`).

    Returns:
        - A pandas DataFrame with expanded sequences; saves to `output_file` if specified.
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` must contain an 'ID' column (or 'DegenerateID' from `compress` output).
        - Useful as a sanity-check that `compress` output covers exactly (and only) the original sequences.
        - `mapping_file` must contain `ID` and `DegenerateID` columns (the `mapping_df` output from `compress`).
            When provided, output includes both `DegenerateID` and original `ID` columns.
        - If `mapping_file` is a string path/basename, `.oligopool.compress.mapping.csv` is appended if missing.
        - Expansion can be exponential (e.g., 10 N's = 4^10 sequences); use `expansion_limit`
            as a safety cap when working with highly degenerate sequences.
    '''

    # Argument Aliasing
    indata   = input_data
    seqcol   = sequence_column
    mapfile  = mapping_file
    outfile  = output_file
    limit    = expansion_limit
    verbose  = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Expansion Verbiage Print
    liner.send('\n[Oligopool Calculator: Degenerate Mode - Expand]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # Parse and validate input data
    (indf,
     indata_valid) = vp.get_parsed_expand_data_info(
        data=indata,
        data_field='     Input Data  ',
        sequence_column=seqcol,
        sequence_column_field='  Sequence Column',
        liner=liner)

    input_rows = len(indf.index) if isinstance(indf, pd.DataFrame) else 0

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Validate mapping file (optional; auto-suffixed to `.oligopool.compress.mapping.csv` when provided as a basename)
    (mapdf,
     mapfile_valid) = vp.get_parsed_compress_mapping_info(
        mapping_data=mapfile,
        mapping_field='   Mapping File  ',
        liner=liner,
        required=False)

    # Validate output file
    outfile_valid = vp.get_outdf_validity(
        outdf=outfile,
        outdf_suffix='.oligopool.expand.csv',
        outdf_field='    Output File  ',
        liner=liner)

    # Validate expansion_limit (optional)
    limit_valid = vp.get_optional_numeric_validity(
        numeric=limit,
        numeric_field=' Expansion Limit ',
        numeric_pre_desc=' At most ',
        numeric_post_desc=' Total Expansion(s)',
        minval=1,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # First Pass Validation
    if not all([
        indata_valid,
        mapfile_valid,
        outfile_valid,
        limit_valid]):
        liner.send('\n')
        raise RuntimeError('Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Adjust outfile suffix
    if outfile is not None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.expand.csv')

    # Adjust Numeric Parameters
    if limit is not None:
        limit = round(limit)

    # Setup stats
    stats = {
        'status': False,
        'basis': 'unsolved',
        'step': 0,
        'step_name': 'initializing',
        'vars': {},
        'warns': {},
    }

    # Schedule file cleanup on error
    ofdeletion = None
    if outfile is not None:
        ofdeletion = ae.register(
            ut.remove_file,
            outfile)

    # Step 1: Parse degenerate sequences
    liner.send('\n[Step 1: Parsing Degenerate Sequences]\n')
    stats['step'] = 1
    stats['step_name'] = 'parsing-sequences'

    # Parse degenerate sequences
    (num_input_seqs,
     total_estimated) = cd.get_parsed_degenerate_info(
        indf=indf,
        seqcol=seqcol,
        liner=liner)

    # Check expansion limit
    if limit is not None and total_estimated > limit:

        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 1,
            'step_name': 'parsing-sequences',
            'vars'    : {
                'input_sequences': num_input_seqs,
                'estimated_expansion': int(total_estimated),
                'expansion_limit': int(limit),
            },
            'warns'   : {},
        }

        liner.send(' Verdict: Estimated Expansion Exceeds Limit\n')
        liner.send(' Time Elapsed: {:.2f} sec\n'.format(tt.time() - t0))

        liner.close()

        stats = ut.stamp_stats(
            stats=stats,
            module='expand',
            input_rows=input_rows,
            output_rows=0)
        return (None, stats)

    liner.send(' Time Elapsed: {:.2f} sec\n'.format(tt.time() - t0))

    # Step 2: Expand sequences
    liner.send('\n[Step 2: Expanding Sequences]\n')
    stats['step'] = 2
    stats['step_name'] = 'expanding-sequences'

    # Expand all sequences
    (expanded_rows,
     total_expanded) = cd.get_expanded_dataframe(
        indf=indf,
        seqcol=seqcol,
        liner=liner)

    liner.send(' Time Elapsed: {:.2f} sec\n'.format(tt.time() - t0))

    # Build output DataFrame
    output_df = pd.DataFrame(expanded_rows)

    # Join with mapping to restore original IDs
    if mapdf is not None:
        # Current 'ID' column contains DegenerateID values
        # Rename to 'DegenerateID' before join
        output_df.rename(columns={'ID': 'DegenerateID'}, inplace=True)

        # Check if mapping has Sequence column for exact matching
        if 'Sequence' in mapdf.columns:
            # Exact match: join on DegenerateID AND ExpandedSeq == Sequence
            # This gives 1:1 mapping from expanded sequence to original ID
            output_df = output_df.merge(
                mapdf[['ID', 'Sequence', 'DegenerateID']],
                left_on=['DegenerateID', 'ExpandedSeq'],
                right_on=['DegenerateID', 'Sequence'],
                how='left')
            # Drop the redundant Sequence column from mapping
            output_df.drop(columns=['Sequence'], inplace=True)
        else:
            # Fallback: join only on DegenerateID (may produce multiple rows per expanded seq)
            output_df = output_df.merge(
                mapdf[['ID', 'DegenerateID']],
                on='DegenerateID',
                how='left')

        # Reorder to put ID first, then DegenerateID
        id_cols = ['ID', 'DegenerateID']
        other_cols = [c for c in output_df.columns if c not in id_cols]
        output_df = output_df[id_cols + other_cols]

        liner.send(' Mapped {:,} Expanded Sequence(s) to Original IDs\n'.format(
            len(output_df)))

    # Reorder columns to put ExpandedSeq and OligoLength at end
    cols = [c for c in output_df.columns if c not in ('ExpandedSeq', 'OligoLength')]
    cols.extend(['ExpandedSeq', 'OligoLength'])
    output_df = output_df[cols]

    # Write output
    if outfile is not None:
        ut.write_df_csv(
            df=output_df,
            outfile=outfile,
            sep=',')

    # Compute statistics
    expansion_factor = total_expanded / num_input_seqs \
        if num_input_seqs > 0 else 1.0

    # Build stats dictionary
    stats['status'] = True
    stats['basis'] = 'solved'
    stats['vars'] = {
        'input_sequences': num_input_seqs,
        'expanded_sequences': total_expanded,
        'expansion_factor': round(expansion_factor, 2),
        'expansion_limited': limit is not None,
        'ids_mapped': mapdf is not None,
    }

    # Expansion Statistics
    liner.send('\n[Expansion Statistics]\n')

    plen = ut.get_printlen(
        value=max(stats['vars'][field] for field in (
            'input_sequences',
            'expanded_sequences')))

    liner.send('  Expansion Status     : Successful\n')
    liner.send(' Degenerate Oligo(s)   : {:{},d}\n'.format(
        num_input_seqs, plen))
    liner.send('   Concrete Sequence(s): {:{},d}\n'.format(
        total_expanded, plen))
    liner.send('  Expansion Factor     : {:.1f}x\n'.format(
        expansion_factor))
    liner.send('  Expansion Limited    : {}\n'.format(
        'Yes' if limit is not None else 'No'))
    liner.send(' Time Elapsed: {:.2f} sec\n'.format(tt.time() - t0))

    # Unschedule file deletion
    if ofdeletion is not None:
        ae.unregister(ofdeletion)

    # Close Liner
    liner.close()

    # Stamp stats
    stats = ut.stamp_stats(
        stats=stats,
        module='expand',
        input_rows=input_rows,
        output_rows=len(output_df))

    return (output_df, stats)
