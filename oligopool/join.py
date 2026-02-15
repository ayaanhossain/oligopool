import time as tt

import atexit as ae

import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp

from typing import Tuple


def join(
    input_data:str|pd.DataFrame,
    other_data:str|pd.DataFrame,
    oligo_length_limit:int,
    join_policy:int|str,
    output_file:str|None=None,
    verbose:bool=True) -> Tuple[pd.DataFrame, dict]:
    '''
    Join two oligo tables on `ID` and insert new columns from `other_data` into the `input_data` column order.
    Useful for recombining parallel branch outputs into a single design table.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated oligo pool variants.
        - `other_data` (`str` / `pd.DataFrame`): Path to a second CSV file or DataFrame with the same IDs as `input_data`.
        - `oligo_length_limit` (`int`): Maximum allowed oligo length (>= 4).
        - `join_policy` (`int` / `str`): Ambiguity policy for column insertion:
            `0`/`left` chooses the left-most valid placement, `1`/`right` chooses the right-most.

    Optional Parameters:
        - `output_file` (`str` / `None`): Filename for output DataFrame; required in CLI usage,
            optional in library usage (default: `None`).
        - `verbose` (`bool`): If `True`, logs progress to stdout (default: `True`).

    Returns:
        - A pandas DataFrame of joined elements; saves to `output_file` if specified.
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - Both inputs must contain a unique 'ID' column; all other columns must be non-empty DNA strings.
        - `input_data` and `other_data` must contain exactly the same `ID` set (order may differ).
        - Any shared non-ID column must match exactly between the two inputs; mismatches are treated as infeasible.
        - New columns from `other_data` are inserted near their nearest shared anchors; `join_policy` resolves ambiguities.
    '''

    # Preserve return style when the caller intentionally used ID as index.
    id_from_index = ut.get_id_index_intent(input_data)

    # Argument Aliasing
    indata     = input_data
    otherdata  = other_data
    oligolimit = oligo_length_limit
    joinpolicy = join_policy
    outfile    = output_file
    verbose    = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Joining Verbiage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Join]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_indata_info(
        indata=indata,
        indata_field='  Input Data  ',
        required_fields=('ID',),
        precheck=False,
        liner=liner)
    input_rows = len(indf.index) if isinstance(indf, pd.DataFrame) else 0

    # Second Pass otherdata Parsing and Validation (precheck; report after ID matching)
    (otherdf,
    otherdata_name,
    otherdata_valid) = vp.get_parsed_indata_info(
        indata=otherdata,
        indata_field='  Other Data  ',
        required_fields=('ID',),
        precheck=True,
        liner=liner)

    # ID-set matching is part of Other Data validity/reporting.
    if indata_valid and otherdata_valid:
        other_ids_match = (len(indf.index) == len(otherdf.index)) and \
                          indf.index.difference(otherdf.index).empty and \
                          otherdf.index.difference(indf.index).empty
        if other_ids_match:
            otherdf = otherdf.reindex(indf.index)
            liner.send(
                '  Other Data  : {} w/ {:,} Record(s)\n'.format(
                    otherdata_name,
                    len(otherdf.index)))
        else:
            liner.send(
                '  Other Data  : {} w/ {:,} Record(s) [COLUMN=\'ID\' DOES NOT MATCH INPUT DATA]\n'.format(
                    otherdata_name,
                    len(otherdf.index)))
        otherdata_valid = otherdata_valid and other_ids_match
    elif otherdata_valid:
        liner.send(
            '  Other Data  : {} w/ {:,} Record(s)\n'.format(
                otherdata_name,
                len(otherdf.index)))

    # Full oligolimit Validation
    oligolimit_valid = vp.get_numeric_validity(
        numeric=oligolimit,
        numeric_field='  Oligo Limit ',
        numeric_pre_desc=' At most ',
        numeric_post_desc=' Base Pair(s)',
        minval=4,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full joinpolicy Validation
    (joinpolicy,
    joinpolicy_valid) = vp.get_typed_categorical_validity(
        category=joinpolicy,
        category_field='   Join Policy',
        category_pre_desc=' ',
        category_post_desc=' Ambiguity Resolution',
        type_name='join_policy',
        liner=liner)

    # Full outfile Validation
    outfile_valid = vp.get_outdf_validity(
        outdf=outfile,
        outdf_suffix='.oligopool.join.csv',
        outdf_field=' Output File  ',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.join.csv')

    # First Pass Validation
    if not all([
        indata_valid,
        otherdata_valid,
        oligolimit_valid,
        joinpolicy_valid,
        outfile_valid,
        ]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Schedule outfile deletion
    ofdeletion = ae.register(
        ut.remove_file,
        outfile)

    outdf = None
    stats = None
    warns = {}

    # Adjust Numeric Paramters
    oligolimit = round(oligolimit)

    # Show Update
    liner.send('\n[Step 1: Checking Join Compatibility]\n')
    st0 = tt.time()

    # Identify shared (overlapping) columns (ID is the index and excluded here).
    shared_columns = [
        col for col in indf.columns.to_list()
            if col in otherdf.columns
    ]

    mismatched_columns = []
    for col in shared_columns:
        if not indf[col].equals(otherdf[col]):
            mismatched_columns.append(col)

    col_printlen = ut.get_printlen(
        value=max(
            len(shared_columns),
            len(mismatched_columns)))

    liner.send(
        '     Shared Column(s): {:{},d} Column(s)\n'.format(
            len(shared_columns),
            col_printlen))

    # Compatibility verdict (treat mismatches as infeasible, not an argument error).
    parsestatus = len(mismatched_columns) == 0

    if not parsestatus:
        example_mismatch_column = mismatched_columns[0]
        mismatch_mask = indf[example_mismatch_column] != otherdf[example_mismatch_column]
        examples = ut.get_row_examples(
            df=indf,
            invalid_mask=mismatch_mask,
            id_col='ID',
            limit=5)
        example_note = ut.format_row_examples(
            examples,
            label='Mismatching ID examples')

        liner.send(
            ' Mismatched Column(s): {:{},d} Column(s) [INFEASIBLE] (example: \'{}\'){}\n'.format(
                len(mismatched_columns),
                col_printlen,
                example_mismatch_column,
                example_note))
        liner.send(
            ' Time Elapsed: {:.2f} sec\n'.format(
                tt.time()-st0))
        liner.send(
            ' Verdict: Join Infeasible due to Shared Column Mismatch\n')

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 1,
            'step_name': 'checking-join-compatibility',
            'vars'    : {
                            'oligo_limit': oligolimit,
                            'join_policy': joinpolicy,
                         'shared_columns': tuple(shared_columns),
                     'mismatched_columns': tuple(mismatched_columns),
                'example_mismatch_column': example_mismatch_column,
                            'example_ids': tuple(examples),
                },
            'warns'   : warns}

        # Return results
        liner.close()
        stats = ut.stamp_stats(
            stats=stats,
            module='join',
            input_rows=input_rows,
            output_rows=0)
        return (outdf, stats)

    liner.send(
        ' Mismatched Column(s): {:{},d} Column(s)\n'.format(
            0,
            col_printlen))
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-st0))
    liner.send(
        ' Verdict: Join Possibly Feasible\n')

    liner.send('\n[Step 2: Computing Join]\n')
    st0 = tt.time()

    # Compute New DataFrame
    outdf = indf.copy()
    output_columns = outdf.columns.to_list()
    source_columns = otherdf.columns.to_list()

    inserted_columns = []
    ignored_columns = []
    ambiguous_resolutions = 0

    for col in source_columns:
        if col in output_columns:
            ignored_columns.append(col)
            continue

        (insert_idx,
        was_ambiguous,
        _,
        _) = ut.get_anchor_insert_index(
            column=col,
            source_columns=source_columns,
            target_columns=output_columns,
            policy=joinpolicy)

        outdf.insert(insert_idx, col, otherdf[col])
        output_columns.insert(insert_idx, col)
        inserted_columns.append(col)
        if was_ambiguous:
            ambiguous_resolutions += 1

    # Show Update
    liner.send(' Join Completed\n')
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-st0))

    liner.send('\n[Step 3: Parsing Oligo Limit]\n')
    st0 = tt.time()

    variantlens    = ut.get_variantlens(indf=outdf)
    minoligolen    = int(variantlens.min())
    maxoligolen    = int(variantlens.max())
    overflow_mask  = variantlens > oligolimit
    overflow_count = int(overflow_mask.sum())
    oligo_printlen = ut.get_printlen(
        value=max(
            minoligolen,
            maxoligolen,
            oligolimit))
    overflow_printlen = ut.get_printlen(
        value=overflow_count)

    liner.send(
        ' Maximum Oligo Length: {:{},d} Base Pair(s)\n'.format(
            oligolimit,
            oligo_printlen))

    if minoligolen == maxoligolen:
        liner.send(
            '  Joined Oligo Length: {:{},d} Base Pair(s)\n'.format(
                minoligolen,
                oligo_printlen))
    else:
        liner.send(
            '  Joined Oligo Length: {:{},d} to {:{},d} Base Pair(s)\n'.format(
                minoligolen,
                oligo_printlen,
                maxoligolen,
                oligo_printlen))

    if overflow_count > 0:
        examples = ut.get_row_examples(
            df=outdf,
            invalid_mask=overflow_mask,
            id_col='ID',
            limit=5)
        example_note = ut.format_row_examples(
            examples,
            label='Overflowing ID examples')
        liner.send(
            '   Overflowing Oligos: {:{},d} Variant(s) [INFEASIBLE]{}\n'.format(
                overflow_count,
                overflow_printlen,
                example_note))
    else:
        liner.send(
            '   Overflowing Oligos: {:{},d} Variant(s)\n'.format(
                0,
                overflow_printlen))

    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-st0))

    if overflow_count > 0:
        liner.send(
            ' Verdict: Join Infeasible due to Oligo Limit Constraints\n')

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 3,
            'step_name': 'parsing-oligo-limit',
            'vars'    : {
                          'oligo_limit': oligolimit,
                       'limit_overflow': True,
                        'min_oligo_len': minoligolen,
                        'max_oligo_len': maxoligolen,
                       'overflow_count': overflow_count,
                          'example_ids': tuple(examples),
                       'shared_columns': tuple(shared_columns),
                       'columns_joined': inserted_columns,
                      'columns_ignored': ignored_columns,
                'ambiguous_resolutions': ambiguous_resolutions,
                          'join_policy': joinpolicy,
                },
            'warns'   : warns}

        # Return results
        liner.close()
        stats = ut.stamp_stats(
            stats=stats,
            module='join',
            input_rows=input_rows,
            output_rows=0)
        return (None, stats)

    liner.send(
        ' Verdict: Join Possibly Feasible\n')

    # Write outdf to file
    if not outfile is None:
        ut.write_df_csv(
            df=outdf,
            outfile=outfile,
            sep=',')

    # Build Stats Dictionary
    stats = {
        'status'   : True,
        'basis'    : 'solved',
        'step'     : 3,
        'step_name': 'parsing-oligo-limit',
        'vars'     : {
                      'oligo_limit': oligolimit,
                   'limit_overflow': False,
                    'min_oligo_len': minoligolen,
                    'max_oligo_len': maxoligolen,
                   'overflow_count': 0,
                   'shared_columns': tuple(shared_columns),
                   'columns_joined': inserted_columns,
                  'columns_ignored': ignored_columns,
            'ambiguous_resolutions': ambiguous_resolutions,
                      'join_policy': joinpolicy,
            },
        'warns'    : warns}

    # Joining Statistics
    liner.send('\n[Joining Statistics]\n')
    liner.send(
        '     Joining Status  : Successful\n')
    liner.send(
        '     Columns Joined  : {:,} Column(s)\n'.format(
            len(stats['vars']['columns_joined'])))
    liner.send(
        '     Columns Ignored : {:,} Column(s)\n'.format(
            len(stats['vars']['columns_ignored'])))
    liner.send(
        ' Ambiguities Resolved: {:,} Columns(s)\n'.format(
            stats['vars']['ambiguous_resolutions']))

    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Unschedule outfile deletion
    ae.unregister(ofdeletion)

    # Close Liner
    liner.close()

    # Return Solution and Statistics
    stats = ut.stamp_stats(
        stats=stats,
        module='join',
        input_rows=input_rows,
        output_rows=len(outdf.index) if outdf is not None else 0)
    outdf_return = outdf
    if (outdf is not None) and (not id_from_index):
        outdf_return = ut.get_df_with_id_column(outdf)
    return (outdf_return, stats)
