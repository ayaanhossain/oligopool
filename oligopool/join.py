import time as tt

import atexit as ae

import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp

from typing import Tuple


def _resolve_insert_index(
    column,
    source_columns,
    target_columns,
    join_policy):
    '''
    Resolve insertion index for a new column from source into target.
    Internal use only.
    '''
    src_idx = source_columns.index(column)

    left_anchor = None
    for idx in range(src_idx - 1, -1, -1):
        candidate = source_columns[idx]
        if candidate in target_columns:
            left_anchor = candidate
            break

    right_anchor = None
    for idx in range(src_idx + 1, len(source_columns)):
        candidate = source_columns[idx]
        if candidate in target_columns:
            right_anchor = candidate
            break

    left_idx = 0 if left_anchor is None else target_columns.index(left_anchor) + 1
    right_idx = len(target_columns) if right_anchor is None else target_columns.index(right_anchor)

    if left_idx > right_idx:
        raise RuntimeError(
            "Unable to place column '{}' due to inconsistent anchor ordering.".format(
                column))

    if left_idx == right_idx:
        return (left_idx, False, left_anchor, right_anchor)

    if join_policy == 0:
        return (left_idx, True, left_anchor, right_anchor)
    return (right_idx, True, left_anchor, right_anchor)

def join(
    input_data:str|pd.DataFrame,
    other_data:str|pd.DataFrame,
    join_policy:int|str,
    output_file:str|None=None,
    verbose:bool=True) -> Tuple[pd.DataFrame, dict]:
    '''
    Join two oligo tables on ID and add new columns from the second table into the first.
    Useful for recombining branch outputs into a single design table.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated oligo pool variants.
        - `other_data` (`str` / `pd.DataFrame`): Path to a second CSV file or DataFrame with the same IDs as `input_data`.
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
        - `input_data` and `other_data` must contain exactly the same IDs (order may differ).
        - Overlapping non-ID column names from `other_data` are ignored (backbone from `input_data` is preserved).
        - Only non-overlapping columns from `other_data` are inserted into `input_data` order.
    '''

    # Preserve return style when the caller intentionally used ID as index.
    id_from_index = ut.get_id_index_intent(input_data)

    # Argument Aliasing
    indata     = input_data
    otherdata  = other_data
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
        other_ids_match = len(indf.index) == len(otherdf.index) and \
                          sorted(indf.index) == sorted(otherdf.index)
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

    # Show Update
    liner.send('\n[Step 1: Computing Join]\n')

    # Compute New DataFrame
    outdf = indf.copy()
    output_columns = outdf.columns.to_list()
    source_columns = otherdf.columns.to_list()
    liner.send(
        '  Input Column(s): {:,} Column(s)\n'.format(
            len(output_columns)))
    for col in output_columns:
        liner.send('    - {}\n'.format(col))

    liner.send(
        '  Other Column(s): {:,} Column(s)\n'.format(
            len(source_columns)))
    for col in source_columns:
        liner.send('    - {}\n'.format(col))

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
        _) = _resolve_insert_index(
            column=col,
            source_columns=source_columns,
            target_columns=output_columns,
            join_policy=joinpolicy)

        outdf.insert(insert_idx, col, otherdf[col])
        output_columns.insert(insert_idx, col)
        inserted_columns.append(col)
        if was_ambiguous:
            ambiguous_resolutions += 1

    liner.send(
        ' Output Column(s): {:,} Column(s)\n'.format(
            len(output_columns)))
    for col in output_columns:
        liner.send('    - {}\n'.format(col))

    # Show Update
    liner.send(' Join Completed\n')
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Write outdf to file
    if not outfile is None:
        ut.write_df_csv(
            df=outdf,
            outfile=outfile,
            sep=',')

    # Build Stats Dictionary
    stats = {
        'status'  : True,
        'basis'   : 'solved',
        'step'    : 1,
        'step_name': 'computing-join',
        'vars'    : {
            'columns_joined': inserted_columns,
            'columns_ignored': ignored_columns,
            'ambiguous_resolutions': ambiguous_resolutions,
            'join_policy': joinpolicy,
            },
        'warns'   : {}}

    # Joining Statistics
    liner.send('\n[Joining Statistics]\n')
    liner.send(
        f'        Join Status  : Successful\n')
    liner.send(
        f'     Columns Joined  : {len(inserted_columns)}\n')
    liner.send(
        f'     Columns Ignored : {len(ignored_columns)}\n')
    liner.send(
        f' Ambiguities Resolved: {ambiguous_resolutions}\n')
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
