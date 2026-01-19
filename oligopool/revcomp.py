import time as tt

import atexit as ae

import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp

from typing import Tuple


def revcomp(
    input_data:str|pd.DataFrame,
    output_file:str|None=None,
    left_context_column:str|None=None,
    right_context_column:str|None=None,
    verbose:bool=True) -> Tuple[pd.DataFrame, dict]:
    '''
    Flip the orientation of a column range by reverse-complementing sequences and reversing column order.
    Useful for mid-pipeline maneuvers where a region must be switched between readout and synthesis orientation.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated oligopool variants.

    Optional Parameters:
        - `output_file` (`str`): Filename for output DataFrame; required in CLI usage,
            optional in library usage (default: `None`).
        - `left_context_column` (`str`): Column for left DNA context (default: `None`).
        - `right_context_column` (`str`): Column for right DNA context (default: `None`).
        - `verbose` (`bool`): If `True`, logs updates to stdout (default: `True`).

    Returns:
        - A pandas DataFrame of revcomped elements; saves to `output_file` if specified.
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` must contain a unique 'ID' column, all other columns must be non-empty DNA strings.
        - `revcomp` module does not require `left_context_column` and `right_context_column` to be adjacent.
        - If `left_context_column` is unspecified, then the first column is considered.
        - Similarly, the last column is considered as the right context column, if unspecified.
        - Useful for inspecting or correcting split-fragment orientation (see `split` output ordering).
        - Useful for mid-pipeline maneuvers where a region is designed in one orientation (assay/readout)
          but must be installed/synthesized in the opposite orientation.
    '''

    # Preserve return style when the caller intentionally used ID as index.
    id_from_index = ut.get_id_index_intent(input_data)

    # Argument Aliasing
    indata       = input_data
    leftcontext  = left_context_column
    rightcontext = right_context_column
    outfile      = output_file
    verbose      = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Reverse Complementing Verbiage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Revcomp]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_indata_info(
        indata=indata,
        indata_field='  Input Data   ',
        required_fields=('ID',),
        precheck=False,
        liner=liner)
    input_rows = len(indf.index) if isinstance(indf, pd.DataFrame) else 0

    # Full outfile Validation
    outfile_valid = vp.get_outdf_validity(
        outdf=outfile,
        outdf_suffix='.oligopool.revcomp.csv',
        outdf_field=' Output File   ',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.revcomp.csv')

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full leftcontext Parsing and Validation
    (_,
    leftcontext_valid) = vp.get_parsed_column_info(
        col=leftcontext,
        df=indf,
        col_field='   Left Context',
        col_desc='Input from Column',
        col_type=0,
        adjcol=None,
        adjval=None,
        iscontext=True,
        typecontext=0,
        liner=liner)

    # Full leftcontext Parsing and Validation
    (_,
    rightcontext_valid) = vp.get_parsed_column_info(
        col=rightcontext,
        df=indf,
        col_field='  Right Context',
        col_desc='Input from Column',
        col_type=0,
        adjcol=None,
        adjval=None,
        iscontext=True,
        typecontext=1,
        liner=liner)

    # First Pass Validation
    if not all([
        indata_valid,
        outfile_valid,
        leftcontext_valid,
        rightcontext_valid,
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
    liner.send('\n[Step 1: Revcomping Elements]\n')

    # Compute New DataFrame
    outdf = indf.copy()

    # Decide on Columns
    columns = outdf.columns.to_list()
    if not left_context_column:
        left_context_column = columns[0]
    if not right_context_column:
        right_context_column = columns[-1]
    min_idx, max_idx = sorted([
        columns.index(left_context_column),
        columns.index(right_context_column)])
    left_context_column  = columns[min_idx]
    right_context_column = columns[max_idx]
    left_chunk = columns[:min_idx]
    elementsrevcomped = columns[min_idx:max_idx+1]
    right_chunk = columns[max_idx+1:]

    # Show Update
    liner.send(f' Starting Column: {elementsrevcomped[+0]}\n')
    liner.send(f'   Ending Column: {elementsrevcomped[-1]}\n')

    # Revcomp Columns
    outdf = outdf[left_chunk + elementsrevcomped[::-1] + right_chunk]
    for col in elementsrevcomped:
        outdf[col] = outdf[col].transform(lambda x: ut.get_revcomp(seq=x))

    # Show Update
    liner.send(' Merging Completed\n')
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
        'step_name': 'revcomping-elements',
        'vars'    : {
            'elements_revcomped': elementsrevcomped,
            },
        'warns'   : {}}

    # Revcomp Statistics
    liner.send('\n[Revcomp Statistics]\n')
    liner.send(
        '  Revcomp Status   : Successful\n')
    liner.send(
        f' Elements Revcomped: {len(elementsrevcomped)}\n')
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
        module='revcomp',
        input_rows=input_rows,
        output_rows=len(outdf.index) if outdf is not None else 0)
    outdf_return = outdf
    if (outdf is not None) and (not id_from_index):
        outdf_return = ut.get_df_with_id_column(outdf)
    return (outdf_return, stats)
