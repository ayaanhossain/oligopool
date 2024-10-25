import time as tt

import atexit as ae

import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp

from typing import Tuple


def merge(
    input_data:str|pd.DataFrame,
    merge_column:str,
    output_file:str|None=None,
    left_context_column:str|None=None,
    right_context_column:str|None=None,
    verbose:bool=True) -> Tuple[pd.DataFrame, dict]:
    '''
    Merges all elements from `left_context` to `right_context` into `merge_column` and removes the
    original ones from the DataFrame. The modified DataFrame can be saved to a CSV file.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated oligopool variants.
        - `merge_column` (`str`): Column name for the merged DNA.

    Optional Parameters:
        - `output_file` (`str`): Filename for output DataFrame (default: `None`).
        - `left_context_column` (`str`): Column for left DNA context (default: `None`).
        - `right_context_column` (`str`): Column for right DNA context (default: `None`).
        - `verbose` (`bool`): If `True`, logs updates to stdout (default: `True`).

    Returns:
        - A pandas DataFrame of merged elements; saves to `output_file` if specified.
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` must contain a unique 'ID' column, all other columns must be non-empty DNA strings.
        - `merge` module does not require the `left_context` or `right_context` to be adjacent.
        - If the left context column is unspecified, then the first column is considered.
        - Similarly, the the last column is considered as the right context column, if unspecified.
    '''

    # Argument Aliasing
    indata       = input_data
    leftcontext  = left_context_column
    rightcontext = right_context_column
    mergecol     = merge_column
    outfile      = output_file
    verbose      = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Merging Verbage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Merge]\n')

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

    # Full outcol Validation
    mergecol_valid = vp.get_parsed_column_info(
        col=mergecol,
        df=indf,
        col_field=' Merged Column ',
        col_desc='Output in Column',
        col_type=1,
        adjcol=None,
        adjval=None,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full outfile Validation
    outfile_valid = vp.get_outdf_validity(
        outdf=outfile,
        outdf_suffix='.oligopool.merge.csv',
        outdf_field=' Output File   ',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.merge.csv')

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
        mergecol_valid,
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
    liner.send('\n[Step 1: Merging Elements]\n')

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
    elementsmerged = columns[min_idx:max_idx+1]

    # Show Update
    liner.send(f' Starting Column: {elementsmerged[+0]}\n')
    liner.send(f'   Ending Column: {elementsmerged[-1]}\n')

    # Build Merged Column
    outdf.insert(min_idx, mergecol, outdf[elementsmerged].astype(str).sum(axis=1))

    # Drop Merged Columns
    outdf = outdf.drop(columns=elementsmerged)

    # Show Update
    liner.send(' Merging Completed\n')
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Write outdf to file
    if not outfile is None:
        outdf.to_csv(
            path_or_buf=outfile,
            sep=',')

    # Build Stats Dictionary
    stats = {
        'status'  : True,
        'basis'   : 'solved',
        'step'    : 1,
        'step_name': 'merging-elements',
        'vars'    : {
            'elements_merged': elementsmerged,
            },
        'warns'   : {}}

    # Merging Statistics
    liner.send('\n[Merging Statistics]\n')
    liner.send(
        '    Merge Status: Successful\n')
    liner.send(
        f' Elements Merged: {len(elementsmerged)}\n')
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Unschedule outfile deletion
    ae.unregister(ofdeletion)

    # Close Liner
    liner.close()

    # Return Solution and Statistics
    return (outdf, stats)
