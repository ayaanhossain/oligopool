import time as tt

import atexit as ae

import numpy  as np
import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp

from typing import Tuple


def final(
    input_data:str|pd.DataFrame,
    output_file:str|None=None,
    verbose:bool=True) -> Tuple[pd.DataFrame, dict]:
    '''
    Concatenates all columns in input_data and returns the final library DataFrame.
    If a path is provided, a CSV file containing the library is written out.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated oligopool variants.

    Optional Parameters:
        - `output_file` (`str`): Filename for output DataFrame (default: `None`).
        - `verbose` (`bool`): If `True`, logs updates to stdout (default: `True`).

    Returns:
        - A pandas DataFrame of final library; saves to `output_file` if specified.
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` must contain a unique 'ID' column, all other columns must be non-empty DNA strings.
        - All annotations are lost in this step, so new elements can only be added as left or right context.
    '''

    # Argument Aliasing
    indata  = input_data
    outfile = output_file
    verbose = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Finalization Verbage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Final]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_indata_info(
        indata=indata,
        indata_field='  Input Data',
        required_fields=('ID',),
        precheck=False,
        liner=liner)

    # Full outfile Validation
    outfile_valid = vp.get_outdf_validity(
        outdf=outfile,
        outdf_suffix='.oligopool.final.csv',
        outdf_field=' Output File',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.final.csv')

    # First Pass Validation
    if not all([
        indata_valid,
        outfile_valid]):
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
    liner.send('\n[Step 1: Finalizing Oligopool]\n')

    # Compute Final DataFrame
    outdf = pd.DataFrame(index=indf.index)
    outdf['CompleteOligo'] = ut.get_df_concat(df=indf)
    outdf['OligoLength']   = list(map(len, outdf['CompleteOligo']))

    # Show Update
    liner.send(' Finalization Completed\n')
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Write outdf to file
    if not outfile is None:
        outdf.to_csv(
            path_or_buf=outfile,
            sep=',')

    # Compute Stats
    minoligolen = np.min(outdf.OligoLength)
    maxoligolen = np.max(outdf.OligoLength)

    # Build Stats Dictionary
    stats = {
        'status'  : True,
        'basis'   : 'solved',
        'step'    : 1,
        'step_name': 'finalizing-oligopool',
        'vars'    : {
            'min_oligo_len': minoligolen,
            'max_oligo_len': maxoligolen},
        'warns'   : {}}

    # Finalization Statistics
    liner.send('\n[Finalization Statistics]\n')

    plen = ut.get_printlen(
        value=max(stats['vars'][field] for field in (
            'min_oligo_len',
            'max_oligo_len')))

    liner.send(
        ' Final Status: Successful\n')

    if minoligolen == maxoligolen:
        liner.send(
            ' Oligo Length: {:{},d} Base Pair(s)\n'.format(
                stats['vars']['min_oligo_len'],
                plen))
    else:
        liner.send(
            ' Oligo Length: {:{},d} to {:{},d} Base Pair(s)\n'.format(
                stats['vars']['min_oligo_len'],
                plen,
                stats['vars']['max_oligo_len'],
                plen))

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Unschedule outfile deletion
    ae.unregister(ofdeletion)

    # Close Liner
    liner.close()

    # Return Solution and Statistics
    return (outdf, stats)
