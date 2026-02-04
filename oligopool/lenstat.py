import time as tt

import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp
from .base import core_lenstat as cl


def lenstat(
    input_data:str|pd.DataFrame,
    oligo_length_limit:int,
    verbose:bool=True) -> dict:
    '''
    Compute element-wise and overall oligo length statistics under an `oligo_length_limit`.
    Useful as a "ruler" during design; returns only a stats dictionary.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated oligo pool variants.
        - `oligo_length_limit` (`int`): Maximum allowed oligo length (≥ 4).

    Optional Parameters:
        - `verbose` (`bool`): If `True`, logs progress to stdout (default: `True`).

    Returns:
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` must contain a unique 'ID' column; all other columns must be non-empty DNA strings.
        - `lenstat` does not modify `input_data` and does not require an `output_file`.
        - Use `lenstat` to track free space under `oligo_length_limit` (e.g., before using `spacer`).
        - `lenstat` assumes all non-ID columns are DNA strings; if your DataFrame contains annotation columns
            or degenerate/IUPAC bases, use `verify` for a more general QC summary.
    '''

    # Argument Aliasing
    indata     = input_data
    oligolimit = oligo_length_limit
    verbose    = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Length Statistics Verbiage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Length Statistics]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_indata_info(
        indata=indata,
        indata_field=' Input Data ',
        required_fields=('ID',),
        precheck=False,
        liner=liner)
    input_rows = len(indf.index) if isinstance(indf, pd.DataFrame) else 0

    # Full oligolimit Validation
    oligolimit_valid = vp.get_numeric_validity(
        numeric=oligolimit,
        numeric_field=' Oligo Limit',
        numeric_pre_desc=' At most ',
        numeric_post_desc=' Base Pair(s)',
        minval=4,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # First Pass Validation
    if not all([
        indata_valid,
        oligolimit_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Paramters
    oligolimit = round(oligolimit)

    # Length Stats Book-keeping
    stats = None
    warns = {}

    # Launching Length Stats Computation
    liner.send('\n[Step 1: Computing Length Statistics]\n')

    # Compute Length Stats
    (intstats,
    minspaceavail,
    maxspaceavail) = cl.lenstat_engine(
        indf=indf,
        oligolimit=oligolimit,
        liner=liner)

    # Length Statistics Results
    liner.send('\n[Length Statistics]\n')

    # Build Stats Print String
    statsprint = ut.get_lenstat_statsprint(
        intstats=intstats)

    # Show Stats String
    if verbose and statsprint:
        print('\n{}'.format(statsprint))

    # Restructure Length Stats Dictionary
    stats = {
        'status'  : True,
        'basis'   : 'solved',
        'step'    : 1,
        'step_name': 'computing-length-statistics',
        'vars'    : {
            'len_stat'     : ut.get_lenstat_dict(
                intstats=intstats),             # Store Element-wise Length Stats
            'oligo_limit'    : oligolimit,            # Specified Oligo Limit
            'min_space_avail': minspaceavail,         # Minimum Space Available
            'max_space_avail': maxspaceavail},        # Maximum Space Available
        'warns'   : warns}

    # Show Free Space Available
    if minspaceavail == maxspaceavail:
        liner.send(
            '\n Free Space Available: {:,} Base Pair(s)\n'.format(
                minspaceavail))
    else:
        liner.send(
            '\n Free Space Available: {:,} to {:,} Base Pair(s)\n'.format(
                minspaceavail,
                maxspaceavail))

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Close Liner
    liner.close()

    # Return Results
    stats = ut.stamp_stats(
        stats=stats,
        module='lenstat',
        input_rows=input_rows,
        output_rows=0)
    return stats
