import time as tt

import collections as cx

import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp
from .base import core_lenstat as cl


def lenstat(
    input_data:str|pd.DataFrame,
    oligo_length_limit:int,
    verbose:bool=True) -> dict:
    '''
    Computes the length statistics of the elements and the resulting oligos in `input_data`.
    Returns only the stats dictionary.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated oligopool variants.
        - `oligo_length_limit` (`int`): Maximum allowed oligo length (≥ 4).

    Optional Parameters:
        - `verbose` (`bool`): If `True`, logs updates to stdout (default: `True`).

    Returns:
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` must contain a unique 'ID' column, all other columns must be non-empty DNA strings.
        - `lenstat` module does not add anything to `input_data`, so there is no `output_file` required.
    '''

    # Argument Aliasing
    indata     = input_data
    oligolimit = oligo_length_limit
    verbose    = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Length Statistics Verbage Print
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
    outdf = None
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
    statsprint = '\n'.join(
        ' ' + line for line in pd.DataFrame.from_dict(
            {k: [u[0]] + [str(v) + ' bp' for v in u[1:-1]] + [u[-1]] \
                for k,u in intstats.items()},
            orient='index',
            columns=pd.MultiIndex.from_arrays(
                [('', '    Min', '    Max', '   Min', '   Max', '   Oligo'),
                 ('', 'Element', 'Element', ' Oligo', ' Oligo', '   Limit'),
                 ('', ' Length', ' Length', 'Length', 'Length', 'Overflow')]
                )).to_string(index=False).split('\n'))

    # Show Stats String
    print('\n{}'.format(
        statsprint))

    # Restructure Length Stats Dictionary
    stats = {
        'status'  : True,
        'basis'   : 'solved',
        'step'    : 1,
        'step_name': 'computing-length-statistics',
        'vars'    : {
            'len_stat'     : cx.OrderedDict({v[0]: { # Store Element-wise Length Stats
                'min_element_len':  v[1],
                'max_element_len':  v[2],
                'min_oligo_len':    v[3],
                'max_oligo_len':    v[4],
                'limit_overflow':   v[5] == 'Yes'} for k,v in intstats.items()}),
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
    return stats
