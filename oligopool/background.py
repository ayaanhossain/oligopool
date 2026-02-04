import time as tt

import atexit as ae

import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp
from .base import core_background as cb


def background(
    input_data:list|str|pd.DataFrame,
    maximum_repeat_length:int,
    output_directory:str,
    verbose:bool=True) -> dict:
    '''
    Build a background k-mer database for screening designed elements against off-target repeat matches.
    Use the resulting directory with `primer`, `barcode`, `motif`, `spacer`, or `verify` modules.

    Required Parameters:
        - `input_data` (`list` / `str` / `pd.DataFrame`): Background sequences to screen against.
        - `maximum_repeat_length` (`int`): Max shared repeat length between designed elements and the background
            reference (between 6 and 20).
        - `output_directory` (`str`): Directory to store the generated background k-mer database.
            A `.oligopool.background` suffix is added if missing.

    Optional Parameters:
        - `verbose` (`bool`): If `True`, logs progress to stdout (default: `True`).

    Returns:
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` can be a list of DNA strings, a CSV/DataFrame with a 'Sequence' column,
            or a FASTA file (.fa/.fasta/.fna, optionally gzipped).
        - If `input_data` is a CSV or DataFrame, it must have a 'Sequence' column with DNA strings.
        - `maximum_repeat_length` here controls non-repetitiveness against `background` only.
        - For advanced manipulation, use `vectorDB` (see `help(oligopool.vectorDB)`).
    '''

    # Argument Aliasing
    indata    = input_data
    maxreplen = maximum_repeat_length
    outdir    = output_directory
    verbose   = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Background Verbiage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Background]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass indata Parsing and Validation
    (background,
    background_valid) = vp.get_parsed_exseqs_info(
        exseqs=indata,
        exseqs_field=' Background Data     ',
        exseqs_desc='Unique Sequence(s)',
        df_field='Sequence',
        required=True,
        liner=liner)
    input_rows = len(background) if background is not None else 0

    # Full maxreplen Validation
    maxreplen_valid = vp.get_numeric_validity(
        numeric=maxreplen,
        numeric_field='    Maximum Repeat   ',
        numeric_pre_desc=' Up to ',
        numeric_post_desc=' Base Pair(s) Background Repeats',
        minval=6,
        maxval=20,
        precheck=False,
        liner=liner)

    # Full outdir Validation
    outdir_valid = vp.get_outdir_validity(
        outdir=outdir,
        outdir_suffix='.oligopool.background',
        outdir_field='     Output Directory',
        liner=liner)

    # First Pass Validation
    if not all([
        background_valid,
        maxreplen_valid,
        outdir_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Adjust Numeric Paramters
    maxreplen = round(maxreplen)

    # Adjust outdir Suffix
    outdir = ut.get_adjusted_path(
        path=outdir,
        suffix='.oligopool.background')

    # Schedule outdir deletion
    oddeletion = ae.register(
        ut.remove_directory,
        outdir)

    # Launching Background Extraction
    liner.send('\n[Computing Background]\n')

    # Define Background Stats
    stats = {
        'status'  : False,
        'basis'   : 'infeasible',
        'step'    : 1,
        'step_name': 'computing-background',
        'vars'    : {
            'kmer_space': 0,  # kmer Space
            'fill_count': 0,  # kmer Fill Count
            'left_count': 0}, # kmer Left Count
        'warns'   : {}}

    # Start Timer
    t0 = tt.time()

    # Extract Background
    stats = cb.background_engine(
        background=background,
        maxreplen=maxreplen,
        outdir=outdir,
        stats=stats,
        liner=liner)

    # Counting Status
    if stats['status']:
        backgroundstatus = 'Successful'
    else:
        backgroundstatus = 'Failed'

    # Background Statistics
    liner.send('\n[Background Statistics]\n')

    plen = ut.get_printlen(
        value=max(stats['vars'][field] for field in (
            'kmer_space',
            'fill_count',
            'left_count')))

    sntn = 'e' if plen > 15 else 'd'

    liner.send(
        ' Background Status: {}\n'.format(
            backgroundstatus))
    liner.send(
        '      k-mer Space : {:{},{}} Unique {:,}-mers\n'.format(
            stats['vars']['kmer_space'],
            plen,
            sntn,
            maxreplen+1))
    liner.send(
        '       Fill Count : {:{},{}} Unique {:,}-mers ({:6.2f} %)\n'.format(
            stats['vars']['fill_count'],
            plen,
            sntn,
            maxreplen+1,
            ut.safediv(
                A=stats['vars']['fill_count'] * 100.,
                B=stats['vars']['kmer_space'])))
    liner.send(
        '       Left Count : {:{},{}} Unique {:,}-mers ({:6.2f} %)\n'.format(
            stats['vars']['left_count'],
            plen,
            sntn,
            maxreplen+1,
            ut.safediv(
                A=stats['vars']['left_count'] * 100.,
                B=stats['vars']['kmer_space'])))

    liner.send(' Time Elapsed: {:.2f} sec\n'.format(tt.time()-t0))

    # Unschedule outfile deletion
    if backgroundstatus == 'Successful':
        ae.unregister(oddeletion)

    # Close Liner
    liner.close()

    # Return Statistics
    stats = ut.stamp_stats(
        stats=stats,
        module='background',
        input_rows=input_rows,
        output_rows=0)
    return stats
