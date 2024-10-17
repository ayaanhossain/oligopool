import time  as tt

import atexit as ae

import pandas as pd

from .base import utils as ut
from .base import validation_parsing  as vp
from .base import core_split as cs

from typing import Tuple


def split(
    input_data:str|pd.DataFrame,
    split_length_limit:int,
    minimum_melting_temperature:float,
    minimum_hamming_distance:int,
    minimum_overlap_length:int,
    maximum_overlap_length:int,
    output_file:str|None=None,
    verbose:bool=True) -> Tuple[pd.DataFrame, dict]:
    '''
    Splits longer oligos into shorter overlapping fragments with all overlaps at particular
    coordinates having specified minimum pairwise Hamming distances and minimum melting temperatures.
    Overlap lengths can be controlled depending on downstream assembly strategy. Produces a derived
    new DataFrame with 'Split' columns which can be saved to a specified CSV file.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated oligopool variants.
        - `split_length_limit` (`int`): Maximum allowed length for split oligos (≥ 4).
        - `minimum_melting_temperature` (`float`): Minimum overlap region Tm (≥ 4°C).
        - `minimum_hamming_distance` (`int`): Minimum overlap region pairwise Hamming distance (≥ 1).
        - `minimum_overlap_length` (`int`): Minimum overlap region length (≥ 15).
        - `maximum_overlap_length` (`int`): Maximum overlap region length (≥ 15).

    Optional Parameters:
        - `output_file` (`str`): Filename for output DataFrame (default: `None`).
        - `verbose` (`bool`): If `True`, logs updates to stdout (default: `True`).

    Returns:
        - A pandas DataFrame with split oligos; saves to `output_file` if specified.
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` must contain a unique 'ID' column, all other columns must be non-empty DNA strings.
        - `minimum_overlap_length` should always be larger than `minimum_hamming_distance`.
        - Total number of fragments is auto determined, and can be variable per oligo, depending on length.
        - Returned DataFrame contains split oligos, annotation from `input_data` is lost.
    '''

    # Argument Aliasing
    indata     = input_data
    splitlimit = split_length_limit
    mintmelt   = minimum_melting_temperature
    minhdist   = minimum_hamming_distance
    minoverlap = minimum_overlap_length
    maxoverlap = maximum_overlap_length
    outfile    = output_file
    verbose    = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Splitting Verbage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Split]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_indata_info(
        indata=indata,
        indata_field='   Input Data       ',
        required_fields=('ID',),
        precheck=False,
        liner=liner)

    # Full splitlimit Validation
    splitlimit_valid = vp.get_numeric_validity(
        numeric=splitlimit,
        numeric_field='   Split Limit      ',
        numeric_pre_desc=' Split Fragments at most ',
        numeric_post_desc=' Base Pair(s) Each',
        minval=4,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full mintmelt Validation
    tmelt_valid = vp.get_numeric_validity(
        numeric=mintmelt,
        numeric_field=' Melting Temperature',
        numeric_pre_desc=' At least ',
        numeric_post_desc=' °C b/w On-Target Overlaps',
        minval=4,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full minhdist Validation
    minhdist_valid = vp.get_numeric_validity(
        numeric=minhdist,
        numeric_field=' Hamming Distance   ',
        numeric_pre_desc=' At least ',
        numeric_post_desc=' Mismatch(es) per Off-Target Overlap Pair ',
        minval=1,
        maxval=splitlimit if splitlimit_valid else float('inf'),
        precheck=False,
        liner=liner)

    # Full minoverlap and maxoverlap Validation
    (minoverlap,
    maxoverlap,
    overlap_valid) = vp.get_parsed_range_info(
        minval=minoverlap,
        maxval=maxoverlap,
        range_field=' Overlap Length     ',
        range_unit='Base Pair(s) Fragment Overlap(s)',
        range_min=15 if not minhdist_valid else max(15, minhdist),
        range_max=splitlimit if splitlimit_valid else float('inf'),
        liner=liner)

    # Full outfile Validation
    outfile_valid = vp.get_outdf_validity(
        outdf=outfile,
        outdf_suffix='.oligopool.split.csv',
        outdf_field='  Output File       ',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.split.csv')

    # First Pass Validation
    if not all([
        indata_valid,
        splitlimit_valid,
        tmelt_valid,
        minhdist_valid,
        overlap_valid,
        outfile_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Paramters
    splitlimit = round(splitlimit)
    minhdist   = round(minhdist)
    minoverlap = round(minoverlap)
    maxoverlap = round(maxoverlap)

    # Primer Design Book-keeping
    outdf = None
    stats = None
    warns = {}

    # Parse Oligopool Split Limit Feasibility
    liner.send('\n[Step 1: Parsing Split Limit]\n')

    # Parse splitlimit
    (parsestatus,
    seqlist,
    oligounderflow,
    unevensplit,
    minvariantlen,
    maxvariantlen,
    minsplitcount,
    maxsplitcount) = cs.get_parsed_splitlimit(
        indf=indf,
        splitlimit=splitlimit,
        liner=liner)

    # splitlimit infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 1,
            'step_name': 'parsing-split-limit',
            'vars'    : {
                    'split_limit': splitlimit,
                'oligo_underflow': oligounderflow,
                   'uneven_split': unevensplit,
                'min_variant_len': minvariantlen,
                'max_variant_len': maxvariantlen,
                'min_split_count': minsplitcount,
                'max_split_count': maxsplitcount,},
            'warns'   : warns}

        # Return results
        liner.close()
        return (outdf, stats)

    # Define spanlen
    # Note: spanlen is the minimum span
    #       for a split region
    spanlen = max(minhdist, minoverlap)

    # Compute Sequence Matrix
    liner.send('\n[Step 2: Computing Sequence Matrix]\n')

    # Compute padvec and seqmat
    (padvec,
    seqmat) = cs.get_seqmat_padvec(
        seqlist=seqlist,
        maxvariantlen=maxvariantlen,
        liner=liner)

    # Compute Sequence Matrix
    liner.send('\n[Step 3: Computing Entropy Vector]\n')

    # Compute padvec and seqmat
    entvec = cs.get_entvec(
        seqmat=seqmat,
        maxvariantlen=maxvariantlen,
        liner=liner)

    # Parse Oligopool Limit Feasibility
    liner.send('\n[Step 4: Parsing Variable Contigs]\n')

    # Parse splitlimit
    (parsestatus,
    varcont,
    varcontcount,
    mergedcontcount,
    filtercontcount) = cs.get_varcont(
        entvec=entvec,
        minhdist=minhdist,
        spanlen=spanlen,
        liner=liner)

    # splitlimit infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 4,
            'step_name': 'parsing-variable-contig',
            'vars'    : {
                'variable_contig_count': varcontcount,
                  'merged_contig_count': mergedcontcount,
                  'filter_contig_count': filtercontcount},
            'warns'   : warns}

        # Return results
        liner.close()
        return (outdf, stats)

    # Launching Split Design
    liner.send('\n[Step 5: Computing Split]\n')

    # Define Split Design Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 8,
        'step_name': 'computing-split',
        'vars'    : {
                      'num_splits': 0,   # Total Number of Splits
                      'split_lens': [],  # Split Oligo Lengths
                    'overlap_lens': [],  # Split Overlap Lengths
                  'mean_Tm_distro': [],  # Mean Tm of Each Split
            'mean_distance_distro': [],  # Mean HDist of Each Split
              'infeasible_contigs': False,  # Infeasible Contigs Flag
                   'uneven_splits': False}, # Uneven Splits Flag
        'warns'   : warns}

    # Schedule outfile deletion
    ofdeletion = ae.register(
        ut.remove_file,
        outfile)

    # Design Split
    (split,
    overlap,
    stats) = cs.split_engine(
        seqlist=seqlist,
        splitlimit=splitlimit,
        mintmelt=mintmelt,
        minhdist=minhdist,
        maxoverlap=maxoverlap,
        minvariantlen=minvariantlen,
        maxvariantlen=maxvariantlen,
        spanlen=spanlen,
        seqmat=seqmat,
        varcont=varcont,
        stats=stats,
        liner=liner)

    # Success Relevant Stats
    if stats['status']:

        # Launching Stats Aggregation
        liner.send('\n[Step 6: Aggregating Stats]\n')

        # Compute Tm and HDist Distribution
        # and Finalize Split Sequences
        (splitstore,
        stats) = cs.aggregate_stats(
            seqlist=seqlist,
            seqmat=seqmat,
            split=split,
            overlap=overlap,
            stats=stats,
            liner=liner)

    # Split Status
    if stats['status']:
        splitstatus = 'Successful'
    else:
        splitstatus = 'Failed'

    # Insert split into outdf
    if stats['status']:

        # Prepare outdf
        outdf = pd.DataFrame(
            index=indf.index)

        # Insert Splits into Columns
        for sidx in range(len(split)):
            outdf['Split{}'.format(sidx+1)] = splitstore[sidx]

        # Write outdf to file
        if not outfile is None:
            outdf.to_csv(
                path_or_buf=outfile,
                sep=',')

    # Split Design Statistics
    liner.send('\n[Split Design Statistics]\n')

    liner.send(
        '     Design Status  : {}\n'.format(
            splitstatus))

    # Success Relevant Stats
    if stats['status']:

        maxval = max(max(min(stats['vars'][field]) for field in (
            'split_lens',
            'overlap_lens',
            'mean_Tm_distro',
            'mean_distance_distro')),
            stats['vars']['num_splits'])

        sntn, plen = ut.get_notelen(
            printlen=ut.get_printlen(
                value=maxval))

        liner.send(
            '     No. of Splits  : {:{},d} Fragments per Variant\n'.format(
                stats['vars']['num_splits'],
                plen))
        liner.send(
            '      Split Length  : {:{},d} {}Base Pair(s)\n'.format(
                min(stats['vars']['split_lens']),
                plen,
                ['', 'to {:,d} '.format(max(stats['vars']['split_lens']))][
                    max(stats['vars']['split_lens']) != min(stats['vars']['split_lens'])
                ]))
        liner.send(
            '    Overlap Length  : {:{},d} {}Base Pair(s)\n'.format(
                min(stats['vars']['overlap_lens']),
                plen,
                ['', 'to {:,d} '.format(max(stats['vars']['overlap_lens']))][
                    max(stats['vars']['overlap_lens']) != min(stats['vars']['overlap_lens'])
                ]))
        liner.send(
            '    Overlap Tm      : {:{},d} {}°C\n'.format(
                min(stats['vars']['mean_Tm_distro']),
                plen,
                ['', 'to {:,d} '.format(max(stats['vars']['mean_Tm_distro']))][
                    max(stats['vars']['mean_Tm_distro']) != min(stats['vars']['mean_Tm_distro'])
                ]))
        liner.send(
            '    Overlap Distance: {:{},d} {}Mismatch(es)\n'.format(
               min(stats['vars']['mean_distance_distro']),
                plen,
                ['', 'to {:,d} '.format(max(stats['vars']['mean_distance_distro']))][
                    max(stats['vars']['mean_distance_distro']) != min(stats['vars']['mean_distance_distro'])
                ]))

    # Failure Relavant Stats
    else:
        liner.send(
            ' Infeasible Contigs : {}\n'.format(
                ['No', 'Yes'][stats['vars']['infeasible_contigs']]))
        liner.send(
            '      Unven Split   : {}\n'.format(
                ['No', 'Yes'][stats['vars']['uneven_splits']]))

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Unschedule outfile deletion
    if splitstatus == 'Successful':
        ae.unregister(ofdeletion)

    # Close Liner
    liner.close()

    # Return Solution and Statistics
    return (outdf, stats)
