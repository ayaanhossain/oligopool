import time  as tt

import collections as cx
import atexit as ae

import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp
from .base import core_barcode as cb

from typing import Tuple


def barcode(
    input_data:str|pd.DataFrame,
    oligo_length_limit:int,
    barcode_length:int,
    minimum_hamming_distance:int,
    maximum_repeat_length:int,
    barcode_column:str,
    output_file:str|None=None,
    barcode_type:int=0,
    left_context_column:str|None=None,
    right_context_column:str|None=None,
    excluded_motifs:list|str|pd.DataFrame|None=None,
    verbose:bool=True) -> Tuple[pd.DataFrame, dict]:
    '''
    Generates constrained barcodes, ensuring a minimum Hamming distance between each pair
    and excluding specified motifs, even when flanked by context sequences. The output is a
    DataFrame of designed barcodes, which can be saved as a CSV file if specified.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated oligopool variants.
        - `oligo_length_limit` (`int`): Maximum allowed oligo length (≥ 4).
        - `barcode_length` (`int`): Length of the designed barcodes (≥ 4).
        - `minimum_hamming_distance` (`int`): Minimum pairwise Hamming distance (≥ 1).
        - `maximum_repeat_length` (`int`): Max shared repeat length with oligos (≥ 4).
        - `barcode_column` (`str`): Column name for the designed barcodes.

    Optional Parameters:
        - `output_file` (`str`): Filename for output DataFrame (default: `None`).
        - `barcode_type` (`int`): Barcode design type
            0 for fast terminus optimized,
            1 for slow spectrum optimized.
            (default: 0)
        - `left_context_column` (`str`): Column for left DNA context (default: `None`).
        - `right_context_column` (`str`): Column for right DNA context (default: `None`).
        - `excluded_motifs` (`list` / `str` / `pd.DataFrame`): Motifs to exclude;
            can be a CSV path or DataFrame (default: `None`).
        - `verbose` (`bool`): If `True`, logs updates to stdout (default: `True`).

    Returns:
        - A pandas DataFrame of generated barcodes; saves to `output_file` if specified.
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` must contain a unique 'ID' column, all other columns must be non-empty DNA strings.
        - Column names in `input_data` must be unique, and exclude `barcode_column`.
        - At least one of `left_context_column` or `right_context_column` must be specified.
        - If `excluded_motifs` is a CSV or DataFrame, it must have 'ID' and 'Exmotif' columns.
        - If barcode design is challenging, consider
            * altering `barcode_length`, or
            * reducing `minimum_hamming_distance`, or
            * switching to terminus optimized barcodes, or
            * increasing `maximum_repeat_length`, or
            * reducing `excluded_motifs` to relax the constraints.
    '''

    # Argument Aliasing
    indata       = input_data
    oligolimit   = oligo_length_limit
    barcodelen   = barcode_length
    minhdist     = minimum_hamming_distance
    maxreplen    = maximum_repeat_length
    barcodecol   = barcode_column
    outfile      = output_file
    barcodetype  = barcode_type
    leftcontext  = left_context_column
    rightcontext = right_context_column
    exmotifs     = excluded_motifs
    verbose      = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Barcoding Verbage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Barcode]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_indata_info(
        indata=indata,
        indata_field='    Input Data    ',
        required_fields=('ID',),
        precheck=False,
        liner=liner)

    # Full oligolimit Validation
    oligolimit_valid = vp.get_numeric_validity(
        numeric=oligolimit,
        numeric_field='    Oligo Limit   ',
        numeric_pre_desc=' At most ',
        numeric_post_desc=' Base Pair(s)',
        minval=4,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full barcodelen Validation
    barcodelen_valid = vp.get_numeric_validity(
        numeric=barcodelen,
        numeric_field='  Barcode Length  ',
        numeric_pre_desc=' Exactly ',
        numeric_post_desc=' Base Pair(s)',
        minval=4,
        maxval=float('inf') if not oligolimit_valid else oligolimit,
        precheck=False,
        liner=liner)

    # Full minhdist Validation
    minhdist_valid = vp.get_numeric_validity(
        numeric=minhdist,
        numeric_field='  Hamming Distance',
        numeric_pre_desc=' At least ',
        numeric_post_desc=' Mismatch(es) per Barcode Pair',
        minval=1,
        maxval=barcodelen if barcodelen_valid else float('inf'),
        precheck=False,
        liner=liner)

    # Full maxreplen Validation
    maxreplen_valid = vp.get_numeric_validity(
        numeric=maxreplen,
        numeric_field='   Repeat Length  ',
        numeric_pre_desc=' Up to ',
        numeric_post_desc=' Base Pair(s) Oligopool Repeats',
        minval=4,
        maxval=barcodelen if barcodelen_valid else float('inf'),
        precheck=False,
        liner=liner)

    # Full outcol Validation
    barcodecol_valid = vp.get_parsed_column_info(
        col=barcodecol,
        df=indf,
        col_field='  Barcode Column  ',
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
        outdf_suffix='.oligopool.barcode.csv',
        outdf_field='   Output File    ',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.barcode.csv')

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full barcodetype Validation
    barcodetype_valid = vp.get_categorical_validity(
        category=barcodetype,
        category_field='  Barcode Type    ',
        category_pre_desc=' ',
        category_post_desc=' Barcodes',
        category_dict={
            0: 'Terminus Optimized',
            1: 'Spectrum Optimized'},
        liner=liner)

    # Store Context Names
    leftcontextname  = leftcontext
    rightcontextname = rightcontext

    # Full leftcontext Parsing and Validation
    (leftcontext,
    leftcontext_valid) = vp.get_parsed_column_info(
        col=leftcontext,
        df=indf,
        col_field='     Left Context ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=rightcontextname,
        adjval=+1,
        iscontext=True,
        typecontext=0,
        liner=liner)

    # Full leftcontext Parsing and Validation
    (rightcontext,
    rightcontext_valid) = vp.get_parsed_column_info(
        col=rightcontext,
        df=indf,
        col_field='    Right Context ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=leftcontextname,
        adjval=-1,
        iscontext=True,
        typecontext=1,
        liner=liner)

    # Full exmotifs Parsing and Validation
    (exmotifs,
    exmotifs_valid) = vp.get_parsed_exseqs_info(
        exseqs=exmotifs,
        exseqs_field=' Excluded Motifs  ',
        exseqs_desc='Unique Motif(s)',
        df_field='Exmotif',
        required=False,
        liner=liner)

    # First Pass Validation
    if not all([
        indata_valid,
        oligolimit_valid,
        barcodelen_valid,
        minhdist_valid,
        maxreplen_valid,
        barcodecol_valid,
        outfile_valid,
        barcodetype_valid,
        leftcontext_valid,
        rightcontext_valid,
        exmotifs_valid,]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Paramters
    oligolimit = round(oligolimit)
    barcodelen = round(barcodelen)
    minhdist   = round(minhdist)
    maxreplen  = round(maxreplen)

    # Define Edge Effect Length
    edgeeffectlength = None

    # Barcode Design Book-keeping
    has_context = False
    outdf = None
    stats = None
    warns = {}

    # Parse Oligopool Limit Feasibility
    liner.send('\n[Step 1: Parsing Oligo Limit]\n')

    # Parse oligolimit
    (parsestatus,
    minvariantlen,
    maxvariantlen,
    minelementlen,
    maxelementlen,
    minspaceavail,
    maxspaceavail) = ut.get_parsed_oligolimit(
        indf=indf,
        variantlens=None,
        oligolimit=oligolimit,
        minelementlen=barcodelen,
        maxelementlen=barcodelen,
        element='Barcode',
        liner=liner)

    # oligolimit infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 1,
            'step_name': 'parsing-oligo-limit',
            'vars'    : {
                   'oligo_limit': oligolimit,
                'limit_overflow': True,
                'min_variant_len': minvariantlen,
                'max_variant_len': maxvariantlen,
                'min_element_len': minelementlen,
                'max_element_len': maxelementlen,
                'min_space_avail': minspaceavail,
                'max_space_avail': maxspaceavail},
            'warns'   : warns}

        # Return results
        liner.close()
        return (outdf, stats)

    # Parse Barcode Length Feasibility
    liner.send('\n[Step 2: Parsing Barcode Length]\n')

    # Parse barcodelen
    (parsestatus,
    designspace,
    targetcount) = cb.get_parsed_barcode_length(
        barcodelen=barcodelen,
        indf=indf,
        liner=liner)

    # barcodelen infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 2,
            'step_name': 'parsing-barcode-length',
            'vars'    : {
                 'barcode_len': barcodelen,
                'design_space': designspace,
                'target_count': targetcount},
            'warns'   : warns}

        # Return results
        liner.close()
        return (outdf, stats)

    # Parse Excluded Motifs
    if not exmotifs is None:

        # Show update
        liner.send('\n[Step 3: Parsing Excluded Motifs]\n')

        # Update Step 3 Warning
        warns[3] = {
            'warn_count': 0,
            'step_name' : 'parsing-excluded-motifs',
            'vars': None}

        # Parse exmotifs
        (parsestatus,
        exmotifs,
        problens,
        _,
        __) = ut.get_parsed_exmotifs(
            exmotifs=exmotifs,
            typer=tuple,
            element='Barcode',
            leftcontext=leftcontext,
            rightcontext=rightcontext,
            warn=warns[3],
            liner=liner)

        # Remove Step 3 Warning
        if not warns[3]['warn_count']:
            warns.pop(3)

        # exmotifs infeasible
        if not parsestatus:

            # Prepare stats
            stats = {
                'status'  : False,
                'basis'   : 'infeasible',
                'step'    : 3,
                'step_name': 'parsing-excluded-motifs',
                'vars'    : {
                     'prob_lens': problens,
                    'prob_count': tuple(list(
                        4**pl for pl in problens))},
                'warns'   : warns}

            # Return results
            liner.close()
            return (outdf, stats)

        # Update Edge-Effect Length
        edgeeffectlength = ut.get_edgeeffectlength(
            exmotifs=exmotifs)

        # Extract Left and Right Context
        if (not leftcontext  is None) or \
           (not rightcontext is None):

           # Set Context Flag
            has_context = True

            # Show update
            liner.send('\n[Step 4: Extracting Context Sequences]\n')

            # Extract Both Contexts
            (leftcontext,
            rightcontext) = ut.get_extracted_context(
                leftcontext=leftcontext,
                rightcontext=rightcontext,
                edgeeffectlength=edgeeffectlength,
                reduce=False,
                liner=liner)

    # Finalize Context
    if not has_context:
        (leftcontext,
         rightcontext) = None, None

    # Parse Oligopool Repeats
    liner.send('\n[Step 5: Parsing Oligopool Repeats]\n')

    # Parse Repeats from indf
    (parsestatus,
    sourcecontext,
    kmerspace,
    fillcount,
    freecount,
    oligorepeats) = ut.get_parsed_oligopool_repeats(
        df=indf,
        maxreplen=maxreplen,
        element='Barcode',
        merge=False,
        liner=liner)

    # Repeat Length infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status': False,
            'basis' : 'infeasible',
            'step'  : 5,
            'vars'  : {
                'source_context': sourcecontext,
                'kmer_space'    : kmerspace,
                'fill_count'    : fillcount,
                'free_count'    : freecount},
            'warns' : warns}

        # Return results
        liner.close()
        return (outdf, stats)

    # Launching Barcode Design
    liner.send('\n[Step 6: Computing Barcodes]\n')

    # Define Barcode Design Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 6,
        'step_name': 'computing-barcodes',
        'vars'    : {
               'target_count': targetcount,   # Required Number of Barcodes
              'barcode_count': 0,             # Barcode Design Count
                  'type_fail': 0,             # Barcode Tyoe Failure Count
              'distance_fail': 0,             # Hamming Distance Fail Count
                'repeat_fail': 0,             # Repeat Fail Count
               'exmotif_fail': 0,             # Exmotif Elimination Fail Count
                  'edge_fail': 0,             # Edge Effect Fail Count
            'distance_distro': None,          # Hamming Distance Distribution
            'exmotif_counter': cx.Counter()}, # Exmotif Encounter Counter
        'warns'   : warns}

    # Schedule outfile deletion
    ofdeletion = ae.register(
        ut.remove_file,
        outfile)

    # Design Barcodes
    (codes,
    store,
    stats) = cb.barcode_engine(
        barcodelen=barcodelen,
        minhdist=minhdist,
        maxreplen=maxreplen,
        barcodetype=barcodetype,
        oligorepeats=oligorepeats,
        leftcontext=leftcontext,
        rightcontext=rightcontext,
        exmotifs=exmotifs,
        targetcount=targetcount,
        stats=stats,
        liner=liner)

    # Success Relevant Stats
    if not store is None:

        # Launching Distance Distribution Analysis
        liner.send('\n[Step 7: Computing Distance Distribution]\n')

        # Compute Hamming Distance Distribution
        stats['vars']['distance_distro'] = cb.get_distro(
            store=store,
            liner=liner)

    # Barcode Status
    if stats['status']:
        barcodestatus = 'Successful'
    else:
        barcodestatus = 'Failed'

    # Insert codes into indf
    if stats['status']:

        # Update indf
        ut.update_df(
            indf=indf,
            lcname=leftcontextname,
            rcname=rightcontextname,
            out=codes,
            outcol=barcodecol)

        # Prepare outdf
        outdf = indf

        # Write outdf to file
        if not outfile is None:
            outdf.to_csv(
                path_or_buf=outfile,
                sep=',')

    # Barcoding Statistics
    liner.send('\n[Barcode Design Statistics]\n')

    plen = ut.get_printlen(
        value=max(stats['vars'][field] for field in (
            'target_count',
            'barcode_count')))

    liner.send(
        '   Design Status   : {}\n'.format(
            barcodestatus))
    liner.send(
        '   Target Count    : {:{},d} Barcode(s)\n'.format(
            stats['vars']['target_count'],
            plen))
    liner.send(
        '  Barcode Count    : {:{},d} Barcode(s) ({:6.2f} %)\n'.format(
            stats['vars']['barcode_count'],
            plen,
            ut.safediv(
                A=stats['vars']['barcode_count'] * 100.,
                B=targetcount)))

    # Success Relevant Stats
    if stats['status']:
        if stats['vars']['distance_distro']:

            dlen = ut.get_printlen(
                value=max(map(
                    lambda x: x[1],
                    stats['vars']['distance_distro'])))

            liner.send('   Pair-wise Distance Distribution\n')

            for percentage,distance in stats['vars']['distance_distro']:
                liner.send(
                    '     - {:6.2f} % Barcode(s) w/ Distance ≥ {:{},d} Mismatch(es)\n'.format(
                        percentage,
                        distance,
                        dlen))

    # Failure Relavant Stats
    else:
        maxval = max(stats['vars'][field] for field in (
            'distance_fail',
            'repeat_fail',
            'exmotif_fail',
            'edge_fail'))

        sntn, plen = ut.get_notelen(
            printlen=ut.get_printlen(
                value=maxval))

        total_conflicts = stats['vars']['distance_fail'] + \
                          stats['vars']['repeat_fail']   + \
                          stats['vars']['exmotif_fail']  + \
                          stats['vars']['edge_fail']
        liner.send(
            ' Distance Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['distance_fail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['distance_fail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '   Repeat Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['repeat_fail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['repeat_fail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '  Exmotif Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['exmotif_fail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['exmotif_fail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '     Edge Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['edge_fail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['edge_fail'] * 100.,
                    B=total_conflicts)))

        # Enumerate Motif-wise Fail Counts
        if stats['vars']['exmotif_counter']:

            qlen = max(len(motif) \
                for motif in stats['vars']['exmotif_counter'].keys()) + 2

            sntn, vlen = ut.get_notelen(
                printlen=ut.get_printlen(
                    value=max(
                        stats['vars']['exmotif_counter'].values())))

            liner.send('   Exmotif-wise Conflict Distribution\n')

            for exmotif,count in stats['vars']['exmotif_counter'].most_common():
                exmotif = '\'{}\''.format(exmotif)
                liner.send(
                    '     - Motif {:>{}} Triggered {:{},{}} Event(s)\n'.format(
                        exmotif, qlen, count, vlen, sntn))

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Unschedule outfile deletion
    if barcodestatus == 'Successful':
        ae.unregister(ofdeletion)

    # Close Liner
    liner.close()

    # Return Solution and Statistics
    return (outdf, stats)
