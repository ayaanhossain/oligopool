import time  as tt

import collections as cx
import atexit as ae

import numpy  as np
import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp
from .base import core_motif as cm

from typing import Tuple


def spacer(
    input_data:str|pd.DataFrame,
    oligo_length_limit:int,
    spacer_column:str,
    output_file:str|None=None,
    spacer_length:int|None=None,
    left_context_column:str|None=None,
    right_context_column:str|None=None,
    excluded_motifs:list|str|pd.DataFrame|None=None,
    verbose:bool=True) -> Tuple[pd.DataFrame, dict]:
    '''
    Inserts a fixed or variable length neutral spacer free from excluded motifs. If the spacer length
    is not specified, then the spacer per oligo is such that the final length of the oligo matches
    specified `oligo_length_limit`. Returned DataFrame can be written as a CSV file.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated oligopool variants.
        - `oligo_length_limit` (`int`): Maximum allowed oligo length (≥ 4).
        - `spacer_column` (`str`): Column name for inserting the designed spacers.

    Optional Parameters:
        - `output_file` (`str`): Filename for output DataFrame (default: `None`).
        - `spacer_length` (`int`): Length of the inserted spacers,
            if `None` automatically determine this length per oligo (default: `None`).
        - `left_context_column` (`str`): Column for left DNA context (default: `None`).
        - `right_context_column` (`str`): Column for right DNA context (default: `None`).
        - `excluded_motifs` (`list` / `str` / `pd.DataFrame`): Motifs to exclude;
            can be a CSV path or DataFrame (default: `None`).
        - `verbose` (`bool`): If `True`, logs updates to stdout (default: `True`).

    Returns:
        - A pandas DataFrame of inserted spacers; saves to `output_file` if specified.
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` must contain a unique 'ID' column, all other columns must be non-empty DNA strings.
        - Column names in `input_data` must be unique, and exclude `spacer_column`.
        - At least one of `left_context_column` or `right_context_column` must be specified.
        - If `excluded_motifs` is a CSV or DataFrame, it must have 'ID' and 'Exmotif' columns.
        - Oligo rows already summing to or exceeding `oligo_length_limit` have a `'-'` (dash) as spacer.
    '''

    # Argument Aliasing
    indata       = input_data
    oligolimit   = oligo_length_limit
    spacercol    = spacer_column
    outfile      = output_file
    spacerlen    = spacer_length
    leftcontext  = left_context_column
    rightcontext = right_context_column
    exmotifs     = excluded_motifs
    verbose      = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Spacer Verbage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Spacer]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_indata_info(
        indata=indata,
        indata_field='    Input Data   ',
        required_fields=('ID',),
        precheck=False,
        liner=liner)

    # Full oligolimit Validation
    oligolimit_valid = vp.get_numeric_validity(
        numeric=oligolimit,
        numeric_field='    Oligo Limit  ',
        numeric_pre_desc=' At most ',
        numeric_post_desc=' Base Pair(s)',
        minval=4,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full spacercol Validation
    spacercol_valid = vp.get_parsed_column_info(
        col=spacercol,
        df=indf,
        col_field='   Spacer Column ',
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
        outdf_suffix='.oligopool.spacer.csv',
        outdf_field='   Output File   ',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.spacer.csv')

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full spacerlen Validation
    (spacerlen,
    spacerlen_valid) = vp.get_parsed_spacerlen_info(
        spacerlen=spacerlen,
        spacerlen_field='   Spacer Length ',
        df_field='Length',
        oligolimit=oligolimit,
        oligolimit_valid=oligolimit_valid,
        indf=indf,
        indata_valid=indata_valid,
        liner=liner)

    # Store Context Names
    leftcontextname  = leftcontext
    rightcontextname = rightcontext

    # Full leftcontext Parsing and Validation
    (leftcontext,
    leftcontext_valid) = vp.get_parsed_column_info(
        col=leftcontext,
        df=indf,
        col_field='     Left Context',
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
        col_field='    Right Context',
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
        exseqs_field=' Excluded Motifs ',
        exseqs_desc='Unique Motif(s)',
        df_field='Exmotif',
        required=False,
        liner=liner)

    # First Pass Validation
    if not all([
        indata_valid,
        oligolimit_valid,
        spacercol_valid,
        outfile_valid,
        spacerlen_valid,
        leftcontext_valid,
        rightcontext_valid,
        exmotifs_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Paramters
    oligolimit = round(oligolimit)

    # Define Edge Effect Length
    edgeeffectlength = None

    # Spacer Design Book-keeping
    targetcount = len(indf.index)
    has_context = False
    variantlens = None
    outdf = None
    stats = None
    warns = {}

    # Extract spacerlen (Auto-Inference)
    if spacerlen is None:
        liner.send('\n[Step 1: Extracting Spacer Length]\n')
        (spacerlen,
        variantlens) = cm.get_extracted_spacerlen(
            indf=indf,
            oligolimit=oligolimit,
            liner=liner)

    # Parse Oligopool Limit Feasibility
    liner.send('\n[Step 2: Parsing Oligo Limit]\n')

    # Parse oligolimit
    (parsestatus,
    minvariantlen,
    maxvariantlen,
    minelementlen,
    maxelementlen,
    minspaceavail,
    maxspaceavail) = ut.get_parsed_oligolimit(
        indf=indf,
        variantlens=variantlens,
        oligolimit=oligolimit,
        minelementlen=np.min(spacerlen),
        maxelementlen=np.max(spacerlen),
        element='Spacer',
        liner=liner)

    # oligolimit infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 2,
            'step_name': 'parsing-oligo-limit',
            'vars'    : {
                   'oligo_limit': oligolimit,
                'limit_overflow': True,
                'min_varian_tlen': minvariantlen,
                'max_variant_len': maxvariantlen,
                'min_element_len': minelementlen,
                'max_element_len': maxelementlen,
                'min_space_avail': minspaceavail,
                'max_space_avail': maxspaceavail},
            'warns'   : warns}

        # Return results
        liner.close()
        return (outdf, stats)

    # Extract Spacer Length Groups
    liner.send('\n[Step 3: Extracting Spacer Groups]\n')

    # Group spacerlen
    spacergroup = cm.get_grouped_spacerlen(
        spacerlen=spacerlen,
        liner=liner)

    # Parse Excluded Motifs
    if not exmotifs is None:

        # Show update
        liner.send('\n[Step 4: Parsing Excluded Motifs]\n')

        # Update Step 4 Warning
        warns[4] = {
            'warn_count': 0,
            'step_name' : 'parsing-excluded-motifs',
            'vars': None}

        # Parse exmotifs
        (parsestatus,
        exmotifs,
        problens,
        leftpartition,
        rightpartition) = ut.get_parsed_exmotifs(
            exmotifs=exmotifs,
            typer=tuple,
            element='Spacer',
            leftcontext=leftcontext,
            rightcontext=rightcontext,
            warn=warns[4],
            liner=liner)

        # Remove Step 4 Warning
        if not warns[4]['warn_count']:
            warns.pop(4)

        # exmotifs infeasible
        if not parsestatus:

            # Prepare stats
            stats = {
                'status'  : False,
                'basis'   : 'infeasible',
                'step'    : 4,
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

        # Parse Edge Effects
        if ((not leftcontext  is None) or \
            (not rightcontext is None)):

            # Set Context Flag
            has_context = True

            # Show Update
            liner.send('\n[Step 5: Extracting Context Sequences]\n')

            # Extract Both Contexts
            (leftcontext,
            rightcontext) = ut.get_extracted_context(
                leftcontext=leftcontext,
                rightcontext=rightcontext,
                edgeeffectlength=edgeeffectlength,
                reduce=False,
                liner=liner)

            # Show update
            liner.send('\n[Step 6: Parsing Edge Effects]\n')

            # Update Step 6 Warning
            warns[6] = {
                'warn_count': 0,
                'step_name' : 'parsing-edge-effects',
                'vars': None}

            # Compute Forbidden Prefixes and Suffixes
            (prefixdict,
            suffixdict) = cm.get_parsed_edgeeffects(
                motifseq='NNNN',
                leftcontext=leftcontext,
                rightcontext=rightcontext,
                leftpartition=leftpartition,
                rightpartition=rightpartition,
                exmotifs=exmotifs,
                warn=warns[6],
                liner=liner)

            # Remove Step 6 Warning
            if not warns[6]['warn_count']:
                warns.pop(6)

    # Finalize Context
    if not has_context:
        (leftcontext,
        rightcontext,
        prefixdict,
        suffixdict) = (None, None, None, None)

    # Launching Motif Design
    liner.send('\n[Step 7: Computing Spacers]\n')

    # Define Motif Design Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 7,
        'step_name': 'computing-spacers',
        'vars'    : {
               'target_count': targetcount,   # Required Number of Spacers
               'spacer_count': 0,             # Spacer Design Count
               'exmotif_fail': 0,             # Exmotif Elimination Fail Count
                  'edge_fail': 0,             # Edge Effect Fail Count
            'exmotif_counter': cx.Counter()}, # Exmotif Encounter Counter
        'warns'   : warns}

    # Schedule outfile deletion
    ofdeletion = ae.register(
        ut.remove_file,
        outfile)

    # Design Spacers
    (spacers,
    stats) = cm.spacer_engine(
        spacergroup=spacergroup,
        leftcontext=leftcontext,
        rightcontext=rightcontext,
        exmotifs=exmotifs,
        edgeeffectlength=edgeeffectlength,
        prefixdict=prefixdict,
        suffixdict=suffixdict,
        targetcount=targetcount,
        stats=stats,
        liner=liner)

    # Spacer Status
    if stats['status']:
        spacerstatus = 'Successful'
    else:
        spacerstatus = 'Failed'

    # Insert spacer into indf
    if stats['status']:

        # Update indf
        ut.update_df(
            indf=indf,
            lcname=leftcontextname,
            rcname=rightcontextname,
            out=spacers,
            outcol=spacercol)

        # Prepare outdf
        outdf = indf

        # Write outdf to file
        if not outfile is None:
            outdf.to_csv(
                path_or_buf=outfile,
                sep=',')

    # Spacer Design Statistics
    liner.send('\n[Spacer Design Statistics]\n')

    plen = ut.get_printlen(
        value=max(stats['vars'][field] for field in (
            'target_count',
            'spacer_count')))

    liner.send(
        '  Design Status   : {}\n'.format(
            spacerstatus))
    liner.send(
        '  Target Count    : {:{},d} Spacer(s)\n'.format(
            stats['vars']['target_count'],
            plen))
    liner.send(
        '  Spacer Count    : {:{},d} Spacer(s) ({:6.2f} %)\n'.format(
            stats['vars']['spacer_count'],
            plen,
            ut.safediv(
                A=stats['vars']['spacer_count'] * 100.,
                B=targetcount)))

    # Failure Relavant Stats
    if not stats['status']:
        maxval = max(stats['vars'][field] for field in (
            'exmotif_fail',
            'edge_fail'))

        sntn, plen = ut.get_notelen(
            printlen=ut.get_printlen(
                value=maxval))

        total_conflicts = stats['vars']['exmotif_fail']  + \
                          stats['vars']['edge_fail']
        liner.send(
            ' Exmotif Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['exmotif_fail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['exmotif_fail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '    Edge Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
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
    if spacerstatus == 'Successful':
        ae.unregister(ofdeletion)

    # Close Liner
    liner.close()

    # Return Solution and Statistics
    return (outdf, stats)
