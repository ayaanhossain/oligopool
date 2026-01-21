import time as tt

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
    maximum_repeat_length:int,
    spacer_column:str,
    output_file:str|None=None,
    spacer_length:int|list|str|pd.DataFrame|None=None,
    left_context_column:str|None=None,
    right_context_column:str|None=None,
    patch_mode:bool=False,
    excluded_motifs:list|str|pd.DataFrame|None=None,
    verbose:bool=True,
    random_seed:int|None=None) -> Tuple[pd.DataFrame, dict]:
    '''
    Insert neutral spacer DNA to meet length targets under repeat and excluded-motif constraints.
    Supports fixed-length, per-variant, or auto-sized spacers (to match `oligo_length_limit`).

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated oligopool variants.
        - `oligo_length_limit` (`int`): Maximum allowed oligo length (≥ 4).
        - `maximum_repeat_length` (`int`): Max shared repeat length with oligos (≥ 4).
        - `spacer_column` (`str`): Column name for inserting the designed spacers.

    Optional Parameters:
        - `output_file` (`str`): Filename for output DataFrame; required in CLI usage,
            optional in library usage (default: `None`).
        - `spacer_length` (`int` / `list` / `str` / `pd.DataFrame`): Length of the inserted spacers,
            can be defined per oligo in a list or a DataFrame; if `None` the spacer length per oligo
            is determined automatically to match `oligo_length_limit` (default: `None`).
        - `left_context_column` (`str`): Column for left DNA context (default: `None`).
        - `right_context_column` (`str`): Column for right DNA context (default: `None`).
        - `patch_mode` (`bool`): If `True`, fill only missing values in an existing spacer column
            (does not overwrite existing spacers). (Default: `False`).
        - `excluded_motifs` (`list` / `str` / `pd.DataFrame`): Motifs to exclude;
            can be a CSV path or DataFrame (default: `None`).
        - `verbose` (`bool`): If `True`, logs updates to stdout (default: `True`).
        - `random_seed` (`int` / `None`): Seed for local RNG (default: `None`).

    Returns:
        - A pandas DataFrame of inserted spacers; saves to `output_file` if specified.
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` must contain a unique 'ID' column, all other columns must be non-empty DNA strings.
        - Column names in `input_data` must be unique, and exclude `spacer_column`.
        - At least one of `left_context_column` or `right_context_column` must be specified.
        - When `spacer_length` is a CSV or DataFrame, it must have 'ID' and 'Length' columns.
        - If `excluded_motifs` is a CSV or DataFrame, it must have an 'Exmotif' column.
        - Oligo rows already summing to or exceeding `oligo_length_limit` have a `'-'` (dash) as spacer.
        - Patch mode (`patch_mode=True`) supports incremental pool extension: existing values in
          `spacer_column` are preserved and only missing values (e.g., `None`/NaN/empty/`'-'`) are
          filled (some rows may still get `'-'` if no spacer can fit under `oligo_length_limit`).
    '''

    # Preserve return style when the caller intentionally used ID as index.
    id_from_index = ut.get_id_index_intent(input_data)

    # Argument Aliasing
    indata       = input_data
    oligolimit   = oligo_length_limit
    maxreplen    = maximum_repeat_length
    spacercol    = spacer_column
    outfile      = output_file
    spacerlen    = spacer_length
    leftcontext  = left_context_column
    rightcontext = right_context_column
    patch_mode   = patch_mode
    exmotifs     = excluded_motifs
    verbose      = verbose
    random_seed  = random_seed

    # Local RNG
    rng = np.random.default_rng(random_seed)

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Spacer Verbiage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Spacer]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # Patch Mode (pre-parse for output column handling)
    patch_mode_valid = isinstance(patch_mode, (bool, int, np.bool_)) and \
        (patch_mode in (0, 1, True, False))
    patch_mode_on = bool(patch_mode) if patch_mode_valid else False
    allow_missing_cols = None
    if patch_mode_on and isinstance(spacercol, str):
        allow_missing_cols = set((spacercol,))

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_indata_info(
        indata=indata,
        indata_field='      Input Data   ',
        required_fields=('ID',),
        precheck=False,
        liner=liner,
        allow_missing_cols=allow_missing_cols)
    input_rows = len(indf.index) if isinstance(indf, pd.DataFrame) else 0

    # Full oligolimit Validation
    oligolimit_valid = vp.get_numeric_validity(
        numeric=oligolimit,
        numeric_field='      Oligo Limit  ',
        numeric_pre_desc=' At most ',
        numeric_post_desc=' Base Pair(s)',
        minval=4,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full maxreplen Validation
    maxreplen_valid = vp.get_numeric_validity(
        numeric=maxreplen,
        numeric_field='     Repeat Length ',
        numeric_pre_desc=' Up to ',
        numeric_post_desc=' Base Pair(s) Oligopool Repeats',
        minval=4,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full spacercol Validation
    spacercol_valid = vp.get_parsed_column_info(
        col=spacercol,
        df=indf,
        col_field='     Spacer Column ',
        col_desc='Output in Column',
        col_type=1,
        adjcol=None,
        adjval=None,
        iscontext=False,
        typecontext=None,
        liner=liner,
        allow_existing=patch_mode_on)

    # Full outfile Validation
    outfile_valid = vp.get_outdf_validity(
        outdf=outfile,
        outdf_suffix='.oligopool.spacer.csv',
        outdf_field='     Output File   ',
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
        spacerlen_field='     Spacer Length ',
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
        col_field='       Left Context',
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
        col_field='      Right Context',
        col_desc='Input from Column',
        col_type=0,
        adjcol=leftcontextname,
        adjval=-1,
        iscontext=True,
        typecontext=1,
        liner=liner)

    # Patch mode parsing (printed after context args; evaluated earlier for indata handling)
    if patch_mode_valid:
        liner.send('{}: {}\n'.format(
            '      Patch Mode   ',
            ['Disabled (Design All Spacers)', 'Enabled (Fill Missing Spacers)'][patch_mode_on]))
    else:
        liner.send('{}: {} [INPUT TYPE IS INVALID]\n'.format(
            '      Patch Mode   ',
            patch_mode))

    # Full exmotifs Parsing and Validation
    (exmotifs,
    exmotifs_valid) = vp.get_parsed_exseqs_info(
        exseqs=exmotifs,
        exseqs_field='   Excluded Motifs ',
        exseqs_desc='Unique Motif(s)',
        df_field='Exmotif',
        required=False,
        liner=liner)

    # First Pass Validation
    if not all([
        indata_valid,
        oligolimit_valid,
        maxreplen_valid,
        spacercol_valid,
        outfile_valid,
        spacerlen_valid,
        patch_mode_valid,
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

    # Patch mode bookkeeping
    spacercol_exists = False
    missing_mask = None
    existing_mask = None
    if patch_mode_on and isinstance(indf, pd.DataFrame):
        (spacercol_exists,
        _) = ut.get_col_exist_idx(
            col=spacercol,
            df=indf)
        if spacercol_exists:
            missing_mask = ut.get_missing_mask(
                series=indf[spacercol],
                allow_dash=True)
            indf[spacercol] = ut.fill_missing_values(
                series=indf[spacercol],
                missing_mask=missing_mask,
                fill='-')
            existing_mask = ~missing_mask
            targetcount = int(missing_mask.sum())
        else:
            missing_mask = np.ones(
                len(indf.index),
                dtype=bool)
            existing_mask = ~missing_mask

        # Validate existing spacers (DNA only)
        if spacercol_exists and existing_mask.any():
            invalid_mask = np.zeros(
                len(indf.index),
                dtype=bool)
            for idx in np.where(existing_mask)[0]:
                value = indf[spacercol].iat[idx]
                if not ut.is_DNA(
                        seq=value,
                        dna_alpha=set('ATGC')):
                    invalid_mask[idx] = True
            if invalid_mask.any():
                examples = ut.get_row_examples(
                    df=indf,
                    invalid_mask=invalid_mask,
                    id_col='ID',
                    limit=5)
                example_note = ut.format_row_examples(examples)
                liner.send(
                    '     Spacer Column : Output in Column \'{}\' [INVALID EXISTING VALUE]{}\n'.format(
                        spacercol,
                        example_note))
                liner.send('\n')
                raise RuntimeError(
                    'Invalid Argument Input(s).')

    # No missing spacers in patch mode
    if patch_mode_on and \
       missing_mask is not None and \
       not missing_mask.any():

        stats = {
            'status'  : True,
            'basis'   : 'complete',
            'step'    : 0,
            'step_name': 'no-missing-spacers',
            'vars'    : {
               'target_count': 0,
                'spacer_count': 0,
               'orphan_oligo': [],
                'repeat_fail': 0,
               'exmotif_fail': 0,
                  'edge_fail': 0,
            'exmotif_counter': cx.Counter()},
            'warns'   : warns}
        stats['random_seed'] = random_seed

        outdf = indf
        if not outfile is None:
            ut.write_df_csv(
                df=outdf,
                outfile=outfile,
                sep=',')

        liner.send('\n[Spacer Design Statistics]\n')
        liner.send(
            '  Design Status   : Successful\n')
        liner.send(
            '  Target Count    : 0 Spacer(s)\n')
        liner.send(
            '  Spacer Count    : 0 Spacer(s) (  0.00 %)\n')
        liner.send(
            '  Orphan Oligo    : 0 Entries\n')
        liner.send(
            ' Time Elapsed: {:.2f} sec\n'.format(
                tt.time()-t0))

        liner.close()
        stats = ut.stamp_stats(
            stats=stats,
            module='spacer',
            input_rows=input_rows,
            output_rows=len(outdf.index))
        outdf_return = outdf
        if not id_from_index:
            outdf_return = ut.get_df_with_id_column(outdf)
        return (outdf_return, stats)

    # Extract spacerlen (Auto-Inference)
    if spacerlen is None:
        liner.send('\n[Step 1: Extracting Spacer Length]\n')
        (spacerlen,
        variantlens) = cm.get_extracted_spacerlen(
            indf=indf,
            oligolimit=oligolimit,
            liner=liner)

    # Target spacer lengths (patch mode)
    spacerlen_target = spacerlen
    variantlens_target = variantlens
    if missing_mask is not None:
        spacerlen_target = spacerlen[missing_mask]
        if not variantlens is None:
            variantlens_target = variantlens[missing_mask]

    # Parse Oligopool Limit Feasibility
    liner.send('\n[Step 2: Parsing Oligo Limit]\n')

    # Parse oligolimit
    (parsestatus,
    minoligolen,
    maxoligolen,
    minelementlen,
    maxelementlen,
    minspaceavail,
    maxspaceavail) = ut.get_parsed_oligolimit(
        indf=indf,
        variantlens=variantlens_target,
        oligolimit=oligolimit,
        minelementlen=np.min(spacerlen_target),
        maxelementlen=np.max(spacerlen_target),
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
                  'min_oligo_len': minoligolen,
                  'max_oligo_len': maxoligolen,
                'min_element_len': minelementlen,
                'max_element_len': maxelementlen,
                'min_space_avail': minspaceavail,
                'max_space_avail': maxspaceavail},
            'warns'   : warns}
        stats['random_seed'] = random_seed
        stats = ut.stamp_stats(
            stats=stats,
            module='spacer',
            input_rows=input_rows,
            output_rows=0)

        # Return results
        liner.close()
        return (outdf, stats)

    # Extract Spacer Length Groups
    liner.send('\n[Step 3: Extracting Spacer Groups]\n')

    # Group spacerlen
    spacergroup = cm.get_grouped_spacerlen(
        spacerlen=spacerlen_target,
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
            stats['random_seed'] = random_seed
            stats = ut.stamp_stats(
                stats=stats,
                module='spacer',
                input_rows=input_rows,
                output_rows=0)

            # Return results
            liner.close()
            return (outdf, stats)

        # Update Edge-Effect Length
        edgeeffectlength = ut.get_edgeeffectlength(
            maxreplen=maxreplen,
            exmotifs=exmotifs)

        # Parse Edge Effects
        if ((not leftcontext  is None) or \
            (not rightcontext is None)):

            # Set Context Flag
            has_context = True

            # Show Update
            liner.send('\n[Step 5: Extracting Context Sequences]\n')

            # Extract Both Contexts
            # In patch mode, only extract context for rows being designed.
            lcontext = leftcontext
            rcontext = rightcontext
            if missing_mask is not None:
                if not lcontext is None:
                    lcontext = lcontext[missing_mask]
                if not rcontext is None:
                    rcontext = rcontext[missing_mask]

            (leftcontext,
            rightcontext) = ut.get_extracted_context(
                leftcontext=lcontext,
                rightcontext=rcontext,
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
                motiftype=0,
                leftcontext=leftcontext,
                rightcontext=rightcontext,
                leftpartition=leftpartition,
                rightpartition=rightpartition,
                exmotifs=exmotifs,
                element='Spacer',
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

    # Parse Oligopool Repeats
    liner.send('\n[Step 7: Parsing Oligopool Repeats]\n')

    # Parse Repeats from indf
    (parsestatus,
    sourcecontext,
    kmerspace,
    fillcount,
    freecount,
    oligorepeats) = ut.get_parsed_oligopool_repeats(
        df=indf if missing_mask is None else indf.loc[missing_mask],
        maxreplen=maxreplen,
        element='Spacer',
        merge=False,
        liner=liner)

    # Repeat Length infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status': False,
            'basis' : 'infeasible',
            'step'  : 7,
            'step_name': 'parsing-oligopool-repeats',
            'vars'  : {
                'source_context': sourcecontext,
                'kmer_space'    : kmerspace,
                'fill_count'    : fillcount,
                'free_count'    : freecount},
            'warns' : warns}
        stats['random_seed'] = random_seed
        stats = ut.stamp_stats(
            stats=stats,
            module='spacer',
            input_rows=input_rows,
            output_rows=0)

        # Return results
        liner.close()
        return (outdf, stats)

    # Launching Motif Design
    liner.send('\n[Step 8: Computing Spacers]\n')

    # Define Motif Design Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 8,
        'step_name': 'computing-spacers',
        'vars'    : {
               'target_count': targetcount,   # Required Number of Spacers
                'spacer_count': 0,             # Spacer Design Count
               'orphan_oligo': None,          # Orphan Oligo Indexes
                'repeat_fail': 0,             # Repeat Fail Count
               'exmotif_fail': 0,             # Exmotif Elimination Fail Count
                  'edge_fail': 0,             # Edge Effect Fail Count
            'exmotif_counter': cx.Counter()}, # Exmotif Encounter Counter
        'warns'   : warns}
    stats['random_seed'] = random_seed

    # Schedule outfile deletion
    ofdeletion = ae.register(
        ut.remove_file,
        outfile)

    # Design Spacers
    (spacers,
    stats) = cm.spacer_engine(
        spacergroup=spacergroup,
        maxreplen=maxreplen,
        oligorepeats=oligorepeats,
        leftcontext=leftcontext,
        rightcontext=rightcontext,
        exmotifs=exmotifs,
        edgeeffectlength=edgeeffectlength,
        prefixdict=prefixdict,
        suffixdict=suffixdict,
        targetcount=targetcount,
        stats=stats,
        liner=liner,
        rng=rng)

    # Spacer Status
    if stats['status']:
        spacerstatus = 'Successful'
    else:
        spacerstatus = 'Failed'

    # Insert spacer into indf
    if stats['status']:

        # Update indf
        if spacercol_exists:
            indf.loc[missing_mask, spacercol] = spacers
        else:
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
            ut.write_df_csv(
                df=outdf,
                outfile=outfile,
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
    liner.send(
        '  Orphan Oligo    : {:{},d} Entries\n'.format(
            len(stats['vars']['orphan_oligo']),
            plen))

    # Failure Relavant Stats
    if not stats['status']:
        maxval = max(stats['vars'][field] for field in (
            'repeat_fail',
            'exmotif_fail',
            'edge_fail'))

        sntn, plen = ut.get_notelen(
            printlen=ut.get_printlen(
                value=maxval))

        total_conflicts = stats['vars']['repeat_fail']  + \
                          stats['vars']['exmotif_fail'] + \
                          stats['vars']['edge_fail']

        liner.send(
            '  Repeat Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['repeat_fail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['repeat_fail'] * 100.,
                    B=total_conflicts)))
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
    stats = ut.stamp_stats(
        stats=stats,
        module='spacer',
        input_rows=input_rows,
        output_rows=len(outdf.index) if outdf is not None else 0)
    outdf_return = outdf
    if (outdf is not None) and (not id_from_index):
        outdf_return = ut.get_df_with_id_column(outdf)
    return (outdf_return, stats)
