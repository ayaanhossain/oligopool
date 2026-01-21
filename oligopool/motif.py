import time as tt

import collections as cx
import atexit as ae

import numpy  as np
import pandas  as pd

from .base import utils as ut
from .base import validation_parsing as vp
from .base import core_motif as cm

from typing import Tuple


def motif(
    input_data:str|pd.DataFrame,
    oligo_length_limit:int,
    motif_sequence_constraint:str,
    maximum_repeat_length:int,
    motif_column:str,
    output_file:str|None=None,
    motif_type:int=0,
    left_context_column:str|None=None,
    right_context_column:str|None=None,
    patch_mode:bool=False,
    excluded_motifs:list|str|pd.DataFrame|None=None,
    verbose:bool=True,
    random_seed:int|None=None) -> Tuple[pd.DataFrame, dict]:
    '''
    Add or design constrained motifs under a sequence constraint with repeat/excluded-motif screening
    (including edge effects with context). Supports per-variant motifs and constant motif anchors for
    building indexable architectures.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated oligopool variants.
        - `oligo_length_limit` (`int`): Maximum allowed oligo length (≥ 4).
        - `motif_sequence_constraint` (`str`): IUPAC degenerate sequence constraint, or a constant.
        - `maximum_repeat_length` (`int`): Max shared repeat length with oligos (≥ 4).
        - `motif_column` (`str`): Column name for inserting the designed motifs.

    Optional Parameters:
        - `output_file` (`str`): Filename for output DataFrame; required in CLI usage,
            optional in library usage (default: `None`).
        - `motif_type` (`int`): Motif type to design
            (0 for per-variant motifs, 1 for a single constant motif shared by all variants; default: 0).
        - `left_context_column` (`str`): Column for left DNA context (default: `None`).
        - `right_context_column` (`str`): Column for right DNA context (default: `None`).
        - `patch_mode` (`bool`): If `True`, fill only missing values in an existing motif/anchor
            column (does not overwrite existing motifs). (Default: `False`).
        - `excluded_motifs` (`list` / `str` / `pd.DataFrame`): Motifs to exclude;
            can be a CSV path or DataFrame (default: `None`).
        - `verbose` (`bool`): If `True`, logs updates to stdout (default: `True`).
        - `random_seed` (`int` / `None`): Seed for local RNG (default: `None`).

    Returns:
        - A pandas DataFrame of added motifs; saves to `output_file` if specified.
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` must contain a unique 'ID' column, all other columns must be non-empty DNA strings.
        - Column names in `input_data` must be unique, and exclude `motif_column`.
        - At least one of `left_context_column` or `right_context_column` must be specified.
        - If `excluded_motifs` is a CSV or DataFrame, it must have an 'Exmotif' column.
        - Constant bases in sequence constraint may lead to `excluded_motifs` and be impossible to solve.
        - Use `motif_type=1` to design constant motif anchors (e.g., barcode prefix/suffix anchors for indexing).
        - For anchors, tune `maximum_repeat_length` to control how distinct the anchor is from the surrounding oligos.
        - For anchors, `motif_sequence_constraint` can be an IUPAC pattern (e.g., 'NNNNNNNNNN') or include fixed bases.
        - Anchors should typically be designed prior to barcode generation.
        - Patch mode (`patch_mode=True`) supports incremental pool extension: existing values in
          `motif_column` are preserved and only missing values (e.g., `None`/NaN/empty/`'-'`) are
          filled. For `motif_type=1`, an existing compatible constant anchor is reused for new rows.
    '''

    # Preserve return style when the caller intentionally used ID as index.
    id_from_index = ut.get_id_index_intent(input_data)

    # Argument Aliasing
    indata       = input_data
    oligolimit   = oligo_length_limit
    motifseq     = motif_sequence_constraint
    maxreplen    = maximum_repeat_length
    motifcol     = motif_column
    outfile      = output_file
    motiftype    = motif_type
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

    # Motif Verbiage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Motif]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # Patch Mode (pre-parse for output column handling)
    patch_mode_valid = isinstance(patch_mode, (bool, int, np.bool_)) and \
        (patch_mode in (0, 1, True, False))
    patch_mode_on = bool(patch_mode) if patch_mode_valid else False
    allow_missing_cols = None
    if patch_mode_on and isinstance(motifcol, str):
        allow_missing_cols = set((motifcol,))

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_indata_info(
        indata=indata,
        indata_field='      Input Data    ',
        required_fields=('ID',),
        precheck=False,
        liner=liner,
        allow_missing_cols=allow_missing_cols)
    input_rows = len(indf.index) if isinstance(indf, pd.DataFrame) else 0

    # Full oligolimit Validation
    oligolimit_valid = vp.get_numeric_validity(
        numeric=oligolimit,
        numeric_field='      Oligo Limit   ',
        numeric_pre_desc=' At most ',
        numeric_post_desc=' Base Pair(s)',
        minval=4,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # First Pass motifseq Validation
    motifseq_valid = vp.get_seqconstr_validity(
        seqconstr=motifseq,
        seqconstr_field='      Motif Sequence',
        minlenval=1,
        element='MOTIF',
        liner=liner)

    # Full maxreplen Validation
    maxreplen_valid = vp.get_numeric_validity(
        numeric=maxreplen,
        numeric_field='     Repeat Length  ',
        numeric_pre_desc=' Up to ',
        numeric_post_desc=' Base Pair(s) Oligopool Repeats',
        minval=4,
        maxval=len(motifseq) if motifseq_valid else float('inf'),
        precheck=False,
        liner=liner)

    # Full motifcol Validation
    motifcol_valid = vp.get_parsed_column_info(
        col=motifcol,
        df=indf,
        col_field='      Motif Column  ',
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
        outdf_suffix='.oligopool.motif.csv',
        outdf_field='     Output File    ',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.motif.csv')

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full motiftype Validation
    motiftype_valid = vp.get_categorical_validity(
        category=motiftype,
        category_field='      Motif Type    ',
        category_pre_desc=' ',
        category_post_desc=' Motifs',
        category_dict={
            0: 'Non-Constant',
            1: 'Constant'},
        liner=liner)

    # Store Context Names
    leftcontextname  = leftcontext
    rightcontextname = rightcontext

    # Full leftcontext Parsing and Validation
    (leftcontext,
    leftcontext_valid) = vp.get_parsed_column_info(
        col=leftcontext,
        df=indf,
        col_field='       Left Context ',
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
        col_field='      Right Context ',
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
            '      Patch Mode    ',
            ['Disabled (Design All Motifs)', 'Enabled (Fill Missing Motifs)'][patch_mode_on]))
    else:
        liner.send('{}: {} [INPUT TYPE IS INVALID]\n'.format(
            '      Patch Mode    ',
            patch_mode))

    # Full exmotifs Parsing and Validation
    (exmotifs,
    exmotifs_valid) = vp.get_parsed_exseqs_info(
        exseqs=exmotifs,
        exseqs_field='   Excluded Motifs  ',
        exseqs_desc='Unique Motif(s)',
        df_field='Exmotif',
        required=False,
        liner=liner)

    # First Pass Validation
    if not all([
        indata_valid,
        oligolimit_valid,
        motifseq_valid,
        maxreplen_valid,
        motifcol_valid,
        outfile_valid,
        motiftype_valid,
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

    # Motif Design Book-keeping
    targetcount = len(indf.index)
    has_context = False
    outdf = None
    stats = None
    warns = {}

    # Patch mode bookkeeping
    motifcol_exists = False
    missing_mask = None
    existing_mask = None
    if patch_mode_on and isinstance(indf, pd.DataFrame):
        (motifcol_exists,
        _) = ut.get_col_exist_idx(
            col=motifcol,
            df=indf)
        if motifcol_exists:
            missing_mask = ut.get_missing_mask(
                series=indf[motifcol],
                allow_dash=True)
            indf[motifcol] = ut.fill_missing_values(
                series=indf[motifcol],
                missing_mask=missing_mask,
                fill='-')
            existing_mask = ~missing_mask
            targetcount = int(missing_mask.sum())
        else:
            missing_mask = np.ones(
                len(indf.index),
                dtype=bool)
            existing_mask = ~missing_mask

        # Validate existing motifs (length + DNA)
        if motifcol_exists and existing_mask.any():
            invalid_mask = np.zeros(
                len(indf.index),
                dtype=bool)
            for idx in np.where(existing_mask)[0]:
                value = indf[motifcol].iat[idx]
                if (not ut.is_DNA(
                        seq=value,
                        dna_alpha=set('ATGC'))) or \
                   (len(value) != len(motifseq)):
                    invalid_mask[idx] = True
            if invalid_mask.any():
                examples = ut.get_row_examples(
                    df=indf,
                    invalid_mask=invalid_mask,
                    id_col='ID',
                    limit=5)
                example_note = ut.format_row_examples(examples)
                liner.send(
                    '      Motif Column  : Output in Column \'{}\' [INVALID EXISTING VALUE]{}\n'.format(
                        motifcol,
                        example_note))
                liner.send('\n')
                raise RuntimeError(
                    'Invalid Argument Input(s).')

        # Constant motif anchors: reuse existing constant motif if present
        if motiftype == 1 and \
           motifcol_exists and \
           existing_mask is not None and \
           existing_mask.any() and \
           missing_mask is not None and \
           missing_mask.any():

            uniques = ut.get_uniques(
                iterable=indf[motifcol].iloc[existing_mask],
                typer=tuple)

            # Must already be constant (patch mode does not overwrite)
            if len(uniques) != 1:
                liner.send(
                    '      Motif Column  : Output in Column \'{}\' [NON-UNIQUE EXISTING VALUES]\n'.format(
                        motifcol))
                liner.send('\n')
                raise RuntimeError(
                    'Invalid Argument Input(s).')

            motif_constant = str(uniques[0]).upper()
            motifseq_upper = str(motifseq).upper()
            compatible = True
            for idx, nt in enumerate(motif_constant):
                if not nt in ut.ddna_space[motifseq_upper[idx]]:
                    compatible = False
                    break
            if not compatible:
                liner.send(
                    '      Motif Column  : Output in Column \'{}\' [EXISTING VALUE VIOLATES SEQUENCE CONSTRAINT]\n'.format(
                        motifcol))
                liner.send('\n')
                raise RuntimeError(
                    'Invalid Argument Input(s).')

            indf.loc[missing_mask, motifcol] = motif_constant
            missing_mask = np.zeros(
                len(indf.index),
                dtype=bool)
            existing_mask = ~missing_mask
            targetcount = 0

    # No missing motifs in patch mode
    if patch_mode_on and \
       missing_mask is not None and \
       not missing_mask.any():

        stats = {
            'status'  : True,
            'basis'   : 'complete',
            'step'    : 0,
            'step_name': 'no-missing-motifs',
            'vars'    : {
               'target_count': 0,
                'motif_count': 0,
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

        liner.send('\n[Motif Design Statistics]\n')
        liner.send(
            '  Design Status   : Successful\n')
        liner.send(
            '  Target Count    : 0 Motif(s)\n')
        liner.send(
            '   Motif Count    : 0 Motif(s) (  0.00 %)\n')
        liner.send(
            '  Orphan Oligo    : 0 Entries\n')
        liner.send(
            ' Time Elapsed: {:.2f} sec\n'.format(
                tt.time()-t0))

        liner.close()
        stats = ut.stamp_stats(
            stats=stats,
            module='motif',
            input_rows=input_rows,
            output_rows=len(outdf.index))
        outdf_return = outdf
        if not id_from_index:
            outdf_return = ut.get_df_with_id_column(outdf)
        return (outdf_return, stats)

    # Parse Oligopool Limit Feasibility
    liner.send('\n[Step 1: Parsing Oligo Limit]\n')

    # Parse oligolimit
    (parsestatus,
    minoligolen,
    maxoligolen,
    minelementlen,
    maxelementlen,
    minspaceavail,
    maxspaceavail) = ut.get_parsed_oligolimit(
        indf=indf,
        variantlens=None if (missing_mask is None) else ut.get_variantlens(
            indf=indf.loc[missing_mask]),
        oligolimit=oligolimit,
        minelementlen=len(motifseq),
        maxelementlen=len(motifseq),
        element='Motif',
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
            module='motif',
            input_rows=input_rows,
            output_rows=0)

        # Return results
        liner.close()
        return (outdf, stats)

    # Parse Excluded Motifs
    if not exmotifs is None:

        # Show update
        liner.send('\n[Step 2: Parsing Excluded Motifs]\n')

        # Update Step 2 Warning
        warns[2] = {
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
            element='Motif',
            leftcontext=leftcontext,
            rightcontext=rightcontext,
            warn=warns[2],
            liner=liner)

        # Remove Step 2 Warning
        if not warns[2]['warn_count']:
            warns.pop(2)

        # exmotifs infeasible
        if not parsestatus:

            # Prepare stats
            stats = {
                'status'  : False,
                'basis'   : 'infeasible',
                'step'    : 2,
                'step_name': 'parsing-excluded-motifs',
                'vars'    : {
                     'prob_lens': problens,
                    'prob_count': tuple(list(
                        4**pl for pl in problens))},
                'warns'   : warns}
            stats['random_seed'] = random_seed
            stats = ut.stamp_stats(
                stats=stats,
                module='motif',
                input_rows=input_rows,
                output_rows=0)

            # Return results
            liner.close()
            return (outdf, stats)

        # Update Edge-Effect Length
        edgeeffectlength = ut.get_edgeeffectlength(
            maxreplen=maxreplen,
            exmotifs=exmotifs)

    # Parsing Sequence Constraint Feasibility
    liner.send('\n[Step 3: Parsing Motif Sequence]\n')

    # Update Step 2 Warning
    warns[3] = {
        'warn_count': 0,
        'step_name' : 'parsing-motif-sequence',
        'vars': None}

    # Parse primerseq
    (optrequired,
    homology,
    fixedbaseindex,
    exmotifindex) = cm.get_parsed_sequence_constraint(
        motifseq=motifseq,
        exmotifs=exmotifs,
        warn=warns[3],
        liner=liner)

    # Remove Step 3 Warning
    if not warns[3]['warn_count']:
        warns.pop(3)

    # Parse Edge Effects
    if ((not optrequired is False) and \
        (not exmotifs    is None)) and \
       ((not leftcontext  is None) or \
        (not rightcontext is None)):

        # Set Context Flag
        has_context = True

        # Show Update
        liner.send('\n[Step 4: Extracting Context Sequences]\n')

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
        liner.send('\n[Step 5: Parsing Edge Effects]\n')

        # Update Step 5 Warning
        warns[5] = {
            'warn_count': 0,
            'step_name' : 'parsing-edge-effects',
            'vars': None}

        # Compute Forbidden Prefixes and Suffixes
        (prefixdict,
        suffixdict) = cm.get_parsed_edgeeffects(
            motifseq=motifseq,
            motiftype=motiftype,
            leftcontext=leftcontext,
            rightcontext=rightcontext,
            leftpartition=leftpartition,
            rightpartition=rightpartition,
            exmotifs=exmotifs,
            element='Motif',
            warn=warns[5],
            liner=liner)

        # Remove Step 5 Warning
        if not warns[5]['warn_count']:
            warns.pop(5)

    # Finalize Context
    if not has_context:
        (leftcontext,
        rightcontext,
        prefixdict,
        suffixdict) = (None, None, None, None)

    # Parse Oligopool Repeats
    liner.send('\n[Step 6: Parsing Oligopool Repeats]\n')

    # Parse Repeats from indf
    (parsestatus,
    sourcecontext,
    kmerspace,
    fillcount,
    freecount,
    oligorepeats) = ut.get_parsed_oligopool_repeats(
        df=indf if (missing_mask is None or motiftype == 1) else indf.loc[missing_mask],
        maxreplen=maxreplen,
        element='Motif',
        merge=motiftype == 1,
        liner=liner)

    # Repeat Length infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status': False,
            'basis' : 'infeasible',
            'step'  : 6,
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
            module='motif',
            input_rows=input_rows,
            output_rows=0)

        # Return results
        liner.close()
        return (outdf, stats)

    # Launching Motif Design
    liner.send('\n[Step 7: Computing Motifs]\n')

    # Define Motif Design Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 7,
        'step_name': 'computing-motifs',
        'vars'    : {
               'target_count': targetcount,   # Required Number of Motifs
                'motif_count': 0,             # Motif Design Count
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

    # Design Motifs
    (motifs,
    stats) = cm.motif_engine(
        motifseq=motifseq,
        motiftype=motiftype,
        homology=homology,
        optrequired=optrequired,
        fixedbaseindex=fixedbaseindex,
        maxreplen=maxreplen,
        oligorepeats=oligorepeats,
        leftcontext=leftcontext,
        rightcontext=rightcontext,
        exmotifs=exmotifs,
        exmotifindex=exmotifindex,
        edgeeffectlength=edgeeffectlength,
        prefixdict=prefixdict,
        suffixdict=suffixdict,
        targetcount=targetcount,
        stats=stats,
        liner=liner,
        rng=rng)

    # Motif Status
    if stats['status']:
        motifstatus = 'Successful'
    else:
        motifstatus = 'Failed'

    # Insert motif into indf
    if stats['status']:

        # Update indf
        if motifcol_exists:
            indf.loc[missing_mask, motifcol] = motifs
        else:
            ut.update_df(
                indf=indf,
                lcname=leftcontextname,
                rcname=rightcontextname,
                out=motifs,
                outcol=motifcol)

        # Prepare outdf
        outdf = indf

        # Write outdf to file
        if not outfile is None:
            ut.write_df_csv(
                df=outdf,
                outfile=outfile,
                sep=',')

    # Motif Design Statistics
    liner.send('\n[Motif Design Statistics]\n')

    plen = ut.get_printlen(
        value=max(stats['vars'][field] for field in (
            'target_count',
            'motif_count')))

    liner.send(
        '  Design Status   : {}\n'.format(
            motifstatus))
    liner.send(
        '  Target Count    : {:{},d} Motif(s)\n'.format(
            stats['vars']['target_count'],
            plen))
    liner.send(
        '   Motif Count    : {:{},d} Motif(s) ({:6.2f} %)\n'.format(
            stats['vars']['motif_count'],
            plen,
            ut.safediv(
                A=stats['vars']['motif_count'] * 100.,
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
    if motifstatus == 'Successful':
        ae.unregister(ofdeletion)

    # Close Liner
    liner.close()

    # Return Solution and Statistics
    stats = ut.stamp_stats(
        stats=stats,
        module='motif',
        input_rows=input_rows,
        output_rows=len(outdf.index) if outdf is not None else 0)
    outdf_return = outdf
    if (outdf is not None) and (not id_from_index):
        outdf_return = ut.get_df_with_id_column(outdf)
    return (outdf_return, stats)
