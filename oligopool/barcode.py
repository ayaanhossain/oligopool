import time as tt

import collections as cx
import atexit as ae

import numpy  as np
import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp
from .base import core_barcode as cb
from .base import vectordb as db

from typing import Tuple


def barcode(
    input_data:str|pd.DataFrame,
    oligo_length_limit:int,
    barcode_length:int,
    minimum_hamming_distance:int,
    maximum_repeat_length:int,
    barcode_column:str,
    output_file:str|None=None,
    barcode_type:int|str=0,
    left_context_column:str|None=None,
    right_context_column:str|None=None,
    patch_mode:bool=False,
    cross_barcode_columns:list|None=None,
    minimum_cross_distance:int|None=None,
    excluded_motifs:list|str|dict|pd.DataFrame|None=None,
    background_directory:str|list|None=None,
    random_seed:int|None=None,
    verbose:bool=True) -> Tuple[pd.DataFrame, dict]:
    '''
    Design per-variant barcodes under Hamming-distance, repeat, and excluded-motif constraints, with
    optional context screening and cross-set separation against existing barcode columns.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated oligo pool variants.
        - `oligo_length_limit` (`int`): Maximum allowed oligo length (>= 4).
        - `barcode_length` (`int`): Length of the designed barcodes (>= 4).
        - `minimum_hamming_distance` (`int`): Minimum pairwise Hamming distance (>= 1).
        - `maximum_repeat_length` (`int`): Max shared repeat length with oligos (>= 4).
        - `barcode_column` (`str`): Column name for the designed barcodes.

    Optional Parameters:
        - `output_file` (`str` / `None`): Filename for output DataFrame; required in CLI usage,
            optional in library usage (default: `None`).
        - `barcode_type` (`int` / `str`): Barcode design mode (default: 0). See Notes.
        - `left_context_column` (`str` / `None`): Column for left DNA context (default: `None`).
        - `right_context_column` (`str` / `None`): Column for right DNA context (default: `None`).
        - `patch_mode` (`bool`): If `True`, fill only missing values in an existing barcode column
            (does not overwrite existing barcodes). (default: `False`).
        - `cross_barcode_columns` (`list` / `None`): Existing barcode column(s) used as a cross-set constraint (default: `None`).
        - `minimum_cross_distance` (`int` / `None`): Minimum cross-set Hamming distance (default: `None`).
        - `excluded_motifs` (`list` / `str` / `dict` / `pd.DataFrame` / `None`): Excluded motif source(s) (default: `None`).
        - `background_directory` (`str` / `list` / `None`): Background k-mer DB(s) from `background()` (default: `None`).
        - `random_seed` (`int` / `None`): Seed for local RNG (default: `None`).
        - `verbose` (`bool`): If `True`, logs progress to stdout (default: `True`).

    Returns:
        - A pandas DataFrame of generated barcodes; saves to `output_file` if specified.
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` must contain a unique 'ID' column; all other columns must be non-empty DNA strings.
            Column names must be unique and exclude `barcode_column`.
        - `barcode_type`:
            0 or 'terminus' = fast terminus optimized, 1 or 'spectrum' = slow spectrum optimized
            (aliases: 'term', 'spec').
        - Terminus optimization targets distinctive 5'/3' ends; spectrum optimization targets k-mer saturation.
        - At least one of `left_context_column` or `right_context_column` must be specified.
        - `excluded_motifs`: one or more sources (list/dict), merged; strict ATGC only; CSV/DataFrame requires 'Exmotif';
            FASTA sources are supported.
        - `background_directory` screens against one or more background k-mer DB(s); designs avoid k-mers in ALL specified DBs.
        - If design is challenging, adjust `barcode_length`, `minimum_hamming_distance`, `maximum_repeat_length`,
            and/or `excluded_motifs`, or switch `barcode_type`.
        - Constant anchors (e.g., index prefix/suffix) are typically designed first (see `motif`, `motif_type=1`).
        - Cross-set mode is global (not per-row): `cross_barcode_columns` and `minimum_cross_distance` must be set
            together; each new barcode must be >= `minimum_cross_distance` mismatches away from every strict-ATGC
            barcode in the union of sequences across `cross_barcode_columns` (length must equal `barcode_length`).
        - Patch mode (`patch_mode=True`) preserves existing values in `barcode_column` and fills only missing values
            (`None`/NaN/empty/`'-'`); existing values must already be strict ATGC of length `barcode_length`.
    '''

    # Preserve return style when the caller intentionally used ID as index.
    id_from_index = ut.get_id_index_intent(input_data)

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
    patch_mode   = patch_mode
    cross_cols   = cross_barcode_columns
    cross_mind   = minimum_cross_distance
    exmotifs     = excluded_motifs
    background   = background_directory
    random_seed  = random_seed
    verbose      = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Barcoding Verbiage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Barcode]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # Patch Mode (pre-parse for output column handling)
    (patch_mode_on,
    patch_mode_valid) = vp.get_parsed_flag_info(
        flag=patch_mode,
        flag_field='      Patch Mode    ',
        flag_desc_off='Disabled (Design All Barcodes)',
        flag_desc_on='Enabled (Fill Missing Barcodes)',
        liner=liner,
        precheck=True)
    allow_missing_cols = None
    if patch_mode_on and isinstance(barcodecol, str):
        allow_missing_cols = set((barcodecol,))

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

    # Full barcodelen Validation
    barcodelen_valid = vp.get_numeric_validity(
        numeric=barcodelen,
        numeric_field='    Barcode Length  ',
        numeric_pre_desc=' Exactly ',
        numeric_post_desc=' Base Pair(s)',
        minval=4,
        maxval=float('inf') if not oligolimit_valid else oligolimit,
        precheck=False,
        liner=liner)

    # Full minhdist Validation
    minhdist_valid = vp.get_numeric_validity(
        numeric=minhdist,
        numeric_field='    Hamming Distance',
        numeric_pre_desc=' At least ',
        numeric_post_desc=' Mismatch(es) per Barcode Pair',
        minval=1,
        maxval=barcodelen if barcodelen_valid else float('inf'),
        precheck=False,
        liner=liner)

    # Full maxreplen Validation
    maxreplen_valid = vp.get_numeric_validity(
        numeric=maxreplen,
        numeric_field='     Repeat Length  ',
        numeric_pre_desc=' Up to ',
        numeric_post_desc=' Base Pair(s) Oligo Pool Repeats',
        minval=4,
        maxval=barcodelen if barcodelen_valid else float('inf'),
        precheck=False,
        liner=liner)

    # Full outcol Validation
    barcodecol_valid = vp.get_parsed_column_info(
        col=barcodecol,
        df=indf,
        col_field='    Barcode Column  ',
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
        outdf_suffix='.oligopool.barcode.csv',
        outdf_field='     Output File    ',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.barcode.csv')

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full barcodetype Validation
    (barcodetype,
    barcodetype_valid) = vp.get_typed_categorical_validity(
        category=barcodetype,
        category_field='    Barcode Type    ',
        category_pre_desc=' ',
        category_post_desc=' Barcodes',
        type_name='barcode_type',
        liner=liner)

    # Store Context Names
    leftcontextname  = leftcontext
    rightcontextname = rightcontext

    # Normalize crosscols input early so skip_cols handling treats
    # string and iterable forms equivalently during context parsing.
    if isinstance(cross_cols, str):
        cross_cols = [cross_cols]

    # Build skip_cols for adjacency: skip over the barcode column (and any
    # cross-barcode columns) when checking whether left/right contexts are
    # adjacent — these designed columns may sit between the contexts in the
    # user's DataFrame without breaking the physical oligo layout.
    _skip = set()
    if barcodecol and isinstance(barcodecol, str) and barcodecol in indf.columns:
        _skip.add(barcodecol)
    if hasattr(cross_cols, '__iter__') and not isinstance(cross_cols, str):
        for _cc in cross_cols:
            if isinstance(_cc, str) and _cc in indf.columns:
                _skip.add(_cc)
    _skip_cols = list(_skip) if _skip else None

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
        liner=liner,
        skip_cols=_skip_cols)

    # Full rightcontext Parsing and Validation
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
        liner=liner,
        skip_cols=_skip_cols)

    # Patch mode parsing (printed after context args; evaluated earlier for indata handling)
    vp.get_parsed_flag_info(
        flag=patch_mode,
        flag_field='      Patch Mode    ',
        flag_desc_off='Disabled (Design All Barcodes)',
        flag_desc_on='Enabled (Fill Missing Barcodes)',
        liner=liner,
        flag_on=patch_mode_on,
        flag_valid=patch_mode_valid)

    # Full crosscols Parsing and Validation
    (cross_cols,
    cross_mind,
    cross_valid) = vp.get_parsed_cross_barcode_info(
        crosscols=cross_cols,
        crosscols_field='      Cross Barcodes',
        mindist=cross_mind,
        mindist_field='      Cross Distance',
        barcodelen=barcodelen if barcodelen_valid else 0,
        outcol=barcodecol,
        df=indf,
        liner=liner)

    # Full exmotifs Parsing and Validation
    (exmotifs,
    exmotifs_valid,
    exmotif_inputs) = vp.get_parsed_exmotifs(
        exmotifs=exmotifs,
        exmotifs_field='   Excluded Motifs  ',
        liner=liner)

    # Full background Parsing and Validation
    (background_valid,
    backgrounds_info,
    max_bg_K) = vp.get_parsed_backgrounds(
        backgrounds=background,
        backgrounds_field=' Background Database',
        liner=liner)

    # Validate random_seed (do not auto-generate)
    (random_seed,
    seed_valid) = vp.get_parsed_random_seed_info(
        random_seed=random_seed,
        random_seed_field='     Random Seed    ',
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
        patch_mode_valid,
        cross_valid,
        exmotifs_valid,
        background_valid,
        seed_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Local RNG
    rng = np.random.default_rng(random_seed)

    # Open Backgrounds
    backgrounds = []
    opened_backgrounds = []  # Track for cleanup
    if backgrounds_info:
        for bg_ref, bg_type, _ in backgrounds_info:
            if bg_type == 'path':
                bg_path = ut.get_adjusted_path(
                    path=bg_ref,
                    suffix='.oligopool.background')
                bg = db.vectorDB(
                    path=bg_path,
                    maximum_repeat_length=None)
                opened_backgrounds.append(bg)
            else:  # 'instance'
                bg = bg_ref
            backgrounds.append(bg)

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Paramters
    oligolimit = round(oligolimit)
    barcodelen = round(barcodelen)
    minhdist   = round(minhdist)
    maxreplen  = round(maxreplen)
    if not cross_mind is None:
        cross_mind = round(cross_mind)

    # Patch mode bookkeeping
    barcodecol_exists = False
    missing_mask = None
    existing_mask = None
    targetcount = len(indf.index)
    if patch_mode_on and isinstance(indf, pd.DataFrame):
        (barcodecol_exists,
        _) = ut.get_col_exist_idx(
            col=barcodecol,
            df=indf)
        if barcodecol_exists:
            missing_mask = ut.get_missing_mask(
                series=indf[barcodecol],
                allow_dash=True)
            indf[barcodecol] = ut.fill_missing_values(
                series=indf[barcodecol],
                missing_mask=missing_mask,
                fill='-')
            existing_mask = ~missing_mask
            targetcount = int(missing_mask.sum())
        else:
            missing_mask = np.ones(
                len(indf.index),
                dtype=bool)
            existing_mask = ~missing_mask

        # Validate existing barcodes (length + DNA)
        if barcodecol_exists and existing_mask.any():
            invalid_mask = np.zeros(
                len(indf.index),
                dtype=bool)
            for idx in np.where(existing_mask)[0]:
                value = indf[barcodecol].iat[idx]
                if (not ut.is_DNA(
                        seq=value,
                        dna_alpha=set('ATGC'))) or \
                   (len(value) != barcodelen):
                    invalid_mask[idx] = True
            if invalid_mask.any():
                examples = ut.get_row_examples(
                    df=indf,
                    invalid_mask=invalid_mask,
                    id_col='ID',
                    limit=5)
                example_note = ut.format_row_examples(examples)
                liner.send(
                    '    Barcode Column  : Output in Column \'{}\' [INVALID EXISTING VALUE]{}\n'.format(
                        barcodecol,
                        example_note))
                liner.send('\n')
                raise RuntimeError(
                    'Invalid Argument Input(s).')

    # Define Edge Effect Length
    edgeeffectlength = None

    # Barcode Design Book-keeping
    has_context = False
    outdf = None
    stats = None
    warns = {}

    # No missing barcodes in patch mode
    if patch_mode_on and \
       missing_mask is not None and \
       not missing_mask.any():

        stats = {
            'status'  : True,
            'basis'   : 'complete',
            'step'    : 0,
            'step_name': 'no-missing-barcodes',
            'vars'    : {
                  'target_count': 0,
                 'barcode_count': 0,
                 'orphan_oligo': [],
                     'type_fail': 0,
                 'distance_fail': 0,
           'existing_distance_fail': 0,
           'cross_distance_fail': 0,
                   'repeat_fail': 0,
              'background_fail': 0,
                  'exmotif_fail': 0,
                     'edge_fail': 0,
               'distance_distro': None,
               'exmotif_counter': cx.Counter(),
               'space_exhausted': False,
               'trial_exhausted': False,
            },
            'warns'   : warns}
        stats['random_seed'] = random_seed
        if not cross_cols is None:
            stats['vars']['cross_barcode_columns'] = tuple(cross_cols)
            stats['vars']['minimum_cross_distance'] = cross_mind
            stats['vars']['cross_set_size'] = 0

        outdf = indf
        if not outfile is None:
            ut.write_df_csv(
                df=outdf,
                outfile=outfile,
                sep=',')

        liner.send('\n[Barcode Design Statistics]\n')
        liner.send(
            '   Design Status   : Successful\n')
        liner.send(
            '   Target Count    : 0 Barcode(s)\n')
        liner.send(
            '  Barcode Count    : 0 Barcode(s) (  0.00 %)\n')
        liner.send(
            '   Orphan Oligo    : 0 Entries\n')
        liner.send(
            ' Time Elapsed: {:.2f} sec\n'.format(
                tt.time()-t0))

        liner.close()
        stats = ut.stamp_stats(
            stats=stats,
            module='barcode',
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
            module='barcode',
            input_rows=input_rows,
            output_rows=0)

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
        indf=indf if missing_mask is None else indf.loc[missing_mask],
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
        stats['random_seed'] = random_seed
        stats = ut.stamp_stats(
            stats=stats,
            module='barcode',
            input_rows=input_rows,
            output_rows=0)

        # Return results
        liner.close()
        return (outdf, stats)

    # Parse Excluded Motifs
    if ((not exmotifs is None) and
        ((not leftcontext  is None) or \
         (not rightcontext is None))):

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
            stats['random_seed'] = random_seed
            stats = ut.stamp_stats(
                stats=stats,
                module='barcode',
                input_rows=input_rows,
                output_rows=0)

            # Return results
            liner.close()
            return (outdf, stats)

        # Update Edge-Effect Length
        edgeeffectlength = ut.get_edgeeffectlength(
            maxreplen=maxreplen,
            exmotifs=exmotifs)

    # Re-calculate Edge-Effect Length
    if edgeeffectlength is None:
        edgeeffectlength = maxreplen
    else:
        edgeeffectlength = max(
            edgeeffectlength,
            maxreplen)

    # Update Edge-Effect Length for Background
    if backgrounds:
        edgeeffectlength = max(edgeeffectlength or 0, max_bg_K)

    # Extract Left and Right Context
    if ((not leftcontext  is None) or \
        (not rightcontext is None)):

        # Set Context Flag
        has_context = True

        # Show update
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

    # Finalize Context
    if not has_context:
        leftcontext,rightcontext = None, None
    elif missing_mask is not None:
        if not leftcontext is None:
            if len(leftcontext) == len(missing_mask):
                leftcontext = [
                    seq for seq, keep in zip(leftcontext, missing_mask)
                    if keep]
        if not rightcontext is None:
            if len(rightcontext) == len(missing_mask):
                rightcontext = [
                    seq for seq, keep in zip(rightcontext, missing_mask)
                    if keep]

    # Parse Oligo Pool Repeats
    liner.send('\n[Step 5: Parsing Oligo Pool Repeats]\n')

    # Parse Repeats from indf
    (parsestatus,
    sourcecontext,
    kmerspace,
    fillcount,
    freecount,
    oligorepeats) = ut.get_parsed_oligopool_repeats(
        df=indf if missing_mask is None else indf.loc[missing_mask],
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
            module='barcode',
            input_rows=input_rows,
            output_rows=0)

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
                     'target_count': targetcount,  # Required Number of Barcodes
                    'barcode_count': 0,            # Barcode Design Count
                     'orphan_oligo': None,         # Orphan Oligo Indexes
                        'type_fail': 0,            # Barcode Type Failure Count
                    'distance_fail': 0,            # Hamming Distance Fail Count
           'existing_distance_fail': 0,        # Existing Barcode Distance Fail Count
              'cross_distance_fail': 0,        # Cross Distance Fail Count
                      'repeat_fail': 0,            # Repeat Fail Count
                  'background_fail': 0,            # Background Fail Count
                     'exmotif_fail': 0,            # Exmotif Elimination Fail Count
                        'edge_fail': 0,            # Edge Effect Fail Count
                  'distance_distro': None,         # Hamming Distance Distribution
                  'exmotif_counter': cx.Counter(), # Exmotif Encounter Counter
                  'space_exhausted': False,        # Space Exhausted Bool
                  'trial_exhausted': False,        # Trial Exhausted Bool
            },
        'warns'   : warns}
    stats['random_seed'] = random_seed
    if not cross_cols is None:
        stats['vars']['cross_barcode_columns'] = tuple(cross_cols)
        stats['vars']['minimum_cross_distance'] = cross_mind

    # Schedule outfile deletion
    ofdeletion = ae.register(
        ut.remove_file,
        outfile)

    # Existing barcode constraints setup (patch mode; enforce distance vs existing)
    existing_store_plus = None
    existing_contigsize = None
    existing_coocache = None
    existing_set_size = 0
    if patch_mode_on and \
       barcodecol_exists and \
       existing_mask is not None and \
       existing_mask.any():
        liner.send(' Preparing Existing Barcodes ...\r')
        (existing_store_plus,
        existing_set_size,
        existing_contigsize,
        existing_coocache) = ut.get_cross_barcode_store(
            df=indf.loc[existing_mask, [barcodecol]],
            crosscols=[barcodecol],
            barcodelen=barcodelen,
            minhdist=minhdist)
        liner.send(
            ' Existing Set Size: {:,} Unique Barcode(s)\n'.format(
                existing_set_size))

    # Cross-set constraints setup (if enabled)
    cross_store_plus = None
    cross_contigsize = None
    cross_coocache = None
    cross_set_size = 0

    if not cross_cols is None:
        liner.send(' Preparing Cross Barcodes ...\r')
        (cross_store_plus,
        cross_set_size,
        cross_contigsize,
        cross_coocache) = ut.get_cross_barcode_store(
            df=indf,
            crosscols=cross_cols,
            barcodelen=barcodelen,
            minhdist=cross_mind)
        stats['vars']['cross_set_size'] = cross_set_size
        liner.send(
            '    Cross Set Size: {:,} Unique Barcode(s)\n'.format(
                cross_set_size))

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
        cross_store_plus=cross_store_plus,
        cross_set_size=cross_set_size,
        cross_contigsize=cross_contigsize,
        cross_coocache=cross_coocache,
        minimum_cross_distance=cross_mind,
        existing_store_plus=existing_store_plus,
        existing_set_size=existing_set_size,
        existing_contigsize=existing_contigsize,
        existing_coocache=existing_coocache,
        existing_min_distance=minhdist,
        exmotifs=exmotifs,
        background=backgrounds if backgrounds else None,
        targetcount=targetcount,
        stats=stats,
        liner=liner,
        rng=rng,)

    # Success Relevant Stats
    if not store is None:

        # Launching Distance Distribution Analysis
        liner.send('\n[Step 7: Computing Distance Distribution]\n')

        # Compute Hamming Distance Distribution
        stats['vars']['distance_distro'] = cb.get_distro(
            store=store,
            liner=liner,
            rng=rng)

    # Barcode Status
    if stats['status']:
        barcodestatus = 'Successful'
    else:
        barcodestatus = 'Failed'

    # Insert codes into indf
    if stats['status']:

        # Update indf
        if barcodecol_exists:
            indf.loc[missing_mask, barcodecol] = codes
        else:
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
            ut.write_df_csv(
                df=outdf,
                outfile=outfile,
                sep=',')

    # Barcoding Statistics
    liner.send('\n[Barcode Design Statistics]\n')

    plen = ut.get_printlen(
        value=max(stats['vars'][field] for field in (
            'target_count',
            'barcode_count')))

    liner.send(
        '     Design Status   : {}\n'.format(
            barcodestatus))
    liner.send(
        '     Target Count    : {:{},d} Barcode(s)\n'.format(
            stats['vars']['target_count'],
            plen))
    liner.send(
        '    Barcode Count    : {:{},d} Barcode(s) ({:6.2f} %)\n'.format(
            stats['vars']['barcode_count'],
            plen,
            ut.safediv(
                A=stats['vars']['barcode_count'] * 100.,
                B=targetcount)))
    liner.send(
        '     Orphan Oligo    : {:{},d} Entries\n'.format(
            len(stats['vars']['orphan_oligo']),
            plen))

    # Success Relevant Stats
    if stats['status']:
        if stats['vars']['distance_distro']:

            dlen = ut.get_printlen(
                value=max(map(
                    lambda x: x[1],
                    stats['vars']['distance_distro'])))

            liner.send(' Pair-wise Distance Distribution\n')

            for percentage,distance in stats['vars']['distance_distro']:
                liner.send(
                    '   - {:6.2f} % Barcode(s) w/ Distance >= {:{},d} Mismatch(es)\n'.format(
                        percentage,
                        distance,
                        dlen))

    # Failure Relavant Stats
    else:
        maxval = max(stats['vars'][field] for field in (
            'distance_fail',
            'existing_distance_fail',
            'cross_distance_fail',
            'repeat_fail',
            'background_fail',
            'exmotif_fail',
            'edge_fail'))

        sntn, plen = ut.get_notelen(
            printlen=ut.get_printlen(
                value=maxval))

        total_conflicts = stats['vars']['distance_fail']       + \
                          stats['vars']['existing_distance_fail'] + \
                          stats['vars']['cross_distance_fail'] + \
                          stats['vars']['repeat_fail']         + \
                          stats['vars']['background_fail']     + \
                          stats['vars']['exmotif_fail']        + \
                          stats['vars']['edge_fail']
        liner.send(
            '   Distance Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['distance_fail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['distance_fail'] * 100.,
                    B=total_conflicts)))
        if patch_mode_on and \
           existing_mask is not None and \
           existing_mask.any():
            liner.send(
                '   Existing Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                    stats['vars']['existing_distance_fail'],
                    plen,
                    sntn,
                    ut.safediv(
                        A=stats['vars']['existing_distance_fail'] * 100.,
                        B=total_conflicts)))
        if not cross_cols is None:
            liner.send(
                '      Cross Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                    stats['vars']['cross_distance_fail'],
                    plen,
                    sntn,
                    ut.safediv(
                        A=stats['vars']['cross_distance_fail'] * 100.,
                        B=total_conflicts)))
        liner.send(
            '     Repeat Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['repeat_fail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['repeat_fail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            ' Background Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['background_fail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['background_fail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '    Exmotif Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['exmotif_fail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['exmotif_fail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '       Edge Conflicts: {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['edge_fail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['edge_fail'] * 100.,
                    B=total_conflicts)))
        liner.send(f'      Space Exhausted: {stats["vars"]["space_exhausted"]}\n')
        liner.send(f'      Trial Exhausted: {stats["vars"]["trial_exhausted"]}\n')

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

    # Close Backgrounds
    for bg in opened_backgrounds:
        bg.close()

    # Excluded-motif attribution for stats
    ut.stamp_exmotif_stats(
        stats=stats,
        exmotif_inputs=exmotif_inputs)

    # Close Liner
    liner.close()

    # Return Solution and Statistics
    stats = ut.stamp_stats(
        stats=stats,
        module='barcode',
        input_rows=input_rows,
        output_rows=len(outdf.index) if outdf is not None else 0)
    outdf_return = outdf
    if (outdf is not None) and (not id_from_index):
        outdf_return = ut.get_df_with_id_column(outdf)
    return (outdf_return, stats)
