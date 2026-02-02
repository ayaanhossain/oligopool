import time as tt

import collections as cx
import atexit as ae

import numpy as np
import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp
from .base import core_primer as cp
from .base import vectordb as db

from typing import Tuple


def primer(
    input_data:str|pd.DataFrame,
    oligo_length_limit:int,
    primer_sequence_constraint:str,
    primer_type:int|str,
    minimum_melting_temperature:float,
    maximum_melting_temperature:float,
    maximum_repeat_length:int,
    primer_column:str,
    output_file:str|None=None,
    left_context_column:str|None=None,
    right_context_column:str|None=None,
    patch_mode:bool=False,
    oligo_sets:list|str|pd.DataFrame|None=None,
    paired_primer_column:str|None=None,
    excluded_motifs:list|str|pd.DataFrame|None=None,
    background_directory:str|None=None,
    random_seed:int|None=None,
    verbose:bool=True) -> Tuple[pd.DataFrame, dict]:
    '''
    Design constrained primers under an IUPAC sequence constraint with Tm/repeat/dimer screening.
    Supports background screening, chained primer design via Tm matching, multiplexed per-set primers
    via `oligo_sets`, and patch mode for incremental pool extension.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated oligo pool variants.
        - `oligo_length_limit` (`int`): Maximum allowed oligo length (≥ 4).
        - `primer_sequence_constraint` (`str`): IUPAC degenerate sequence constraint.
        - `primer_type` (`int` / `str`): Primer direction/type (default: required). See Notes.
        - `minimum_melting_temperature` (`float`): Minimum primer Tm (≥ 25°C).
        - `maximum_melting_temperature` (`float`): Maximum primer Tm (≤ 95°C).
        - `maximum_repeat_length` (`int`): Max shared repeat length with oligos (≥ 6).
        - `primer_column` (`str`): Column name for the designed primer.

    Optional Parameters:
        - `output_file` (`str`): Filename for output DataFrame; required in CLI usage,
            optional in library usage (default: `None`).
        - `left_context_column` (`str`): Column for left DNA context (default: `None`).
        - `right_context_column` (`str`): Column for right DNA context (default: `None`).
        - `patch_mode` (`bool`): If `True`, fill only missing values in an existing primer column
            (does not overwrite existing primers). (default: `False`).
        - `oligo_sets` (`list` / `str` / `pd.DataFrame`): Per-oligo grouping labels to design
            set-specific primers; can be a list aligned to `input_data` or a CSV/DataFrame
            with 'ID' and 'OligoSet' columns (default: `None`).
        - `paired_primer_column` (`str`): Column for paired primer sequence (default: `None`).
        - `excluded_motifs` (`list` / `str` / `pd.DataFrame`): Motifs to exclude;
            can be a list, CSV, DataFrame, or FASTA file (default: `None`).
        - `background_directory` (`str`): Directory for background k-mer sequences (default: `None`).
        - `random_seed` (`int` / `None`): Seed for local RNG (default: `None`).
        - `verbose` (`bool`): If `True`, logs updates to stdout (default: `True`).

    Returns:
        - A pandas DataFrame of designed primers; saves to `output_file` if specified.
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` must contain a unique 'ID' column; all other columns must be non-empty DNA strings.
            Column names must be unique and exclude `primer_column`.
        - At least one of `left_context_column` or `right_context_column` must be specified.
        - `primer_type`:
            0 or 'forward' for forward, 1 or 'reverse' for reverse (aliases: 'fwd', 'f', 'rev', 'r').
        - The paired primer type is inferred based on the current primer type.
        - If `paired_primer_column` is specified, Tm of the designed primer is optimized within 1°C of it.
        - `maximum_repeat_length` controls non-repetitiveness against `input_data` only; to screen against a
            background, build a background DB with `background(...)` and pass it via `background_directory`.
        - If `excluded_motifs` is a CSV or DataFrame, it must have an 'Exmotif' column.
        - `excluded_motifs` can be a list, CSV, DataFrame, or FASTA file.
        - Constant motifs in sequence constraint may lead to sub-optimal primers.
        - Chained primer design: design one primer, then its partner via `paired_primer_column`
            (pairing is inferred from `primer_type`).
        - `oligo_sets` can be any group labels. Primers are designed per set and screened for cross-set
            compatibility; if `paired_primer_column` is provided, it must be constant within each set and
            pairing/Tm matching is applied per set.
        - Patch mode (`patch_mode=True`) preserves existing values in `primer_column` and fills only missing values
            (`None`/NaN/empty/`'-'`). With `oligo_sets`, existing per-set primers are reused and missing-only sets
            trigger new per-set primer design.
    '''

    # Preserve return style when the caller intentionally used ID as index.
    id_from_index = ut.get_id_index_intent(input_data)

    # Argument Aliasing
    indata       = input_data
    oligolimit   = oligo_length_limit
    primerseq    = primer_sequence_constraint
    primertype   = primer_type
    mintmelt     = minimum_melting_temperature
    maxtmelt     = maximum_melting_temperature
    maxreplen    = maximum_repeat_length
    primercol    = primer_column
    outfile      = output_file
    pairedcol    = paired_primer_column
    leftcontext  = left_context_column
    rightcontext = right_context_column
    patch_mode   = patch_mode
    oligosets    = oligo_sets
    exmotifs     = excluded_motifs
    background   = background_directory
    random_seed  = random_seed
    verbose      = verbose

    # Local RNG
    rng = np.random.default_rng(random_seed)

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Primer Verbiage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Primer]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # Patch Mode (pre-parse for output column handling)
    (patch_mode_on,
    patch_mode_valid) = vp.get_parsed_flag_info(
        flag=patch_mode,
        flag_field='      Patch Mode       ',
        flag_desc_off='Disabled (Design All Primers)',
        flag_desc_on='Enabled (Fill Missing Primers)',
        liner=liner,
        precheck=True)
    allow_missing_cols = None
    if patch_mode_on and isinstance(primercol, str):
        allow_missing_cols = set((primercol,))

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_indata_info(
        indata=indata,
        indata_field='      Input Data       ',
        required_fields=('ID',),
        precheck=False,
        liner=liner,
        allow_missing_cols=allow_missing_cols)
    input_rows = len(indf.index) if isinstance(indf, pd.DataFrame) else 0

    # Full oligolimit Validation
    oligolimit_valid = vp.get_numeric_validity(
        numeric=oligolimit,
        numeric_field='      Oligo Limit      ',
        numeric_pre_desc=' At most ',
        numeric_post_desc=' Base Pair(s)',
        minval=4,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # First Pass primerseq Validation
    primerseq_valid = vp.get_seqconstr_validity(
        seqconstr=primerseq,
        seqconstr_field='     Primer Sequence   ',
        minlenval=15,
        element='PRIMER',
        liner=liner)

    # Full primertype Validation
    (primertype,
    primertype_valid) = vp.get_typed_categorical_validity(
        category=primertype,
        category_field='     Primer Type       ',
        category_pre_desc=' ',
        category_post_desc=' Primer Design',
        type_name='primer_type',
        liner=liner)

    # Full mintmelt and maxtmelt Validation
    (mintmelt,
    maxtmelt,
    tmelt_valid) = vp.get_parsed_range_info(
        minval=mintmelt,
        maxval=maxtmelt,
        range_field='    Melting Temperature',
        range_unit='°C',
        range_min=25,
        range_max=95,
        liner=liner)

    # Full maxreplen Validation
    maxreplen_valid = vp.get_numeric_validity(
        numeric=maxreplen,
        numeric_field='     Repeat Length     ',
        numeric_pre_desc=' Up to ',
        numeric_post_desc=' Base Pair(s) Oligo Pool Repeats',
        minval=6,
        maxval=len(primerseq) if primerseq_valid else float('inf'),
        precheck=False,
        liner=liner)

    # Full primercol Validation
    primercol_valid = vp.get_parsed_column_info(
        col=primercol,
        df=indf,
        col_field='     Primer Column     ',
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
        outdf_suffix='.oligopool.primer.csv',
        outdf_field='     Output File       ',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.primer.csv')

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Store Context Names
    leftcontextname  = leftcontext
    rightcontextname = rightcontext

    # Full leftcontext Parsing and Validation
    (leftcontext,
    leftcontext_valid) = vp.get_parsed_column_info(
        col=leftcontext,
        df=indf,
        col_field='       Left Context    ',
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
        col_field='      Right Context    ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=leftcontextname,
        adjval=-1,
        iscontext=True,
        typecontext=1,
        liner=liner)

    # Patch mode parsing (printed after context args; evaluated earlier for indata handling)
    vp.get_parsed_flag_info(
        flag=patch_mode,
        flag_field='      Patch Mode       ',
        flag_desc_off='Disabled (Design All Primers)',
        flag_desc_on='Enabled (Fill Missing Primers)',
        liner=liner,
        flag_on=patch_mode_on,
        flag_valid=patch_mode_valid)

    # Full oligosets Parsing and Validation
    (oligosets,
    oligosets_valid) = vp.get_parsed_oligosets_info(
        oligosets=oligosets,
        oligosets_field='      Oligo Sets       ',
        df_field='OligoSet',
        indf=indf,
        indata_valid=indata_valid,
        liner=liner)

    # Full pairedprimer Validation
    (pairedprimer,
    pairedprimer_valid) = vp.get_parsed_pairedprimer_info(
        pairedcol=pairedcol,
        pairedcol_field='     Paired Primer     ',
        df=indf,
        oligosets=oligosets,
        allow_missing=patch_mode_on,
        liner=liner)

    # Full exmotifs Parsing and Validation
    (exmotifs,
    exmotifs_valid) = vp.get_parsed_exseqs_info(
        exseqs=exmotifs,
        exseqs_field='   Excluded Motifs     ',
        exseqs_desc='Unique Motif(s)',
        df_field='Exmotif',
        required=False,
        liner=liner)

    # Full background Parsing and Validation
    (background_valid,
    background_type) = vp.get_parsed_background(
        background=background,
        background_field=' Background Database   ',
        liner=liner)

    # First Pass Validation
    if not all([
        indata_valid,
        oligolimit_valid,
        primerseq_valid,
        primertype_valid,
        tmelt_valid,
        maxreplen_valid,
        primercol_valid,
        outfile_valid,
        patch_mode_valid,
        oligosets_valid,
        pairedprimer_valid,
        leftcontext_valid,
        rightcontext_valid,
        exmotifs_valid,
        background_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Open Background
    if background_type == 'path':
        background = ut.get_adjusted_path(
            path=background,
            suffix='.oligopool.background')
        background = db.vectorDB(
            path=background,
            maximum_repeat_length=None)

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Paramters
    oligolimit = round(oligolimit)
    maxreplen  = round(maxreplen)

    # Patch mode bookkeeping
    primercol_exists = False
    missing_mask = None
    existing_mask = None
    design_mask = None
    targetcount = len(indf.index)

    if patch_mode_on and isinstance(indf, pd.DataFrame):
        (primercol_exists,
        _) = ut.get_col_exist_idx(
            col=primercol,
            df=indf)
        if primercol_exists:
            missing_mask = ut.get_missing_mask(
                series=indf[primercol],
                allow_dash=True)
            indf[primercol] = ut.fill_missing_values(
                series=indf[primercol],
                missing_mask=missing_mask,
                fill='-')
            existing_mask = ~missing_mask
            targetcount = int(missing_mask.sum())
        else:
            missing_mask = np.ones(
                len(indf.index),
                dtype=bool)
            existing_mask = ~missing_mask

        # Validate existing primers (length + DNA)
        if primercol_exists and existing_mask.any():
            invalid_mask = np.zeros(
                len(indf.index),
                dtype=bool)
            for idx in np.where(existing_mask)[0]:
                value = indf[primercol].iat[idx]
                if (not ut.is_DNA(
                        seq=value,
                        dna_alpha=set('ATGC'))) or \
                   (len(value) != len(primerseq)):
                    invalid_mask[idx] = True
            if invalid_mask.any():
                examples = ut.get_row_examples(
                    df=indf,
                    invalid_mask=invalid_mask,
                    id_col='ID',
                    limit=5)
                example_note = ut.format_row_examples(examples)
                liner.send(
                    '     Primer Column     : Output in Column \'{}\' [INVALID EXISTING VALUE]{}\n'.format(
                        primercol,
                        example_note))
                liner.send('\n')
                raise RuntimeError(
                    'Invalid Argument Input(s).')

    # Define Edge Effect Length
    edgeeffectlength = None

    # Primer Design Book-keeping
    outdf = None
    stats = None
    warns = {}
    pairedprimer_map = None
    set_labels = None
    set_groups = None
    set_sizes = None
    target_indf = indf
    target_oligosets = oligosets
    target_labels = None
    target_groups = None
    target_sizes = None
    existing_primers = {}

    # Define oligo set groups (if provided)
    if not oligosets is None:
        (set_labels,
        set_groups,
        set_sizes) = ut.get_oligoset_groups(
            oligosets=oligosets)
        target_labels = set_labels
        target_groups = set_groups
        target_sizes = set_sizes

    # Patch mode set handling
    if patch_mode_on and missing_mask is not None:

        # No oligo sets: reuse a constant primer if present
        if oligosets is None and primercol_exists and existing_mask.any():
            uniques = ut.get_uniques(
                iterable=indf[primercol].iloc[existing_mask],
                typer=tuple)
            if len(uniques) != 1:
                liner.send(
                    '     Primer Column     : Output in Column \'{}\' [NON-UNIQUE EXISTING VALUES]\n'.format(
                        primercol))
                liner.send('\n')
                raise RuntimeError(
                    'Invalid Argument Input(s).')

            indf.loc[missing_mask, primercol] = uniques[0]
            design_mask = np.zeros(
                len(indf.index),
                dtype=bool)
        elif oligosets is None:
            design_mask = missing_mask
            if design_mask.any():
                target_indf = indf.loc[design_mask]

        # Oligo sets: reuse existing primers per set, design for empty sets
        else:
            design_labels = []
            for label in set_labels:
                idx = np.array(set_groups[label])
                set_missing = missing_mask[idx]
                set_existing = ~set_missing

                # Existing primers in this set
                if set_existing.any():
                    values = indf[primercol].iloc[idx][set_existing]
                    uniques = ut.get_uniques(
                        iterable=values,
                        typer=tuple)
                    if len(uniques) != 1:
                        liner.send(
                            '     Primer Column     : Output in Column \'{}\' [NON-UNIQUE WITHIN SET=\'{}\']\n'.format(
                                primercol,
                                label))
                        liner.send('\n')
                        raise RuntimeError(
                            'Invalid Argument Input(s).')

                    existing_primers[label] = uniques[0]
                    if set_missing.any():
                        indf.loc[indf.index[idx][set_missing], primercol] = uniques[0]

                # Missing-only set
                else:
                    existing_primers[label] = None
                    design_labels.append(label)

            design_mask = np.zeros(
                len(indf.index),
                dtype=bool)
            for label in design_labels:
                design_mask[set_groups[label]] = True

            if design_labels:
                target_indf = indf.loc[design_mask]
                target_oligosets = oligosets[design_mask]
                (target_labels,
                target_groups,
                target_sizes) = ut.get_oligoset_groups(
                    oligosets=target_oligosets)

    # Patch mode with no design targets
    if patch_mode_on and \
       design_mask is not None and \
       not design_mask.any():

        stats = {
            'status'  : True,
            'basis'   : 'complete',
            'step'    : 0,
            'step_name': 'no-missing-primers',
            'vars'    : {
                   'primer_Tm': None,          # Primer Melting Temperature
                   'primer_GC': None,          # Primer GC Content
                 'hairpin_MFE': None,          # Primer Hairpin Free Energy
               'homodimer_MFE': None,          # Homodimer Free Energy
             'heterodimer_MFE': None,          # Heterodimer Free Energy
                     'Tm_fail': 0,             # Melting Temperature Fail Count
                 'repeat_fail': 0,             # Repeat Fail Count
             'background_fail': 0,             # Background Fail Count
              'homodimer_fail': 0,             # Homodimer Fail Count
            'heterodimer_fail': 0,             # Heterodimer Fail Count
             'crossdimer_fail': 0,             # Cross-primer Dimer Fail Count
                'exmotif_fail': 0,             # Exmotif Elimination Fail Count
                   'edge_fail': 0,             # Edge Effect Fail Count
             'exmotif_counter': cx.Counter()}, # Exmotif Encounter Counter
            'warns'   : warns}
        stats['random_seed'] = random_seed
        if not oligosets is None:
            stats['vars']['oligo_set_count'] = len(set_labels)
            stats['vars']['oligo_set_sizes'] = [
                {'oligo_set': label, 'count': set_sizes[label]}
                for label in set_labels]
            stats['vars']['primer_sets'] = []

        outdf = indf
        if not outfile is None:
            ut.write_df_csv(
                df=outdf,
                outfile=outfile,
                sep=',')

        liner.close()
        stats = ut.stamp_stats(
            stats=stats,
            module='primer',
            input_rows=input_rows,
            output_rows=len(outdf.index))
        outdf_return = outdf
        if not id_from_index:
            outdf_return = ut.get_df_with_id_column(outdf)
        return (outdf_return, stats)

    # Store Full Contexts
    leftcontext_full = leftcontext
    rightcontext_full = rightcontext

    # Track per-set paired primer usage
    if isinstance(pairedprimer, dict):
        pairedprimer_map = pairedprimer
        pairedprimer = None

    # Parse Oligopool Limit Feasibility
    liner.send('\n[Step 1: Parsing Oligo Limit]\n')

    # Parse oligolimit
    variantlens = None
    if patch_mode_on and design_mask is not None:
        variantlens = ut.get_variantlens(
            indf=target_indf)

    (parsestatus,
    minoligolen,
    maxoligolen,
    minelementlen,
    maxelementlen,
    minspaceavail,
    maxspaceavail) = ut.get_parsed_oligolimit(
        indf=indf,
        variantlens=variantlens,
        oligolimit=oligolimit,
        minelementlen=len(primerseq),
        maxelementlen=len(primerseq),
        element='Primer',
        liner=liner)

    # oligolimit infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'   : False,
            'basis'    : 'infeasible',
            'step'     : 1,
            'step_name': 'parsing-oligo-limit',
            'vars'     : {
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
            module='primer',
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
            'step_name': 'parsing-excluded-motifs',
            'vars': None}

        # Parse exmotifs
        (parsestatus,
        exmotifs,
        problens,
        leftpartition,
        rightpartition) = ut.get_parsed_exmotifs(
            exmotifs=exmotifs,
            typer=tuple,
            element='Primer',
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
                module='primer',
                input_rows=input_rows,
                output_rows=0)

            # Return results
            liner.close()
            return (outdf, stats)

        # Update Edge-Effect Length
        edgeeffectlength = ut.get_edgeeffectlength(
            maxreplen=maxreplen,
            exmotifs=exmotifs)

    # Update Edge-Effect Length for Background
    if background_type is not None:
        edgeeffectlength = max(edgeeffectlength or 0, background.K)

    # Parsing Sequence Constraint Feasibility
    liner.send('\n[Step 3: Parsing Primer Sequence]\n')

    # Update Step 3 Warning
    warns[3] = {
        'warn_count': 0,
        'step_name' : 'parsing-primer-sequence',
        'vars': None}

    # Parse primerseq
    (parsestatus,
    primerseq,
    homology,
    fixedbaseindex,
    exmotifindex,
    designspace,
    internalrepeats,
    pairedrepeats) = cp.get_parsed_sequence_constraint(
        primerseq=primerseq,
        primertype=primertype,
        exmotifs=exmotifs,
        pairedprimer=pairedprimer,
        warn=warns[3],
        liner=liner)

    # Remove Step 3 Warning
    if not warns[3]['warn_count']:
        warns.pop(3)

    # primerseq infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 3,
            'step_name': 'parsing-primer-sequence',
            'vars'    : {
                    'design_space': designspace,
                'internal_repeats': internalrepeats,
                  'paired_repeats': pairedrepeats},
            'warns'   : warns}
        stats['random_seed'] = random_seed
        stats = ut.stamp_stats(
            stats=stats,
            module='primer',
            input_rows=input_rows,
            output_rows=0)

        # Return results
        liner.close()
        return (outdf, stats)

    # Define pairedrepeats
    pairedrepeats = None
    if pairedprimer_map is None and \
       not pairedprimer is None:
        pairedrepeats = set(ut.stream_canon_spectrum(
            seq=pairedprimer,
            k=6))

    # Parse Melting Temperature
    liner.send('\n[Step 4: Parsing Melting Temperature]\n')

    tmelt_bounds = None
    tmelt_by_set = None

    # Paired primer is global or absent
    if pairedprimer_map is None:

        # Parse mintmelt and maxtmelt
        (parsestatus,
        estimatedminTm,
        estimatedmaxTm,
        higherminTm,
        lowermaxTm,
        mintmelt,
        maxtmelt) = cp.get_parsed_primer_tmelt_constraint(
            primerseq=primerseq,
            pairedprimer=pairedprimer,
            mintmelt=mintmelt,
            maxtmelt=maxtmelt,
            element='Primer',
            liner=liner,
            rng=rng)

        # mintmelt and maxtmelt infeasible
        if not parsestatus:

            # Prepare stats
            stats = {
                'status'  : False,
                'basis'   : 'infeasible',
                'step'    : 4,
                'step_name': 'parsing-melting-temperature',
                'vars'    : {
                    'estimated_min_Tm': estimatedminTm,
                    'estimated_max_Tm': estimatedmaxTm,
                       'higher_min_Tm': higherminTm,
                        'lower_max_Tm': lowermaxTm},
                'warns'   : warns}
            stats['random_seed'] = random_seed
            stats = ut.stamp_stats(
                stats=stats,
                module='primer',
                input_rows=input_rows,
                output_rows=0)

            # Return results
            liner.close()
            return (outdf, stats)

    # Paired primer is set-specific
    else:

        # Compute feasible Tm bounds once
        (posminTm,
        posmaxTm) = cp.get_primer_tmelt_bounds(
            primerseq=primerseq,
            liner=liner,
            rng=rng)
        tmelt_bounds = (posminTm, posmaxTm)
        tmelt_by_set = {}

        # Parse paired Tm per set
        for idx,label in enumerate(target_labels):
            if idx > 0:
                liner.send('\n')
            liner.send(
                ' Set {}: {:,} Record(s)\n'.format(
                    label,
                    target_sizes[label]))

            (parsestatus,
            estimatedminTm,
            estimatedmaxTm,
            higherminTm,
            lowermaxTm,
            set_mintmelt,
            set_maxtmelt) = cp.get_parsed_primer_tmelt_constraint(
                primerseq=primerseq,
                pairedprimer=pairedprimer_map[label],
                mintmelt=mintmelt,
                maxtmelt=maxtmelt,
                element='Primer',
                liner=liner,
                rng=rng,
                tmelt_bounds=tmelt_bounds)

            # mintmelt and maxtmelt infeasible for this set
            if not parsestatus:

                # Prepare stats
                stats = {
                    'status'  : False,
                    'basis'   : 'infeasible',
                    'step'    : 4,
                    'step_name': 'parsing-melting-temperature',
                    'vars'    : {
                        'estimated_min_Tm': estimatedminTm,
                        'estimated_max_Tm': estimatedmaxTm,
                           'higher_min_Tm': higherminTm,
                            'lower_max_Tm': lowermaxTm,
                             'failed_set': label},
                    'warns'   : warns}
                stats['random_seed'] = random_seed
                stats = ut.stamp_stats(
                    stats=stats,
                    module='primer',
                    input_rows=input_rows,
                    output_rows=0)

                # Return results
                liner.close()
                return (outdf, stats)

            tmelt_by_set[label] = (set_mintmelt, set_maxtmelt)

    # Parse Edge Effects
    leftedges = None
    rightedges = None
    prefixdicts = None
    suffixdicts = None
    prefixdict = None
    suffixdict = None

    if not exmotifs is None and \
       ((not leftcontext_full  is None) or \
        (not rightcontext_full is None)):

        # Show update
        liner.send('\n[Step 5: Extracting Context Sequences]\n')

        # Extract Primer Contexts
        if oligosets is None:
            leftcontext_target = leftcontext_full
            rightcontext_target = rightcontext_full
            if design_mask is not None:
                if not leftcontext_full is None:
                    leftcontext_target = leftcontext_full[design_mask]
                if not rightcontext_full is None:
                    rightcontext_target = rightcontext_full[design_mask]

            (leftcontext,
            rightcontext) = ut.get_extracted_context(
                leftcontext=leftcontext_target,
                rightcontext=rightcontext_target,
                edgeeffectlength=edgeeffectlength,
                reduce=True,
                liner=liner)
        else:
            leftcontext_target = leftcontext_full
            rightcontext_target = rightcontext_full
            if design_mask is not None:
                if not leftcontext_full is None:
                    leftcontext_target = leftcontext_full[design_mask]
                if not rightcontext_full is None:
                    rightcontext_target = rightcontext_full[design_mask]

            (leftedges,
            rightedges) = ut.get_extracted_context(
                leftcontext=leftcontext_target,
                rightcontext=rightcontext_target,
                edgeeffectlength=edgeeffectlength,
                reduce=False,
                liner=liner)

        # Show update
        if oligosets is None:
            liner.send('\n[Step 6: Parsing Edge Effects]\n')

            # Update Step 6 Warning
            warns[6] = {
                'warn_count': 0,
                'step_name' : 'parsing-edge-effects',
                'vars': None}

            # Compute Forbidden Prefixes and Suffixes
            (prefixdict,
            suffixdict) = cp.get_parsed_edgeeffects(
                primerseq=primerseq,
                leftcontext=leftcontext,
                rightcontext=rightcontext,
                leftpartition=leftpartition,
                rightpartition=rightpartition,
                exmotifs=exmotifs,
                element='Primer',
                warn=warns[6],
                liner=liner)

            # Remove Step 6 Warning
            if not warns[6]['warn_count']:
                warns.pop(6)

        else:
            liner.send('\n[Step 6: Parsing Edge Effects (Per-Set)]\n')

            # Update Step 6 Warning
            warns[6] = {
                'warn_count': 0,
                'step_name' : 'parsing-edge-effects',
                'vars': None}

            # Compute Forbidden Prefixes and Suffixes per set
            prefixdicts = {}
            suffixdicts = {}
            for label in target_labels:

                # Build per-set contexts
                set_left = None
                set_right = None
                if not leftedges is None:
                    set_left = list(set(
                        leftedges[idx] for idx in target_groups[label]))
                if not rightedges is None:
                    set_right = list(set(
                        rightedges[idx] for idx in target_groups[label]))

                liner.send(
                    ' Set {}: {:,} Record(s)\n'.format(
                        label,
                        target_sizes[label]))

                setwarn = {
                    'warn_count': 0,
                    'vars': None}

                (prefixdict,
                suffixdict) = cp.get_parsed_edgeeffects(
                    primerseq=primerseq,
                    leftcontext=set_left,
                    rightcontext=set_right,
                    leftpartition=leftpartition,
                    rightpartition=rightpartition,
                    exmotifs=exmotifs,
                    element='Primer',
                    warn=setwarn,
                    liner=liner)

                prefixdicts[label] = prefixdict
                suffixdicts[label] = suffixdict
                warns[6]['warn_count'] += setwarn['warn_count']

            # Remove Step 6 Warning
            if not warns[6]['warn_count']:
                warns.pop(6)

    # Parse Oligo Pool Repeats
    liner.send('\n[Step 7: Parsing Oligo Pool Repeats]\n')

    # Parse Repeats from indf
    (parsestatus,
    sourcecontext,
    kmerspace,
    fillcount,
    freecount,
    oligorepeats) = ut.get_parsed_oligopool_repeats(
        df=indf,
        maxreplen=maxreplen,
        element='Primer',
        merge=True,
        liner=liner)

    # Repeat Length infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 7,
            'step_name': 'parsing-oligopool-repeats',
            'vars'    : {
                'source_context': sourcecontext,
                'kmer_space'    : kmerspace,
                'fill_count'    : fillcount,
                'free_count'    : freecount},
            'warns'   : warns}
        stats['random_seed'] = random_seed
        stats = ut.stamp_stats(
            stats=stats,
            module='primer',
            input_rows=input_rows,
            output_rows=0)

        # Return results
        liner.close()
        return (outdf, stats)

    # Launching Primer Design
    liner.send('\n[Step 8: Computing Primer]\n')

    # Define Primer Design Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 8,
        'step_name': 'computing-primer',
        'vars'    : {
                   'primer_Tm': None,          # Primer Melting Temperature
                   'primer_GC': None,          # Primer GC Content
                 'hairpin_MFE': None,          # Primer Hairpin Free Energy
               'homodimer_MFE': None,          # Homodimer Free Energy
             'heterodimer_MFE': None,          # Heterodimer Free Energy
                     'Tm_fail': 0,             # Melting Temperature Fail Count
                 'repeat_fail': 0,             # Repeat Fail Count
             'background_fail': 0,             # Background Fail Count
              'homodimer_fail': 0,             # Homodimer Fail Count
            'heterodimer_fail': 0,             # Heterodimer Fail Count
             'crossdimer_fail': 0,             # Cross-primer Dimer Fail Count
                'exmotif_fail': 0,             # Exmotif Elimination Fail Count
                   'edge_fail': 0,             # Edge Effect Fail Count
             'exmotif_counter': cx.Counter()}, # Exmotif Encounter Counter
        'warns'   : warns}
    stats['random_seed'] = random_seed
    if not oligosets is None:
        stats['vars']['oligo_set_count'] = len(set_labels)
        stats['vars']['oligo_set_sizes'] = [
            {'oligo_set': label, 'count': set_sizes[label]}
            for label in set_labels]
        stats['vars']['primer_sets'] = []

    # Schedule outfile deletion
    ofdeletion = ae.register(
        ut.remove_file,
        outfile)

    # Design Primer(s)
    primer = None
    primer_by_set = None

    if oligosets is None:
        (primer,
        stats) = cp.primer_engine(
            primerseq=primerseq,
            primerspan=None,
            homology=homology,
            primertype=primertype,
            fixedbaseindex=fixedbaseindex,
            mintmelt=mintmelt,
            maxtmelt=maxtmelt,
            maxreplen=maxreplen,
            oligorepeats=oligorepeats,
            pairedprimer=pairedprimer,
            pairedspan=None,
            pairedrepeats=pairedrepeats,
            crossprimers=None,
            exmotifs=exmotifs,
            exmotifindex=exmotifindex,
            edgeeffectlength=edgeeffectlength,
            prefixdict=prefixdict,
            suffixdict=suffixdict,
            background=background,
            stats=stats,
            liner=liner,
            rng=rng)

    else:
        primer_by_set = {}
        crossprimers = []
        if patch_mode_on and existing_primers:
            crossprimers.extend(
                [primer for primer in existing_primers.values()
                 if not primer is None])

        liner.send('\n')

        for idx,label in enumerate(target_labels):
            if idx > 0:
                liner.send('\n')
            liner.send(
                ' Set {}: {:,} Record(s)\n'.format(
                    label,
                    target_sizes[label]))

            # Per-set paired primer and Tm constraints
            set_pairedprimer = pairedprimer
            if not pairedprimer_map is None:
                set_pairedprimer = pairedprimer_map[label]

            set_pairedrepeats = None
            if not set_pairedprimer is None:
                set_pairedrepeats = set(ut.stream_canon_spectrum(
                    seq=set_pairedprimer,
                    k=6))

            set_mintmelt = mintmelt
            set_maxtmelt = maxtmelt
            if not tmelt_by_set is None:
                set_mintmelt, set_maxtmelt = tmelt_by_set[label]

            # Per-set edge effects
            set_prefixdict = prefixdict
            set_suffixdict = suffixdict
            if not prefixdicts is None:
                set_prefixdict = prefixdicts[label]
                set_suffixdict = suffixdicts[label]

            # Per-set stats
            setstats = {
                'status'  : False,
                'basis'   : 'unsolved',
                'step'    : 8,
                'step_name': 'computing-primer',
                'vars'    : {
                       'primer_Tm': None,          # Primer Melting Temperature
                       'primer_GC': None,          # Primer GC Content
                     'hairpin_MFE': None,          # Primer Hairpin Free Energy
                   'homodimer_MFE': None,          # Homodimer Free Energy
                 'heterodimer_MFE': None,          # Heterodimer Free Energy
                         'Tm_fail': 0,             # Melting Temperature Fail Count
                     'repeat_fail': 0,             # Repeat Fail Count
                 'background_fail': 0,             # Background Fail Count
                  'homodimer_fail': 0,             # Homodimer Fail Count
                'heterodimer_fail': 0,             # Heterodimer Fail Count
                 'crossdimer_fail': 0,             # Cross-primer Dimer Fail Count
                    'exmotif_fail': 0,             # Exmotif Elimination Fail Count
                       'edge_fail': 0,             # Edge Effect Fail Count
                 'exmotif_counter': cx.Counter()}, # Exmotif Encounter Counter
                'warns'   : warns}
            setstats['random_seed'] = random_seed

            # Design primer for the set
            (setprimer,
            setstats) = cp.primer_engine(
                primerseq=primerseq,
                primerspan=None,
                homology=homology,
                primertype=primertype,
                fixedbaseindex=fixedbaseindex,
                mintmelt=set_mintmelt,
                maxtmelt=set_maxtmelt,
                maxreplen=maxreplen,
                oligorepeats=oligorepeats,
                pairedprimer=set_pairedprimer,
                pairedspan=None,
                pairedrepeats=set_pairedrepeats,
                crossprimers=crossprimers,
                exmotifs=exmotifs,
                exmotifindex=exmotifindex,
                edgeeffectlength=edgeeffectlength,
                prefixdict=set_prefixdict,
                suffixdict=set_suffixdict,
                background=background,
                stats=setstats,
                liner=liner,
                rng=rng)

            # Aggregate stats
            stats['vars']['Tm_fail'] += setstats['vars']['Tm_fail']
            stats['vars']['repeat_fail'] += setstats['vars']['repeat_fail']
            stats['vars']['background_fail'] += setstats['vars']['background_fail']
            stats['vars']['homodimer_fail'] += setstats['vars']['homodimer_fail']
            stats['vars']['heterodimer_fail'] += setstats['vars']['heterodimer_fail']
            stats['vars']['crossdimer_fail'] += setstats['vars']['crossdimer_fail']
            stats['vars']['exmotif_fail'] += setstats['vars']['exmotif_fail']
            stats['vars']['edge_fail'] += setstats['vars']['edge_fail']
            stats['vars']['exmotif_counter'] += setstats['vars']['exmotif_counter']

            # Design status
            if setstats['status']:
                primer_by_set[label] = setprimer
                crossprimers.append(setprimer)

                # Extend oligopool repeats with designed primer
                if not oligorepeats is None:
                    oligorepeats.update(ut.stream_canon_spectrum(
                        seq=setprimer,
                        k=maxreplen+1))

                # Record per-set metrics
                stats['vars']['primer_sets'].append({
                    'oligo_set'      : label,
                    'primer'         : setprimer,
                    'primer_Tm'      : setstats['vars']['primer_Tm'],
                    'primer_GC'      : setstats['vars']['primer_GC'],
                    'hairpin_MFE'    : setstats['vars']['hairpin_MFE'],
                    'homodimer_MFE'  : setstats['vars']['homodimer_MFE'],
                    'heterodimer_MFE': setstats['vars']['heterodimer_MFE'],
                })

            else:
                stats['vars']['failed_set'] = label
                stats['vars']['failed_set_size'] = set_sizes[label]
                stats['step'] = setstats['step']
                stats['step_name'] = setstats['step_name']
                stats['status'] = False
                stats['basis'] = 'unsolved'
                break

        if len(primer_by_set) == len(target_labels):
            stats['status'] = True
            stats['basis'] = 'solved'

    # Primer Status
    if stats['status']:
        primerstatus = 'Successful'
    else:
        primerstatus = 'Failed'

    # Resolve per-row primer output
    primer_out = primer
    primer_out_mask = None
    if stats['status'] and not oligosets is None:
        if patch_mode_on and design_mask is not None:
            primer_out_mask = design_mask
            primer_out = [primer_by_set[label] for label in target_oligosets]
        else:
            primer_out = [primer_by_set[label] for label in oligosets]

    # Insert primer into indf
    if stats['status']:

        # Update indf
        if patch_mode_on and primercol_exists:
            if primer_out_mask is None:
                primer_out_mask = missing_mask
            indf.loc[primer_out_mask, primercol] = primer_out
        else:
            ut.update_df(
                indf=indf,
                lcname=leftcontextname,
                rcname=rightcontextname,
                out=primer_out,
                outcol=primercol)

        # Prepare outdf
        outdf = indf

        # Write outdf to file
        if not outfile is None:
            ut.write_df_csv(
                df=outdf,
                outfile=outfile,
                sep=',')

    # Primer Design Statistics
    liner.send('\n[Primer Design Statistics]\n')

    liner.send(
        '      Design Status     : {}\n'.format(
            primerstatus))

    # Success Relevant Stats
    if stats['status']:
        if oligosets is None:
            liner.send(
                '     Melting Temperature: {:6.2f} °C\n'.format(
                    stats['vars']['primer_Tm']))
            liner.send(
                '          GC Content    : {:6.2f} %\n'.format(
                    stats['vars']['primer_GC']))
            liner.send(
                '     Hairpin MFE        : {:6.2f} kcal/mol\n'.format(
                    stats['vars']['hairpin_MFE']))
            liner.send(
                '   Homodimer MFE        : {:6.2f} kcal/mol\n'.format(
                    stats['vars']['homodimer_MFE']))

            if not pairedprimer is None:
                liner.send(
                    ' Heterodimer MFE        : {:6.2f} kcal/mol\n'.format(
                        stats['vars']['heterodimer_MFE']))
        else:
            liner.send(
                '      Primer Sets       : {:,} Set(s)\n'.format(
                    len(stats['vars']['primer_sets'])))
            liner.send(
                '   Primer Set Summary\n')

            for entry in stats['vars']['primer_sets']:
                line = (
                    '     - Set {}: Tm {:6.2f} °C, GC {:6.2f} %, '
                    'Hairpin {:6.2f} kcal/mol, Homodimer {:6.2f} kcal/mol'
                ).format(
                    entry['oligo_set'],
                    entry['primer_Tm'],
                    entry['primer_GC'],
                    entry['hairpin_MFE'],
                    entry['homodimer_MFE'])
                if not pairedprimer_map is None:
                    line += ', Heterodimer {:6.2f} kcal/mol'.format(
                        entry['heterodimer_MFE'])
                liner.send('{}\n'.format(line))

    # Failure Relavant Stats
    else:
        maxval = max(stats['vars'][field] for field in (
            'Tm_fail',
            'repeat_fail',
            'background_fail',
            'homodimer_fail',
            'heterodimer_fail',
            'crossdimer_fail',
            'exmotif_fail',
            'edge_fail'))

        sntn, plen = ut.get_notelen(
            printlen=ut.get_printlen(
                value=maxval))

        total_conflicts = stats['vars']['Tm_fail']          + \
                          stats['vars']['repeat_fail']      + \
                          stats['vars']['background_fail']  + \
                          stats['vars']['homodimer_fail']   + \
                          stats['vars']['heterodimer_fail'] + \
                          stats['vars']['crossdimer_fail']  + \
                          stats['vars']['exmotif_fail']     + \
                          stats['vars']['edge_fail']

        liner.send(
            ' Melt. Temp. Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['Tm_fail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['Tm_fail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '      Repeat Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['repeat_fail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['repeat_fail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '  Background Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['background_fail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['background_fail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '   Homodimer Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['homodimer_fail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['homodimer_fail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            ' Heterodimer Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['heterodimer_fail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['heterodimer_fail'] * 100.,
                    B=total_conflicts)))
        if not oligosets is None:
            liner.send(
                ' Crossdimer Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                    stats['vars']['crossdimer_fail'],
                    plen,
                    sntn,
                    ut.safediv(
                        A=stats['vars']['crossdimer_fail'] * 100.,
                        B=total_conflicts)))
        liner.send(
            '     Exmotif Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['exmotif_fail'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['exmotif_fail'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '        Edge Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
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
                    '     - Exmotif {:>{}} Triggered {:{},{}} Event(s)\n'.format(
                        exmotif, qlen, count, vlen, sntn))

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Unschedule outfile deletion
    if primerstatus == 'Successful':
        ae.unregister(ofdeletion)

    # Close Background
    if background_type == 'path':
        background.close()

    # Close Liner
    liner.close()

    # Return Solution and Statistics
    stats = ut.stamp_stats(
        stats=stats,
        module='primer',
        input_rows=input_rows,
        output_rows=len(outdf.index) if outdf is not None else 0)
    outdf_return = outdf
    if (outdf is not None) and (not id_from_index):
        outdf_return = ut.get_df_with_id_column(outdf)
    return (outdf_return, stats)
