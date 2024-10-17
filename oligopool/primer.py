import time as tt

import collections as cx
import atexit as ae

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
    primer_type:int,
    minimum_melting_temperature:float,
    maximum_melting_temperature:float,
    maximum_repeat_length:int,
    primer_column:str,
    output_file:str|None=None,
    paired_primer_column:str|None=None,
    left_context_column:str|None=None,
    right_context_column:str|None=None,
    excluded_motifs:list|str|pd.DataFrame|None=None,
    background_directory:str|None=None,
    verbose:bool=True) -> Tuple[pd.DataFrame, dict]:
    '''
    Designs constrained primers with specified melting temperature and non-repetitiveness
    for all variants in the oligopool. Ensures compatibility with paired primers and minimizes
    dimer formation. Returns a DataFrame of designed primers, optionally saving to `output_file`
    in CSV format.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated oligopool variants.
        - `oligo_length_limit` (`int`): Maximum allowed oligo length (≥ 4).
        - `primer_sequence_constraint` (`int`): IUPAC degenerate sequence constraint.
        - `primer_type` (`int`): Primer type (0 for forward, 1 for reverse).
        - `minimum_melting_temperature` (`float`): Minimum primer Tm (≥ 25°C).
        - `maximum_melting_temperature` (`float`): Maximum primer Tm (≤ 95°C).
        - `maximum_repeat_length` (`int`): Max shared repeat length with oligos (≥ 6).
        - `primer_column` (`str`): Column name for the designed primer.

    Optional Parameters:
        - `output_file` (`str`): Filename for output DataFrame (default: `None`).
        - `paired_primer_column` (`str`): Column for paired primer sequence (default: `None`).
        - `left_context_column` (`str`): Column for left DNA context (default: `None`).
        - `right_context_column` (`str`): Column for right DNA context (default: `None`).
        - `excluded_motifs` (`list` / `str` / `pd.DataFrame`): Motifs to exclude;
            can be a CSV path or DataFrame (default: `None`).
        - `background_directory` (`str`): Directory for background k-mer sequences (default: `None`).
        - `verbose` (`bool`): If `True`, logs updates to stdout (default: `True`).

    Returns:
        - A pandas DataFrame of designed primers; saves to `output_file` if specified.
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` must contain a unique 'ID' column, all other columns must be non-empty DNA strings.
        - Column names in `input_data` must be unique, and exclude `primer_column`.
        - At least one of `left_context_column` or `right_context_column` must be specified.
        - The paired primer type is inferred based on the current primer type.
        - If a paired primer is specified, Tm of the designed primer is optimized within 1°C of it.
        - `maximum_repeat_length` controls non-repetitiveness against `input_data` only, not `background`.
        - If `excluded_motifs` is a CSV or DataFrame, it must have 'ID' and 'Exmotif' columns.
        - Constant motifs in sequence constraint may lead to sub-optimal primers.
    '''

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
    exmotifs     = excluded_motifs
    background   = background_directory
    verbose      = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Primer Verbage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Primer]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # First Pass indata Parsing and Validation
    (indf,
    indata_valid) = vp.get_parsed_indata_info(
        indata=indata,
        indata_field='      Input Data       ',
        required_fields=('ID',),
        precheck=False,
        liner=liner)

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
    primertype_valid = vp.get_categorical_validity(
        category=primertype,
        category_field='     Primer Type       ',
        category_pre_desc=' ',
        category_post_desc=' Primer Design',
        category_dict={
            0: 'Forward',
            1: 'Reverse'},
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
        numeric_post_desc=' Base Pair(s) Oligopool Repeats',
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
        liner=liner)

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

    # Full pairedprimer Validation
    (pairedprimer,
    pairedprimer_valid) = vp.get_constantcol_validity(
        constantcol=pairedcol,
        constantcol_field='     Paired Primer     ',
        df=indf,
        liner=liner)

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

    # Define Edge Effect Length
    edgeeffectlength = None

    # Primer Design Book-keeping
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

            # Return results
            liner.close()
            return (outdf, stats)

        # Update Edge-Effect Length
        edgeeffectlength = ut.get_edgeeffectlength(
            exmotifs=exmotifs)

    # Parsing Sequence Constraint Feasibility
    liner.send('\n[Step 3: Parsing Primer Sequence]\n')

    # Update Step 3 Warning
    warns[3] = {
        'warn_count': 0,
        'step_name' : 'parsing-primer-sequence',
        'vars': None}

    # Parse primerseq
    (parsestatus,
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

        # Return results
        liner.close()
        return (outdf, stats)

    # Define pairedrepeats
    if not pairedprimer is None:
        pairedrepeats = set(ut.stream_canon_spectrum(
            seq=pairedprimer,
            k=6))
    else:
        pairedrepeats = None

    # Parse Melting Temperature
    liner.send('\n[Step 4: Parsing Melting Temperature]\n')

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
        liner=liner)

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

        # Return results
        liner.close()
        return (outdf, stats)

    # Parse Edge Effects
    if not exmotifs is None and \
       ((not leftcontext  is None) or \
        (not rightcontext is None)):

        # Show update
        liner.send('\n[Step 5: Extracting Context Sequences]\n')

        # Extract Primer Contexts
        (leftcontext,
        rightcontext) = ut.get_extracted_context(
            leftcontext=leftcontext,
            rightcontext=rightcontext,
            edgeeffectlength=edgeeffectlength,
            reduce=True,
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
        (prefixdict,
        suffixdict) = None, None

    # Parse Oligopool Repeats
    liner.send('\n[Step 7: Parsing Oligopool Repeats]\n')

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
              'homodimer_fail': 0,             # Homodimer Fail Count
            'heterodimer_fail': 0,             # Heterodimer Fail Count
                'exmotif_fail': 0,             # Exmotif Elimination Fail Count
                   'edge_fail': 0,             # Edge Effect Fail Count
             'exmotif_counter': cx.Counter()}, # Exmotif Encounter Counter
        'warns'   : warns}

    # Schedule outfile deletion
    ofdeletion = ae.register(
        ut.remove_file,
        outfile)

    # Design Primer
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
        exmotifs=exmotifs,
        exmotifindex=exmotifindex,
        edgeeffectlength=edgeeffectlength,
        prefixdict=prefixdict,
        suffixdict=suffixdict,
        background=background,
        stats=stats,
        liner=liner)

    # Primer Status
    if stats['status']:
        primerstatus = 'Successful'
    else:
        primerstatus = 'Failed'

    # Insert primer into indf
    if stats['status']:

        # Update indf
        ut.update_df(
            indf=indf,
            lcname=leftcontextname,
            rcname=rightcontextname,
            out=primer,
            outcol=primercol)

        # Prepare outdf
        outdf = indf

        # Write outdf to file
        if not outfile is None:
            outdf.to_csv(
                path_or_buf=outfile,
                sep=',')

    # Primer Design Statistics
    liner.send('\n[Primer Design Statistics]\n')

    liner.send(
        '      Design Status     : {}\n'.format(
            primerstatus))

    # Success Relevant Stats
    if stats['status']:

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

    # Failure Relavant Stats
    else:
        maxval = max(stats['vars'][field] for field in (
            'Tm_fail',
            'repeat_fail',
            'homodimer_fail',
            'heterodimer_fail',
            'exmotif_fail',
            'edge_fail'))

        sntn, plen = ut.get_notelen(
            printlen=ut.get_printlen(
                value=maxval))

        total_conflicts = stats['vars']['Tm_fail']          + \
                          stats['vars']['repeat_fail']      + \
                          stats['vars']['homodimer_fail']   + \
                          stats['vars']['heterodimer_fail'] + \
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
    return (outdf, stats)
