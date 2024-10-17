import time  as tt

import collections as cx
import atexit as ae

import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp
from .base import core_primer as cp

from typing import Tuple


def pad(
    input_data:str|pd.DataFrame,
    oligo_length_limit:int,
    split_column:str,
    typeIIS_system:str,
    minimum_melting_temperature:float,
    maximum_melting_temperature:float,
    maximum_repeat_length:int,
    output_file:str|None=None,
    verbose:bool=True) -> Tuple[pd.DataFrame, dict]:
    '''
    Pads split oligos with optimized paired primers with a 3prime TypeIIS restriction site of choice
    and adds optional flanking spacers to reach oligo length limit. Returned DataFrame adds new columns
    to `input_data`, is compatible with `final` module, and is optionally written to a CSV file.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated oligopool variants.
        - `oligo_length_limit` (`int`): Maximum allowed padded oligo length (≥ 60).
        - `split_column` (`str`): Column name containing split fragments.
        - `typeIIS_system` (`str`): Type IIS restriction enzyme to be used for pad excision. See notes.
        - `minimum_melting_temperature` (`float`): Minimum padding primer Tm (≥ 25°C).
        - `maximum_melting_temperature` (`float`): Maximum padding primer Tm (≤ 95°C).
        - `maximum_repeat_length` (`int`): Max shared repeat length b/w padding primers & oligos (between 6 and 20).

    Optional Parameters:
        - `output_file` (`str`): Filename for output DataFrame (default: `None`).
        - `verbose` (`bool`): If `True`, logs updates to stdout (default: `True`).

    Returns:
        - A pandas DataFrame with padded oligos; saves to `output_file` if specified.
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` must contain a unique 'ID' column, all other columns must be non-empty DNA strings.
        - Column names in `input_data` must be unique, and exclude `primer_column`.
        - Oligo rows already summing to or exceeding `oligo_length_limit` have a `'-'` (dash) as spacer.
        - The following 34 unique Type IIS systems are supported to faciliate scarless padding primer excision.
            * `AcuI`,  `AlwI`,  `BbsI`,  `BccI`,   `BceAI`,    `BciVI`,
              `BcoDI`, `BmrI`,  `BpuEI`, `BsaI`,   `BseRI`,    `BsmAI`,
              `BsmBI`, `BsmFI`, `BsmI`,  `BspCNI`, `BspQI`,    `BsrDI`,
              `BsrI`,  `BtgZI`, `BtsCI`, `BtsI`,   `BtsIMutI`, `EarI`,
              `EciI`,  `Esp3I`, `FauI`,  `HgaI`,   `HphI`,     `HpyAV`,
              `MlyI`,  `MnlI`,  `SapI`,  `SfaNI`
    '''

    # Argument Aliasing
    indata     = input_data
    oligolimit = oligo_length_limit
    splitcol   = split_column
    typeIIS    = typeIIS_system
    mintmelt   = minimum_melting_temperature
    maxtmelt   = maximum_melting_temperature
    maxreplen  = maximum_repeat_length
    outfile    = output_file
    verbose    = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Barcoding Verbage Print
    liner.send('\n[Oligopool Calculator: Design Mode - Pad]\n')

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

    # Full oligolimit Validation
    oligolimit_valid = vp.get_numeric_validity(
        numeric=oligolimit,
        numeric_field='   Oligo Limit      ',
        numeric_pre_desc=' Design ',
        numeric_post_desc=' Base Pair(s) Padded Oligos',
        minval=60,
        maxval=float('+inf'),
        precheck=False,
        liner=liner)

    # Store Split Column Name
    splitcolname = splitcol

    # Full splitcol Validation
    (splitcol,
    splitcol_valid) = vp.get_parsed_column_info(
        col=splitcol,
        df=indf,
        col_field='   Split Column     ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=None,
        adjval=None,
        iscontext=False,
        typecontext=None,
        liner=liner)

    # Full typeIIS Validation
    (typeIIS,
    typeIISname,
    typeIIS_valid) = vp.get_parsed_typeIIS_info(
        typeIIS=typeIIS,
        typeIIS_field=' TypeIIS System     ',
        liner=liner)

    # Full mintmelt and maxtmelt Validation
    (mintmelt,
    maxtmelt,
    tmelt_valid) = vp.get_parsed_range_info(
        minval=mintmelt,
        maxval=maxtmelt,
        range_field=' Melting Temperature',
        range_unit='°C',
        range_min=25,
        range_max=95,
        liner=liner)

    # Full maxreplen Validation
    maxreplen_valid = vp.get_numeric_validity(
        numeric=maxreplen,
        numeric_field='  Repeat Length     ',
        numeric_pre_desc=' Up to ',
        numeric_post_desc=' Base Pair(s) Oligopool Repeats',
        minval=6,
        maxval=20,
        precheck=False,
        liner=liner)

    # Full outfile Validation
    outfile_valid = vp.get_outdf_validity(
        outdf=outfile,
        outdf_suffix='.oligopool.pad.csv',
        outdf_field='  Output File       ',
        liner=liner)

    # Adjust outfile Suffix
    if not outfile is None:
        outfile = ut.get_adjusted_path(
            path=outfile,
            suffix='.oligopool.pad.csv')

    # First Pass Validation
    if not all([
        indata_valid,
        splitcol_valid,
        typeIIS_valid,
        oligolimit_valid,
        tmelt_valid,
        maxreplen_valid,
        outfile_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Paramters
    oligolimit = round(oligolimit)
    maxreplen  = round(maxreplen)

    # Define Additional Variables
    typeIISmotif = typeIIS.replace('N', '')
    homology     = len(typeIISmotif) + 1
    background   = None

    # Primer Design Book-keeping
    outdf = None
    stats = None
    warns = {}

    # Parse Split Column
    liner.send('\n[Step 1: Parsing Split Column]\n')

    # Parse splitcol
    (parsestatus,
    minfragmentlen,
    maxfragmentlen,
    maxallowedlen,
    paddingbalance) = cp.get_parsed_splitcol(
        splitcol=splitcol,
        oligolimit=oligolimit,
        liner=liner)

    # splitcol infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 1,
            'step_name': 'parsing-split-column',
            'vars'    : {
                'maxfragmentlen': maxfragmentlen,
                 'maxallowedlen': maxallowedlen},
            'warns'   : warns}

        # Return results
        liner.close()
        return (outdf, stats)

    # Parse TypeIIS Constraint
    liner.send('\n[Step 2: Parsing TypeIIS System]\n')

    # Parse typeIIS
    (parsestatus,
    fwdcore,
    revcore,
    fwdseq,
    revseq,
    minpadlen,
    maxpadlen,
    typeIISfree) = cp.get_parsed_typeIIS_constraint(
        typeIIS=typeIIS,
        typeIISname=typeIISname,
        minfragmentlen=minfragmentlen,
        maxfragmentlen=maxfragmentlen,
        oligolimit=oligolimit,
        liner=liner)

    # typeIIS infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 2,
            'step_name': 'parsing-typeIIS-system',
            'vars'    : {
                  'minpadlen': minpadlen,
                  'maxpadlen': maxpadlen,
                'typeIISfree': typeIISfree},
            'warns'   : warns}

        # Return results
        liner.close()
        return (outdf, stats)

    # Parse Melting Temperature
    liner.send('\n[Step 3: Parsing Melting Temperature]\n')

    # Parse mintmelt and maxtmelt
    (parsestatus,
    estimatedminTm,
    estimatedmaxTm,
    higherminTm,
    lowermaxTm,
    mintmelt,
    maxtmelt) = cp.get_parsed_primer_tmelt_constraint(
        primerseq=revseq[:revcore],
        pairedprimer=None,
        mintmelt=mintmelt,
        maxtmelt=maxtmelt,
        element='Pad',
        liner=liner)

    # mintmelt and maxtmelt infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 3,
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

    # Parse Excluded Motifs
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
        exmotifs=(typeIISmotif,ut.get_revcomp(typeIISmotif)),
        typer=tuple,
        element='Pad',
        leftcontext=splitcol,
        rightcontext=splitcol,
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

    # Show update
    liner.send('\n[Step 5: Extracting Context Sequences]\n')

    # Extract Pad Contexts
    (leftcontext,
    rightcontext) = ut.get_extracted_context(
        leftcontext=splitcol,
        rightcontext=splitcol,
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
        primerseq=fwdseq[-fwdcore:],
        leftcontext=leftcontext,
        rightcontext=rightcontext,
        leftpartition=leftpartition,
        rightpartition=rightpartition,
        exmotifs=exmotifs,
        element='Pad',
        warn=warns[6],
        liner=liner)

    # Remove Step 6 Warning
    if not warns[6]['warn_count']:
        warns.pop(6)

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
        element='Pad',
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

    # Define Pad Design Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 8,
        'step_name': 'computing-primer',
        'vars'    : {
                'fwd_pad_primer_Tm': None,          # Forward Pad Melting Temperature
                'rev_pad_primer_Tm': None,          # Reverse Pad Melting Temperature
                'fwd_pad_primer_GC': None,          # Forward Pad GC Content
                'rev_pad_primer_GC': None,          # Reverse Pad GC Content
              'fwd_pad_hairpin_MFE': None,          # Forward Pad Hairpin Free Energy
              'rev_pad_hairpin_MFE': None,          # Reverse Pad Hairpin Free Energy
            'fwd_pad_homodimer_MFE': None,          # Forward Pad Homodimer Free Energy
            'rev_pad_homodimer_MFE': None,          # Reverse Pad Homodimer Free Energy
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

    # Define Forward Primer-Pad Design Stats
    fwdstats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 8,
        'step_name': 'computing-forward-pad',
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

    # Define Forward Primer-Pad Attributes
    pairedrepeats = set()
    exmotifindex  = set([fwdseq.index(
        typeIISmotif) + len(typeIISmotif)])

    # Launching Forward Primer-Pad Design
    liner.send('\n[Step 8: Computing Forward Pad]\n')

    # Design Forward Primer-Pad
    (fwdpad,
    fwdstats) = cp.primer_engine(
        primerseq=fwdseq,
        primerspan=fwdcore,
        homology=homology,
        primertype=0,
        fixedbaseindex=cp.get_fixedbaseindex(seq=fwdseq),
        mintmelt=mintmelt,
        maxtmelt=maxtmelt,
        maxreplen=maxreplen,
        oligorepeats=oligorepeats,
        pairedprimer=None,
        pairedspan=None,
        pairedrepeats=pairedrepeats,
        exmotifs=exmotifs,
        exmotifindex=exmotifindex,
        edgeeffectlength=edgeeffectlength,
        prefixdict=None,
        suffixdict=suffixdict,
        background=background,
        stats=fwdstats,
        liner=liner)

    # Define Reverse Primer-Pad Design Stats
    revstats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 9,
        'step_name': 'computing-reverse-pad',
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

    # Do we Continue?
    if fwdstats['status']:

        # Define Reverse Primer-Pad Attributes
        pairedrepeats = set(ut.stream_canon_spectrum(
            seq=fwdpad[-fwdcore:],
            k=len(typeIIS)))
        exmotifindex  = set([revseq.index(
            ut.get_revcomp(typeIISmotif)) + len(typeIISmotif)])

        # Launching Reverse Primer-Pad Design
        liner.send('\n[Step 9: Computing Reverse Pad]\n')

        # Design Reverse Primer-Pad
        (revpad,
        revstats) = cp.primer_engine(
            primerseq=revseq,
            primerspan=revcore,
            homology=homology,
            primertype=1,
            fixedbaseindex=cp.get_fixedbaseindex(seq=revseq),
            mintmelt=fwdstats['vars']['primer_Tm']-1,
            maxtmelt=fwdstats['vars']['primer_Tm']+1,
            maxreplen=maxreplen,
            oligorepeats=oligorepeats,
            pairedprimer=fwdpad,
            pairedspan=fwdcore,
            pairedrepeats=pairedrepeats,
            exmotifs=exmotifs,
            exmotifindex=exmotifindex,
            edgeeffectlength=edgeeffectlength,
            prefixdict=prefixdict,
            suffixdict=None,
            background=background,
            stats=revstats,
            liner=liner)

    # Meta Merge
    stats['status']    = fwdstats['status'] and \
                         revstats['status']
    stats['basis']     = 'solved' if stats['status'] else 'unsolved'
    stats['step']      = fwdstats['step'] if not revstats['status'] \
                                          else   revstats['step']
    stats['step_name'] = fwdstats['step_name'] if not revstats['status'] \
                                               else   revstats['step_name']

    # Forward Stats Merge
    stats['vars']['fwd_pad_primer_Tm']     = fwdstats['vars']['primer_Tm']
    stats['vars']['fwd_pad_primer_GC']     = fwdstats['vars']['primer_GC']
    stats['vars']['fwd_pad_hairpin_MFE']   = fwdstats['vars']['hairpin_MFE']
    stats['vars']['fwd_pad_homodimer_MFE'] = fwdstats['vars']['homodimer_MFE']
    stats['vars']['Tm_fail']             = fwdstats['vars']['Tm_fail']
    stats['vars']['repeat_fail']         = fwdstats['vars']['repeat_fail']
    stats['vars']['homodimer_fail']      = fwdstats['vars']['homodimer_fail']
    stats['vars']['heterodimer_fail']      = fwdstats['vars']['heterodimer_fail']
    stats['vars']['exmotif_fail']          = fwdstats['vars']['exmotif_fail']
    stats['vars']['edge_fail']             = fwdstats['vars']['edge_fail']
    stats['vars']['exmotif_counter']       = fwdstats['vars']['exmotif_counter']

    # Reverse Stats Merge
    stats['vars']['rev_pad_primer_Tm']     = revstats['vars']['primer_Tm']
    stats['vars']['rev_pad_primer_GC']     = revstats['vars']['primer_GC']
    stats['vars']['rev_pad_hairpin_MFE']   = revstats['vars']['hairpin_MFE']
    stats['vars']['rev_pad_homodimer_MFE'] = revstats['vars']['homodimer_MFE']
    stats['vars']['heterodimer_MFE']       = revstats['vars']['heterodimer_MFE']
    stats['vars']['Tm_fail']              += revstats['vars']['Tm_fail']
    stats['vars']['repeat_fail']          += revstats['vars']['repeat_fail']
    stats['vars']['homodimer_fail']       += revstats['vars']['homodimer_fail']
    stats['vars']['heterodimer_fail']     += revstats['vars']['heterodimer_fail']
    stats['vars']['exmotif_fail']         += revstats['vars']['exmotif_fail']
    stats['vars']['edge_fail']            += revstats['vars']['edge_fail']
    stats['vars']['exmotif_counter']      += revstats['vars']['exmotif_counter']

    # Primer Status
    if stats['status']:
        padstatus = 'Successful'
    else:
        padstatus = 'Failed'

    # Insert primer into indf
    if stats['status']:

        # Extract indf
        indf = indf[[splitcolname]]

        # Compute columns
        LeftSpacer    = []
        ForwardPrimer = []
        RightSpacer   = []
        ReversePrimer = []

        # Decompose Balance
        for balance in paddingbalance:

            # Get the Current Balance
            p,q = balance

            # Left Pad Extration from Balance
            xfwdpad = fwdpad[-p:]
            s = p - fwdcore

            # Do we have Left Spacer?
            if s > 0:
                leftspacer = xfwdpad[:s]
            else:
                leftspacer = '-'
            fwdprimer = xfwdpad[-fwdcore:]

            # Right Pad Extration from Balance
            xrevpad = revpad[:q]
            s = q - revcore

            # Do we have Right Spacer?
            if s > 0:
                rightspacer = xrevpad[-s:]
            else:
                rightspacer = '-'
            revprimer = xrevpad[:+revcore]

            # Add Elements to Columns
            LeftSpacer.append(leftspacer)
            RightSpacer.append(rightspacer)
            ForwardPrimer.append(fwdprimer)
            ReversePrimer.append(revprimer)

        # Add columns
        indf['5primeSpacer']     = LeftSpacer
        indf['3primeSpacer']     = RightSpacer
        indf['ForwardPrimer'] = ForwardPrimer
        indf['ReversePrimer'] = ReversePrimer

        # Prepare outdf
        outdf = indf
        outdf = outdf[['5primeSpacer',
            'ForwardPrimer',
            splitcolname,
            'ReversePrimer',
            '3primeSpacer']]

        # Write outdf to file
        if not outfile is None:
            outdf.to_csv(
                path_or_buf=outfile,
                sep=',')

    # Primer Design Statistics
    liner.send('\n[Pad Design Statistics]\n')

    liner.send(
        '       Design Status    : {}\n'.format(
            padstatus))

    # Success Relevant Stats
    if stats['status']:

        liner.send(
            '     Melting Temperature: {:6.2f} °C and {:6.2f} °C\n'.format(
                stats['vars']['fwd_pad_primer_Tm'],
                stats['vars']['rev_pad_primer_Tm']))
        liner.send(
            '          GC Content    : {:6.2f} %  and {:6.2f} %\n'.format(
                stats['vars']['fwd_pad_primer_GC'],
                stats['vars']['rev_pad_primer_GC']))
        liner.send(
            '     Hairpin MFE        : {:6.2f} kcal/mol and {:6.2f} kcal/mol\n'.format(
                stats['vars']['fwd_pad_hairpin_MFE'],
                stats['vars']['rev_pad_hairpin_MFE']))
        liner.send(
            '   Homodimer MFE        : {:6.2f} kcal/mol and {:6.2f} kcal/mol\n'.format(
                stats['vars']['fwd_pad_homodimer_MFE'],
                stats['vars']['rev_pad_homodimer_MFE']))
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

        total_conflicts = stats['vars']['Tm_fail']        + \
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
    if padstatus == 'Successful':
        ae.unregister(ofdeletion)

    # Close Liner
    liner.close()

    # Return Solution and Statistics
    return (outdf, stats)