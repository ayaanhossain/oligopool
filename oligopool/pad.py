import time as tt

import collections as cx
import atexit as ae

import numpy as np
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
    random_seed:int|None=None,
    verbose:bool=True) -> Tuple[pd.DataFrame, dict]:
    '''
    Pad split fragments with paired primers (embedding a chosen 3' Type IIS site) plus optional
    spacers to meet `oligo_length_limit`. Produces an assembly-ready pad layout suitable for `final`.

    Required Parameters:
        - `input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with annotated oligo pool variants.
        - `oligo_length_limit` (`int`): Maximum allowed padded oligo length (>= 60).
        - `split_column` (`str`): Column name containing split fragments.
        - `typeIIS_system` (`str`): Type IIS restriction enzyme to be used for pad excision. See Notes.
        - `minimum_melting_temperature` (`float`): Minimum padding primer Tm (>= 25 °C).
        - `maximum_melting_temperature` (`float`): Maximum padding primer Tm (<= 95 °C).
        - `maximum_repeat_length` (`int`): Max shared repeat length b/w padding primers & oligos (between 6 and 20).

    Optional Parameters:
        - `output_file` (`str` / `None`): Filename for output DataFrame; required in CLI usage,
            optional in library usage (default: `None`).
        - `random_seed` (`int` / `None`): Seed for local RNG (default: `None`).
        - `verbose` (`bool`): If `True`, logs progress to stdout (default: `True`).

    Returns:
        - A pandas DataFrame with padded oligos; saves to `output_file` if specified.
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - `input_data` must contain a unique 'ID' column; all other columns must be non-empty DNA strings.
        - `pad` expects `split_column` to contain DNA fragments (typically output from `split`).
            Other columns in `input_data` are not preserved in the output.
        - Run `pad` separately for each split fragment column (e.g., `Split1`, `Split2`, ...).
        - Output columns are: `5primeSpacer`, `ForwardPrimer`, `<split_column>`, `ReversePrimer`, `3primeSpacer`.
        - Oligo rows already summing to or exceeding `oligo_length_limit` have a `'-'` (dash) as spacer.
        - Supported Type IIS systems (34): AcuI, AlwI, BbsI, BccI, BceAI, BciVI, BcoDI, BmrI, BpuEI,
            BsaI, BseRI, BsmAI, BsmBI, BsmFI, BsmI, BspCNI, BspQI, BsrDI, BsrI, BtgZI, BtsCI, BtsI,
            BtsIMutI, EarI, EciI, Esp3I, FauI, HgaI, HphI, HpyAV, MlyI, MnlI, SapI, SfaNI.
        - `pad` supports Type IIS systems modeled as motif + 3' cut offset into adjacent `N` bases; enzymes
            that cut upstream (or require more complex cut models) are not included in the built-in list.
        - `pad` checks that the chosen Type IIS recognition site is absent from `split_column` (both
            forward and reverse-complement) and treats it as an excluded motif to avoid reintroducing
            the site at pad junctions.
        - Type IIS is for pad removal; the overlaps from `split` drive downstream assembly.
    '''

    # Preserve return style when the caller intentionally used ID as index.
    id_from_index = ut.get_id_index_intent(input_data)

    # Argument Aliasing
    indata      = input_data
    oligolimit  = oligo_length_limit
    splitcol    = split_column
    typeIIS     = typeIIS_system
    mintmelt    = minimum_melting_temperature
    maxtmelt    = maximum_melting_temperature
    maxreplen   = maximum_repeat_length
    outfile     = output_file
    random_seed = random_seed
    verbose     = verbose

    # Start Liner
    liner = ut.liner_engine(verbose)

    # Padding Verbiage Print
    liner.send('\n[Oligopool Calculator: Assembly Mode - Pad]\n')

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
    input_rows = len(indf.index) if isinstance(indf, pd.DataFrame) else 0

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
        numeric_post_desc=' Base Pair(s) Oligo Pool Repeats',
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

    # Validate random_seed (do not auto-generate)
    (random_seed,
    seed_valid) = vp.get_parsed_random_seed_info(
        random_seed=random_seed,
        random_seed_field='  Random Seed      ',
        liner=liner)

    # First Pass Validation
    if not all([
        indata_valid,
        splitcol_valid,
        typeIIS_valid,
        oligolimit_valid,
        tmelt_valid,
        maxreplen_valid,
        outfile_valid,
        seed_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Local RNG
    rng = np.random.default_rng(random_seed)

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
             'max_fragment_length': maxfragmentlen,
              'max_allowed_length': maxallowedlen},
            'warns'   : warns}
        stats['random_seed'] = random_seed
        stats = ut.stamp_stats(
            stats=stats,
            module='pad',
            input_rows=input_rows,
            output_rows=0)

        # Return results
        liner.close()
        return (outdf, stats)

    # Check TypeIIS Compatibility
    liner.send('\n[Step 2: Checking TypeIIS Compatibility]\n')

    # Verify Type IIS recognition site is absent from split fragments
    (parsestatus,
    conflict_count,
    typeIIS_rc) = cp.get_typeIIS_compatibility(
        splitcol=splitcol,
        typeIISmotif=typeIISmotif,
        indf=indf,
        liner=liner)

    # typeIIS compatibility infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 2,
            'step_name': 'checking-typeIIS-compatibility',
            'vars'    : {
                'type_iis_motif': typeIISmotif,
                'type_iis_reverse_complement_motif': typeIIS_rc,
                'internal_type_iis_sites': conflict_count},
            'warns'   : warns}
        stats['random_seed'] = random_seed
        stats = ut.stamp_stats(
            stats=stats,
            module='pad',
            input_rows=input_rows,
            output_rows=0)

        # Return results
        liner.close()
        return (outdf, stats)

    # Parse TypeIIS Constraint
    liner.send('\n[Step 3: Parsing TypeIIS Constraint]\n')

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
            'step'    : 3,
            'step_name': 'parsing-typeIIS-constraint',
            'vars'    : {
              'min_pad_length': minpadlen,
              'max_pad_length': maxpadlen,
                'type_iis_free': typeIISfree},
            'warns'   : warns}
        stats['random_seed'] = random_seed
        stats = ut.stamp_stats(
            stats=stats,
            module='pad',
            input_rows=input_rows,
            output_rows=0)

        # Return results
        liner.close()
        return (outdf, stats)

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
        primerseq=revseq[:revcore],
        pairedprimer=None,
        mintmelt=mintmelt,
        maxtmelt=maxtmelt,
        element='Pad',
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
                'estimated_minimum_melting_temperature': estimatedminTm,
                'estimated_maximum_melting_temperature': estimatedmaxTm,
                'higher_minimum_melting_temperature'   : higherminTm,
                'lower_maximum_melting_temperature'    : lowermaxTm},
            'warns'   : warns}
        stats['random_seed'] = random_seed
        stats = ut.stamp_stats(
            stats=stats,
            module='pad',
            input_rows=input_rows,
            output_rows=0)

        # Return results
        liner.close()
        return (outdf, stats)

    # Parse Excluded Motifs
    liner.send('\n[Step 5: Parsing Excluded Motifs]\n')

    # Update Step 5 Warning
    warns[5] = {
        'warn_count': 0,
        'step_name' : 'parsing-excluded-motifs',
        'vars': {}}

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
        warn=warns[5],
        liner=liner)

    # Remove Step 5 Warning
    if not warns[5]['warn_count']:
        warns.pop(5)

    # exmotifs infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 5,
            'step_name': 'parsing-excluded-motifs',
            'vars'    : {
                 'problem_lengths': problens,
                'problematic_sequence_counts': tuple(list(
                    4**pl for pl in problens))},
            'warns'   : warns}
        stats['random_seed'] = random_seed
        stats = ut.stamp_stats(
            stats=stats,
            module='pad',
            input_rows=input_rows,
            output_rows=0)

        # Return results
        liner.close()
        return (outdf, stats)

    # Update Edge-Effect Length
    edgeeffectlength = ut.get_edgeeffectlength(
        maxreplen=maxreplen,
        exmotifs=exmotifs)

    # Show update
    liner.send('\n[Step 6: Extracting Context Sequences]\n')

    # Extract Pad Contexts
    (leftcontext,
    rightcontext) = ut.get_extracted_context(
        leftcontext=splitcol,
        rightcontext=splitcol,
        edgeeffectlength=edgeeffectlength,
        reduce=True,
        liner=liner)

    # Show update
    liner.send('\n[Step 7: Parsing Forward Pad Edge Effects]\n')

    # Update Step 7 Warning
    warns[7] = {
        'warn_count': 0,
        'step_name' : 'parsing-forward-pad-edge-effects',
        'vars': {}}

    # Compute Forbidden Prefixes and Suffixes
    (_,
    suffixdict) = cp.get_parsed_edgeeffects(
        primerseq=fwdseq[-fwdcore:],
        leftcontext=None,
        rightcontext=rightcontext,
        leftpartition=None,
        rightpartition=rightpartition,
        exmotifs=exmotifs,
        element='Forwad Pad',
        warn=warns[7],
        liner=liner)

    # Remove Step 7 Warning
    if not warns[7]['warn_count']:
        warns.pop(7)

    # Show update
    liner.send('\n[Step 8: Parsing Reverse Pad Edge Effects]\n')

    # Update Step 8 Warning
    warns[8] = {
        'warn_count': 0,
        'step_name' : 'parsing-reverse-pad-edge-effects',
        'vars': {}}

    # Compute Forbidden Prefixes and Suffixes
    (prefixdict,
    _) = cp.get_parsed_edgeeffects(
        primerseq=revseq[:revcore],
        leftcontext=leftcontext,
        rightcontext=None,
        leftpartition=leftpartition,
        rightpartition=None,
        exmotifs=exmotifs,
        element='Reverse Pad',
        warn=warns[8],
        liner=liner)

    # Remove Step 8 Warning
    if not warns[8]['warn_count']:
        warns.pop(8)

    # Parse Oligo Pool Repeats
    liner.send('\n[Step 9: Parsing Oligo Pool Repeats]\n')

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
            'step'    : 9,
            'step_name': 'parsing-oligopool-repeats',
            'vars'    : {
                'source_context': sourcecontext,
                'kmer_space'    : kmerspace,
                'filled_kmer_count'    : fillcount,
                'remaining_kmer_count' : freecount},
            'warns'   : warns}
        stats['random_seed'] = random_seed
        stats = ut.stamp_stats(
            stats=stats,
            module='pad',
            input_rows=input_rows,
            output_rows=0)

        # Return results
        liner.close()
        return (outdf, stats)

    # Define Pad Design Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 10,
        'step_name': 'computing-pad',
        'vars'    : {
            'forward_pad_primer_melting_temperature'     : None, # Forward Pad Melting Temperature
            'reverse_pad_primer_melting_temperature'     : None, # Reverse Pad Melting Temperature
            'forward_pad_primer_guanine_cytosine_content': None, # Forward Pad GC Content
            'reverse_pad_primer_guanine_cytosine_content': None, # Reverse Pad GC Content
            'forward_pad_hairpin_minimum_free_energy'    : None, # Forward Pad Hairpin Free Energy
            'reverse_pad_hairpin_minimum_free_energy'    : None, # Reverse Pad Hairpin Free Energy
            'forward_pad_homodimer_minimum_free_energy'  : None, # Forward Pad Homodimer Free Energy
            'reverse_pad_homodimer_minimum_free_energy'  : None, # Reverse Pad Homodimer Free Energy
            'heterodimer_minimum_free_energy'            : None, # Heterodimer Free Energy
            'melting_temperature_failure_count'          : 0,    # Melting Temperature Fail Count
            'repeat_failure_count'                       : 0,    # Repeat Fail Count
            'homodimer_failure_count'                    : 0,    # Homodimer Fail Count
            'heterodimer_failure_count'                  : 0,    # Heterodimer Fail Count
            'excluded_motif_failure_count'               : 0,    # Exmotif Elimination Fail Count
            'edge_effect_failure_count'                  : 0,    # Edge Effect Fail Count
            'excluded_motif_encounter_counter': cx.Counter()},   # Exmotif Encounter Counter
        'warns'   : warns}
    stats['random_seed'] = random_seed

    # Schedule outfile deletion
    ofdeletion = ae.register(
        ut.remove_file,
        outfile)

    # Define Forward Primer-Pad Design Stats
    fwdstats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 10,
        'step_name': 'computing-forward-pad',
        'vars'    : {
            'primer_melting_temperature'       : None,         # Primer Melting Temperature
            'primer_guanine_cytosine_content'  : None,         # Primer GC Content
            'hairpin_minimum_free_energy'      : None,         # Primer Hairpin Free Energy
            'homodimer_minimum_free_energy'    : None,         # Homodimer Free Energy
            'heterodimer_minimum_free_energy'  : None,         # Heterodimer Free Energy
            'melting_temperature_failure_count': 0,            # Melting Temperature Fail Count
            'repeat_failure_count'             : 0,            # Repeat Fail Count
            'homodimer_failure_count'          : 0,            # Homodimer Fail Count
            'heterodimer_failure_count'        : 0,            # Heterodimer Fail Count
            'cross_dimer_failure_count'        : 0,            # Cross-primer Dimer Fail Count
            'excluded_motif_failure_count'     : 0,            # Exmotif Elimination Fail Count
            'edge_effect_failure_count'        : 0,            # Edge Effect Fail Count
            'excluded_motif_encounter_counter': cx.Counter()}, # Exmotif Encounter Counter
        'warns'   : warns}
    fwdstats['random_seed'] = random_seed

    # Define Forward Primer-Pad Attributes
    pairedrepeats = set()
    exmotifindex  = set([fwdseq.index(
        typeIISmotif) + len(typeIISmotif)])

    # Launching Forward Primer-Pad Design
    liner.send('\n[Step 10: Computing Forward Pad]\n')

    # Design Forward Primer-Pad
    (fwdpad,
    fwdstats) = cp.primer_engine(
        primerseq=fwdseq,
        primerspan=fwdcore,
        homology=homology,
        primertype=0,
        fixedbaseindex=ut.get_fixed_base_index(seqconstr=fwdseq),
        mintmelt=mintmelt,
        maxtmelt=maxtmelt,
        maxreplen=maxreplen,
        oligorepeats=oligorepeats,
        pairedprimer=None,
        pairedspan=None,
        pairedrepeats=pairedrepeats,
        crossprimers=None,
        exmotifs=exmotifs,
        exmotifindex=exmotifindex,
        edgeeffectlength=edgeeffectlength,
        prefixdict=None,
        suffixdict=suffixdict,
        background=background,
        stats=fwdstats,
        liner=liner,
        rng=rng)

    # Define Reverse Primer-Pad Design Stats
    revstats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 11,
        'step_name': 'computing-reverse-pad',
        'vars'    : {
            'primer_melting_temperature'       : None,          # Primer Melting Temperature
            'primer_guanine_cytosine_content'  : None,          # Primer GC Content
            'hairpin_minimum_free_energy'      : None,          # Primer Hairpin Free Energy
            'homodimer_minimum_free_energy'    : None,          # Homodimer Free Energy
            'heterodimer_minimum_free_energy'  : None,          # Heterodimer Free Energy
            'melting_temperature_failure_count': 0,             # Melting Temperature Fail Count
            'repeat_failure_count'             : 0,             # Repeat Fail Count
            'homodimer_failure_count'          : 0,             # Homodimer Fail Count
            'heterodimer_failure_count'        : 0,             # Heterodimer Fail Count
            'cross_dimer_failure_count'        : 0,             # Cross-primer Dimer Fail Count
            'excluded_motif_failure_count'     : 0,             # Exmotif Elimination Fail Count
            'edge_effect_failure_count'        : 0,             # Edge Effect Fail Count
            'excluded_motif_encounter_counter': cx.Counter()},  # Exmotif Encounter Counter
        'warns'   : warns}
    revstats['random_seed'] = random_seed

    # Do we Continue?
    if fwdstats['status']:

        # Define Reverse Primer-Pad Attributes
        pairedrepeats = set(ut.stream_canon_spectrum(
            seq=fwdpad[-fwdcore:],
            k=len(typeIIS)))
        exmotifindex  = set([revseq.index(
            ut.get_revcomp(typeIISmotif)) + len(typeIISmotif)])

        # Launching Reverse Primer-Pad Design
        liner.send('\n[Step 11: Computing Reverse Pad]\n')

        # Design Reverse Primer-Pad
        (revpad,
        revstats) = cp.primer_engine(
            primerseq=revseq,
            primerspan=revcore,
            homology=homology,
            primertype=1,
            fixedbaseindex=ut.get_fixed_base_index(seqconstr=revseq),
            mintmelt=fwdstats['vars']['primer_melting_temperature']-1,
            maxtmelt=fwdstats['vars']['primer_melting_temperature']+1,
            maxreplen=maxreplen,
            oligorepeats=oligorepeats,
            pairedprimer=fwdpad,
            pairedspan=fwdcore,
            pairedrepeats=pairedrepeats,
            crossprimers=None,
            exmotifs=exmotifs,
            exmotifindex=exmotifindex,
            edgeeffectlength=edgeeffectlength,
            prefixdict=prefixdict,
            suffixdict=None,
            background=background,
            stats=revstats,
            liner=liner,
            rng=rng)

    # Meta Merge
    stats['status']    = fwdstats['status'] and \
                         revstats['status']
    stats['basis']     = 'solved' if stats['status'] else 'unsolved'
    stats['step']      = fwdstats['step'] if not revstats['status'] \
                                          else   revstats['step']
    stats['step_name'] = fwdstats['step_name'] if not revstats['status'] \
                                               else   revstats['step_name']

    # Forward Stats Merge
    stats['vars']['forward_pad_primer_melting_temperature'] = fwdstats['vars']['primer_melting_temperature']
    stats['vars']['forward_pad_primer_guanine_cytosine_content'] = fwdstats['vars']['primer_guanine_cytosine_content']
    stats['vars']['forward_pad_hairpin_minimum_free_energy'] = fwdstats['vars']['hairpin_minimum_free_energy']
    stats['vars']['forward_pad_homodimer_minimum_free_energy'] = fwdstats['vars']['homodimer_minimum_free_energy']
    stats['vars']['melting_temperature_failure_count'] = fwdstats['vars']['melting_temperature_failure_count']
    stats['vars']['repeat_failure_count'] = fwdstats['vars']['repeat_failure_count']
    stats['vars']['homodimer_failure_count'] = fwdstats['vars']['homodimer_failure_count']
    stats['vars']['heterodimer_failure_count'] = fwdstats['vars']['heterodimer_failure_count']
    stats['vars']['excluded_motif_failure_count'] = fwdstats['vars']['excluded_motif_failure_count']
    stats['vars']['edge_effect_failure_count'] = fwdstats['vars']['edge_effect_failure_count']
    stats['vars']['excluded_motif_encounter_counter'] = fwdstats['vars']['excluded_motif_encounter_counter']

    # Reverse Stats Merge
    stats['vars']['reverse_pad_primer_melting_temperature'] = revstats['vars']['primer_melting_temperature']
    stats['vars']['reverse_pad_primer_guanine_cytosine_content'] = revstats['vars']['primer_guanine_cytosine_content']
    stats['vars']['reverse_pad_hairpin_minimum_free_energy'] = revstats['vars']['hairpin_minimum_free_energy']
    stats['vars']['reverse_pad_homodimer_minimum_free_energy'] = revstats['vars']['homodimer_minimum_free_energy']
    stats['vars']['heterodimer_minimum_free_energy'] = revstats['vars']['heterodimer_minimum_free_energy']
    stats['vars']['melting_temperature_failure_count']              += revstats['vars']['melting_temperature_failure_count']
    stats['vars']['repeat_failure_count'] += revstats['vars']['repeat_failure_count']
    stats['vars']['homodimer_failure_count'] += revstats['vars']['homodimer_failure_count']
    stats['vars']['heterodimer_failure_count'] += revstats['vars']['heterodimer_failure_count']
    stats['vars']['excluded_motif_failure_count'] += revstats['vars']['excluded_motif_failure_count']
    stats['vars']['edge_effect_failure_count'] += revstats['vars']['edge_effect_failure_count']
    stats['vars']['excluded_motif_encounter_counter'] += revstats['vars']['excluded_motif_encounter_counter']

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
        indf['5primeSpacer']  = LeftSpacer
        indf['3primeSpacer']  = RightSpacer
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
            ut.write_df_csv(
                df=outdf,
                outfile=outfile,
                sep=',')

    # Primer Design Statistics
    liner.send('\n[Pad Design Statistics]\n')

    liner.send(
        '      Design Status     : {}\n'.format(
            padstatus))

    # Success Relevant Stats
    if stats['status']:

        liner.send(
            '     Melting Temperature: {:6.2f} °C and {:6.2f} °C\n'.format(
                stats['vars']['forward_pad_primer_melting_temperature'],
                stats['vars']['reverse_pad_primer_melting_temperature']))
        liner.send(
            '          GC Content    : {:6.2f} %  and {:6.2f} %\n'.format(
                stats['vars']['forward_pad_primer_guanine_cytosine_content'],
                stats['vars']['reverse_pad_primer_guanine_cytosine_content']))
        liner.send(
            '     Hairpin MFE        : {:6.2f} kcal/mol and {:6.2f} kcal/mol\n'.format(
                stats['vars']['forward_pad_hairpin_minimum_free_energy'],
                stats['vars']['reverse_pad_hairpin_minimum_free_energy']))
        liner.send(
            '   Homodimer MFE        : {:6.2f} kcal/mol and {:6.2f} kcal/mol\n'.format(
                stats['vars']['forward_pad_homodimer_minimum_free_energy'],
                stats['vars']['reverse_pad_homodimer_minimum_free_energy']))
        liner.send(
            ' Heterodimer MFE        : {:6.2f} kcal/mol\n'.format(
                stats['vars']['heterodimer_minimum_free_energy']))

    # Failure Relavant Stats
    else:
        maxval = max(stats['vars'][field] for field in (
            'melting_temperature_failure_count',
            'repeat_failure_count',
            'homodimer_failure_count',
            'heterodimer_failure_count',
            'excluded_motif_failure_count',
            'edge_effect_failure_count'))

        sntn, plen = ut.get_notelen(
            printlen=ut.get_printlen(
                value=maxval))

        total_conflicts = stats['vars']['melting_temperature_failure_count']        + \
                          stats['vars']['repeat_failure_count']      + \
                          stats['vars']['homodimer_failure_count']   + \
                          stats['vars']['heterodimer_failure_count'] + \
                          stats['vars']['excluded_motif_failure_count']     + \
                          stats['vars']['edge_effect_failure_count']

        liner.send(
            ' Melt. Temp. Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['melting_temperature_failure_count'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['melting_temperature_failure_count'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '      Repeat Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['repeat_failure_count'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['repeat_failure_count'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '   Homodimer Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['homodimer_failure_count'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['homodimer_failure_count'] * 100.,
                    B=total_conflicts)))
        liner.send(
            ' Heterodimer Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['heterodimer_failure_count'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['heterodimer_failure_count'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '     Exmotif Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['excluded_motif_failure_count'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['excluded_motif_failure_count'] * 100.,
                    B=total_conflicts)))
        liner.send(
            '        Edge Conflicts  : {:{},{}} Event(s) ({:6.2f} %)\n'.format(
                stats['vars']['edge_effect_failure_count'],
                plen,
                sntn,
                ut.safediv(
                    A=stats['vars']['edge_effect_failure_count'] * 100.,
                    B=total_conflicts)))

        # Enumerate Motif-wise Fail Counts
        if stats['vars']['excluded_motif_encounter_counter']:

            qlen = max(len(motif) \
                for motif in stats['vars']['excluded_motif_encounter_counter'].keys()) + 2

            sntn, vlen = ut.get_notelen(
                printlen=ut.get_printlen(
                    value=max(
                        stats['vars']['excluded_motif_encounter_counter'].values())))

            liner.send('   Exmotif-wise Conflict Distribution\n')

            for exmotif,count in stats['vars']['excluded_motif_encounter_counter'].most_common():
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
    stats = ut.stamp_stats(
        stats=stats,
        module='pad',
        input_rows=input_rows,
        output_rows=len(outdf.index) if outdf is not None else 0)
    outdf_return = outdf
    if (outdf is not None) and (not id_from_index):
        outdf_return = ut.get_df_with_id_column(outdf)
    return (outdf_return, stats)
