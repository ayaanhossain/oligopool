import time  as tt

import collections as cx
import atexit      as ae

import nrpcalc     as nr

import utils       as ut
import valparse    as vp
import coreprimer  as cp
import primer      as pr


# Gamplan
'''
1. Figure out Max and Min Span
2. Build out longest left and right filler constraints Ns + Motif or Motif + Ns
3. Local Model Functions:
   - Homodimer   on the longest  span should be feasible
   - Heterodimer on the shortest span should be feasible
   - Tm on the shortest span should be feasible
   - NO exmotifs
   - NO Edge Constraints
4. Build Maker Object
5. Build Background of Fragments
6. Run Maker to Build Forward Pad
7. Run Maker to Build Reverse Pad
'''

def pad_engine(
    fo,
    homology,
    primertype,
    fixedbaseindex,
    mintmelt,
    maxtmelt,
    maxreplen,
    oligorepeats,
    pairedprimer,
    pairedrepeats,
    exmotifs,
    exmotifindex,
    edgeeffectlength,
    prefixdict,
    suffixdict,
    background,
    stats,
    liner):
    '''
    TBW
    '''

    # Book-keeping
    t0 = tt.time() # Start Timer
    # Flexible Structure Constraint
    primerstruct = ''.join('x.'[idx in fixedbaseindex] \
        for idx in fixedbaseindex)

    # Update Paired Primer Orientation
    # Since Forward Primer Design,
    # interpret Paired Primer as
    # Reverse Primer specified in
    # terms of Forward Strand
    if primertype == 0:
        if not pairedprimer is None:
            pairedprimer = ut.get_revcomp(
                seq=pairedprimer)
    # Correct Free Base Index
    else:
        fixedbaseindex = set(len(primerseq)-1-idx \
            for idx in fixedbaseindex)

    # Optimize exmotifs
    if not exmotifs is None:
        exmotifs = ut.get_grouped_sequences(
            sequences=exmotifs)

    # Define Objective Function
    objectivefunction = lambda primer: cp.primer_objectives(
        primer=primer,
        primerlen=len(primerseq),
        primertype=primertype,
        fixedbaseindex=fixedbaseindex,
        mintmelt=mintmelt,
        maxtmelt=maxtmelt,
        maxreplen=maxreplen,
        oligorepeats=oligorepeats,
        pairedprimer=pairedprimer,
        pairedrepeats=pairedrepeats,
        exmotifs=exmotifs,
        exmotifindex=exmotifindex,
        edgeeffectlength=edgeeffectlength,
        prefixforbidden=prefixdict,
        suffixforbidden=suffixdict,
        background=background,
        inittime=t0,
        stats=stats,
        liner=liner)

    # Define Maker Instance
    maker = nr.base.maker.NRPMaker(
        part_type='DNA',
        seed=None)

    # Design Primer via Maker
    primer = maker.nrp_maker(
        homology=min(len(primerseq), homology),
        seq_constr=primerseq,
        struct_constr=primerstruct,
        target_size=1,
        background=None,
        struct_type='both',
        synth_opt=False,
        local_model_fn=objectivefunction,
        jump_count=1000,
        fail_count=100000,
        output_file=None,
        verbose=False,
        abortion=True,
        allow_internal_repeat=False,
        check_constraints=False)

    # Check Status and Return Solution
    if len(primer) > 0:

        # Final Update
        cp.show_update(
            primer=primer[0],
            element='Primer',
            optstatus=2,
            optstate=0,
            inittime=t0,
            terminal=True,
            liner=liner)

        # Extract Primer
        primer = primer[0]

        # We solved this!
        stats['status'] = True
        stats['basis']  = 'solved'

        # Update Tm
        stats['vars']['primerTm'] = ut.get_tmelt(
            seq=primer)

        # Update GC Percentage
        stats['vars']['primerGC'] = (primer.count('G') + \
                                     primer.count('C')) \
                                        / (len(primer) * 0.01)

        # Update Hairpin Free Energy
        stats['vars']['hairpinMFE'] = cp.folder.evaluate_mfe(
            seq=primer,
            dg=True)[-1]

        # Correct Primer Orientation
        cprimer = ut.get_revcomp(
            seq=primer) if primertype == 1 else primer

        # Update Heterodimer Free Energy
        stats['vars']['homodimerMFE'] = cp.folder.evaluate_mfe_dimer(
            seq1=cprimer,
            seq2=cprimer)[-1]

        # Update Homodimer Free Energy
        if not pairedprimer is None:
            stats['vars']['heterodimerMFE'] = cp.folder.evaluate_mfe_dimer(
                seq1=cprimer,
                seq2=pairedprimer)[-1]

        # Return Results
        return (primer, stats)

    # Design Unsuccessful
    else:

        # Final Update
        liner.send('\* Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - t0))

        # Return Results
        return (None, stats)

def get_background(Lmax):

    # Build the Background Instance
    _, bkgpath = ut.setup_workspace(
            outfile='',
            outfile_suffix=None)
    bkg = bk.background(
        path='./{}'.format(bkgpath),
        Lmax=Lmax,
        verbose=True)
    bkg.clear()

    # Return Background
    return bkg

def get_frag_background(
    fragments,
    Lmax,
    liner):

    # Get Fragment Background
    fragbkg = get_background(Lmax=Lmax)

    # Populate Background Instance
    fragbkg.multiadd(fragments)

    # Return Background
    return fragbkg

def is_pad_local_feasible(
    pad,
    padtype,
    padlen,
    evalspan,
    mintmelt,
    maxtmelt,
    pairedpad,
    pairedevalspan,
    background,
    liner):

    # Book-keeping
    allstate = False
    traceloc = len(pad)

    # Background Non-Repetitiveness
    bcond, bfail = cp.is_background_feasible(
        primer=pad,
        background=background)
    if not bcond:
        liner.send(
            ' Candidate: {} Rejected due to Background Infeasibility'.format(
                pad))
        return False, bfail

    # Melting Temperature
    tcond, tfail = cp.is_tmelt_feasible(
        primer=pad,
        primerlen=padlen,
        mintmelt=mintmelt,
        maxtmelt=maxtmelt)
    if not tcond:
        liner.send(
            ' Candidate: {} Rejected due to Tm Infeasibility'.format(
                pad))
        return False, cp.get_tmelt_traceback(
            primer=pad,
            failtype=tfail)


    # Homodimer Feasibility
    hcond, htloc = cp.is_dimer_feasible(
        primer=pad,
        primertype=padtype,
        primerlen=padlen,
        primerspan=evalspan,
        pairedprimer=None,
        pairedspan=evalspan,
        dimertype=0)
    if not hcond:
        liner.send(
            ' Candidate: {} Rejected due to Homodimer Infeasibility'.format(
                pad))
        return False, htloc

    # Heterodimer Feasibility
    qcond, qtloc = cp.is_dimer_feasible(
        primer=pad,
        primertype=padtype,
        primerlen=padlen,
        primerspan=evalspan,
        pairedprimer=pairedpad,
        pairedspan=pairedevalspan,
        dimertype=1)
    if not qcond:
        liner.send(
            ' Candidate: {} Rejected due to Heterodimer Infeasibility'.format(
                pad))
        return False, qtloc

    # Everything OK!
    if len(pad) < padlen:
        liner.send(
            ' Candidate: {} Partially Accepted'.format(
                pad))
    else:
        liner.send(
            ' Candidate: {} Completely Accepted'.format(
                pad))

    return True

def padding_engine(
    typeIIS,
    fragments,
    finallength,
    mintmelt,
    maxtmelt,
    deltatmelt,
    Lmax,
    liner):

    # Book-keeping
    fwdpads = None
    revpads = None
    typeIIS = ut.typeIIS_dict[typeIIS.lower()]

    # Estimate the Padding Span Length
    liner.send('\n[Computing Padding Lengths]\n')
    (lconsspan,
     rconsspan,
     levalspan,
     revalspan) = get_padding_lengths(
        typeIIS=typeIIS,
        fragments=fragments,
        finallength=finallength,
        liner=liner)

    # Build Forward and Reverse Padding Constraint
    fwdpadseq = ('N' * lconsspan) + typeIIS
    fwdstruct = 'x' * len(fwdpadseq)
    revpadseq = ut.get_revcomp(typeIIS) + ('N' * rconsspan)
    revstruct = 'x' * len(revpadseq)

    # Determine TypeIIS Feasibility
    liner.send('\n[Checking Type IIS Motif Feasibility]\n')
    (typeIIS,
     fwdavail,
     revavail) = cp.evaluate_typeIIS_constraint(
         typeIIS=typeIIS,
         levalspan=levalspan,
         revalspan=revalspan,
         liner=liner)

    if not typeIIS:
        liner.send(' Verdict: Pad Design Infeasible due to Type IIS Motif Constraint\n')
        return (False, None, None)
    else:
        liner.send(' Verdict: Pad Design Possibly Feasible\n')

    # Determine Pad Tm Feasibility
    for locseq,loc in ((fwdpadseq, 'Forward'), (revpadseq, 'Reverse')):

        # Determine Forward / Reverse Pad Tm Feasibility
        liner.send(
            '\n[Checking {} Pad Tm Feasibility]\n'.format(
                loc))
        (tmeltstatus,
         emintmelt,
         emaxtmelt) = cp.evaluate_tmelt_constraint(
            primerseq=locseq,
            mintmelt=mintmelt,
            maxtmelt=maxtmelt,
            liner=liner)

        if not tmeltstatus:
            liner.send(' Verdict: {} Pad Infeasible due to Tm Constraint\n'.format(loc))
            return (False, None, None)
        else:
            liner.send(' Verdict: {} Pad Design Possibly Feasible\n'.format(loc))

    # Build Fragment Background
    fragbkg = get_frag_background(
        fragments=fragments,
        Lmax=Lmax,
        liner=liner)

    # Define Maker Instance
    maker = nr.base.maker.NRPMaker(
        part_type='DNA',
        seed=None)

    # Define Padding Background
    padbkg = get_background(Lmax=7)

    # Show Update
    liner.send('\n[Computing Forward Pad]\n')

    # Forward Pad Local Objective
    fwd_local_model_fn = lambda pad: is_pad_local_feasible(
        pad=pad,
        padtype=0,
        padlen=len(fwdpadseq),
        evalspan=levalspan,
        mintmelt=mintmelt,
        maxtmelt=maxtmelt,
        pairedpad=None,
        pairedevalspan=None,
        background=fragbkg,
        liner=liner)

    # Construct Forward Pad
    fwdpad = maker.nrp_maker(
        homology=8,
        seq_constr=fwdpadseq,
        struct_constr=fwdstruct,
        target_size=1,
        background=padbkg,
        struct_type='both',
        synth_opt=True,
        local_model_fn=fwd_local_model_fn,
        jump_count=1000,
        fail_count=100000,
        output_file=None,
        verbose=False,
        abortion=True,
        allow_internal_repeat=False,
        check_constraints=False)

    # Check Status and Return Solution
    status = not (not fwdpad)
    if not status:
        liner.send('\* Forward Pad Constructed? No\n')
        return (status, None, None)
    else:
        fwdpad = fwdpad[0]
        liner.send('\* Forward Pad Constructed? Yes\n')

    # fwdpad = 'TATGTAGTTATTCATGGTCTCTTCCC'

    # Compute Requirements for Reverse Primer
    fwdtmelt = ut.get_tmelt(seq=fwdpad)
    padbkg.add(fwdpad)

    # Show Update
    liner.send('\n[Computing Reverse Pad]\n')

    # Forward Pad Local Objective
    rev_local_model_fn = lambda pad: is_pad_local_feasible(
        pad=pad,
        padtype=1,
        padlen=len(revpadseq),
        evalspan=revalspan,
        mintmelt=fwdtmelt-deltatmelt,
        maxtmelt=fwdtmelt+deltatmelt,
        pairedpad=fwdpad,
        pairedevalspan=levalspan,
        background=fragbkg,
        liner=liner)

    # Construct Forward Pad
    revpad = maker.nrp_maker(
        homology=6,
        seq_constr=revpadseq,
        struct_constr=revstruct,
        target_size=1,
        background=padbkg,
        struct_type='both',
        synth_opt=True,
        local_model_fn=rev_local_model_fn,
        jump_count=1000,
        fail_count=100000,
        output_file=None,
        verbose=False,
        abortion=True,
        allow_internal_repeat=False,
        check_constraints=False)

    # Check Status and Return Solution
    status = not (not revpad)
    if not status:
        liner.send('\* Reverse Pad Constructed? No\n')
        return (status, None, None)
    else:
        revpad = revpad[0]
        liner.send('\* Reverse Pad Constructed? Yes\n')

    # Drop Backgrounds
    fragbkg.drop()
    padbkg.drop()

    return (True, fwdpad, revpad)

def pad(
    indata,
    splitcol,
    typeIIS,
    oligolimit,
    mintmelt,
    maxtmelt,
    maxreplen,
    outfile=None,
    verbose=True):
    '''

    Note 2. 34 unique TypeIIS systems available for padding
            AcuI,  AlwI,  BbsI,  BccI,   BceAI,    BciVI,
            BcoDI, BmrI,  BpuEI, BsaI,   BseRI,    BsmAI,
            BsmBI, BsmFI, BsmI,  BspCNI, BspQI,    BsrDI,
            BsrI,  BtgZI, BtsCI, BtsI,   BtsIMutI, EarI,
            EciI,  Esp3I, FauI,  HgaI,   HphI,     HpyAV,
            MlyI,  MnlI,  SapI,  SfaNI
    '''

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

    # Full maxreplen Validation
    oligolimit_valid = vp.get_numeric_validity(
        numeric=oligolimit,
        numeric_field='   Oligo Limit      ',
        numeric_pre_desc=' Design ',
        numeric_post_desc=' Base Pair(s) Padded Oligos',
        minval=60,
        maxval=float('+inf'),
        precheck=False,
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

    # Define homology
    homology = len(typeIIS.replace('N', ''))

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
            'stepname': 'parsing-split-column',
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
            'stepname': 'parsing-typeIIS-system',
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
            'stepname': 'parsing-melting-temperature',
            'vars'    : {
                'estimatedminTm': estimatedminTm,
                'estimatedmaxTm': estimatedmaxTm,
                   'higherminTm': higherminTm,
                    'lowermaxTm': lowermaxTm},
            'warns'   : warns}

        # Return results
        liner.close()
        return (outdf, stats)

    # Parse Excluded Motifs
    liner.send('\n[Step 4: Parsing Excluded Motifs]\n')

    # Update Step 4 Warning
    warns[4] = {
        'warncount': 0,
        'stepname' : 'parsing-excluded-motifs',
        'vars': None}

    # Parse exmotifs
    (parsestatus,
    exmotifs,
    problens,
    leftpartition,
    rightpartition) = ut.get_parsed_exmotifs(
        exmotifs=(typeIIS.replace('N', ''), ),
        typer=tuple,
        element='Pad',
        leftcontext=splitcol,
        rightcontext=splitcol,
        warn=warns[4],
        liner=liner)

    # Remove Step 4 Warning
    if not warns[4]['warncount']:
        warns.pop(4)

    # exmotifs infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status'  : False,
            'basis'   : 'infeasible',
            'step'    : 4,
            'stepname': 'parsing-excluded-motifs',
            'vars'    : {
                 'problens': problens,
                'probcount': tuple(list(
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
        'warncount': 0,
        'stepname' : 'parsing-edge-effects',
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
    if not warns[6]['warncount']:
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
            'stepname': 'parsing-oligopool-repeats',
            'vars'    : {
                'sourcecontext': sourcecontext,
                'kmerspace'    : kmerspace,
                'fillcount'    : fillcount,
                'freecount'    : freecount},
            'warns'   : warns}

        # Return results
        liner.close()
        return (outdf, stats)

    # Define Pad Design Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 8,
        'stepname': 'computing-primer',
        'vars'    : {
                'fwdpadprimerTm': None,          # Forward Pad Melting Temperature
                'revpadprimerTm': None,          # Reverse Pad Melting Temperature
                'fwdpadprimerGC': None,          # Forward Pad GC Content
                'revpadprimerGC': None,          # Reverse Pad GC Content
              'fwdpadhairpinMFE': None,          # Forward Pad Hairpin Free Energy
              'revpadhairpinMFE': None,          # Reverse Pad Hairpin Free Energy
            'fwdpadhomodimerMFE': None,          # Forward Pad Homodimer Free Energy
            'revpadhomodimerMFE': None,          # Reverse Pad Homodimer Free Energy
                'heterodimerMFE': None,          # Heterodimer Free Energy
                        'Tmfail': 0,             # Melting Temperature Fail Count
                    'repeatfail': 0,             # Repeat Fail Count
                 'homodimerfail': 0,             # Homodimer Fail Count
               'heterodimerfail': 0,             # Heterodimer Fail Count
                   'exmotiffail': 0,             # Exmotif Elimination Fail Count
                      'edgefail': 0,             # Edge Effect Fail Count
                'exmotifcounter': cx.Counter()}, # Exmotif Encounter Counter
        'warns'   : warns}

    # Schedule outfile deletion
    ofdeletion = ae.register(
        ut.remove_file,
        outfile)

    # Launching Forward Primer Design
    liner.send('\n[Step 8: Computing Forward Primer]\n')

    # Define Forward Primer Design Stats
    fwdstats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 8,
        'stepname': 'computing-primer',
        'vars'    : {
                   'primerTm': None,          # Primer Melting Temperature
                   'primerGC': None,          # Primer GC Content
                 'hairpinMFE': None,          # Primer Hairpin Free Energy
               'homodimerMFE': None,          # Homodimer Free Energy
             'heterodimerMFE': None,          # Heterodimer Free Energy
                     'Tmfail': 0,             # Melting Temperature Fail Count
                 'repeatfail': 0,             # Repeat Fail Count
              'homodimerfail': 0,             # Homodimer Fail Count
            'heterodimerfail': 0,             # Heterodimer Fail Count
                'exmotiffail': 0,             # Exmotif Elimination Fail Count
                   'edgefail': 0,             # Edge Effect Fail Count
             'exmotifcounter': cx.Counter()}, # Exmotif Encounter Counter
        'warns'   : warns}

    # Design Primer
    (fwdprimer,
    stats) = pr.primer_engine(
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
        prefixdict=prefixdict,
        suffixdict=suffixdict,
        background=background,
        stats=stats,
        liner=liner)




    '''
    --1. Parse Split Column--
    2. Parse Oligo Length
    3. Extract Sequence Constraint (both)
    3. Parse TypeIIS Feasibility
    4. Parse Melting Temperature
    5. Call Engine
    6. Parse Output
    '''