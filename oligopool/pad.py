import time  as tt

import collections as cx
import atexit      as ae

import nrpcalc     as nr

import utils       as ut
import valparse    as vp
import coreprimer  as cp


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
    oligolength,
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

    # Full splitcol Validation
    splitcol_valid = vp.get_parsed_column_info(
        col=splitcol,
        df=indf,
        col_field='   Split Column     ',
        col_desc='Input from Column',
        col_type=0,
        adjcol=None,
        adjval=None,
        liner=liner)

    # Full typeIIS Validation
    (typeIIS,
    typeIISname,
    typeIIS_valid) = vp.get_parsed_typeIIS_info(
        typeIIS=typeIIS,
        typeIIS_field=' TypeIIS System     ',
        liner=liner)

    # Full maxreplen Validation
    oligolength_valid = vp.get_numeric_validity(
        numeric=oligolength,
        numeric_field='   Oligo Length     ',
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
        oligolength_valid,
        tmelt_valid,
        maxreplen_valid,
        outfile_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Paramters
    oligolength = round(oligolength)
    maxreplen   = round(maxreplen)

    # Primer Design Book-keeping
    outdf = None
    stats = None

    # Parse Split Column
    liner.send('\n[Parsing Split Column]\n')

    # Parse splitcol
    (parsestatus,
    minfragmentlen,
    maxfragmentlen,
    maxallowedlen) = cp.get_parsed_splitcol(
        indf=indf,
        splitcol=splitcol,
        oligolength=oligolength,
        liner=liner)

    # splitcol infeasible
    if not parsestatus:

        # Prepare stats
        stats = {
            'status': False,
            'basis' : 'infeasible',
            'step'  : 1,
            'vars'  : {
                'maxfragmentlen': maxfragmentlen,
                 'maxallowedlen': maxallowedlen}}

        # Return results
        liner.close()
        return (outdf, stats)

    # Parse TypeIIS Constraint
    liner.send('\n[Parsing TypeIIS System]\n')

    # Parse typeIIS
    # (parsestatus,
    # forwardpadlen,
    # reversepadlen,
    # forwardfillerlen,
    # reversefillerlen) =

    cp.get_parsed_typeIIS_constraint(
        typeIIS=typeIIS,
        typeIISname=typeIISname,
        minfragmentlen=minfragmentlen,
        maxfragmentlen=maxfragmentlen,
        oligolength=oligolength,
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