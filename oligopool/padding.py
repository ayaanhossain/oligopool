import time  as tt

import collections as cx

import nrpcalc    as nr

import background as bk
import coreprimer as cp
import utils      as ut


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

def get_padding_lengths(
    typeIIS,
    fragments,
    finallength,
    liner):

    # Time-keeping
    t0 = tt.time()

    # Figure out Extreme Lengths
    cutlen = len(typeIIS)
    minlen = float('+inf')
    maxlen = float('-inf')

    for idx,frag in enumerate(fragments):
        liner.send(' Analyzing Fragment: {}'.format(idx))
        if len(frag) < minlen:
            minlen = len(frag)
        if len(frag) > maxlen:
            maxlen = len(frag)
    
    liner.send(
        '    Shortest  Pre-Padding Length: {} bp\n'.format(
            minlen))
    liner.send(
        '     Longest  Pre-Padding Length: {} bp\n'.format(
            maxlen))

    # Compute Spans
    fminspan  = finallength - maxlen
    fmaxspan  = finallength - minlen
    
    lconsspan = fmaxspan // 2
    rconsspan = fmaxspan - lconsspan

    levalspan = fminspan // 2
    revalspan = fminspan - levalspan

    shpadspan = lconsspan + minlen + rconsspan
    lgpadspan = levalspan + maxlen + revalspan

    # Show Final Updates
    liner.send(
        ' Forward Pad Construction Length: {} bp\n'.format(
            lconsspan))
    liner.send(
        ' Reverse Pad Construction Length: {} bp\n'.format(
            rconsspan))
    
    liner.send(
        ' Forward Pad   Evaluation Length: {} bp\n'.format(
            levalspan))
    liner.send(
        ' Reverse Pad   Evaluation Length: {} bp\n'.format(
            revalspan))
    
    liner.send(
        '    Shortest Post-Padding Length: {} bp\n'.format(
            shpadspan))
    liner.send(
        '     Longest Post-Padding Length: {} bp\n'.format(
            lgpadspan))
    
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))
    
    # Adjust Construction Based on TypeIIS Cutsite
    lconsspan -= cutlen
    rconsspan -= cutlen

    # Return Results
    return (lconsspan,
        rconsspan,
        levalspan,
        revalspan,)

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

import numpy as np

def get_context():
    with open('promoters.txt') as infile:
        cntx = [x.strip() for x in infile.readlines()]
    cntx = cntx[:-1]  
    rng = np.random.default_rng(7)  
    fcntx = []
    for c in cntx:
        tl = rng.integers(0, 6)
        tr = rng.integers(0, 6)
        c = c[tl:len(c)-tr]
        fcntx.append(c)
    return fcntx

def main():

    liner = ut.liner_engine()

    t0 = tt.time()

    context = get_context()

    fragments = context

    finallength = 120

    (status,
     fwdpad,
     revpad) = padding_engine(
        typeIIS='bsai',
        fragments=fragments,
        finallength=finallength,
        mintmelt=52,
        maxtmelt=54,
        deltatmelt=0.2,
        Lmax=9,
        liner=liner)
    
    if status:
        print()
        print(fwdpad)
        print(revpad)

    print()
    print('Time Elapsed = {:.2f} sec'.format(tt.time()-t0))

if __name__ == '__main__':
    main()