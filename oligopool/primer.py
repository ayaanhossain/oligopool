import time  as tt

import collections as cx

import nrpcalc    as nr

import background as bk
import coreprimer as cp
import utils      as ut


def is_primer_local_feasible(
    primer,
    primertype,
    primerlen,
    mintmelt,
    maxtmelt,
    exmotifs,
    prefixgroup,
    suffixgroup,
    pairedprimer,
    background,
    liner):

    # Book-keeping
    allstate = False
    traceloc = len(primer)

    # Background Non-Repetitiveness
    bcond, bfail = cp.is_background_feasible(
        primer=primer,
        background=background)
    if not bcond:
        liner.send(
            ' Candidate: {} Rejected due to Background Infeasibility'.format(
                primer))
        return False, bfail

    # Melting Temperature
    tcond, tfail = cp.is_tmelt_feasible(
        primer=primer,
        primerlen=primerlen,
        mintmelt=mintmelt,
        maxtmelt=maxtmelt)
    if not tcond:
        liner.send(
            ' Candidate: {} Rejected due to Tm Infeasibility'.format(
                primer))
        return False, cp.get_tmelt_traceback(
            primer=primer,
            failtype=tfail)

    # Motif Embedding
    mcond, motif = cp.is_motif_feasible(
        primer=primer,
        exmotifs=exmotifs)
    if not mcond:
        liner.send(
            ' Candidate: {} Rejected due to Motif Infeasibility'.format(
                primer))
        return False, max(0, len(primer)-len(motif))

    # Context Feasibility
    ccond, ctloc = cp.is_context_feasible(
        primer=primer,
        primerlen=primerlen,
        prefixgroup=prefixgroup,
        suffixgroup=suffixgroup)
    if not ccond:
        liner.send(
            ' Candidate: {} Rejected due to Context Infeasibility'.format(
                primer))
        return False, ctloc

    # Homodimer Feasibility
    hcond, htloc = cp.is_dimer_feasible(
        primer=primer,
        primertype=primertype,
        primerlen=primerlen,
        primerspan=None,
        pairedprimer=None,
        pairedspan=None,
        dimertype=0)
    if not hcond:
        liner.send(
            ' Candidate: {} Rejected due to Homodimer Infeasibility'.format(
                primer))
        return False, htloc

    # Heterodimer Feasibility
    qcond, qtloc = cp.is_dimer_feasible(
        primer=primer,
        primertype=primertype,
        primerlen=primerlen,
        primerspan=None,
        pairedprimer=pairedprimer,
        pairedspan=None,
        dimertype=1)
    if not qcond:
        liner.send(
            ' Candidate: {} Rejected due to Heterodimer Infeasibility'.format(
                primer))
        return False, qtloc

    # Everything OK!
    if len(primer) < primerlen:
        liner.send(
            ' Candidate: {} Partially Accepted'.format(
                primer))
    else:
        liner.send(
            ' Candidate: {} Completely Accepted'.format(
                primer))

    return True

def primer_engine(
    primerseq,
    primertype,
    mintmelt,
    maxtmelt,
    Lmax,
    pairedprimer,
    pairedprimertype,
    deltatmelt,
    exmotifs,
    leftcontext,
    rightcontext,
    background,
    liner):
    
    # Book-keeping
    primerstruct = 'x'*len(primerseq)    # Primer Structure Constraint
    mfails       = cx.Counter()          # Motif Fail Counter
    emmotifs     = None                  # Sequence Constraint Motif Conflicts
    lcf          = True                  # Left  Context Feasibility
    rcf          = True                  # Right Context Feasibility
    sxfails      = set()                 # Problematic Suffix Fails
    pxfails      = set()                 # Problematic Predix Fails
    emintmelt    = None
    emaxtmelt    = None

    # Paired Primer Based Adjustment
    if pairedprimer:

        # Paired Primer Sequence
        if pairedprimertype == 1:
            pairedprimer = ut.get_revcomp(
                seq=pairedprimer)

        # Pairing based Tm
        pairedtmelt = ut.get_tmelt(
            seq=pairedprimer)
        mintmelt = pairedtmelt - deltatmelt
        maxtmelt = pairedtmelt + deltatmelt

    # Context Setup
    lcnum = ut.get_context_num(          # Total Number of Left Context
        context=leftcontext)
    rcnum = ut.get_context_num(          # Total Number of Right Context
        context=rightcontext)
    lcifn = ut.get_context_inference_fn( # Left  Context Selector
        context=leftcontext)
    rcifn = ut.get_context_inference_fn( # Right Context Selector
        context=rightcontext)

    # Setup Exmotifs
    cntxlen = 0
    if exmotifs:
        liner.send('\n[Preprocessing Excluded Motifs]\n')
        exmotifs = ut.prep_exmotifs(
            exmotifs=exmotifs,
            packing=tuple,
            liner=liner)
        cntxlen = ut.get_context_len(
            exmotifs=exmotifs)
    
    # Determine Sequence Constraint Feasibility
    liner.send('\n[Checking Primer Sequence Feasibility]\n')
    (seqstatus,
     emmotifs) = cp.evaluate_seq_constraint(
        primerseq=primerseq,
        exmotifs=exmotifs,
        liner=liner)
    
    if not seqstatus:
        liner.send(' Verdict: Primer Design Infeasible due to Sequence Constraint\n')
        return (False, None)
    else:
        liner.send(' Verdict: Primer Design Possibly Feasible\n')
    

    # Determine Tm Feasibility
    liner.send('\n[Checking Primer Tm Feasibility]\n')
    (tmeltstatus,
     emintmelt,
     emaxtmelt) = cp.evaluate_tmelt_constraint(
        primerseq=primerseq,
        mintmelt=mintmelt,
        maxtmelt=maxtmelt,
        liner=liner)

    if not tmeltstatus:
        liner.send(' Verdict: Primer Design Infeasible due to Tm Constraints\n')
        return (False, None)
    else:
        liner.send(' Verdict: Primer Design Possibly Feasible\n')

    # Extract Edge Effect Constraints
    liner.send('\n[Extracting Primer Edge Constraints]\n')
    (edgestatus,
     prefixgroup,
     suffixgroup,
     prefixdict,
     suffixdict) = cp.get_edge_contraints(
        lcnum=lcnum,
        rcnum=rcnum,
        lcifn=lcifn,
        rcifn=rcifn,
        exmotifs=exmotifs,
        liner=liner)

    # Evaluate Edge Effect Constraints
    if edgestatus:
        liner.send('\n[Checking Primer Edge Feasibility]\n')
        
        if lcnum:
            liner.send(' Infeasible Left Context Sequences:\n')
            lcf, sxfails = cp.evaluate_edge_constraint(
                edgedict=prefixgroup,
                culpdict=suffixdict,
                contexttype=0,
                liner=liner)

        if rcnum:
            liner.send(' Infeasible Right Context Sequences:\n')
            rcf, pxfails = cp.evaluate_edge_constraint(
                edgedict=suffixgroup,
                culpdict=prefixdict,
                contexttype=1,
                liner=liner)

    # Determine Edge Feasibility
    feasible = lcf and rcf
    if not feasible:
        liner.send(' Verdict: Primer Design Infeasible due to Context Sequences\n')
        return (False, None)
    else:
        liner.send(' Verdict: Primer Design Possibly Feasible\n')

    # Define Maker Instance
    maker = nr.base.maker.NRPMaker(
        part_type='DNA',
        seed=None)

    # Show Update
    liner.send('\n[Computing Primer]\n')

    # Local and Global Objectives
    local_model_fn = lambda primer: is_primer_local_feasible(
        primer=primer,
        primertype=primertype,
        primerlen=len(primerseq),
        mintmelt=mintmelt,
        maxtmelt=maxtmelt,
        exmotifs=exmotifs,
        prefixgroup=prefixgroup,
        suffixgroup=suffixgroup,
        pairedprimer=pairedprimer,
        background=background,
        liner=liner)

    # Setup Primer Background
    primerbkg = None
    if pairedprimer:
        _, primerbkgdir = ut.setup_workspace(
            outfile='',
            outfile_suffix=None)
        primerbkg = bk.background(
            path='./{}'.format(primerbkgdir),
            Lmax=Lmax,
            verbose=True)
        primerbkg.clear()
        primerbkg.add(pairedprimer)

    # Run Maker
    primer = maker.nrp_maker(
        homology=Lmax+1,
        seq_constr=primerseq,
        struct_constr=primerstruct,
        target_size=1,
        background=primerbkg,
        struct_type='both',
        synth_opt=True,
        local_model_fn=local_model_fn,
        jump_count=1000,
        fail_count=100000,
        output_file=None,
        verbose=False,
        abortion=True,
        allow_internal_repeat=False,
        check_constraints=False)

    # Drop Primer Background
    if primerbkg:
        primerbkg.drop()

    # Check Status and Return Solution
    status = not (not primer)
    if not status:
        liner.send('\* Primer Constructed? No\n')
        return (status, None)
    else:
        primer = primer[0]
        liner.send('\* Primer Constructed? Yes\n')
        return (status, primer)

def get_context():
    with open('promoters.txt') as infile:
        cntx = [x.strip() for x in infile.readlines()]
    return cntx[:-1]

def main():

    # print(folder.evaluate_mfe_dimer(
    #     seq1='GTGCCTCGTCCGAAACCA',
    #     seq2=ut.get_revcomp('CTGCTAACACACGCCCCA')))

    # return

    liner = ut.liner_engine()

    t0 = tt.time()

    # exmotifs = ['AAAA', 'GGGG', 'CCCC', 'TTTT', 'GGATCC', 'TCTAGA', 'GAATTC'][::-1]
    exmotifs = ['GGATCC', 'TCTAGA', 'CCCCC']

    context = get_context()

    p1 = 'GCAGCAATCGTGACAGGGATCC'
    p2 = 'TTCATACCGACGAGAGATCAGC'
    p3 = 'TCTAGACGCCGATGACACGCAA'
    fwdprimers = []
    revprimers = []

    # Populate background
    background = None
    _, backgrounddir = ut.setup_workspace(
        outfile='',
        outfile_suffix=None)
    background = bk.background(
        path='./{}'.format(backgrounddir),
        Lmax=10,
        verbose=True)
    background.clear()
    background.multiadd(context)

    # Build out P1
    primerseq = 'NSNSNSNSNSNSNSNSNSNSNS'
    primerseq = 'NNNNNNNNNNNNNNNNNNNNNN'
    primerseq = 'NNNNNNNNNNNNNNNCCM'
    pairedseq = 'CTGCTAACACACGCCCCA'
    # pairedseq = None

    status, primer = primer_engine(
        primerseq=primerseq,
        primertype=0,
        mintmelt=53,
        maxtmelt=57,
        Lmax=5,
        pairedprimer=pairedseq,
        pairedprimertype=1,
        deltatmelt=0.2,
        exmotifs=exmotifs,
        leftcontext=None,
        rightcontext=context,
        background=background,
        liner=liner)

    if status:
        print()
        print(primer)
        print('Primer Tm = {}'.format(ut.get_tmelt(primer)))
    
    print()
    print('Time Elapsed = {:.2f} sec'.format(tt.time()-t0))


if __name__ == '__main__':
    main()