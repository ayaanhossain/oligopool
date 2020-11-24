import time  as tt

import collections as cx

import numpy as np

import nrpcalc    as nr

import background as bk
import utils      as ut

# NRPCalc Fold Object
folder = nr.base.utils.Fold(
    temp=37.0,
    dangles=2,
    part_type='DNA')

def is_background_feasible(
    primer,
    background):
    '''
    Local
    '''
    
    if (background is None) or \
       (len(primer) < background.K):
        return True, None

    if primer[-background.K:] in background:
        return False, len(primer)-1

    return True, None

def is_tmelt_feasible(
    primer,
    primerlen,
    mintmelt,
    maxtmelt):
    '''
    Can be Local
    '''

    if len(primer) < primerlen:
        return True, None

    tmelt = ut.get_tmelt(
        seq=primer)

    if tmelt < mintmelt:
        return False, 0

    if tmelt > maxtmelt:
        return False, 1

    return True, None

def get_tmelt_traceback(
    primer,
    failtype):

    tidx = len(primer)
    
    if failtype == 0:
        lastweak = max(
            primer.rfind('A'),
            primer.rfind('T'))
        return min(tidx, lastweak)

    if failtype == 1:
        laststrong = max(
            primer.rfind('G'),
            primer.rfind('C'))
        return min(tidx, laststrong)

def is_deltatmelt_feasible(
    primer,
    pairedtmelt,
    deltatment):
    
    '''
    Uses tmelt feasibility
    '''
    
    pass

def is_motif_feasible(
    primer,
    exmotifs,
    partial):

    '''
    Local
    '''

    return ut.get_motif_conflict(
        seq=primer,
        seqlen=len(primer),
        exmotifs=exmotifs,
        partial=partial)

def is_context_feasible(
    primer,
    primerlen,
    prefixgroup,
    suffixgroup):    
    '''
    Local
    '''

    # Resolve Context on Left
    if len(primer) in prefixgroup:
        if primer in prefixgroup[len(primer)]:
            return False, len(primer) - 1

    # Resolve Context on Right
    if len(primer) == primerlen:
        for mlen in suffixgroup:
            primersuffix = primer[-mlen:]
            if primersuffix in suffixgroup[mlen]:
                return False, len(primer)-mlen

    # All OK
    return True, None

def is_structure_feasible(
    struct1,
    struct2,
    energy):
    '''
    Local
    '''

    minDG  = -8 # kcal/mol
    free3  = 4 # last free bases
    status = True

    if energy <= minDG:
        status = False

    stloc = struct1.find('(')

    if ')' in struct2[-free3:]:
        status = False

    if '(' in struct1[-free3:]:
        status = False

    return status, stloc

def is_dimer_feasible(
    primer,
    primertype,
    primerlen,
    pairedprimer,
    dimertype):
    '''
    Cofold business
    '''
    
    if len(primer) < primerlen:
        return True, None

    if (dimertype == 1) and \
       (pairedprimer is None):
        return True, None

    revert = False
    if primertype == 1:
        primer = ut.get_revcomp(
            seq=primer)
        revert = True

    if dimertype == 0:
        pairedprimer = primer

    s1, s2, eg = folder.evaluate_mfe_dimer(
        seq1=primer,
        seq2=pairedprimer)

    status, stloc = is_structure_feasible(
        struct1=s1,
        struct2=s2,
        energy=eg)

    if revert:
        stloc = primerlen - 1 - stloc

    if not status:
        # if dimertype == 1:
        print()
        print((s1, s2, eg, stloc))
        print('(\'{}\''.format(primer))

    return status, stloc

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
    bcond, bfail = is_background_feasible(
        primer=primer,
        background=background)
    if not bcond:
        liner.send(
            ' Candidate: {} Rejected due to Background Infeasibility'.format(
                primer))
        return False, bfail

    # Melting Temperature
    tcond, tfail = is_tmelt_feasible(
        primer=primer,
        primerlen=primerlen,
        mintmelt=mintmelt,
        maxtmelt=maxtmelt)
    if not tcond:
        liner.send(
            ' Candidate: {} Rejected due to Tm Infeasibility'.format(
                primer))
        return False, get_tmelt_traceback(
            primer=primer,
            failtype=tfail)

    # Motif Embedding
    mcond, motif = is_motif_feasible(
        primer=primer,
        exmotifs=exmotifs,
        partial=True)
    if not mcond:
        liner.send(
            ' Candidate: {} Rejected due to Motif Infeasibility'.format(
                primer))
        return False, max(0, len(primer)-len(motif))

    # Context Feasibility
    ccond, ctloc = is_context_feasible(
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
    hcond, htloc = is_dimer_feasible(
        primer=primer,
        primertype=primertype,
        primerlen=primerlen,
        pairedprimer=None,
        dimertype=0)
    if not hcond:
        liner.send(
            ' Candidate: {} Rejected due to Homodimer Infeasibility'.format(
                primer))
        return False, htloc

    # Heterodimer Feasibility
    qcond, qtloc = is_dimer_feasible(
        primer=primer,
        primertype=primertype,
        primerlen=primerlen,
        pairedprimer=pairedprimer,
        dimertype=1)
    if not qcond:
        liner.send(
            ' Candidate: {} Rejected due to Heterodimer Infeasibility\n'.format(
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

def is_primer_global_feasible(
    primer,
    exmotifs,
    lcifn,
    rcifn,
    cntxlen,
    liner):
    
    # Motif Embedding
    mcond, motif = is_motif_feasible(
        primer=primer,
        exmotifs=exmotifs,
        partial=True)
    if not mcond:
        liner.send(
            ' Candidate: {} Rejected due to Motif Infeasibility'.format(
                primer))
        return False

    # Primer Assignment


    return True

def stream_motif_splits(motif):
    return (
        (motif[:i], motif[i:]) \
        for i in range(1, len(motif)-1))

def get_exmotif_partition(exmotifs):
    
    partition = cx.defaultdict(list)

    for motif in exmotifs:
        for u,v in stream_motif_splits(motif):
            partition[u].append(v)

    return partition

def get_inverted_exmotif_partition(partition):

    inv_partition = cx.defaultdict(list)

    for u,v in partition.items():
        for vi in v:
            inv_partition[vi].append(u)

    return inv_partition

def get_grouped_edge_constraints(
    edgeset):

    edgedict = cx.defaultdict(set)
    for edge in sorted(edgeset, key=len):
        edgedict[len(edge)].add(edge)
    return edgedict

def get_edge_contraints(
    lcnum,
    rcnum,
    lcifn,
    rcifn,
    exmotifs,
    liner):

    # Book-keeping
    px  = set() # Prefixes to avoid
    sx  = set() # Suffixes to avoid
    ps  = False # Prefix Status
    ss  = False # Suffix Status
    lcp = {}    # Left  Context Partition
    rcp = {}    # Right Context Parition

    # Time-keeping
    t0 = tt.time()

    # Do we have any context?
    if not (lcnum or rcnum):
        liner.send(' Context Specified? No\n')
        liner.send(' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))
        return (False, {}, {}, lcp, rcp) # We have no constraints!
    else:
        liner.send(' Context Specified? Yes\n')

    # Build Left Context Partition
    lcp = get_exmotif_partition(
        exmotifs=exmotifs)
    
    # Left Context Checks
    if lcnum:

        # Show Updates
        liner.send(' Checking Left Context ...')

        # Compute Conficts
        for i in range(lcnum):
            lcseq = lcifn(i)
            for suffix in lcp:
                if lcseq.endswith(suffix):
                    # print((lcseq, suffix, lcp[suffix]))
                    px.update(lcp[suffix])
                    ps = True # Prefix Constraints Exist!

        # Show Updates
        liner.send(' Found {} Primer Prefix Constraints\n'.format(
            len(px)))

    # Build Right Context Parition
    rcp = get_inverted_exmotif_partition(
        partition=lcp)

    # Right Context Checks
    if rcnum:

        # Show Updates
        liner.send(' Checking Right Context ...')

        # Compute Conflicts
        for i in range(rcnum):
            rcseq = rcifn(i)
            for prefix in rcp:
                if rcseq.startswith(prefix):
                    # print((rcseq, prefix, rcp[prefix]))
                    sx.update(rcp[prefix])
                    ss = True # Suffix Constraints Exist!

        # Show Updates
        liner.send(' Found {} Primer Suffix Constraints\n'.format(
            len(sx)))

    # Final Update
    liner.send(' Time Elapsed: {:.2f} sec\n'.format(
        tt.time()-t0))

    # Group Edge Constraints by Length
    if px:
        px = get_grouped_edge_constraints(
            edgeset=px)
    if sx:
        sx = get_grouped_edge_constraints(
            edgeset=sx)
    
    # Return Results
    return (ps or ss, px, sx, lcp, rcp)

def evaluate_edge_constraint(
    edgedict,
    culpdict,
    contexttype,
    liner):

    # Book-keeping
    infedges = set()
    state    = True

    for edgelen in edgedict:
        if len(edgedict[edgelen]) == 4**edgelen:
            for motif in edgedict[edgelen]:
                for counterpart in culpdict[motif]:
                    liner.send(
                        '  Sequences {} {}\n'.format(
                            ['Ending in', 'Starting with'][contexttype],
                            counterpart))
                infedges.add(counterpart)
                state = False

    if state:
        liner.send('  No Infeasible Sequences Found\n')

    return (state, infedges)

def get_primer_extreme(
    primerseq,
    exttype):
    
    extbases = []

    for ib in primerseq:
        space = list(ut.dna_space[ib])
        extrm = list(ut.dna_space[
            ['W', 'S'][exttype]].intersection(space))
        if extrm:
            extbases.append(np.random.choice(extrm))
        else:
            extbases.append(np.random.choice(space))
    
    return ''.join(extbases)

def evaluate_tmelt_constraint(
    primerseq,
    mintmelt,
    maxtmelt,
    liner):
    
    # Book-keeping
    extmintmelt = float('inf')  # Minimum Feasible Tm
    extmaxtmelt = float('-inf') # Maximum Feasible Tm
    status = True               # Feasibility Status

    # Time-keeping
    t0 = tt.time()
    
    # Estimate Minimum Feasible Tm
    for i in range(100 * len(primerseq)):
        minprimer = get_primer_extreme(
            primerseq=primerseq,
            exttype=0)
        extmintmelt = min(extmintmelt, ut.get_tmelt(
            seq=minprimer))
    liner.send(
        ' Estimated Minimum Tm = {:.2f} C\n'.format(
            extmintmelt))

    # Estimate Maximum Feasible Tm
    for i in range(100 * len(primerseq)):
        maxprimer = get_primer_extreme(
            primerseq=primerseq,
            exttype=1)
        extmaxtmelt = max(extmaxtmelt, ut.get_tmelt(
            seq=maxprimer))
    liner.send(
        ' Estimated Maximum Tm = {:.2f} C\n'.format(
            extmaxtmelt))

    # Minimum Feasibility
    if extmaxtmelt < mintmelt:
        liner.send(
            '  Required Minimum Tm of {:.2f} C > Estimated Maximum Tm Infeasible\n'.format(
                mintmelt))
        status = status and False
    else:
        liner.send(
            '  Required Minimum Tm of {:.2f} C Feasible\n'.format(
                mintmelt))

    # Maximum Feasibility
    if extmintmelt > maxtmelt:
        liner.send(
            '  Required Maximum Tm of {:.2f} C > Estimated Minimum Tm Infeasible\n'.format(
                maxtmelt))
        status = status and False
    else:
        liner.send(
            '  Required Maximum Tm of {:.2f} C Feasible\n'.format(
                maxtmelt))

    # Show Final Update
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))
    
    # Return Results
    return (status, extmintmelt, extmaxtmelt)

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

    # Determine Tm Feasibility
    liner.send('\n[Checking Primer Tm Feasibility]\n')
    (tmeltstatus,
     emintmelt,
     emaxtmelt) = evaluate_tmelt_constraint(
        primerseq=primerseq,
        mintmelt=mintmelt,
        maxtmelt=maxtmelt,
        liner=liner)

    if not tmeltstatus:
        liner.send(' Verdict: Primer Design Infeasible due to Tm Constraints\n')
        return (False, None)
    else:
        liner.send(' Verdict: Primer Design Possibly Feasible\n')

    # Context Setup
    lcnum = ut.get_context_num(          # Total Number of Left Context
        context=leftcontext)
    rcnum = ut.get_context_num(          # Total Number of Right Context
        context=rightcontext)
    # tcnum = max(lcnum, rcnum)            # Total Number of Unique Contexts
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

    # Extract Edge Effect Constraints
    liner.send('\n[Extracting Edge Constraints]\n')
    (edgestatus,
     prefixgroup,
     suffixgroup,
     prefixdict,
     suffixdict) = get_edge_contraints(
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
            lcf, sxfails = evaluate_edge_constraint(
                edgedict=prefixgroup,
                culpdict=suffixdict,
                contexttype=0,
                liner=liner)

        if rcnum:
            liner.send(' Infeasible Right Context Sequences:\n')
            rcf, pxfails = evaluate_edge_constraint(
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

    # global_model_fn = lambda primer: is_primer_global_feasible(
    #     primer=primer,
    #     exmotifs=exmotifs,
    #     lcifn=lcifn,
    #     rcifn=rcifn,
    #     cntxlen=cntxlen,
    #     liner=liner)
    global_model_fn = lambda primer: True

    # Define Maker Instance
    maker = nr.base.maker.NRPMaker(
        part_type='DNA',
        seed=None)

    # Show Update
    liner.send('\n[Computing Primer]\n')

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
        global_model_fn=global_model_fn,
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
    print(folder.evaluate_mfe_dimer(
        seq1='GTGCCTCGTCCGAAACCA',
        seq2=ut.get_revcomp('CTGCTAACACACGCCCCA')))

    return


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

if __name__ == '__main__':
    main()