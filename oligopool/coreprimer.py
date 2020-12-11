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
    exmotifs):

    '''
    Local
    '''

    return ut.get_motif_conflict(
        seq=primer,
        seqlen=len(primer),
        exmotifs=exmotifs,
        partial=True,
        checkall=False)

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
    primerspan,
    pairedprimer,
    pairedspan,
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
    
    # For Padding
    if primerspan:
        primer = primer[-primerspan:]

    if dimertype == 0:
        pairedprimer = primer
    
    # For Padding
    if pairedspan:
        pairedprimer = pairedprimer[-pairedspan:]

    s1, s2, eg = folder.evaluate_mfe_dimer(
        seq1=primer,
        seq2=pairedprimer)

    status, stloc = is_structure_feasible(
        struct1=s1,
        struct2=s2,
        energy=eg)

    if revert:
        stloc = primerlen - 1 - stloc

    # if not status:
    #     # if dimertype == 1:
    #     print()
    #     print((s1, s2, eg, stloc))
    #     print('(\'{}\''.format(primer))

    return status, stloc

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

def evaluate_seq_constraint(
    primerseq,
    exmotifs,
    liner):

    t0        = tt.time()
    primerlen = len(primerseq)

    liner.send(' Conflicting Motifs:\n')
    
    status, pmotif = ut.get_motif_conflict(
        seq=primerseq,
        seqlen=primerlen,
        exmotifs=exmotifs,
        partial=False,
        checkall=True)
    
    if not status:
        for motif in pmotif:
            if len(motif) < primerlen:
                liner.send(
                    '  Sequence Contains Motif {}\n'.format(
                        motif))
            else:
                liner.send(
                    '  Sequence Contained in Motif {}\n'.format(
                        motif))
    else:
        liner.send('  No Conflicting Motifs Found\n')
    
    liner.send(' Time Elapsed: {:.2f} sec\n'.format(
        tt.time()-t0))
    
    return (status, pmotif)

def evaluate_typeIIS_constraint(
    typeIIS,
    levalspan,
    revalspan,
    liner):

    t0 = tt.time()
    typeIISlen = len(typeIIS)

    liner.send(
        ' Minimum Forward Padding Length: {} bp\n'.format(
            levalspan))
    liner.send(
        ' Minimum Reverse Padding Length: {} bp\n'.format(
            revalspan))
    liner.send(
        ' Required Type IIS Motif Length: {} bp\n'.format(
            typeIISlen))

    fwdavail      = levalspan - typeIISlen
    revavail      = revalspan - typeIISlen
    typeIISstatus = True

    for locavail,loc in ((fwdavail, 'Forward'), (revavail, 'Reverse')):
        liner.send(
            '  {} Pad 5\' End Gap Length: {} bp\n'.format(
                loc, locavail))
        if locavail >= 6:
            liner.send('   Gap Length ≥ 6 bp? Yes\n')
        else:
            liner.send('   Gap Length ≥ 6 bp? No\n')
            typeIISstatus = typeIISstatus and False
    
    liner.send(' Time Elapsed: {:.2f} sec\n'.format(
        tt.time()-t0))
    
    return typeIISstatus, fwdavail, revavail
