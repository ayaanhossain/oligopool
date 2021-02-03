import time  as tt

import collections as cx

import numpy as np
from numpy.core.fromnumeric import repeat

import nrpcalc    as nr

import background as bk
import vectordb   as db
import utils      as ut


# NRPCalc Fold Object
folder = nr.base.utils.Fold(
    temp=37.0,
    dangles=2,
    part_type='DNA')

# Parser and Setup Functions

def get_constant_regions(primerseq):
    '''
    Extract all constant, non-degenerate
    regions degined in primerseq.
    Internal use only.

    :: primerseq
       type - string
       desc - primer sequence constraint
    '''

    # Make a local copy
    cprimerseq = str(primerseq)

    # Remove dgenerate nucleotides
    for nt in set(ut.ddna_space.keys()) - set('ATGC'):
        if nt in cprimerseq:
            cprimerseq = cprimerseq.replace(nt, '-')

    # Split and extract defined regions
    regions = set(cprimerseq.split('-'))
    if '' in regions:
        regions.remove('')

    # Return results
    return regions

def get_palindromes(seq):
    '''
    Return all palindromic hexamers
    in sequence. Internal use only.

    :: seq
       type - string
       desc - sequence to analyze
    '''

    palindromes = set()
    if len(seq) >= 6:
        for kmer in ut.stream_spectrum(seq=seq, k=6):
            if kmer == ut.get_revcomp(kmer):
                palindromes.add(kmer) # Found a 6-mer core palindrome
    return palindromes

def get_palindrome_conflicts(regions):
    '''
    Check if regions contains hexamer
    palindromes. Internal use only.

    :: regions
       type - list
       desc - list of defined regions
              in a primer sequence
    '''

    # Book-keeping
    palindromes = set()

    # Check for palindromic regions
    for region in regions:
        palindromes.update(get_palindromes(
        seq=region))

    # Return results
    if not palindromes:
        return (True, None) # We're OK
    else:
        return (False, tuple(palindromes))

def get_internal_repeat_conflicts(regions):
    '''
    Check if regions contain internal
    repeats. Internal use only.

    :: regions
       type - list
       desc - list of defined regions
              in a primer sequence
    '''

    # Book-keeping
    internalrepeats = cx.Counter()
    status = True # No Conflict

    # Compute internalrepeats
    for region in regions:
        internalrepeats.update(ut.stream_canon_spectrum(
            seq=region, k=6))

    # Compute results
    if internalrepeats:
        status = max(internalrepeats.values()) == 1

    if not status:
        internalrepeats = tuple([repeat for repeat in internalrepeats \
            if internalrepeats[repeat] > 1])

    # Return results
    if status:
        return (True, None)
    else:
        return (status, internalrepeats)

def get_paired_repeat_conflicts(
    primerseq,
    primertype,
    pairedprimer):
    '''
    Check if primerseq and pairedprimer
    share repeats. Internal use only,

    :: primerseq
       type - string
       desc - primer sequence constraint
    :: primertype
       type - integer
       desc - primer type idenfier
              0 = forward primer
              1 = reverse primer
    :: pairedprimer
       type - string / None
       desc - primer sequence paired with
              current primer design job
    '''

    # Book-keeping
    p = primerseq
    q = pairedprimer
    t = primertype

    # Correct orientation
    if t == 0: # p == FWD, q == REV in FWD Format
        p = p
        q = q
        # Note: A direct repeat is a hetero conflict
    else:      # p == REV, q == FWD in FWD Format
        p = p                 # or p = ut.get_revcomp(p)
        q = ut.get_revcomp(q) # or q = q
        # Note: An inverted repeat is a hetero conflict

    # Dynamic programming computation
    # of Longest Common Substring
    mat = np.zeros((len(p)+1, len(q)+1))
    for i in range(1, len(p)+1):
        for j in range(1, len(q)+1):
            if p[i-1] == q[j-1]:
                mat[i, j] = mat[i-1, j-1] + 1
    lcs = int(np.amax(mat))
    pairedrepeats = tuple(p[e-lcs:e] \
        for e in np.where(mat == np.amax(mat))[0])

    # Return results
    return (lcs <= 5, pairedrepeats)

def get_parsed_sequence_constraint(
    primerseq,
    primertype,
    exmotifs,
    pairedprimer,
    liner):
    '''
    Check primer sequence feasibility.
    Internal use only.

    :: primerseq
       type - string
       desc - primer sequence constraint
    :: primertype
       type - integer
       desc - primer type idenfier
              0 = forward primer
              1 = reverse primer
    :: exmotifs
       type - deque / None
       desc - deque of all motifs
              to be excluded
    :: pairedprimer
       type - string / None
       desc - primer sequence paired with
              current primer design job
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Time-keeping
    t0 = tt.time()

    # Design Space Analysis
    liner.send(' Computing Design Space ...')

    dspace = 1
    for nt in primerseq:
        dspace *= len(ut.ddna_space[nt])
    sntn, plen = ut.get_notelen(
        printlen=ut.get_printlen(
            value=dspace))

    # Show Update
    if dspace == 1:
        dspace_ok = False
        liner.send(
            ' Design Space: 1 Possible Primer [INFEASIBLE] (Degenerate Nucleotides Absent)\n')
    else:
        dspace_ok = True
        liner.send(
            ' Design Space: {:{},{}} Possible Primer(s)\n'.format(
                dspace,
                plen,
                sntn))

    # Exmotifs Analysis
    liner.send(' Computing Motif Conflicts ...')

    motif_ok, excludedmotifs = ut.get_motif_conflict(
        seq=primerseq,
        seqlen=len(primerseq),
        exmotifs=exmotifs,
        partial=False,
        checkall=True)

    # Show Update
    if not motif_ok:
        excludedmotifs = tuple(excludedmotifs)
        liner.send(
            ' Found {:,} Excluded Motif(s)\n'.format(
                len(excludedmotifs)))
        plen = max(map(len, excludedmotifs)) + 2
        for motif in excludedmotifs:
            motif = '\'{}\''.format(motif)
            liner.send(
                '   - Excluded Motif {:>{}} Present [INFEASIBLE]\n'.format(
                    motif,
                    plen))
    else:
        liner.send(
            ' Found 0 Excluded Motif(s)\n')

     # Region Analysis
    liner.send(' Computing Constant Regions ...')

    regions = get_constant_regions(
        primerseq=primerseq)

    # Internal Repeat Analysis
    liner.send(' Computing Internal Repeat Conflicts ...')

    (repeat_ok,
    internalrepeats) = get_internal_repeat_conflicts(
        regions=regions)

    # Show Update
    if not repeat_ok:
        liner.send(
            ' Found {:,} Internal Repeat(s)\n'.format(
                len(internalrepeats)))
        plen = max(map(len, internalrepeats)) + 2
        for repeat in internalrepeats:
            repeat = '\'{}\''.format(repeat)
            liner.send(
                '   - Internal Repeat {:>{}} Present [INFEASIBLE] (Strong Homodimer Motif)\n'.format(
                    repeat,
                    plen))
    else:
        liner.send(
            ' Found 0 Internal Repeat(s)\n')

    # Palindrome Analysis
    liner.send(' Computing Palindrome Conflicts ...')

    palindrome_ok, palindromes = get_palindrome_conflicts(
        regions=regions)

    # Show Update
    if not palindrome_ok:
        liner.send(
            ' Found {:,} Palindrome(s)\n'.format(
                len(palindromes)))
        plen = max(map(len, palindromes)) + 2
        for palindrome in palindromes:
            palindrome = '\'{}\''.format(palindrome)
            liner.send(
                '   - Palindrome {:>{}} Present [INFEASIBLE] (Strong Homodimer Motif)\n'.format(
                    palindrome,
                    plen))
    else:
        liner.send(
            ' Found 0 Palindrome(s)\n')

    # Maximum Paired Repeat Analysis
    if not pairedprimer is None:
        liner.send(' Computing Paired Repeat Conflicts ...')

        (repeat_ok,
        pairedrepeats) = get_paired_repeat_conflicts(
            primerseq=primerseq,
            primertype=primertype,
            pairedprimer=pairedprimer)
    else:
        repeat_ok     = True
        pairedrepeats = None

    # Show Update
    if not repeat_ok:
        liner.send(
            ' Found {:,} Paired Repeat(s)\n'.format(
                len(pairedrepeats)))
        plen = max(map(len, pairedrepeats)) + 2
        for repeat in pairedrepeats:
            repeat = '\'{}\''.format(repeat)
            liner.send(
                '   - Paired Repeat {:>{}} Present [INFEASIBLE] (Strong Heterodimer Motif)\n'.format(
                    repeat,
                    plen))
    else:
        liner.send(
            ' Found 0 Paired Repeat(s)\n')

    # Show Time Elapsed
    liner.send(' Time Elapsed: {:.2f} sec\n'.format(
        tt.time()-t0))

    # Compute primerseq Feasibility
    parsestatus = all([
        dspace_ok,
        motif_ok,
        repeat_ok,
        palindrome_ok])

    if not parsestatus:
        liner.send(
            ' Verdict: Primer Design Infeasible due to Sequence Constraint\n')
    else:
        liner.send(
            ' Verdict: Primer Design Possibly Feasible\n')

    # Return Results
    return (parsestatus,
        dspace,
        excludedmotifs,
        internalrepeats,
        palindromes,
        pairedrepeats)

def get_primer_extreme(
    primerseq,
    exttype):
    '''
    Get candidate sequences with strong
    and weak bases for evaluation.
    Internal use only.

    :: primerseq
       type - string
       desc - primer sequence constraint
    :: exttype
       type - integer
       desc - base group identifier
    '''

    extbases = []
    for ib in primerseq:
        space = list(ut.ddna_space[ib])
        extrm = list(ut.ddna_space[
            ['W', 'S'][exttype]].intersection(space))
        if extrm:
            extbases.append(np.random.choice(extrm))
        else:
            extbases.append(np.random.choice(space))
    return ''.join(extbases)

def get_parsed_tmelt_constraint(
    primerseq,
    pairedprimer,
    mintmelt,
    maxtmelt,
    liner):
    '''
    Determine melting temperature bounds for
    primerseq, subject to pairedprimer and
    mintmelt and maxtmelt. Internal use only.

    :: primerseq
       type - string
       desc - primer sequence constraint
    :: pairedprimer
       type - string
       desc - full paired primer sequence
    :: mintmelt
       type - float
       desc - melting temperature lowerbound
    :: maxtmelt
       type - float
       desc - melting temperature upperbound
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    posminTm    = float('inf')  # Minimum Possible Tm
    posmaxTm    = float('-inf') # Maximum Possible Tm
    minTm       = None          # Updated Minimum Tm
    maxTm       = None          # Updated Maximum Tm
    lowermaxTm  = False         # Maximum Tm Underflow
    higherminTm = False         # Minimum Tm Overflow
    parsestatus = True          # Feasibility parsestatus
    deltatmelt  = 1             # Tm Difference Tolerance

    # Time-keeping
    t0 = tt.time()

    # Required Melting Temperature Range
    if pairedprimer is None:
        # As is ... no update
        liner.send('  Parameter Source: From Input Argument(s)\n')
        liner.send(
            ' Required Tm Range: {:.2f} to {:.2f} °C\n'.format(
                mintmelt,
                maxtmelt))
    else:
        # Updated Tm
        pairedtmelt = ut.get_tmelt(
            seq=pairedprimer)
        mintmelt = pairedtmelt - deltatmelt
        maxtmelt = pairedtmelt + deltatmelt

        liner.send('  Parameter Source: Altered to Match Paired Primer\n')
        liner.send(
            '   Paired Tm Range: {:.2f} to {:.2f} °C\n'.format(
                mintmelt,
                maxtmelt))

    # Estimate Tm Range
    liner.send('Estimating Feasible Melting Temperature Range ...')

    # Estimate Minimum Feasible Tm
    for i in range(100 * len(primerseq)):
        minprimer = get_primer_extreme(
            primerseq=primerseq,
            exttype=0)
        posminTm = min(posminTm, ut.get_tmelt(
            seq=minprimer))
        liner.send(
            ' Estimated Minimum Tm: {:.2f} °C'.format(
                posminTm))

    # Estimate Maximum Feasible Tm
    for i in range(100 * len(primerseq)):
        maxprimer = get_primer_extreme(
            primerseq=primerseq,
            exttype=1)
        posmaxTm = max(posmaxTm, ut.get_tmelt(
            seq=maxprimer))
        liner.send(
            ' Estimated Maximum Tm: {:.2f} °C'.format(
                posmaxTm))

    # Show Update
    liner.send(
        ' Possible Tm Range: {:.2f} to {:.2f} °C\n'.format(
            posminTm,
            posmaxTm))

    # Compute Tm Conflicts
    conflictcount = 0
    if posminTm > maxtmelt:
        conflictcount += 1
    if posmaxTm < mintmelt:
        conflictcount += 1

    # Found Tm Conflicts?
    if conflictcount:

        # Show Update
        liner.send(
            ' Found {:,} Melting Temperature Conflict(s)\n'.format(
                conflictcount))

        # Minimum Feasibility
        if posmaxTm < mintmelt:
            liner.send(
                '   - Required Min Tm of {:.2f} °C is Higher than Possible Max Tm of {:.2f} °C [INFEASIBLE]\n'.format(
                    mintmelt,
                    posmaxTm))
            higherminTm = True
            parsestatus = parsestatus and False

        # Maximum Feasibility
        if posminTm > maxtmelt:
            liner.send(
                '   - Required Max Tm of {:.2f} °C is Lower than Possible Min Tm of {:.2f} °C [INFEASIBLE]\n'.format(
                    maxtmelt,
                    posminTm))
            lowermaxTm  = True
            parsestatus = parsestatus and False

    # No Tm Conflicts
    else:
        minTm = max(mintmelt, posminTm)
        maxTm = min(maxtmelt, posmaxTm)

        # Show Update
        liner.send(
            '  Updated Tm Range: {:.2f} to {:.2f} °C\n'.format(
                minTm,
                maxTm))
        liner.send(
            ' Found 0 Melting Temperature Conflict(s)\n')

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Show Final Update
    if parsestatus:
        liner.send(
            ' Verdict: Primer Design Possibly Feasible\n')
    else:
        liner.send(
            ' Verdict: Primer Design Infeasible due to Melting Temperature Constraints\n')

    # Return Results
    return (parsestatus,
        posminTm,
        posmaxTm,
        higherminTm,
        lowermaxTm,
        minTm,
        maxTm)

def stream_motif_splits(motif):
    '''
    Stream all prefix-suffix splits
    for given motif. This is also
    left partition of motif.
    Internal use only.

    :: motif
       type - string
       desc - motif to split
    '''

    return ((motif[:i], motif[i:]) \
        for i in range(1, len(motif)))

def get_exmotif_partition(exmotifs):
    '''
    Compute the left partition of
    exmotifs. Internal use only.

    :: exmotifs
       type - cx.deque
       desc - deque of all motifs
              to be excluded
    '''

    partition = cx.defaultdict(list)
    for motif in exmotifs:
        for u,v in stream_motif_splits(
            motif=motif):
            partition[u].append(v)
    return partition

def get_inverted_exmotif_partition(partition):
    '''
    Compute the right partition of
    left partitioned exmotifs by
    inverting left partition.
    Internal use only.

    :: partition
       type - dict
       desc - dictionary of left
              partitioned exmotifs,
              prefix -> suffix
    '''

    inv_partition = cx.defaultdict(list)
    for u,v in partition.items():
        for vi in v:
            inv_partition[vi].append(u)
    return inv_partition

def get_parsed_edgeeffects(
    primerseq,
    leftcontext,
    rightcontext,
    exmotifs,
    liner):
    '''
    Check edge effect feasibility.
    Internal use only.

    :: primerseq
       type - string
       desc - primer sequence constraint
    :: leftcontext
       type - tuple
       desc - tuple of all unique left
              context sequences
    :: rightcontext
       type - tuple
       desc - tuple of all unique right
              context sequences
    :: exmotifs
       type - cx.deque
       desc - deque of all motifs
              to be excluded
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    prefixforbidden = None # Prefixes to avoid
    suffixforbidden = None # Suffixes to avoid
    leftpartition   = None # Left  Context Partition
    rightpartition  = None # Right Context Parition
    prefixstatus    = True # No Prefixes Forbidden
    suffixstatus    = True # No Suffixes Forbidden

    probseqprefixes = None # Problematic Sequence Prefixes
    probprefixlens  = None # Problematic Sequence Prefix Lengths

    probseqsuffixes = None # Problematic Sequence suffixes
    probsuffixlens  = None # Problematic Sequence suffix Lengths

    # Time-keeping
    t0 = tt.time()

    # Show Update
    liner.send(' Left Paritioning Excluded Motifs ...')

    # Left Partition exmotifs
    leftpartition = get_exmotif_partition(
        exmotifs=exmotifs)

    # Show Update
    liner.send(' Right Paritioning Excluded Motifs ...')

    # Right Partition exmotifs
    rightpartition = get_inverted_exmotif_partition(
        partition=leftpartition)

    # Assess Left Edge Effects
    if leftcontext:

        # Show Updates
        liner.send(' Parsing Left Context ...')

        # Compute Prefix Conficts
        prefixforbidden = set()
        for lcseq in leftcontext:
            for prefix in leftpartition:
                if lcseq.endswith(prefix):
                    # print((lcseq, prefix, leftpartition[prefix]))
                    prefixforbidden.update(leftpartition[prefix])

        # Show Updates
        liner.send(' Forbidden Prefix: {:,} Excluded Motif Prefix(es)\n'.format(
            len(prefixforbidden)))

        # Show Updates
        liner.send(' Computng Primer Sequence Prefix Conflict ...')

        # Compute Sequence Conflicts
        probseqprefixes = []
        for prefix in prefixforbidden:
            if primerseq.startswith(prefix):
                probseqprefixes.append(prefix)

        # Show Updates
        if probseqprefixes:
            probseqprefixes = tuple(probseqprefixes)
            prefixstatus    = False

            liner.send(
                ' Found {:,} Primer Sequence Prefix Conflict(s)\n'.format(
                    len(probseqprefixes)))

            plen = max(map(len, probseqprefixes)) + 2

            for prefix in probseqprefixes:
                prefix = '\'{}\''.format(prefix)
                liner.send(
                    '   - Primer Starts with {:>{}} [INFEASIBLE] (Forbidden Prefix)\n'.format(
                        prefix,
                        plen))
        else:
            probseqprefixes = None
            liner.send(' Found 0 Primer Sequence Prefix Conflict(s)\n')

        # Show Updates
        liner.send(' Computing Primer Sequence Prefix Length Distribution ...')

        # Group Forbidden Prefixes
        prefixforbidden = ut.get_grouped_sequences(
            sequences=prefixforbidden)

        # Compute Length Feasibility
        liner.send(' Primer Sequence Forbidden Prefix Length Distribution\n')

        probprefixlens = []

        plen = ut.get_printlen(
            value=max(prefixforbidden.keys()))
        qlen = ut.get_printlen(
            value=max(map(
                len, prefixforbidden.values())))

        for prefixlen in prefixforbidden:

            # Compute infeasibility
            parsemsg = ''
            if len(prefixforbidden[prefixlen]) == 4**prefixlen:
                parsemsg     = ' [INFEASIBLE] (All {}-mers Forbidden)'.format(prefixlen)
                prefixstatus = False
                probprefixlens.append(prefixlen)

            # Show update
            liner.send('   - {:{},d} Prefix(es) of Length {:{},d}{}\n'.format(
                len(prefixforbidden[prefixlen]),
                plen,
                prefixlen,
                qlen,
                parsemsg))

        # Finalize probprefixlens
        probprefixlens = None if len(probprefixlens) == 0 else probprefixlens

    # Assess Right Edge Effects
    if rightcontext:

        # Show Updates
        liner.send(' Parsing Right Context ...')

        # Compute Suffix Conficts
        suffixforbidden = set()
        for rcseq in rightcontext:
            for suffix in rightpartition:
                if rcseq.startswith(suffix):
                    # print((rcseq, suffix, rightpartition[suffix]))
                    suffixforbidden.update(rightpartition[suffix])

        # Show Updates
        liner.send(' Forbidden Suffix: {:,} Excluded Motif Suffix(es)\n'.format(
            len(suffixforbidden)))


        # Show Updates
        liner.send(' Computng Primer Sequence Suffix Conflict ...')

        # Compute Sequence Conflicts
        probseqsuffixes = []
        for suffix in suffixforbidden:
            if primerseq.endswith(suffix):
                probseqsuffixes.append(suffix)

        # Show Updates
        if probseqsuffixes:
            probseqsuffixes = tuple(probseqsuffixes)
            suffixstatus    = False

            liner.send(
                ' Found {:,} Primer Sequence Suffix Conflict(s)\n'.format(
                    len(probseqsuffixes)))

            plen = max(map(len, probseqsuffixes)) + 2

            for suffix in probseqsuffixes:
                suffix = '\'{}\''.format(suffix)
                liner.send(
                    '   - Primer Ends with {:>{}} [INFEASIBLE] (Forbidden Suffix)\n'.format(
                        suffix,
                        plen))
        else:
            probseqsuffixes = None
            liner.send(' Found 0 Primer Sequence Suffix Conflict(s)\n')

        # Show Updates
        liner.send(' Computing Primer Sequence Suffix Length Distribution ...')

        # Group Forbidden Suffixes
        suffixforbidden = ut.get_grouped_sequences(
            sequences=suffixforbidden)

        # Compute Length Feasibility
        liner.send(' Primer Sequence Forbidden Suffix Length Distribution\n')

        probsuffixlens = []

        plen = ut.get_printlen(
            value=max(suffixforbidden.keys()))
        qlen = ut.get_printlen(
            value=max(map(
                len, suffixforbidden.values())))

        for suffixlen in suffixforbidden:

            # Compute infeasibility
            parsemsg = ''
            if len(suffixforbidden[suffixlen]) == 4**suffixlen:
                parsemsg     = ' [INFEASIBLE] (All {}-mers Forbidden)'.format(suffixlen)
                suffixstatus = False
                probsuffixlens.append(suffixlen)

            # Show update
            liner.send('   - {:{},d} Suffix(es) of Length {:{},d}{}\n'.format(
                len(suffixforbidden[suffixlen]),
                plen,
                suffixlen,
                qlen,
                parsemsg))

        # Finalize probsuffixlens
        probsuffixlens = None if len(probsuffixlens) == 0 else probsuffixlens

    # Compute Final Results
    parsestatus = prefixstatus and suffixstatus

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Show Final Update
    if parsestatus:
        liner.send(
            ' Verdict: Primer Design Possibly Feasible\n')
    else:
        liner.send(
            ' Verdict: Primer Design Infeasible due to Edge Effect(s)\n')

    # Return Results
    return (parsestatus,
        probseqprefixes,
        probprefixlens,
        probseqsuffixes,
        probsuffixlens,
        prefixforbidden,
        suffixforbidden)

# Engine Objective and Helper Functions

def is_background_feasible(
    primer,
    background):
    '''
    Determine if primer contains a repeat
    with background. Internal use only.

    :: primer
       type - string
       desc - a partially explored primer
              sequence path
    :: background
       type - db.vectorDB
       desc - vectorDB instance storing
              background repeats
    '''

    # No Background or Too Short
    # a Primer Candidate
    if (background is None) or \
       (len(primer) < background.K):
        return True, None # No Conflict

    # Background Repeat Found?
    if primer[-background.K:] in background:
        return False, len(primer)-1 # Conflict

    # No Conflict
    return True, None

def is_oligopool_feasible(
    primer,
    maxreplen,
    oligorepeats):
    '''
    Determine if primer contains a repeat
    with oligopool. Internal use only.

    :: primer
       type - string
       desc - a partially explored primer
              sequence path
    :: oligorepeats
       type - set
       desc - set storing oligopool repeats
    '''

    # Too Short a Primer Candidate
    if len(primer) < maxreplen+1:
        return True, None # No Conflict

    # Oligopool Repeat Found?
    canonkmer = min(
        primer[-maxreplen+1:],
        ut.get_revcomp(
            seq=primer[-maxreplen+1:]))
    if canonkmer in oligorepeats:
        return False, len(primer)-1 # Conflict

    # No Conflict
    return True, None

def is_paired_feasible(
    primer,
    pairedrepeats):
    '''
    Determine if primer contains a repeat
    with paired primer. Internal use only.

    :: primer
       type - string
       desc - a partially explored primer
              sequence path
    :: pairedrepeats
       type - set
       desc - set storing paired primer
              repeats
    '''

    # No Paired Repeats or Too Short
    # a Primer Candidate
    if (pairedrepeats is None) or \
       (len(primer) < 6):
        return True, None # No Conflict

    # Paired Primer Repeat Found?
    canonkmer = min(
        primer[-6:],
        ut.get_revcomp(
            seq=primer[-6:]))
    if canonkmer in pairedrepeats:
        return False, len(primer)-1 # Conflict

    # No Conflict
    return True, None

def is_tmelt_feasible(
    primer,
    primerlen,
    mintmelt,
    maxtmelt):
    '''
    Determine if primer Tm in bounded range.
    Internal use only.

    :: primer
       type - string
       desc - a partially explored primer
              sequence path
    :: primerlen
       type - integer
       desc - full primer sequence length
    :: mintmelt
       type - float
       desc - melting temperature lower bound
    :: maxtmelt
       type - float
       desc - melting temperature upper bound
    '''

    # Too Short a Primer Candidate
    if len(primer) < primerlen:
        return True, None # No Conflict

    # Compute Melting Temperature
    tmelt = ut.get_tmelt(
        seq=primer)

    # Lower Bound Violated
    if tmelt < mintmelt:
        return False, 0   # Conflict Type 0

    # Upper Bound Violated
    if tmelt > maxtmelt:
        return False, 1   # Conflict Type 1

    # Melting Temprature Bounded
    return True, None     # No Conflict

def get_tmelt_traceback(
    primer,
    failtype):
    '''
    Compute cautious traceback location
    to meet melting temperature failure.
    Internal use only.

    :: primer
       type - string
       desc - a partially explored primer
              sequence path
    :: failtype
       type - integer
       desc - melting temperature failure
              type identifier
    '''

    # Define Default Traceback
    defaultidx = len(primer) - 1

    # Lower Bound Mitigation
    if failtype == 0:

        # Locate Latest Weak Base
        lastweakidx = max(
            primer.rfind('A'),
            primer.rfind('T'))

        # No Weak Bases
        if lastweakidx == -1:
            # Abort
            return 0
        else:
            # Traceback
            return min(
                defaultidx,
                lastweakidx)

    # Upper Bound Mitigation
    if failtype == 1:

        # Locate Latest Strong Base
        laststrongidx = max(
            primer.rfind('G'),
            primer.rfind('C'))

        # No Strong Base
        if laststrongidx == -1:
            # Abort
            return 0
        else:
            # Traceback
            return min(
                defaultidx,
                laststrongidx)

def is_motif_feasible(
    primer,
    exmotifs):
    '''
    Determine if primer devoid of exmotifs.
    Internal use only.

    :: primer
       type - string
       desc - a partially explored primer
              sequence path
    :: exmotifs
       type - cx.deque
       desc - deque of all motifs
              to be excluded
    '''

    # Quick Lookup
    for mlen in exmotifs:

        # Path is Partial
        if len(primer) < mlen:
            # No Additional Checks Needed
            break

        # Check Trace
        trace = primer[-mlen:]
        if trace in exmotifs[mlen]:
            # Found Conflict
            return False, trace

    # No Conflict
    return True, None

def is_edge_feasible(
    primer,
    primerlen,
    prefixforbidden,
    suffixforbidden):
    '''
    Determine if primer prefix and suffix
    is forbidden. Internal use only.

    :: primer
       type - string
       desc - a paritally explored primer
              sequence path
    :: primerlen
       type - integer
       desc - full primer sequence length
    :: prefixforbidden
       type - dict
       desc - dictionary of forbidden primer
              prefix sequences
    :: suffixforbidden
       type - dict
       desc - dictionary of forbidden primer
              suffix sequences
    '''

    # Primer Prefix Selected Forbidden?
    if not prefixforbidden is None:
        if len(primer) in prefixforbidden:
            if primer in prefixforbidden[len(primer)]:
                return False, len(primer) - 1

    # Primer Prefix Selected Forbidden?
    if not suffixforbidden is None:
        if len(primer) == primerlen:
            for mlen in suffixforbidden:
                primersuffix = primer[-mlen:]
                if primersuffix in suffixforbidden[mlen]:
                    return False, len(primer) - 1

    # Everything is OK
    return True, None

def is_structure_feasible(
    struct1,
    struct2,
    energy):
    '''
    Determine if structure is feasible.
    Internal use only.

    :: struct1
       type - string
       desc - MFE dimer structure of primer
    :: struct2
       type - string
       desc - MFE dimer structure of paired
              primer
    :: energy
       type - float
       desc - minimum free energy of dimer
              complex
    '''

    # Book-keeping
    minMFE    = -7   # kcal/mol
    free3ends = +4   # Free Bases
    objective = True # No Conflict

    # MFE Very Favorable
    if energy <= minMFE:
        objective = False # Conflict in Energy

    # Locate First Paired Base
    traceloc = struct1.find('(')

    # No Structure Found: No Conflict in Structure
    if traceloc == -1:
        objective = objective and True
        traceloc  = 0 # Abortive

    # Structures Exist: Check 3' Free Bases
    else:
        # Conflict: Paired Primer 3' Paired
        if ')' in struct2[-free3ends:]:
            objective = objective and False

        # Conflict: Primer 3' Paired
        elif '(' in struct1[-free3ends:]:
            objective = objective and False

        # No Conflict: 3' Ends Clear
        else:
            objective = objective and True

    # Return Results
    return objective, traceloc

def is_dimer_feasible(
    primer,
    primertype,
    primerlen,
    primerspan,
    pairedprimer,
    pairedspan,
    dimertype):
    '''
    Determine if primer is devoid of dimers.
    Internal use only.

    :: primer
       type - string
       desc - a paritally explored primer
              sequence path
    :: primerlen
       type - integer
       desc - full primer sequence length
    :: primerlen
       type - integer
       desc - full primer sequence length
    :: primerspan
       type - integer / None
       desc - primer length span to evaluate
              for padding
    :: pairedprimer
       type - string / None
       desc - paired primer sequence
    :: pairedspan
       type - integer / None
       desc - paired primer length span to
              evaluate for padding
    :: dimertype
       type - integer
       desc - dimer type computation identifier
              0 = compute homodimer
              1 = compute heterdimer
    '''

    # Too Short a Primer Candidate
    if len(primer) < primerlen:
        return True, None

    # Unpaired Heterodimer: No Computation
    if (dimertype == 1) and \
       (pairedprimer is None):
        return True, None

    # Update Primer Orientation
    revert = False
    if primertype == 1:
        primer = ut.get_revcomp(
            seq=primer)
        revert = True

    # Primer Span Extraction for Padding
    if not primerspan is None:
        primer = primer[-primerspan:]

    # Homodimer Adjustment
    if dimertype == 0:
        pairedprimer = primer

    # Paired Span Extraction for Padding
    if not pairedspan is None:
        pairedprimer = pairedprimer[-pairedspan:]

    # Compute the Dimer Minimum Free Energy Structures
    x, y, z = folder.evaluate_mfe_dimer(
        seq1=primer,
        seq2=pairedprimer)

    # Compute Result
    objective, traceloc = is_structure_feasible(
        struct1=x,
        struct2=y,
        energy=z)

    # Update Traceback Orientation
    if revert:
        traceloc = primerlen - 1 - traceloc

    # Return Result
    return objective, traceloc

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
