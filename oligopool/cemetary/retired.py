# Code Snippets that are good ideas, but unused

def get_runvec(seqvec):
    '''
    Return GC run vector for melting temp calc.
    Internal use only.

    :: seq
       type - string
       desc - sequence to decompose
    '''

    # Compute run vector
    rv = np.zeros(seqvec.shape[0]+1, dtype=np.float64)
    i  = 0
    gc = set([float(ord('G')), float(ord('C'))])
    while i < seqvec.shape[0]:
        rv[i+1] = rv[i] + (1. if seqvec[i] in gc else 0.)
        i += 1

    # Return result
    return rv

def get_runmat(seqmat):
    '''
    Return the GC run matrix for given seqmat.
    Internal use only.

    :: seqmat
       type - np.array
       desc - numeric sequence matrix
    '''

    j,k = seqmat.shape
    rm  = np.zeros((j, k+1), dtype=np.float64)
    i   = 0
    while i < j:
        rm[i, :] = get_runvec(seqmat[i, :])
        i += 1
    return rm

def get_approx_Tm(rv, i, j, mvc=50):
    '''
    Return approximate melting temperature
    based on the GC run vector. Internal
    use only.

    :: rv
       type - np.array
       desc - GC run vector
    :: i
       type - integer
       desc - sequence starting index (0-based)
    :: j
       type - integer / None
       desc - sequence ending index (0-based)
    :: mvc
       type - float
       desc - monovalent conc in mM
              (default=50 mM)
    '''

    # Calculate terms
    na = mvc / 1000.
    t1 = 81.5
    t2 = 16.6 * np.log10(na)
    gc = rv[j] - rv[i]
    rl = j - i
    t3 = 0.41 * 100. * (gc / rl)
    t4 = 600. / rl

    # Calculate correction
    cr = 0.
    if (j-i) < 60:
        if rv[j]   - rv[j-1] > 0 or \
           rv[i+1] - rv[i]   > 0:
            cr += 1.2
        if rv[j-1] - rv[j-2] > 0 or \
           rv[i+2] - rv[i+1] > 0:
            cr += 0.8
        if cr == 0.:
            cr = -1.0

    # Return result
    return t1 + t2 + t3 - t4 + cr

def is_all_assignable(
    seq,
    exmotifs,
    cifn,
    cntxtype,
    cntxnum,
    cntxlen):

    aidx  = 0

    while aidx < cntxnum:

        if cntxtype == 0:
            lcntx = cifn(aidx)[-cntxlen:]
            incntxseq = lcntx + seq
        else:
            rcntx = cifn(aidx)[:+cntxlen]
            incntxseq = seq + rcntx

        # Determine context feasibility
        mcond, motif = get_motif_conflict(
            seq=incntxseq,
            seqlen=len(incntxseq),
            exmotifs=exmotifs,
            partial=True)

        if not mcond:
            mstart = incntxseq.find(motif)

            # Context on Left
            if cntxtype == 0:
                # print((lcntx, seq, 'left'))
                tloc = len(motif) - len(lcntx) + mstart - 1
            # Context on Right
            else:
                # print((seq, rcntx, 'right'))
                tloc = mstart

            return False, tloc

        aidx += 1

    return True, None

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
    leftpartition = ut.get_exmotif_partition(
        exmotifs=exmotifs)

    # Show Update
    liner.send(' Right Paritioning Excluded Motifs ...')

    # Right Partition exmotifs
    rightpartition = ut.get_inverted_exmotif_partition(
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