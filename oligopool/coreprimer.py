import time  as tt

import collections as cx

import numpy       as np

import nrpcalc     as nr

import background  as bk
import vectordb    as db
import utils       as ut


# NRPCalc Fold Object
folder = nr.base.utils.Fold(
    temp=37.0,
    dangles=2,
    part_type='DNA')

# Parser and Setup Functions

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
    warn,
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
    :: warn
       type - dict
       desc - warning dictionary entry
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    t0 = tt.time()
    exmotifindex = None  # No Conflicts
    suboptimal   = False # Suboptimality Tracker
    homology     = 6     # Initially, for Maker

    # Compute Fixed Base Index
    fixedbaseindex = set()
    for idx in range(len(primerseq)):
        if len(ut.ddna_space[primerseq[idx]]) == 1:
            fixedbaseindex.add(idx)

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
            ' Design Space: 1 Possible Primer(s) [INFEASIBLE] (Degenerate Nucleotides Absent)\n')
    else:
        dspace_ok = True
        liner.send(
            ' Design Space: {:{},{}} Possible Primer(s)\n'.format(
                dspace,
                plen,
                sntn))

    # Exmotifs Analysis
    liner.send(' Computing Motif Conflicts ...')

    motif_ok, excludedmotifs = ut.get_exmotif_conflict(
        seq=primerseq,
        seqlen=len(primerseq),
        exmotifs=exmotifs,
        partial=False,
        checkall=True)

    # Show Update
    if not motif_ok:

        # Update Warning Entry
        warn['vars'] = {'exmotifembedded': set()}

        # Compute Embedded Motif Indices
        # to be Ignored Downstream
        exmotlocdict = ut.get_exmotif_conflict_index(
            seq=primerseq,
            conflicts=excludedmotifs)
        exmotifindex = set()
        for exmotif in exmotlocdict:
            for loc in exmotlocdict[exmotif]:
                exmotifindex.add(loc+len(exmotif))

        # Make exmotifs immutable
        excludedmotifs = tuple(excludedmotifs)

        # Show Updates
        liner.send(
            ' Found {:,} Excluded Motif(s)\n'.format(
                len(excludedmotifs)))

        # Record Warnings
        warn['warncount'] += len(excludedmotifs)
        warn['vars']['exmotifembedded'].update(excludedmotifs)
        suboptimal = True

        # Show Excluded Motifs
        plen = max(map(len, excludedmotifs)) + 2
        for motif in excludedmotifs:
            motif = '\'{}\''.format(motif)
            liner.send(
                '   - Excluded Motif {:>{}} Present [WARNING] (Excluded Motif Embedded)\n'.format(
                    motif,
                    plen))
    else:
        liner.send(
            ' Found 0 Excluded Motif(s)\n')

    # Region Analysis
    liner.send(' Computing Constant Regions ...')

    # Compute Constant Regions
    regions = ut.get_constant_regions(
        seqconstr=primerseq)

    # Finalize homology
    if regions:
        homology = max(homology,
                       max(map(len, regions)) + 1)

    # Internal Repeat Analysis
    liner.send(' Computing Internal Repeat Conflicts ...')

    # Compute Internal Repeats
    (repeat_ok,
    internalrepeats) = get_internal_repeat_conflicts(
        regions=regions)

    # Show Update
    if not repeat_ok:

        # Show Updates
        liner.send(
            ' Found {:,} Internal Repeat(s)\n'.format(
                len(internalrepeats)))

        # Show Internal Repeats
        plen = max(map(len, internalrepeats)) + 2
        for repeat in internalrepeats:
            repeat = '\'{}\''.format(repeat)
            liner.send(
                '   - Internal Repeat {:>{}} Present [INFEASIBLE] (Strong Hairpin Motif)\n'.format(
                    repeat,
                    plen))
    else:
        liner.send(
            ' Found 0 Internal Repeat(s)\n')

    # Palindrome Analysis
    liner.send(' Computing Palindrome Conflicts ...')

    # Update Warning Entry
    warn['vars'] = {'palindromeembedded': set()}

    # Compute Panlindrome
    palindrome_ok, palindromeembedded = get_palindrome_conflicts(
        regions=regions)

    # Show Update
    if not palindrome_ok:

        # Show Updates
        liner.send(
            ' Found {:,} Palindrome(s)\n'.format(
                len(palindromeembedded)))

        # Record Warnings
        warn['warncount'] += len(palindromeembedded)
        warn['vars']['palindromeembedded'].update(palindromeembedded)
        suboptimal = True

        # Show Palindromes
        plen = max(map(len, palindromeembedded)) + 2
        for palindrome in palindromeembedded:
            palindrome = '\'{}\''.format(palindrome)
            liner.send(
                '   - Palindrome {:>{}} Present [WARNING] (Strong Homodimer Motif)\n'.format(
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

        # Show Updates
        liner.send(
            ' Found {:,} Paired Repeat(s)\n'.format(
                len(pairedrepeats)))

        # Show Paired Repeats
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
        repeat_ok])

    if not parsestatus:
        liner.send(
            ' Verdict: Primer Design Infeasible due to Sequence Constraint\n')
    elif suboptimal:
        liner.send(
            ' Verdict: Primer Design Potentially Suboptimal\n')
    else:
        liner.send(
            ' Verdict: Primer Design Possibly Feasible\n')

    # Return Results
    return (parsestatus,
        homology,
        fixedbaseindex,
        exmotifindex,
        dspace,
        internalrepeats,
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

def get_parsed_primer_tmelt_constraint(
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

def get_parsed_edgeeffects(
    primerseq,
    leftcontext,
    rightcontext,
    leftpartition,
    rightpartition,
    exmotifs,
    warn,
    liner):
    '''
    Cluster left and right context sequences
    and record forbidden prefix and suffixes.
    Internal use only.

    :: primerseq
       type - string
       desc - motif sequence constraint
    :: leftcontext
       type - tuple
       desc - tuple of all left context
              sequences
    :: rightcontext
       type - tuple
       desc - tuple of all right context
              sequences
    :: exmotifs
       type - cx.deque
       desc - deque of all motifs
              to be excluded
    :: warn
       type - dict
       desc - warning dictionary entry
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    return ut.get_parsed_edgeeffects(
        sequence=primerseq,
        element='Primer',
        leftcontext=leftcontext,
        rightcontext=rightcontext,
        leftpartition=leftpartition,
        rightpartition=rightpartition,
        exmotifs=exmotifs,
        merge=True,
        warn=warn,
        liner=liner)

def get_parsed_splitcol(
    indf,
    splitcol,
    oligolength,
    liner):
    '''
    TBW
    '''

    # Book-keeping
    t0 = tt.time()
    parsestatus   = False
    maxallowedlen = oligolength - 30 # Ensure minimum padding feasible

    # Show Update
    liner.send(' Analyzing Split Fragments ...')

    # Compute Properties
    fragments = set(indf[splitcol])
    fraglens  = list(map(len, fragments))
    minffragmentlen    = min(fraglens)
    maxfragmentlen    = max(fraglens)

    # Show Update
    plen = ut.get_printlen(
        value=max(len(fragments), maxfragmentlen, minffragmentlen))

    liner.send('   Split Fragments: {:{},d} Unique Sequences(s)\n'.format(
        len(fragments),
        plen))

    liner.send(' Minimum Length   : {:{},d} Base Pair(s)\n'.format(
        minffragmentlen,
        plen))

    # Compute Feasibility
    if maxfragmentlen <= maxallowedlen:
        parsestatus = True
        parsemsg    = ''
    else:
        parsemsg = ' [INFEASIBLE] (Fragment(s) Longer than Oligo Length - 30 = {} Base Pair(s))'.format(
            maxallowedlen)

    # Show Updates
    liner.send(' Maximum Length   : {:{},d} Base Pair(s){}\n'.format(
        maxfragmentlen,
        plen,
        parsemsg))

    liner.send(' Time Elapsed: {:.2f} sec\n'.format(
        tt.time()-t0))

    # Show Verdict
    if not parsestatus:
        liner.send(
            ' Verdict: Pad Design Infeasible due to Split Column Constraints\n')
    else:
        liner.send(
            ' Verdict: Pad Design Possibly Feasible\n')

    # Return Results
    return (parsestatus,
        minffragmentlen,
        maxfragmentlen,
        maxallowedlen)

def get_parsed_typeIIS_constraint(
    typeIIS,
    typeIISname,
    minfragmentlen,
    maxfragmentlen,
    oligolength,
    liner):
    '''
    TBW
    '''

    # Book-keeping
    t0 = tt.time()
    parsestatus = True

    # How much padding space available?
    minpadlen = oligolength - maxfragmentlen
    maxpadlen = oligolength - minfragmentlen

    plen = ut.get_printlen(
        value=max(
            oligolength,
            maxfragmentlen,
            minfragmentlen))

    liner.send(
        ' Minimum Padding: {:{},d} bp due to {:{},d} bp Fragment(s)\n'.format(
            minpadlen,
            plen,
            maxfragmentlen,
            plen))

    liner.send(
        ' Maximum Padding: {:{},d} bp due to {:{},d} bp Fragment(s)\n'.format(
            maxpadlen,
            plen,
            minfragmentlen,
            plen))

    # Compute TypeIIS Padding Lengths
    typeIISlen         = (6 + len(typeIIS)) * 2          # 2x 6 Ns + TypeIIS Motif + Cut
    typeIISfree1       = max(0, 30 - typeIISlen)         # 2x 5' Ns to pad upto 15 bp TypeIIS Core
    typeIISrequiredlen = typeIISfree1 + typeIISlen       # Total Length of TypeIIS Core >= 30 bp
    typeIISfree2       = minpadlen - typeIISrequiredlen  # 2x 5' Ns to pad upto minpadlen

    # Compute feasibility
    if typeIISfree2 < 0:
        parsestatus = False
        parsemsg = ' [INFEASIBLE] (Using \'{}\' Requires {} bp Longer Minimum Padding)'.format(
            typeIISname,
            np.abs(typeIISfree2))
    else:
        parsemsg = ''

    liner.send(
        ' TypeIIS System : Requires {:,} bp Minimum Padding{}\n'.format(
            typeIISrequiredlen,
            parsemsg))

    print(typeIISlen)
    print(typeIISfree1)
    print(typeIISrequiredlen)
    print(typeIISfree2)

    # Finalize Padding Lengths
    typeIIScore = 'N'*(typeIISfree1 // 2) + 'N'*6 + typeIIS

    print(typeIIScore, len(typeIIScore))

    for val in ut.typeIIS_dict.values():
        print(val[0], len(val[1]) + val[-1])

    # # How much space available
    # forwardpadlen    = maxpadlen // 2
    # forwardprimerlen = minpadlen // 2
    # forwardfillerlen = forwardlen - forwardpadlen

    # liner.send(
    #     ' Forward Pad Constructed: {:{},d} to Base Pair(s)\n'.format(
    #         forwardconsspan,
    #         plen))

    # liner.send(
    #     ' Padding Required: {:,} to {:,} Base Pair(s)\n'.format(
    #         minpadlen,
    #         maxpadlen))


    # reverseconsspan = maxpadlen - forwardconsspan
    # reverseevalspan = minpadlen - forwardevalspan

    # Figure out Max and Min Pad Spans (both)

    # Check if after Motif, >= 6 bp left on 5' of Pad (both)

    # Return parsestatus, available forward pad length, required fwd pad len, ... rev pad

def get_padding_lengths(
    typeIIS,
    fragments,
    finallength,
    liner):
    '''
    TBW
    '''

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

def get_parsed_pad_tmelt_constraint(
    typeIIS,
    mintmelt,
    maxtmelt,
    liner):
    pass

    # Same as Primer, but both pads evaluated

    '''
    for pad in [fwd, rev]:
        for padlength in range(len(5' free size)):
            primerseq = build contraint (based on pad type)

            get extreme min and max

            if min > required max:
                break because infeasible further
                previous min,max is the Tm window

            elif max < required min:
                try increasing padding length

            elif min,max found is compatible:
                break because solved

        adjust min,max --> this is pad Tm range
    '''

# Engine Objective and Helper Functions

def show_update(
    element,
    primer,
    optstatus,
    optstate,
    inittime,
    terminal,
    liner):
    '''
    Display the current progress in primer
    generation. Internal use only.

    :: element
       type - string
       desc - primer element name, e.g.
              'Primer' or 'Pad'
    :: primer
       type - string
       desc - a partially explored primer
              sequence path
    :: optstatus
       type - integer
       desc - primer feasibility status
    :: optstate
       type - integer
       desc - feasibility failure state marker
    :: inittime
       type - tt.time
       desc - initial time stamp
    :: terminal
       type - boolean
       desc - if True will terminate update to newline
              otherwise, rewrite previous update
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    if len(primer) >= 30:
        design = primer[:14] + '..' + primer[-14:]
    else:
        design = primer

    liner.send(' Candidate: {} {} is {}{}'.format(
        element,
        design,
        ['Rejected', 'Provisionally Accepted', 'Accepted'][optstatus],
        ['',
        ' due to Background Repeat',
        ' due to Oligopool Repeat',
        ' due to Paired Repeat',
        ' due to Excluded Motif',
        ' due to Edge Effect',
        ' due to Melting Temperature',
        ' due to Homodimer',
        ' due to Heterodimer'][optstate]))

    if terminal:
        liner.send('\* Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - inittime))

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
       type - db.vectorDB / None
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
       type - set / None
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
       type - set / None
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
       desc - primer melting temperature
              lower bound
    :: maxtmelt
       type - float
       desc - primer melting temperature
              upper bound
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
    failtype,
    fixedbaseindex):
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
    :: fixedbaseindex
       type - set
       desc - set of all fixed base indices
    '''

    # Define Default Traceback
    defaultidx = len(primer) - 1

    # Lower Bound Mitigation
    if failtype == 0:

        # Locate All Weak Base Indices
        weakidx = set(idx for idx in range(len(primer)) \
            if primer[idx] in 'AT') - fixedbaseindex

        # Locate Latest Weak Base
        lastweakidx = max(weakidx) if weakidx else -1

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

        # Locate All Strong Base Indices
        strongidx = set(idx for idx in range(len(primer)) \
            if primer[idx] in 'GC') - fixedbaseindex

        # Locate Latest Strong Base
        laststrongidx = max(strongidx) if strongidx else -1

        # No Strong Bases
        if laststrongidx == -1:
            # Abort
            return 0
        else:
            # Traceback
            return min(
                defaultidx,
                laststrongidx)

def is_exmotif_feasible(
    primer,
    exmotifs,
    exmotifindex):
    '''
    Determine if primer devoid of exmotifs.
    Internal use only.

    :: primer
       type - string
       desc - a partially explored primer
              sequence path
    :: exmotifs
       type - set / None
       desc - set of all excluded motifs
    :: exmotifindex
       type - set
       desc - set of constraint embedded
              exmotif ending indices
    '''

    return ut.is_local_exmotif_feasible(
        seq=primer,
        exmotifs=exmotifs,
        exmotifindex=exmotifindex)

def is_edge_feasible(
    primer,
    primerlen,
    edgeeffectlength,
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
    :: edgeeffectlength
       type - integer
       desc - context length for edge effects
    :: prefixforbidden
       type - dict / None
       desc - dictionary of forbidden primer
              prefix sequences
    :: suffixforbidden
       type - dict / None
       desc - dictionary of forbidden primer
              suffix sequences
    '''

    return ut.is_local_edge_feasible(
        seq=primer,
        seqlen=primerlen,
        lcseq=None,
        rcseq=None,
        edgeeffectlength=edgeeffectlength,
        prefixforbidden=prefixforbidden,
        suffixforbidden=suffixforbidden)

def is_structure_feasible(
    struct1,
    struct2,
    energy,
    fixedbaseindex):
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
    :: fixedbaseindex
       type - set
       desc - set of all fixed base indices
    '''

    # Compute Mutable Paired Base Index
    pairedindex = set()
    for idx in range(len(struct1)):
        if struct1[idx] == '(':
            pairedindex.add(idx)
    pairedindex -= fixedbaseindex

    # Book-keeping
    minHFE    = -7   # kcal/mol
    free3ends = +4   # Free Bases
    objective = True # No Conflict

    # HFE Very Favorable
    if (energy <= minHFE) and \
       len(pairedindex) > 0:
        objective = False # Conflict in Energy, but Mutable

    # Locate First Paired Base
    if pairedindex:
        traceloc = min(pairedindex)
    else:
        traceloc = -1

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
    fixedbaseindex,
    dimertype):
    '''
    Determine if primer is devoid of dimers.
    Internal use only.

    :: primer
       type - string
       desc - a paritally explored primer
              sequence path
    :: primertype
       type - integer
       desc - primer design type identifier
              0 = forward primer design
              1 = reverse primer design
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
    :: fixedbaseindex
       type - set
       desc - set of all fixed base indices
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

    # Infer Primer Orientation
    revert = False
    if primertype == 1:
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
        energy=z,
        fixedbaseindex=fixedbaseindex)

    # Update Traceback Orientation
    if revert:
        traceloc = primerlen - 1 - traceloc

    # Return Result
    return objective, traceloc

def primer_objectives(
    primer,
    primerlen,
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
    prefixforbidden,
    suffixforbidden,
    background,
    inittime,
    stats,
    liner):
    '''
    Determine if a primer satisfies all
    local objectives. Internal use only.

    :: primer
       type - string
       desc - a partially explored primer
              sequence path
    :: primerlen
       type - integer
       desc - length of the full primer
              sequence being designed
    :: primertype
       type - integer
       desc - primer design type identifier
              0 = forward primer design
              1 = reverse primer design
    :: fixedbaseindex
       type - set
       desc - set of all fixed base indices
    :: mintmelt
       type - float
       desc - primer melting temperature lower
              bound
    :: maxtmelt
       type - float
       desc - primer melting temperature upper
              bound
    :: oligorepeats
       type - set
       desc - set storing oligopool repeats
    :: pairedprimer
       type - string / None
       desc - paired primer sequence
    :: pairedrepeats
       type - set / None
       desc - set storing paired primer repeats
    :: exmotifs
       type - cx.deque / None
       desc - deque of all excluded motifs
    :: exmotifindex
       type - set / None
       desc - set of constraint embedded
              exmotif ending indices
    :: edgeeffectlength
       type - integer
       desc - context length for edge effects
    :: prefixforbidden
       type - dict / None
       desc - dictionary of forbidden primer
              prefix sequences
    :: suffixforbidden
       type - dict / None
       desc - dictionary of forbidden primer
              suffix sequences
    :: background
       type - db.vectorDB / None
       desc - vectorDB instance storing
              background repeats
    :: inittime
       type - tt.time
       desc - initial time stamp
    :: stats
       type - dict
       desc - primer design stats storage
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Objective 1: Background Non-Repetitiveness
    obj1, traceloc = is_background_feasible(
        primer=primer,
        background=background)

    # Objective 1 Failed
    if not obj1:

        # Show Update
        show_update(
            element='Primer',
            primer=primer,
            optstatus=0,
            optstate=1,
            inittime=inittime,
            terminal=False,
            liner=liner)

        # Update Stats
        stats['vars']['repeatfail'] += 1

        # Return Traceback
        return False, traceloc

    # Objective 2: Oligopool Non-Repetitiveness
    obj2, traceloc = is_oligopool_feasible(
        primer=primer,
        maxreplen=maxreplen,
        oligorepeats=oligorepeats)

    # Objective 2 Failed
    if not obj2:

        # Show Update
        show_update(
            element='Primer',
            primer=primer,
            optstatus=0,
            optstate=2,
            inittime=inittime,
            terminal=False,
            liner=liner)

        # Update Stats
        stats['vars']['repeatfail'] += 1

        # Return Traceback
        return False, traceloc

    # Objective 3: Paired Primer Non-Repetitiveness
    obj3, traceloc = is_paired_feasible(
        primer=primer,
        pairedrepeats=pairedrepeats)

    # Objective 3 Failed
    if not obj3:

        # Show Update
        show_update(
            element='Primer',
            primer=primer,
            optstatus=0,
            optstate=3,
            inittime=inittime,
            terminal=False,
            liner=liner)

        # Update Stats
        stats['vars']['repeatfail'] += 1

        # Return Traceback
        return False, traceloc

    # Objective 4: Motif Embedding
    obj4, exmotif = is_exmotif_feasible(
        primer=primer,
        exmotifs=exmotifs,
        exmotifindex=exmotifindex)

    # Objective 4 Failed
    if not obj4:

        # Show Update
        show_update(
            element='Primer',
            primer=primer,
            optstatus=0,
            optstate=4,
            inittime=inittime,
            terminal=False,
            liner=liner)

        # Update Stats
        stats['vars']['exmotiffail'] += 1
        stats['vars']['exmotifcounter'][exmotif] += 1

        # Return Traceback
        return False, max(
            0,
            len(primer)-1)

    # Objective 5: Edge Feasibility (Edge-Effects)
    obj5, dxmotifs, traceloc = is_edge_feasible(
        primer=primer,
        primerlen=primerlen,
        edgeeffectlength=edgeeffectlength,
        prefixforbidden=prefixforbidden,
        suffixforbidden=suffixforbidden)

    # Objective 5 Failed
    if not obj5:

        # Show Update
        show_update(
            element='Primer',
            primer=primer,
            optstatus=0,
            optstate=5,
            inittime=inittime,
            terminal=False,
            liner=liner)

        # Update Stats
        stats['vars']['edgefail'] += len(dxmotifs)
        stats['vars']['exmotifcounter'].update(dxmotifs)

        # Return Traceback
        return False, traceloc

    # Is Primer Fully Explored?
    if len(primer) == primerlen:

        # Objective 6: Melting Temeperature Bounded
        obj6, failtype = is_tmelt_feasible(
            primer=primer,
            primerlen=primerlen,
            mintmelt=mintmelt,
            maxtmelt=maxtmelt)

        # Objective 6 Failed
        if not obj6:

            # Show Update
            show_update(
                element='Primer',
                primer=primer,
                optstatus=0,
                optstate=6,
                inittime=inittime,
                terminal=False,
                liner=liner)

            # Compute Traceback
            traceloc = get_tmelt_traceback(
                primer=primer,
                failtype=failtype,
                fixedbaseindex=fixedbaseindex)

            # Update Stats
            stats['vars']['Tmfail'] += 1

            # Return Traceback
            return False, traceloc

        # Update Primer Orientation
        cprimer = ut.get_revcomp(
            seq=primer) if primertype == 1 else primer

        # Objective 7: Homodimer Feasibility
        obj7, traceloc = is_dimer_feasible(
            primer=cprimer,
            primertype=primertype,
            primerlen=primerlen,
            primerspan=None,
            pairedprimer=None,
            pairedspan=None,
            fixedbaseindex=fixedbaseindex,
            dimertype=0)

        # Objective 7 Failed
        if not obj7:

            # Show Update
            show_update(
                element='Primer',
                primer=primer,
                optstatus=0,
                optstate=7,
                inittime=inittime,
                terminal=False,
                liner=liner)

            # Update Stats
            stats['vars']['homodimerfail'] += 1

            # Return Traceback
            return False, traceloc

        # Objective 8: Heterodimer Feasibility
        obj8, traceloc = is_dimer_feasible(
            primer=cprimer,
            primertype=primertype,
            primerlen=primerlen,
            primerspan=None,
            pairedprimer=pairedprimer,
            pairedspan=None,
            fixedbaseindex=fixedbaseindex,
            dimertype=1)

        # Objective 8 Failed
        if not obj8:

            # Show Update
            show_update(
                element='Primer',
                primer=primer,
                optstatus=0,
                optstate=8,
                inittime=inittime,
                terminal=False,
                liner=liner)

            # Update Stats
            stats['vars']['heterodimerfail'] += 1

            # Return Traceback
            return False, traceloc

    # Show Update
    show_update(
        primer=primer,
        element='Primer',
        optstatus=1,
        optstate=0,
        inittime=inittime,
        terminal=False,
        liner=liner)

    # All Objectives OK!
    return True

