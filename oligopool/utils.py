import sys
import os
import shutil  as su
import time    as tt

import uuid    as uu
import atexit  as ae

import collections as cx

import numpy   as np
import numba   as nb
import primer3 as p3


# Global Lookups
complement_table = str.maketrans(
    'ACGTURYSWKMBVDHN',
    'TGCATKYWSRMBDHVN')

dna_alpha = set(
    'ATGC')

ddna_space = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'S': {'G', 'C'},
    'W': {'A', 'T'},
    'K': {'G', 'T'},
    'M': {'A', 'C'},
    'B': {'C', 'G', 'T'},
    'V': {'A', 'C', 'G'},
    'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'},
    'N': {'A', 'T', 'G', 'C'}}

ddna_alpha = set(ddna_space.keys())

typeIIS_dict = {
    'acui'    : 'CTGAAG'  + 'N' * 16,
    'alwi'    : 'GGATC'   + 'N' *  5,
    'bbsi'    : 'GAAGAC'  + 'N' *  6,
    'bcci'    : 'CCATC'   + 'N' *  5,
    'bceai'   : 'ACGGC'   + 'N' * 14,
    'bcivi'   : 'GTATCC'  + 'N' *  6,
    'bcodi'   : 'GTCTC'   + 'N' *  5,
    'bmri'    : 'ACTGGG'  + 'N' *  5,
    'bpuei'   : 'CTTGAG'  + 'N' * 16,
    'bsai'    : 'GGTCTC'  + 'N' *  5,
    'bseri'   : 'GAGGAG'  + 'N' * 10,
    'bsmai'   : 'GTCTC'   + 'N' *  5,
    'bsmbi'   : 'CGTCTC'  + 'N' *  5,
    'bsmfi'   : 'GGGAC'   + 'N' * 14,
    'bsmi'    : 'GAATGC'  + 'N' *  1,
    'bspcni'  : 'CTCAG'   + 'N' *  9,
    'bspqi'   : 'GCTCTTC' + 'N' *  4,
    'bsrdi'   : 'GCAATG'  + 'N' *  2,
    'bsri'    : 'ACTGG'   + 'N' *  1,
    'btgzi'   : 'GCGATG'  + 'N' * 14,
    'btsci'   : 'GGATG'   + 'N' *  2,
    'btsi'    : 'GCAGTG'  + 'N' *  2,
    'btsimuti': 'CAGTG'   + 'N' *  2,
    'eari'    : 'CTCTTC'  + 'N' *  4,
    'ecii'    : 'GGCGGA'  + 'N' * 11,
    'esp3i'   : 'CGTCTC'  + 'N' *  5,
    'faui'    : 'CCCGC'   + 'N' *  6,
    'hgai'    : 'GACGC'   + 'N' * 10,
    'hphi'    : 'GGTGA'   + 'N' *  8,
    'hpyav'   : 'CCTTC'   + 'N' *  6,
    'mlyi'    : 'GAGTC'   + 'N' *  5,
    'mnli'    : 'CCTC'    + 'N' *  7,
    'sapi'    : 'GCTCTTC' + 'N' *  4,
    'sfani'   : 'GCATC'   + 'N' *  9}

# Decorators

def coroutine(func):
    '''
    Decorator used to prime a coroutine.
    Internal use only.

    :: func
       type - coroutine
       desc - a coroutine we want online
    '''

    def online(*args, **kwargs):
        _coroutine = func(*args, **kwargs)
        next(_coroutine)
        return _coroutine

    return online

# Printing / Logging

@coroutine
def liner_engine(online=True):
    '''
    Return a coroutine to print stuff
    to stdout and log run information.
    Internal use only.

    :: online
       type - boolean
       desc - if True will print received
              strings to stdout
              (default=True)
    '''

    # Book-keeping
    clrlen = 0

    # Online
    try:
        while True:

            # Receive String
            printstr = (yield)

            # Line Preservation
            if printstr.startswith('\*'):
                sys.stdout.write('\n')
                clrlen   = 0
                printstr = printstr.lstrip('\*')

            # Line String
            if online:
                sys.stdout.write('\r' + ' '*clrlen)
                sys.stdout.write('\r' + printstr)
                clrlen = len(printstr)
                sys.stdout.flush()

    # Closure
    except GeneratorExit:
        if online:
            sys.stdout.write('\n')

def get_printlen(value):
    '''
    Return printing length.
    Internal use only.

    :: value
       type - nu.Real
       desc - value to evaluate
    '''

    svalue = str(value)
    if len(svalue) <= 3:
        return len(svalue)
    else:
        return len(str(value)) + \
            int(safelog10(value) / 3)

def get_notelen(printlen):
    '''
    Return notation length.
    Internal use only.

    :: value
       type - integer
       desc - printlen to evaluate
    '''

    if printlen <= 15:
        return 'd', printlen
    else:
        return 'e', 0

# Numeric Functions

def get_trials(prob):
    '''
    Return the number of required trials
    given a probability of success.
    Internal use only.

    :: prob
       type - float
       desc - probability of a rate event
    '''
    growth_rate = -100.0 # Lower is Rarer
    # print(-10 / np.log(1. - prob))
    return 2 * int(np.ceil(
        growth_rate / np.log(1. - prob)))

def get_prob(success, trials):
    '''
    Return the probability of success given
    the number of successes in some finite
    trials. Internal use only.

    :: success
       type - integer
       desc - total number of successes
    :: trials
       type - integer
       desc - total number of attempts made
    '''
    return (success + 1.) / (trials + 2.)

def safediv(A, B):
    '''
    Return A / B, except when B is zero return
    zero (no more pesky Div by Zero Errors).
    Internal use only.

    :: A
       type - nu.Real
       desc - numerator
    :: B
       type - nu.Real
       desc - denominator
    '''

    return 0. if B == 0. else float(A) / B

def safelog10(A):
    '''
    Return the log10 of A, except when A is less
    than or equal to zero return zero (no more
    pesky Math Domain Errors).
    Internal use only.

    :: A
       type - nu.Real
       desc - quantity
    '''

    return np.log10(A) if A > 0. else 0.

def get_sample(value, lf, uf):
    '''
    Return a random integer between
    (lf*value, uf*value).
    Internal use only.

    :: value
       type - nu.Real
       desc - numeric value
    :: lf
       type - nu.Real
       desc - sample lower bound fraction
              of value
    :: uf
       type - nu.Real
       desc - sample upper bound fraction
              of value
    '''

    return np.random.randint(
        *map(np.round, (value * lf,
                        value * uf)))

# DataFrame Functions

def get_uniques(
    iterable,
    typer):
    '''
    Return the unique elements in iterable
    packed via typer. Internal use only.

    :: iterable
       type - iterable
       desc - iterable to uniquify
    :: typer
       type - function
       desc - factory function to wrap
              uniques extracted
    '''
    uniques = []
    seen    = set()
    for element in iterable:
        if not element in seen:
            uniques.append(element)
            seen.add(element)
    return typer(uniques)

def get_col_exist_idx(
    col,
    df):
    '''
    Determine if col exists as a column in df,
    and return its index. Internal use only.

    :: col
       type - string
       desc - column name to check for existence
    ::df
       type - pd.DataFrame
       desc - DataFrame to check existence in
    '''

    col = col.lower()
    collist = tuple(c.lower() for c in df.columns)
    if col == df.index.name.lower() or \
       col in collist:
        return True, collist.index(col)
    else:
        return False, None

def get_df_concat(df):
    '''
    Concatenate all columns in df.
    Internal use only.

    :: df
       type - pd.DataFrame
       desc - DataFrame to concatenate
              columns from
    '''

    return tuple(df.values.sum(axis=1))

def update_df(
    indf,
    lcname,
    rcname,
    out,
    outcol):
    '''
    Insert out values in outcol between
    lcname and rcname cols if provided.
    Internal use only.

    :: indf
       type - pd.DataFrame
       desc - DataFrame to update
    :: lcname
       type - string / None
       desc - left context column name
    :: rcname
       type - string
       desc - right context column name
    :: out
       type - list / string
       desc - value(s) to insert
    :: outcol
       type - string
       desc - column name for values
    '''

    # Book-keeping
    insidx    = len(indf.columns)
    insstatus = False

    # Do we have a left context?
    if not insstatus and \
       not lcname is None:

        # Determine existence and index
        (lcexists,
        lcidx) = get_col_exist_idx(
            col=lcname,
            df=indf)

        # Update Book-keeping
        if lcexists:
            insidx    = lcidx + 1
            insstatus = True

    # Do we have a right context?
    if not insstatus and \
       not rcname is None:

        # Determine existence and index
        (rcexists,
        rcidx) = get_col_exist_idx(
            col=rcname,
            df=indf)

        # Update Book-keeping
        if rcexists:
            insidx    = max(0, rcidx-1)
            insstatus = True

    # Update indf
    indf.insert(
        loc=insidx,
        column=outcol,
        value=out)

# Oligo Functions

@nb.njit
def get_store_hdist(
    store,
    idx,
    i=0,
    j=None,
    direction=0):
    '''
    Return the minimum pairwise hamming
    distance between store[idx, i:j] and
    store[*: i:j] depending on direction.
    Internal use only.

    :: store
       type - np.array
       desc - numeric sequence array
    :: idx
       type - integer
       desc - location index of sequence
              being compared
    :: i
       type - integer
       desc - starting index of comparison
    :: j
       type - integer
       desc - ending   index of comparison
    :: direction
       type - integer
       desc - direction flag;
              0 = upward   comparison
              1 = downward comparison
              2 = all-pair comparison
    '''

    # Default Result
    hdist = store.shape[1]
    if i >= 0 and not j is None:
        hdist = j-i

    # Something to compare against?
    if idx > 0:

        # Upward / All-Pair Comparison
        if direction == 0 or direction == 2:
            hdist = min(
                hdist,
                (store[:idx, i:j] != store[idx, i:j]).sum(1).min())

    # Something to compare against?
    if idx < store.shape[0] - 1:

        # Downward / All-Pair Comparison
        if direction == 1 or direction == 2:
            hdist = min(
                hdist,
                (store[idx+1:, i:j] != store[idx, i:j]).sum(1).min())

    # Return Result
    return hdist

def get_tmelt(
    seq,
    i=0,
    j=None,
    mvc=50.,
    dvc=0.,
    ntc=0.8,
    olc=50.):
    '''
    Return the melting temperature of seq.
    Internal use only.

    Note: primer3-py API docs

    :: seq
       type - string
       desc - DNA string in context
    :: i
       type - integer
       desc - starting index of comparison
              (default=)
    :: j
       type - integer
       desc - ending index of comparison
    :: mvc
       type - float
       desc - monovalent cation conc. (mM)
              (default=50.0 nM)
    :: dvc
       type - float
       desc - divalent cation conc. (mM)
              (default=0.0 nM)
    :: ntc
       type - float
       desc - dNTP conc. (mM)
              (default=0.8 mM)
    :: olc
       type - float
       desc - oligo conc. (nM)
              (default=50 nM)
    '''
    return p3.bindings.calcTm(
        seq=seq[i:j],
        mv_conc=mvc,
        dv_conc=dvc,
        dntp_conc=ntc,
        dna_conc=olc)

# Motif and Context Functions

def get_parsed_exmotifs(
    exmotifs,
    typer,
    element,
    liner):
    '''
    Check feasibility and return all
    exmotifs sorted by length.
    Internal use only.

    :: exmotifs
       type - list
       desc - list of all motifs
              to be excluded
    :: typer
       type - function
       desc - factory function for
              wrapping exmotifs
    :: element
       type - string
       desc - element being designed
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Time-keeping
    t0 = tt.time()

    # Sort Enque all motifs by length
    liner.send(' Sorting and Enqueing Motif(s) ...')

    exmotifs = sorted(
        exmotifs,
        key=len)
    dq = typer(exmotifs)

    liner.send(' Sorted and Enqued: {:,} Unique Motif(s)\n'.format(
        len(exmotifs)))

    # Check motif feasibility
    liner.send(' Computing Motif Length Distribution ...')

    # Compute length distribution
    cr = cx.Counter(len(m) for m in dq)

    liner.send(' Motif Length Distribution\n')

    klen = get_printlen(
        value=max(cr.keys()))
    vlen = get_printlen(
        value=max(cr.values()))

    # Check for infeasible lengths
    parsestatus = True
    problens    = []

    for mlen in sorted(cr.keys()):

        # Compute infeasibility
        parsemsg = ''
        if cr[mlen] == 4**(mlen):
            parsemsg    = ' [INFEASIBLE] (All {}-mers Excluded)'.format(mlen)
            parsestatus = False
            problens.append(mlen)

        # Show update
        liner.send('   - {:{},d} Motif(s) of Length {:{},d}{}\n'.format(
            cr[mlen],
            vlen,
            mlen,
            klen,
            parsemsg))

    # Finalize problens
    problens = None if len(problens) == 0 else problens

    # Show feasibility verdict
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    if not parsestatus:
        liner.send(
            ' Verdict: {} Design Infeasible due to Excluded Motif Constraints\n'.format(
                element))
    else:
        liner.send(
            ' Verdict: {} Design Possibly Feasible\n'.format(
                element))

    # Return Results
    return (parsestatus,
        dq,
        problens)

def get_motif_conflict(
    seq,
    seqlen,
    exmotifs,
    partial=False,
    checkall=False):
    '''
    Determine if the sequence does not contain or is
    contained in one of the exluded motifs (motif
    feasibility). Internal use only.

    :: seq
       type - string
       desc - sequence for conflict checking
    :: seqlen
       type - integer
       desc - length of sequence
    :: exmotifs
       type - iterable / None
       desc - list of motifs to exclude
              in designed barcodes,
              None otherwise
    :: partial
       type - boolean
       desc - if True will not check for
              conflicts for motifs longer
              than seqlen, otherwise run
              all checks
              (default=False)
    :: checkall
       type - boolean
       desc - if True will check for all
              possible conflicts due to
              exmotifs
              (default=False)
    '''

    # Book-keeping
    status = True # No Conflict
    pmotif = set() if checkall else None

    # Do we have anything to exclude?
    if exmotifs:

        # Loop through all motifs
        for motif in exmotifs:

            # Is sequence partial?
            if partial:
                if len(motif) > seqlen:
                    break # No need to check further

            # Embedding and Embedded Conflict
            if (len(motif) <= seqlen and motif in seq) or \
               (len(motif)  > seqlen and seq   in motif):

                # We got Conflict!
                status = False

                # What is the Conflict?
                if checkall:
                    pmotif.add(motif)
                else:
                    pmotif = motif
                    break

    # Return result
    if status:
        return (True, None)
    else:
        return (False, pmotif)

def get_edgeeffectlength(exmotifs):
    '''
    Get the context length to search in for exmotif
    edge-effects. Internal use only.

    :: exmotifs
       type - iterable / None
       desc - iterable of motifs to exlude,
              None otherwise
    '''

    if not exmotifs is None:
        return len(exmotifs[-1])
    return 0

def get_grouped_sequences(sequences):
    '''
    Group all selected strings by length.
    Internal use only.

    :: sequences
       type - iterable
       desc - set of all sequences / motifs
              to group by length
    '''

    groupdict = cx.defaultdict(set)
    for sequence in sorted(sequences, key=len):
        groupdict[len(sequence)].add(sequence)
    return groupdict

def get_extracted_edge(
    contextseq,
    position,
    edgeeffectlength):
    '''
    Extract context sequence edge based on
    context position. Internal use only.

    :: contextseq
       type - string
       desc - context sequence to evaluate
    :: position
       type - integer
       desc - context position identifier
              0 =  left context
              1 = right context
    :: edgeeffectlength
       type - integer
       desc - length of context sequence to
              extract for edge-effect eval
    '''
    # Do we have edges to consider?
    if not edgeeffectlength is None:

        # Extract Left Context
        if position == 0:
            return contextseq[-edgeeffectlength:]

        # Extract Right Context
        else:
            return contextseq[:+edgeeffectlength]

    # No edges to consider
    return None

def get_extracted_context(
    leftcontext,
    rightcontext,
    edgeeffectlength,
    reduce,
    liner):
    '''
    Return parsed context sequences.
    Internal use only.

    :: leftcontext
       type - pd.Series
       desc - left context sequences
    :: rightcontext
       type - pd.Series
       desc - right context sequences
    :: edgeeffectlength
       type - integer / None
       desc - length of context sequence to
              extract for edge-effect eval
    :: reduce
       type - boolean
       desc - if True will reduce multi-unique
              context sequences to uniques;
              otherwise return full context
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Start Timer
    t0 = tt.time()

    # Setup Contexts
    contexts         = ((leftcontext, 0),  (rightcontext, 1))
    uniquecount      = []
    extractedcontext = []

    # Parse Context
    for context,position in contexts:

        # No Context
        if context is None:
            uniquecount.append(0)
            extractedcontext.append(None)

        # Context Present
        else:
            # Extract Uniques
            uniques = list(context.unique())
            uniquecount.append(
                len(uniques))

            # Context is Constant
            if uniquecount[-1] == 1:

                # Extract Edge Sequence
                edgeseq = get_extracted_edge(
                    contextseq=uniques[-1],
                    position=position,
                    edgeeffectlength=edgeeffectlength)

                # Store Edge Sequence
                extractedcontext.append(
                    edgeseq)

            # Context is Variable
            else:

                # Reduce and Extract Context
                if reduce:
                    extractedcontext.append(
                        uniques)
                else:
                    extractedcontext.append(
                        context.to_list())

                # Compute Ignored Indices
                for idx,seq in enumerate(extractedcontext[-1]):

                    # Extract Edge Sequence
                    edgeseq = get_extracted_edge(
                        contextseq=seq,
                        position=position,
                        edgeeffectlength=edgeeffectlength)

                    # Show Update
                    liner.send(
                        ' Extracting {} Context: Processed {:,} Sequences'.format(
                            ['Left', 'Right'][position],
                            idx+1))

                    # Update Edge Sequence
                    extractedcontext[-1][idx] = edgeseq

    # Compute plen
    plen = get_printlen(
        value=max(uniquecount))

    # Show Updates
    liner.send(
        '  Left Context: {:{},d} Unique Sequence(s)\n'.format(
            uniquecount[0],
            plen))
    liner.send(
        ' Right Context: {:{},d} Unique Sequence(s)\n'.format(
            uniquecount[1],
            plen))
    liner.send('  Time Elapsed: {:.2f} sec\n'.format(tt.time()-t0))

    # Return Parsed Context
    return (extractedcontext[0],
        extractedcontext[1])

def get_parsed_oligo_repeats(
    df,
    maxreplen,
    element,
    liner):
    '''
    Check if oligopool repeats are
    feasible. Internal use only.

    :: df
       type - pd.DataFrame
       desc - DataFrame to extract
              oligopool repeats from
    :: maxreplen
       type - integer
       desc - maximum shared repeat
              length with oligopool
    :: element
       type - string
       desc - element being designed
    :: liner
       type - coroutine
       desc - dyanamic printing
    '''

    # Book-keeping
    oligorepeats = set()
    t0           = tt.time()
    kmerspace    = ((4**(maxreplen+1)) // 2)
    fillcount    = None
    leftcount    = None

    # Verbage Stuff
    plen = get_printlen(
        value=kmerspace)
    if plen > 15:
        sntn = 'e'
        plen = len('{:e}'.format(kmerspace))
    else:
        sntn = 'd'

    # Show Update
    liner.send('  k-mer Space: {:{},{}} Unique {:,}-mers\n'.format(
        kmerspace,
        plen,
        sntn,
        maxreplen+1))
    liner.send(' Extracting Repeats ...')

    # Extract Repeats
    for oligo in get_df_concat(df=df):
        oligorepeats.update(
            stream_canon_spectrum(
                seq=oligo,
                k=maxreplen+1))

    # Compute Feasibility
    fillcount   = min(kmerspace, len(oligorepeats))
    leftcount   = kmerspace - fillcount
    parsestatus = (leftcount * 1.) / kmerspace > .01
    if parsestatus:
        statusmsg = ''
    else:
        statusmsg = ' [INFEASIBLE]'

    # Show Final Updates
    liner.send('   Fill Count: {:{},{}} Unique {:,}-mers ({:6.2f} %)\n'.format(
        fillcount,
        plen,
        sntn,
        maxreplen+1,
        safediv(
            A=fillcount*100.,
            B=kmerspace)))
    liner.send('   Left Count: {:{},{}} Unique {:,}-mers ({:6.2f} %){}\n'.format(
        leftcount,
        plen,
        sntn,
        maxreplen+1,
        safediv(
            A=leftcount*100.,
            B=kmerspace),
        statusmsg))

    if not parsestatus:
        liner.send(
            ' Verdict: {} Design Infeasible due to Repeat Length Constraint\n'.format(
                element))
    else:
        liner.send(
            ' Verdict: {} Design Possibly Feasible\n'.format(
                element))

    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Return Results
    return (parsestatus,
        kmerspace,
        fillcount,
        leftcount,
        oligorepeats)

# Sequence Analysis

def is_DNA(seq, dna_alpha=dna_alpha):
    '''
    Determine if seq is a DNA string.
    Internal use only.

    :: seq
       type - string
       desc - a candidate sequence
    :: dna_alpha
       type - set
       desc - alphabet set considered
              valid for seq
              (defaul={A, T, G, C})
    '''

    if isinstance(seq, str) and \
       set(seq) <= dna_alpha:
        return True
    return False

def get_comp(seq):
    '''
    Return the complement of seq.
    Internal use only.

    :: seq
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    '''

    return seq.translate(complement_table)

def get_revcomp(seq):
    '''
    Return the reverse complement of seq.
    Internal use only.

    :: seq
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    '''

    return get_comp(seq)[::-1]

def stream_spectrum(seq, k):
    '''
    Stream the k-mer spectrum of seq.
    Internal use only.

    :: seq
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: k
       type - integer
       desc - k-mer length / k-value
    '''

    return (seq[i:i+k] for i in range(len(seq)-k+1))

def stream_canon_spectrum(seq, k):
    '''
    Stream the canonical k-mer spectrum
    of seq. Internal use only.

    :: seq
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: k
       type - integer
       desc - k-mer length / k-value
    '''

    return map(
        lambda x: min(x, get_revcomp(x)),
        stream_spectrum(seq=seq, k=k))

def get_hdist(seq1, seq2, max_hd=None):
    '''
    Compute the hamming distance between
    seq1 and seq2, within max_hd bound.
    Internal use only.

    :: seq1
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: seq2
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: max_hd
       type - integer / None
       desc - compute hamming distance
              between seq1 and seq2,
              within max_hd bound,
              return None otherwise
              (default=None)
    '''

    seq_len = min(len(seq1), len(seq2))
    i = 0
    hdist = 0
    while i < seq_len:
        hdist += seq1[i] != seq2[i]
        i += 1
        if not max_hd is None and \
           hdist > max_hd:
            return None
    return hdist

def get_edist(seq1, seq2, max_ed=None, mode='NW'):
    '''
    Compute the mode edit distance between
    seq1 and seq2, within max_ed bound.
    Internal use only.

    :: seq1
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: seq2
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: max_hd
       type - integer / None
       desc - compute edit distance
              between seq1 and seq2,
              within max_ed bound,
              return None otherwise
    :: mode
       type - string
       desc - 'NW' / 'HW' for global
              / infix alignments
    '''

    if len(seq2) > len(seq2):
        seq1, seq2 = seq2, seq1
    i = 0
    if max_ed is None:
        max_ed = -1
    ed_align = ed.align(
        seq2,
        seq1,
        mode=mode,
        task='distance',
        k=max_ed)
    edist = ed_align['editDistance']
    if edist == -1:
        return None
    return edist

# Workspace Functions

def get_adjusted_path(path, suffix):
    '''
    Adjust outfile with suffix.
    Internal use only.

    :: path
       type - string
       desc - some path which must
              have given suffix
    :: suffix
       type - string / None
       desc - suffix to concatenate
              with given path
    '''

    if path.endswith('/'):
        path = path.rstrip('/')
    if not suffix is None and \
       not path.endswith(suffix):
        path += str(suffix)
    return path

def get_path_status(
    path,
    suffix=None,
    readable=True,
    writable=False,
    creatable=False):
    '''
    Arbitrary path evaluation.
    Internal use only.

    :: path
       type - string
       desc - path to store or output
              information
    :: suffix
       type - string / None
       desc - required path suffix
              (default=None)
    :: readable
       type - boolean
       desc - if True check is path is
              readable  (file/dir   read)
              (default=True)
    :: writable
       type - boolean
       desc - if True check if path is
              writable  (file/dir  write)
              (default=False)
    :: creatable
       type - boolean
       desc - if True check if path is
              creatable (file/dir create)
              (default=False)
    '''

    # Status:
    # 0  = Path is invalid      input
    # 1  = Path is unreadable   file
    # 2  = Path is unwritable   file
    # 3  = Path is empty        file
    # 4  = Path is non-empty    file
    # 5  = Path is unreadable   directory
    # 6  = Path is unwritable   directory
    # 7  = Path is empty        directory
    # 8  = Path is non-empty    directory
    # 9  = Path is non-existent
    # 10 = Path is creatable
    # 11 = Path is uncreatable
    # X  = Path is existent     special file

    # Is the path valid?
    if not isinstance(path, str) or not path:
        return 0 # only non-empty strings are paths

    # Adjust Path for Checks
    path = get_adjusted_path(
        path=path,
        suffix=suffix)

    # Is the path in existence?
    if os.path.exists(path):

        # Path points to a file?
        if os.path.isfile(path):

            # Read accessible?
            if  readable and not \
                os.access(path, os.R_OK):
                    return 1 # Unreadable file

            # Is the file empty?
            if not os.stat(path).st_size:

                # Can we write to this file?
                if  writable and not \
                    os.access(path, os.W_OK):
                        return 2 # Unwritable file

                return 3 # Empty file

            else:
                return 4 # Non-empty file

        # Path points to a directory?
        elif os.path.isdir(path):

            # What's the accessibility like?
            if  readable and not \
                os.access(path, os.R_OK):
                    return 5 # Unreadable directory

            # Directory empty?
            if not any(os.scandir(path)):

                # Can we write to this directory?
                if  writable and not \
                    os.access(path, os.W_OK):
                        return 6 # Unwritable directory

                return 7 # Empty Directory

            else:
                return 8 # Non-empty Directory

        # Path points to a special file such as
        # socket, FIFO, device file etc.
        else:
            return 'X' # Special File Identifier

    # Is the path creatable?
    if creatable:

        try:
            # Build parent directory path
            parentpath = os.path.dirname(path) or \
                         os.getcwd()
            # Try building parent path
            os.makedirs(
                name=parentpath,
                exist_ok=True)

        except:
            return 11 # Cannot create path

        else:
            return 10 # Can create path

    # Path is non-existent
    return 9

def setup_directory(dirpath):
    '''
    Build dirpath. Internal use only.

    :: dirpath
       type - string
       desc - a directory built to
              store information
    '''

    dirpath = get_adjusted_path(
        path=dirpath,
        suffix=None)
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)

def setup_workspace(
    outfile,
    outfile_suffix):
    '''
    Adjust given output file and setup
    output directory, and schedule its
    removal at exit. Internal use only.

    :: outfile
       type - string
       desc - output file
    :: outfile_suffix
       type - string
       desc - output file suffix
    '''

    # Adjust outfile with suffix
    outfile = get_adjusted_path(
        path=outfile,
        suffix=outfile_suffix)

    # Define outdir
    outdir = '.'.join([
        outfile,
        str(uu.uuid4()),
        'dir'])

    # Schedule outdir deletion
    ae.register(
        remove_directory,
        outdir)

    # Setup outdir
    setup_directory(
        dirpath=outdir)

    # Return Workspace
    return outfile, outdir

def remove_file(filepath):
    '''
    Delete file stored at filepath.
    Internal use only.

    :: filepath
       type - string
       desc - a file to be removed
    '''

    if not filepath is None:
        filepath = get_adjusted_path(
            path=filepath,
            suffix=None)
        if os.path.exists(filepath):
            os.remove(filepath)

def remove_directory(dirpath):
    '''
    Delete dirpath and its contents.
    Internal use only.

    :: dirpath
       type - string
       desc - a directory (with content)
              to be removed
    '''

    if not dirpath is None:
        dirpath = get_adjusted_path(
            path=dirpath,
            suffix=None)
        if os.path.exists(dirpath):
            su.rmtree(dirpath)
