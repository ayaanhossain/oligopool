import sys

import time    as tt

import collections as cx

import numpy   as np
import numba   as nb
import primer3 as p3


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
        sys.stdout.write('\n')

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

# Oligo Functions

@nb.njit
def get_hdist(
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

# Motif Functions

def prep_exmotifs(exmotifs, packing, liner):
    '''
    Return all exmotifs sorted by
    length. Internal use only.

    :: exmotifs
       type - list
       desc - list of all motifs
              to be excluded
    :: packing
       type - function
       desc - factory function for
              wrapping exmotifs
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Sort motifs by length
    liner.send(' Sorting Motifs ...')
    t0 = tt.time()
    exmotifs.sort(key=len)
    liner.send(' Sorted {} Motifs in {:.2f} sec\n'.format(
        len(exmotifs), tt.time()-t0))

    # Enque all motifs
    liner.send(' Enqueing Motifs ...')
    t0 = tt.time()
    dq = packing(exmotifs)
    liner.send(' Enqued {} Motifs in {:.2f} sec\n'.format(
        len(exmotifs), tt.time()-t0))

    # Return Results
    return dq

def get_context_len(exmotifs):
    '''
    Get the exmotif context length.
    Internal use only.

    :: exmotifs
       type - iterable / None
       desc - iterable of motifs to exlude,
              None otherwise
    '''

    return len(exmotifs[-1]) - 1

def get_context_inference_fn(context):
    '''
    Return selector functions for building
    barcode context. Internal use only.

    :: context
       type - list / string / None
       desc - list of context sequence,
              or a single context sequence,
              or None
    '''

    # All contexts unique
    if isinstance(context, list):
        cifn = lambda x: context[x]

    # All contexts constant
    elif isinstance(context, str):
        cifn = lambda x: context
    
    # No context
    elif context == None:
        cifn = lambda x: ''

    # Return Selector
    return cifn

def get_motif_conflict(
    seq,
    seqlen,
    exmotifs):
    '''
    Determine if the barcode does not contain or is
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
    '''

    # Do we have anything to exclude?
    if exmotifs:

        # Loop through all motifs
        for motif in exmotifs:

            # Embedding Conflict
            if len(motif) <= seqlen and \
               motif in seq:
                return False, motif

            # Embedded Conflict
            if len(motif)  > seqlen and \
               seq in motif:
                return False, motif
    
    # No Motif Conflct
    return (True, None)

def get_assignment_index(
    seq,
    exmotifs,
    lcifn,
    rcifn,
    cntxlen,
    carr):
    '''
    Determine the context assignment index
    for designed sequence. Internal use only.

    :: seq
       type - string
       desc - designed sequence to be
              assigned
    :: exmotifs
       type - list / None
       desc - list of motifs to exclude
              in designed barcodes,
              None otherwise
    :: lcifn
       type - lambda
       desc - selector for the left
              sequence context
    :: rcifn
       type - lambda
       desc - selector for the right
              sequence context
    :: cntxlen
       type - lambda
       desc - total length of context
              to extract from either
              left or right
    :: carr
       type - np.array
       desc - context assignment array
    '''
    
    # Book-keeping
    i      = 0
    cfails = cx.Counter()

    # Loop through all contexts for assignment
    while i < len(carr):

        # Fetch Context
        aidx = carr.popleft()

        # Fetch left and right contexts
        lcntx = lcifn(aidx)[-cntxlen:]
        rcntx = rcifn(aidx)[:+cntxlen]

        # Build out in-context sequence
        incntxseq = lcntx + seq + rcntx

        # Determine context feasibility
        mcond, motif = get_motif_conflict(
            seq=incntxseq,
            seqlen=len(incntxseq),
            exmotifs=exmotifs)

        # We have a winner folks!
        if mcond:
            return (True, aidx, None)

        # We lost an assignment, record cause
        else:
            cfails[motif] += 1

        # Update Iteration
        i += 1
        
        # We will try this context again
        carr.append(aidx)

    # Failed Assignment
    return (False, None, cfails)