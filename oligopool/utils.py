import sys

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

