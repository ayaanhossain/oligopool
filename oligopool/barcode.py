import time  as tt

import collections as cx
import numpy       as np
import numba       as nb

import utils as ut


def get_jumpertype(barcodelength):
    '''
    Return appropriate jumper type.
    Internal use only.

    :: barcodelength
       type - integer
       desc - required barcode length
    '''

    if barcodelength <= 12:
        return 1 # Finite   Random Jumps
    else:
        return 2 # Infinite Random Jumps

def get_finite_jumper(upb):
    '''
    Return finite space jumper.
    Internal use only.

    :: upb
       type - integer
       desc - space uppberbound
    '''
    
    # Setup Permutation
    idx = 0
    qtb = upb // 4
    fld = np.random.permutation(qtb)
    hfl = (upb // 2) - 1
    upl = upb - 1
    
    # Stream Coordinates
    while idx < len(fld):
        coo = fld[idx]
        yield coo
        yield upl - coo
        yield hfl + coo
        yield hfl - coo
        idx += 1

def get_infinite_jumper(upb):
    '''
    Return infinite space jumper.
    Internal use only.

    :: upb
       type - integer
       desc - space upperbound
    '''

    # Setup RNG
    rng = np.random.default_rng()
    qtb = upb // 4
    hfl = (upb // 2) - 1
    upl = upb - 1

    # Stream Coordinates
    while True:
        coo = rng.integers(0, qtb)
        yield coo
        yield upl - coo
        yield hfl + coo
        yield hfl - coo


def stream_barcodes(barcodelength):
    '''
    Randomly jump through sequence space and
    stream barcodes. Internal use only.

    :: barcodelength
       type - integer
       desc - required barcode length
    '''

    # Jumper Setup
    upb = np.power(4, barcodelength)
    jtp = get_jumpertype(
        barcodelength=barcodelength)
    if jtp == 1:
        jumper = get_finite_jumper(upb)
    else:
        jumper = get_infinite_jumper(upb)

    # Transformer Setup
    txf = lambda x: float(ord(x)) - 48.
    
    # Stream Jumps!
    while True:
        # Fetch a Jump Coordinate
        try:
            idx = next(jumper)
        except StopIteration:
            yield None # Finite Iterator Exhaused

        # Stream Coordinate Transformed to Barcode
        yield np.array(
            tuple(map(
                txf,
                np.base_repr(idx, base=4).zfill(
                    barcodelength))),
            dtype=np.float64)

@nb.njit
def is_barcode_feasible(
    store,
    barcode,
    count,
    minhdist):
    '''
    Determine if the minimum distance between
    barcode and store is greater or equal to
    given minhdist (feasibility).
    Internal use only.

    :: store
       type - np.array
       desc - vectorized storage of numeric
              encoding of all previous barcodes
    :: barcode
       type - np.array
       desc - numeric vector encoding barcode
    :: count
       type - integer
       desc - current storage fill count    
    '''

    return (store[:count, :] != barcode).sum(
        1).min() >= minhdist

def show_update(
    count,
    plen,
    barcode,
    cf,
    liner):
    '''
    Display the current progress in barcode
    generation. Internal use only.

    :: count
       type - integer
       desc - current storage fill count
    :: barcode
       type - np.array
       desc - numeric vector encoding barcode
    :: cf
       type - boolean
       desc - barcode feasibility status
    :: plen
       type - integer
       desc - barcode index width
    :: liner
       type - coroutine
       desc - dynamic printing
    '''
    liner.send(' [Slot: {:{},d}] Candidate {} {}'.format(
        count,
        plen,
        barcode,
        # ''.join(map(lambda x: chr(int(x)), barcode)),
        ['Rejected', 'Accepted'][cf]))

def barcode_engine(
    targetsize,
    barcodelength,
    minhdist,
    liner):
    '''
    Stream barcodes fulfilling constraints.
    Internal use only.

    :: targetsize
       type - integer
       desc - required library size
    :: barcodelength
       type - integer
       desc - required barcode length
    :: minhdist
       type - integer
       desc - minimum pairwise hamming
              distance between a pair
              of barcodes
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Determine Jumper
    jtp = get_jumpertype(
        barcodelength=barcodelength)

    # Book-keeping
    count = 0
    store = np.zeros(
        (targetsize, barcodelength),
        dtype=np.float64)
    
    prob  = ut.get_prob(
        success=1,
        trials=np.power(
            4, barcodelength - minhdist))
    trial = ut.get_trials(
        prob=prob)
    sscnt = 1
    flcnt = 0
    excnt = trial
    
    verbage_reach  = 0
    verbage_target = np.random.randint(
        *map(np.round, (targetsize * 0.080,
                        targetsize * 0.120)))
    plen = len(str(targetsize)) + \
        int(np.log10(targetsize) / 3)

    # Generator Setup
    barcodes = stream_barcodes(
        barcodelength=barcodelength)
    barcode  = None
    acccode  = None
    
    # Build Barcodes
    while True:

        # Sample a Barcode in Space
        barcode  = next(barcodes)

        # Space Exhausted?
        if barcode is None:
            show_update(
                count=count,
                plen=plen,
                barcode=acccode,
                cf=True,
                liner=liner)
            return store[:count, :]
        
        # Compare Sample with Stored
        if count:
            cf = is_barcode_feasible(
                store=store,
                barcode=barcode,
                count=count,
                minhdist=minhdist)
        else:
            cf = True

        # Inifinite Jumper Book-keeping Update
        if jtp == 2:
            excnt += 1
        
        # Accept Sample into Store?
        if cf:

            # Show Update on Success
            show_update(
                count=count,
                plen=plen,
                barcode=barcode,
                cf=cf,
                liner=liner)

            # Update Barcode Store
            store[count, :] = barcode
            acccode = barcode

            # Update Store Fill Count
            count += 1
            
            # Inifinite Jumper Book-keeping Update
            if jtp == 2:
                sscnt += 1
                prob  = ut.get_prob(
                    success=sscnt,
                    trials=excnt)
                trial = ut.get_trials(
                    prob=prob)
                flcnt = 0
        else:
            # Inifinite Jumper Book-keeping Update
            if jtp == 2:
                flcnt += 1

        # Verbage Book-keeping
        if verbage_reach >= verbage_target:
            show_update(
                count=count,
                plen=plen,
                barcode=barcode,
                cf=cf,
                liner=liner)
            verbage_reach = -1
        verbage_reach += 1
        
        # Target Reached?
        if count == targetsize:
            show_update(
                count=count,
                plen=plen,
                barcode=barcode,
                cf=cf,
                liner=liner)
            return store[:count, :]

        # Trials Exhausted for Inifinite Jumper?
        if jtp == 2:
            if flcnt == trial:
                show_update(
                    count=count,
                    plen=plen,
                    barcode=acccode,
                    cf=True,
                    liner=liner)
                return store[:count, :]

def get_stats():
    '''
    TBC, later.
    '''
    pass

def barcode():
    '''
    Given a library size of t, generate
    l-bp barcodes such that the minimum
    hamming distance between any two
    barcodes is d.

    Driver, later
    '''
    pass

def main():
    liner = ut.liner_engine()

    t0 = tt.time()
    
    store = barcode_engine(
        targetsize=1000,
        barcodelength=10,
        minhdist=3,
        liner=liner)

    print('\n{}'.format(store.shape))

    print(tt.time()-t0)

if __name__ == '__main__':
    main()