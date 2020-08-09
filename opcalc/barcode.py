import time  as tt

import collections as cx
import numpy       as np

import utils as ut


def get_walkertype(barcodelength):
    '''
    Return appropriate walker type.
    Internal use only.

    :: barcodelength
       type - integer
       desc - required barcode length
    '''

    if barcodelength <= 12:
        return 1 # Finite   Random Walker
    else:
        return 2 # Infinite Random Walker

def get_finite_walker(upb):
    '''
    Return finite space walker.
    Internal use only.

    :: upb
       type - integer
       desc - space uppberbound
    '''
    
    # Setup Basis
    deq = cx.deque(np.random.permutation(upb))
    
    # Stream Coordinates
    while deq:
        yield deq.popleft()

def get_infinite_walker(upb):
    '''
    Return infinite space walker.
    Internal use only.

    :: upb
       type - integer
       desc - space upperbound
    '''

    # Setup Basis
    rng = np.random.default_rng()

    # Stream Coordinates
    while True:
        yield rng.integers(0, upb)


def hyperion(barcodelength):
    '''
    Hyper jump through sequence space and
    stream barcodes. Internal use only.

    :: barcodelength
       type - integer
       desc - required barcode length
    '''

    # Walker Setup
    upb = np.power(4, barcodelength)
    wtp = get_walkertype(
        barcodelength=barcodelength)
    if wtp == 1:
        walker = get_finite_walker(upb)
    else:
        walker = get_infinite_walker(upb)

    # Transformer Setup
    txf = lambda x: float(ord(x))
    
    # Hyperion Jump!
    while True:
        # Fetch a Coordinate
        try:
            idx = next(walker)
        except StopIteration:
            yield None # Finite Iterator Exhaused

        # Stream Coordinate Transformed to Barcode
        yield tuple(map(txf, np.base_repr(
            idx, base=4).zfill(barcodelength)))

def barcode_engine(
    targetsize,
    barcodelength,
    minhdist):
    '''
    Given a library size of t, generate
    l-bp barcodes such that the minimum
    hamming distance between any two
    barcodes is d.

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
    '''

    # Determine Walker
    wtp = get_walkertype(
        barcodelength=barcodelength)

    # Book-keeping
    count = 0
    store = np.zeros((
        targetsize,
        barcodelength))
    prob  = ut.get_prob(
        success=1,
        trials=np.power(
            4, barcodelength - minhdist))
    trial = ut.get_trials(
        prob=prob)
    sscnt = 1
    flcnt = 0
    excnt = trial

    # Generator Setup
    barcodes = hyperion(
        barcodelength=barcodelength)
    
    # Build Barcodes
    while True:

        # Sample a Barcode in Space
        barcode = next(barcodes)

        # Space Exhausted?
        if barcode is None:
            return store[:count, :]
        
        # Compare Sample with Stored
        if count:
            _d = (store[:count, :] != barcode).sum(1).min()
        else:
            _d = minhdist

        if wtp == 2:
            excnt += 1
        
        # Accept Sample into Store?
        if _d >= minhdist:
            store[count, :] = barcode
            count += 1            
            
            # Inifinite Walker Book-keeping Update
            if wtp == 2:
                sscnt += 1
                prob  = ut.get_prob(
                    success=sscnt,
                    trials=excnt)
                trial = ut.get_trials(
                    prob=prob)
                flcnt = 0

            if (count+1) % 10 == 0:
                print(count, _d, barcode)
        else:
            # Inifinite Walker Book-keeping Update
            if wtp == 2:
                flcnt += 1
        
        # Target Reached?
        if count == targetsize:
            return store[:count, :]

        # Trials Exhausted for Inifinite Walker?
        if wtp == 2:
            if flcnt == trial:
                return store[:count, :]

def main():
    t0 = tt.time()
    print(barcode_engine(
        targetsize=80000,
        barcodelength=8,
        minhdist=5).shape)
    print(tt.time()-t0)

if __name__ == '__main__':
    main()