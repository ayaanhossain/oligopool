import time  as tt

import numpy as np


def hyperion(l):
    '''
    Hyper jump through sequence space and
    stream barcodes. Internal use only.

    :; l
       type - integer
       desc - required barcode length
    '''

    # RNG Setup
    rng = np.random.default_rng()
    lwb = 0
    upb = np.power(4, l)
    txf = lambda x: float(ord(x))

    # Iterator Setup
    if l <= 12:
        # Finite Random Walk
        hype = iter(np.random.permutation(upb))
    else:
        # Infinite Random Walk
        hype = iter(lambda: rng.integers(lwb, upb), None)
    
    # Hyper Jump!
    while True:
        # Fetch a Walk Coordinate
        try:
            idx = next(hype)
        except StopIteration:
            yield None # Finite Iterator Exhaused

        # Stream Coordinate Transformed to Barcode
        yield tuple(map(txf, np.base_repr(
            idx, base=4).zfill(l)))

def barcode_engine(
    target_size,
    barcode_length,
    min_hdist):
    '''
    Given a library size of t, generate
    l-bp barcodes such that the minimum
    hamming distance between any two
    barcodes is d.

    :: target_size
       type - integer
       desc - required library size
    :: barcode_length
       type - integer
       desc - required barcode length
    :: min_hdist
       type - integer
       desc - minimum pairwise hamming
              distance between a pair
              of barcodes
    '''

    # Book-keeping
    fill  = 0
    store = np.zeros((
        target_size,
        barcode_length))

    # Generator Setup
    barcodes = hyperion(l=barcode_length)
    
    # Build Barcodes
    while True:

        # Sample a Barcode in Space
        barcode = next(barcodes)

        # Space Exhausted?
        if barcode is None:
            return store[:fill, :]
        
        # Compare Sample with Stored
        if fill:
            _d = (store[:fill, :] != barcode).sum(1).min()
        else:
            _d = min_hdist
        
        # Accept Sample into Store?
        if _d >= min_hdist:
            store[fill, :] = barcode
            fill += 1
            if (fill+1) % 100 == 0:
                print(fill, _d, barcode)
        
        # Target Reached?
        if fill == target_size:
            return store[:fill, :]

def main():
    t0 = tt.time()
    print(barcode_engine(
        target_size=100000,
        barcode_length=10,
        min_hdist=3).shape)
    print(tt.time()-t0)

if __name__ == '__main__':
    main()