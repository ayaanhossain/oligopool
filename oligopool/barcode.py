import time  as tt

import collections as cx
import numpy       as np
import numba       as nb

import utils as ut


decoder = dict(zip((0., 1., 2., 3.), 'AGTC'))


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
    hfb = upb // 2
    qtb = upb // 4
    fld = np.random.permutation(qtb)

    # Split Field Index
    idx = 0

    # Stream Coordinates
    while idx < len(fld):
        coo = fld[idx]
        yield coo
        yield upb - 1 - coo
        yield hfb + coo
        yield hfb - 1 - coo
        idx += 1

def get_infinite_jumper(upb):
    '''
    Return infinite space jumper.
    Internal use only.

    :: upb
       type - integer
       desc - space upperbound
    '''

    # Setup RNG Permutation
    rng = np.random.default_rng()

    # Stream Coordinates
    while True:
        yield rng.integers(0, upb)

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
            # print(idx)
        except StopIteration:
            yield None # Finite Iterator Exhaused

        # Stream Coordinate Transformed to Barcode
        yield np.array(tuple(map(
                txf,
                np.base_repr(
                    idx,
                    base=4).zfill(
                        barcodelength))),
            dtype=np.float64)

def get_barcodeseq(barcode):
    '''
    Return the decoded sequence from
    barcode coordinate array.
    Internal use only.

    :: barcode
       type - np.array
       desc - encoded barcode array
    '''

    return ''.join(map(
        lambda x: decoder[x],
        barcode))

def is_hamming_feasible(store, count, minhdist):
    '''
    Determine if the minimum distance between
    barcode and store is greater or euqal to
    given minhdist (hamming feasibility).
    Internal use only.

    :: store
       type - np.array
       desc - vectorized storage of numeric
              encoding of all previous barcodes
    :: count
       type - integer
       desc - current storage fill count
    :: minhdist
       type - integer
       desc - minimum pairwise hamming
              distance between a pair
              of barcodes
    '''

    return ut.get_store_hdist(
        store=store,
        idx=count,
        direction=0) >= minhdist

def is_motif_feasible(
    barcodeseq,
    barcodelength,
    exmotifs):
    '''
    Determine if the barcode does not contain or is
    contained in one of the exluded motifs (motif
    feasibility). Internal use only.

    :: barcodeseq
       type - string
       desc - decoded barcode sequence
    :: barcodelength
       type - string
       desc - barcode length
    :: exmotifs
       type - list / None
       desc - list of motifs to exclude
              in designed barcodes,
              None otherwise
    '''

    return ut.get_motif_conflict(
        seq=barcodeseq,
        seqlen=barcodelength,
        exmotifs=exmotifs,
        partial=False)

def is_assignment_feasible(
    barcodeseq,
    exmotifs,
    lcifn,
    rcifn,
    cntxlen,
    carr):
    '''
    Determine the context assignment index
    for designed sequence. Internal use only.

    :: barcodeseq
       type - string
       desc - decoded barcode sequence
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
    
    return ut.get_assignment_index(
        seq=barcodeseq,
        exmotifs=exmotifs,
        lcifn=lcifn,
        rcifn=rcifn,
        cntxlen=cntxlen,
        carr=carr)

def is_barcode_feasible(
    store,
    count,
    minhdist,
    barcodeseq,
    barcodelength,
    exmotifs,
    lcifn,
    rcifn,
    cntxlen,
    carr):
    '''
    Determine if the barcode is hamming,
    motif and edge satisfiable.
    Internal use only.

    :: store
       type - np.array
       desc - vectorized storage of numeric
              encoding of all previous barcodes
    :: count
       type - integer
       desc - current storage fill count
    :: minhdist
       type - integer
       desc - minimum pairwise hamming
              distance between a pair
              of barcodes
    :: barcodeseq
       type - string
       desc - decoded barcode sequence
    :: barcodelength
       type - string
       desc - barcode length
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

    # Hamming Condition
    hcond = is_hamming_feasible(
        store=store,
        count=count,
        minhdist=minhdist)
    if not hcond:
        return (False, 1, None, None)       # Hamming Failure

    # Motif Embedding
    mcond, motif = is_motif_feasible(
        barcodeseq=barcodeseq,
        barcodelength=barcodelength,
        exmotifs=exmotifs)
    if not mcond:
        return (False, 2, None, {motif: 1})  # Motif Failure

    # Assignment Feasibility (Edge-Effects)
    acond, aidx, afails = is_assignment_feasible(
        barcodeseq=barcodeseq,
        exmotifs=exmotifs,
        lcifn=lcifn,
        rcifn=rcifn,
        cntxlen=cntxlen,
        carr=carr)
    if not acond:
        return (False, 3, None, afails)      # Edge Failure

    # All conditions met!
    return (True, 0, aidx, None)

def show_update(
    count,
    plen,
    barcodeseq,
    cf,
    stm,
    liner):
    '''
    Display the current progress in barcode
    generation. Internal use only.

    :: count
       type - integer
       desc - current storage fill count
    :: plen
       type - integer
       desc - barcode index width
    :: barcodeseq
       type - np.array
       desc - numeric vector encoding barcode
    :: cf
       type - boolean
       desc - barcode feasibility status
    :: stm
       type - integer
       desc - feasibility failure state marker
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    liner.send(' [Slot: {:{},d}] Barcode {} is {} {}'.format(
        count,
        plen,
        barcodeseq,
        ['Rejected', 'Accepted'][cf],
        ['',
         'due to Hamming Infeasibility',
         'due to Motif Infeasibility',
         'due to Edge Infeasibility'][stm]))

def barcode_engine(
    targetsize,
    barcodelength,
    minhdist,
    exmotifs,
    leftcontext,
    rightcontext,
    liner):
    '''
    Return barcodes fulfilling constraints.
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
    :: exmotifs
       type - list / None
       desc - list of motifs to exclude
              in designed barcodes,
              None otherwise
    :: leftcontext
       type - list / None
       desc - list of sequence to the
              left of barcodes,
              None otherwise
    :: rightcontext
       type - list / None
       desc - list of sequences to the
              right of barcodes,
              None otherwise
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    count  = 0                            # Barcode Count
    store  = np.zeros(                    # Store Encoded Barcodes
        (targetsize, barcodelength),
        dtype=np.float64)
    codes  = []                           # Store Decoded Barcodes
    aarr   = []                           # Assignment Array
    mfails = cx.Counter()                 # Motif Fail Counter

    # Determine Jumper
    jtp = get_jumpertype(
        barcodelength=barcodelength)
    
    # Context Setup
    carr   = cx.deque(range(targetsize))  # Context Array
    lcifn  = ut.get_context_inference_fn( # Left  Context Selector
        context=leftcontext)
    rcifn  = ut.get_context_inference_fn( # Right Context Selector
        context=rightcontext)
    
    # Infinite Jumper Failure
    prob  = ut.get_prob(                  # Probability of Success
        success=1,
        trials=np.power(
            4, barcodelength - minhdist))
    trial = ut.get_trials(                # Trials Required
        prob=prob)
    sscnt = 1                             # Success Count
    flcnt = 0                             # Failure Count
    excnt = trial                         # Trial   Count

    # Error Type Fails
    hammfc  = 0
    motiffc = 0
    edgefc  = 0
    
    # Verbage
    verbage_reach  = 0
    verbage_target = np.random.randint(
        *map(np.round, (targetsize * 0.080,
                        targetsize * 0.120)))
    plen = len(str(targetsize)) + \
        int(np.log10(targetsize) / 3)

    # Setup Exmotifs
    cntxlen = 0
    if exmotifs:
        liner.send('\n[Preprocessing Excluded Motifs]\n')
        exmotifs = ut.prep_exmotifs(
            exmotifs=exmotifs,
            packing=tuple,
            liner=liner)
        cntxlen = ut.get_context_len(
            exmotifs=exmotifs)

    # Show Update
    liner.send('\n[Computing {}-bp Barcodes]\n'.format(
        barcodelength))

    # Generator Setup
    barcodes = stream_barcodes(
        barcodelength=barcodelength)
    barcode  = None # Current Candidate
    accseq   = None # Last Successful Candidate
    
    # Build Barcodes
    while True:

        # Sample a Barcode in Space
        barcode = next(barcodes)

        # Space Exhausted?
        if barcode is None:
            
            # Final Update
            show_update(
                count=count,
                plen=plen,
                barcodeseq=accseq,
                cf=True,
                stm=0,
                liner=liner)

            # Solution Status
            liner.send('\* Targetsize Reached? No\n')

            return (0,
                None,
                hammfc,
                motiffc,
                edgefc,
                mfails) # No Solution

        # Update Barcode Store
        store[count, :] = barcode

        # Decode Barcode Sequence
        barcodeseq = get_barcodeseq(barcode)
        
        # Check Feasibility
        cf, stm, aidx, afails = is_barcode_feasible(
            store=store,
            count=count,
            minhdist=minhdist,
            barcodeseq=barcodeseq,
            barcodelength=barcodelength,
            exmotifs=exmotifs,
            lcifn=lcifn,
            rcifn=rcifn,
            cntxlen=cntxlen,
            carr=carr)

        # Inifinite Jumper Book-keeping Update
        if jtp == 2:
            excnt += 1
        
        # Accept Sample into Store?
        if cf:

            # Update Store Fill Count
            count += 1

            # Show Update on Success
            show_update(
                count=count,
                plen=plen,
                barcodeseq=barcodeseq,
                cf=cf,
                stm=stm,
                liner=liner)

            # Record Assignment Index
            aarr.append(aidx)

            # Update Codes
            codes.append(barcodeseq)

            # Update Last Accounted Barcode
            accseq = barcodeseq
            
            # Inifinite Jumper Book-keeping Update
            if jtp == 2:
                sscnt += 1
                prob  = ut.get_prob(
                    success=sscnt,
                    trials=excnt)
                trial = ut.get_trials(
                    prob=prob)
                flcnt = 0
        
        # Rejected from Store
        else:

            # Update Fail counters
            if stm == 1:
                hammfc  += 1
            if stm == 2:
                motiffc += 1
            if stm == 3:
                edgefc  += 1

            # Inifinite Jumper Book-keeping Update
            if jtp == 2:
                flcnt += 1

            # Record Context Failure Motifs
            if not afails is None:
                mfails += afails

        # Verbage Book-keeping
        if verbage_reach >= verbage_target:
            show_update(
                count=count,
                plen=plen,
                barcodeseq=barcodeseq,
                cf=cf,
                stm=stm,
                liner=liner)
            verbage_reach = -1
        verbage_reach += 1
        
        # Target Reached?
        if count == targetsize:
            
            # Construct the Sorted Barcodes
            codes = [code for aidx,code in sorted(
                zip(aarr, codes))]

            # Solution Status
            liner.send('\* Targetsize Reached? Yes\n')

            # Return Solution
            return (1,
                codes,
                hammfc,
                motiffc,
                edgefc,
                None) # Solution Reached!

        # Trials Exhausted for Inifinite Jumper?
        if jtp == 2:
            if flcnt == trial:

                # Final Update
                show_update(
                    count=count,
                    plen=plen,
                    barcodeseq=accseq,
                    cf=True,
                    stm=0,
                    liner=liner)

                # Solution Status
                liner.send('\* Targetsize Reached? No\n')

                return (0,
                    None,
                    hammfc,
                    motiffc,
                    edgefc,
                    mfails) # No Solution

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

def get_context():
    with open('promoters.txt') as infile:
        cntx = [x.strip() for x in infile.readlines()]
    return cntx

def main():

    liner = ut.liner_engine()

    t0 = tt.time()

    # exmotifs = ['AAAA', 'GGGG', 'CCCC', 'TTTT', 'GGATCC', 'TCTAGA', 'GAATTC'][::-1]
    exmotifs = ['GGATCC', 'TCTAGA']
    # exmotifs = []

    lcntx = get_context()
    
    results = barcode_engine(
        targetsize=len(lcntx),
        barcodelength=10,
        minhdist=3,
        exmotifs=exmotifs,
        leftcontext=lcntx,
        rightcontext=lcntx,
        liner=liner)

    if results[0] == 1:
        print('\n{} Barcodes Stored\n'.format(len(results[1])))
    else:
        if results[-1]:
            print('\n{} Motifs Conflicting\n'.format(len(results[-1])))
            print(results[-1].most_common())
    print('Hamming Conlict: {}'.format(results[-4]))
    print('Motif   Conlict: {}'.format(results[-3]))
    print('Edge    Conlict: {}'.format(results[-2]))

    print(tt.time()-t0)

if __name__ == '__main__':
    main()