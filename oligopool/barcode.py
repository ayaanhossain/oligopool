import time  as tt

import collections as cx
import numpy       as np
import numba       as nb
import bitarray    as ba

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

def is_hamming_feasible(store, count):
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
    '''

    return ut.get_hdist(
        store=store,
        idx=count,
        direction=0)

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
       desc - barcode sequence
    :: exmotifs
       type - list / None
       desc - list of motifs to exclude
              in designed barcodes,
              None otherwise
    '''

    if exmotifs:

        # Loop through all motifs
        for motif in exmotifs:

            # Embedding Conflict
            if len(motif) <= barcodelength and \
               motif in barcodeseq:
                return False, motif

            # Embedded Conflict
            if len(motif)  > barcodelength and \
               barcodeseq in motif:
                return False, motif
    
    # No Motif Conflct
    return (True, None)

def is_barcode_feasible(
    store,
    count,
    minhdist,
    barcodeseq,
    barcodelength,
    exmotifs):
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
    '''

    # Hamming Condition
    hcond = is_hamming_feasible(
        store=store,
        count=count) >= minhdist
    if not hcond:
        return (False, 1) # Hamming Failure

    # Motif Embedding
    mcond, motif = is_motif_feasible(
        barcodeseq=barcodeseq,
        barcodelength=barcodelength,
        exmotifs=exmotifs)
    if not mcond:
        return (False, 2) # Motif Failure

    # All conditions met!
    return (True, 0)

def get_barcodeseq(barcode):

    return ''.join(map(
        lambda x: decoder[x],
        barcode))

def get_context_inference_fn(context):

    if isinstance(context, list):
        linf = lambda x: context[x]
    elif isinstance(context, str):
        linf = lambda x: context
    elif context == None:
        linf = lambda x: ''

    return linf

def get_assignment_index(
    barcodeseq,
    exmotifs,
    linfn,
    rinfn,
    cntxlen,
    midx,
    marr):
    
    aidx = midx

    while aidx < len(marr):

        if marr[aidx] is True:
            # print(marr[aidx])
            aidx += 1
            continue

        lcntx = linfn(aidx)[-cntxlen:]
        rcntx = rinfn(aidx)[:+cntxlen]

        incntxseq = lcntx + barcodeseq + rcntx

        mcond, motif = is_motif_feasible(
            barcodeseq=incntxseq,
            barcodelength=len(incntxseq),
            exmotifs=exmotifs)

        if mcond:
            return (True, aidx)

        else:
            print('{}: {} | {} | {} | {}'.format(aidx, motif, linfn(aidx), incntxseq, barcodeseq))
        
        aidx += 1

    return (False, None)


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

    liner.send(' [Slot: {:{},d}] Barcode {} is {} {}\n'.format(
        count,
        plen,
        barcodeseq,
        ['Rejected', 'Accepted'][cf],
        ['',
         'due to Hamming Infeasibility',
         'due to Motif Infeasibility',
         'due to Context Infeasibility'][stm]))

    # if stm == 3:
    #     tt.sleep(5)

def barcode_engine(
    targetsize,
    barcodelength,
    minhdist,
    exmotifs,
    leftcontext,
    rightcontext,
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

    # Determine Jumper
    jtp = get_jumpertype(
        barcodelength=barcodelength)

    # Book-keeping
    count = 0
    store = np.zeros(
        (targetsize, barcodelength),
        dtype=np.float64)
    
    mstore = np.ones(
        targetsize,
        dtype=np.float64) * -1.
    marr  = ba.bitarray('0'*targetsize)
    linfn = get_context_inference_fn(
        context=leftcontext)
    rinfn = get_context_inference_fn(
        context=rightcontext)
    midx  = 0
    
    codes = []
    
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
            show_update(
                count=count,
                plen=plen,
                barcodeseq=accseq,
                cf=True,
                stm=0,
                liner=liner)
            return store[:count, :]

        # print(barcode)

        # Update Barcode Store
        store[count, :] = barcode

        # Decode Barcode Sequence
        barcodeseq = get_barcodeseq(barcode)
        # print(barcodeseq)
        
        # Check Feasibility
        cf, stm = is_barcode_feasible(
            store=store,
            count=count,
            minhdist=minhdist,
            barcodeseq=barcodeseq,
            barcodelength=barcodelength,
            exmotifs=exmotifs)

        # if cf:
        #     print('Basic OK')

        # Get Barcode Matching / Assignment
        if cf:
            cf, aidx = get_assignment_index(
                barcodeseq=barcodeseq,
                exmotifs=exmotifs,
                linfn=linfn,
                rinfn=rinfn,
                cntxlen=cntxlen,
                midx=midx,
                marr=marr)

            if not cf:
                stm = 3

            # liner.send('Barcode {} -> Context {}\n'.format(midx, aidx))

        # Inifinite Jumper Book-keeping Update
        if jtp == 2:
            excnt += 1
        
        # Accept Sample into Store?
        if cf:

            # Update Matching Index
            if not aidx is None:
                if aidx == midx:
                    midx += 1
                marr[aidx] = True
                mstore[count] = aidx

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

            # Update Codes
            codes.append(barcodeseq)

            # Update Last Accounted Barcode
            accseq = barcode
            
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
                barcodeseq=barcodeseq,
                cf=cf,
                stm=stm,
                liner=liner)
            verbage_reach = -1
        verbage_reach += 1
        
        # Target Reached?
        if count == targetsize:
            return store[:count, :]

        # Trials Exhausted for Inifinite Jumper?
        if jtp == 2:
            if flcnt == trial:
                show_update(
                    count=count,
                    plen=plen,
                    barcodeseq=accseq,
                    cf=True,
                    stm=0,
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

def get_context():
    with open('promoters.txt') as infile:
        cntx = [x.strip() for x in infile.readlines()]
    return cntx

def main():

    liner = ut.liner_engine()

    t0 = tt.time()

    # exmotifs = ['AAAA', 'GGGG', 'CCCC', 'TTTT', 'GGATCC', 'TCTAGA', 'GAATTC'][::-1]
    exmotifs = ['GGATCC', 'TCTAGA']

    lcntx = get_context()
    
    store = barcode_engine(
        targetsize=len(lcntx),
        barcodelength=10,
        minhdist=3,
        exmotifs=exmotifs,
        leftcontext=lcntx,
        rightcontext=None,
        liner=liner)

    print('\n{}'.format(store.shape))

    print(tt.time()-t0)

if __name__ == '__main__':
    main()