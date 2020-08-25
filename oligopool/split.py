import time  as tt

import collections as cx

import numpy as np
import numba as nb

import utils as ut


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

def get_seqvec(seq):
    '''
    Return the numeric vector for seq.
    Internal use only.

    :: seq
       type - string
       desc - a sequence to split
    '''
    
    return np.array(
        tuple(float(ord(nt)) for nt in seq),
        dtype=np.float64)

def get_seqmat(seqlist):
    '''
    Return the numeric representation
    of seqlist. Internal use only.

    :: seqlist
       type - iterable
       desc - list of sequences to split
    '''
    
    seqmat = np.zeros(
        (len(seqlist), len(seqlist[0])),
        dtype=np.float64)
    for idx,seq in enumerate(seqlist):
        seqmat[idx, :] = get_seqvec(seq=seq)
    return np.array(seqmat, dtype=np.float64)

def get_entvec(seqmat):
    '''
    Return entropy of sequence matrix.
    Internal use only.

    :: seqmat
       type - np.array
       desc - numeric sequence matrix
    '''

    # Compute Count Frequency
    d = np.zeros(
        (4, seqmat.shape[1]),
        dtype=np.float64)
    p = np.zeros(4)
    for i in range(seqmat.shape[1]):
        m = np.unique(
            seqmat[:, i],
            return_counts=True)[1]
        p[:m.shape[0]] = m
        d[:, i] += p
        p *= 0.

    # Convert to Normalized Distribution
    d = d / d.sum(0)

    # Return Results
    return cx.deque(np.abs(
        (d*(np.log(d, where=d>0.) / np.log(4))).sum(0)))

def get_varcont(entvec):
    '''
    Return all variable region span indices.
    Internal use only.

    :: entvec
       type - cx.deque
       desc - positional entropy
    '''

    # Setup Parsing
    varcont = cx.deque()
    start   = None
    end     = None
    i       = -1

    # Parse Contigs
    while entvec:

        # Update index
        i += 1

        # Extract entropy
        ent = entvec.popleft()

        # Constant Region
        if ent <= 0.25:

            # A contig built?
            if not start is None and \
               not end   is None:
                varcont.append((start, end))
            
            # Reset for next contig
            start = None
            end   = None

            # Next
            continue

        # Variable Region
        else:

            # Contig continues?
            if start is None:
                start = i
            if end is None:
                end = i
            
            # Update ending index
            end += 1

            # Next
            continue

    # Remnant Update
    if not end is None:
        varcont.append((start, end))

    # Return Results
    return varcont

def get_merged_varcont(varcont, mergefactor):
    '''
    Return a merged varcont, merging variables
    separated by at most mergefactor constants.
    Internal use only.

    :: varcont
       type - cx.deque
       desc - all variable region span indices
    :: mergefactor
       type - integer
       desc - maximum gap length between two
              variable regions to be merged
              into a single contig
    '''

    # Setup
    merged  = cx.deque()
    current = list(varcont.popleft())
    
    # Merge Contigs
    while varcont:
        contig = varcont.popleft()
        if contig[0]-current[1] <= mergefactor:
            current[1] = contig[1]
        else:
            merged.append(tuple(current))
            current = list(contig)
    merged.append(tuple(current))
    
    # Return Result
    return merged

def get_endpoints(varcont, start, pcend, splitlen):
    '''
    Return all end-point variable contigs given
    current oligo starting point and splitlen.
    Internal use only.

    :: varcont
       type - cx.deque
       desc - all variable region span indices
    :: start
       type - integer
       desc - current split starting point
    :: piend
       type - integer
       desc - previous contig interval end
    :: splitlen
       type - integer
       desc - maximum split length
    '''

    # Setup endpoint Queue
    epq = cx.deque()
    end = start + splitlen

    # Absorb Potential Contigs
    for p,q in varcont:

        # Contig Left of/at (Start, PrevContigEnd)
        if p <= start and q <= pcend:
            continue

        # Contig at/before End
        if p <= end:

            # Contig before End
            if q <= end:
                epq.appendleft((
                    max((p, start, pcend)),
                    min(q, end)))
                continue

            # Contig encompasses End
            if end <= q:
                epq.appendleft((
                    max((p, start, pcend)),
                    min(q, end)))
                break # We're done
        
        # Contig beyond End
        if p > end:
            break # We're done!

    # Return Results
    return epq

def is_splittable(varcont, splitlen, minhdist):
    '''
    Determine if the variable contigs are
    within splitlen range, otherwise there
    is no solution to problem instance.
    Internal use only.

    :: varcont
       type - cx.deque
       desc - all variable region span indices
    :: splitlen
       type - integer
       desc - maximum split length
    :: minhdist
       type - integer
       desc - minimum pairwise hamming distance
    '''

    ci = iter(varcont)
    qq = next(ci)[1]
    for p,q in varcont:
        if p - qq + 2*minhdist > splitlen:
            return False
        qq = q
    return True

def is_last_fragment(start, splitlen, seqlen):
    '''
    Determine if the current split is the last
    frament in split. Internal use only.

    :: start
       type - integer
       desc - current split starting point
    :: splitlen
       type - integer
       desc - maximum split length
    :: seqlen
       type - integer
       desc - complete oligo length
    '''
    if start+splitlen >= seqlen:
        return True
    return False

def is_spannable(p, q, minhdist):
    '''
    Determine if a contig span satisfies minhdist.
    Internal use only.

    :: p
       type - integer
       desc - span start
    :: q
       type - integer
       desc - span end
    :: minhdist
       type - integer
       desc - minimum pairwise hamming distance
    '''
    if q - p < minhdist:
        return False
    return True

def split_engine(
    seqlist,
    splitlen=170,
    mintm=50,
    minhdist=10):

    liner = ut.liner_engine()
    
    seqmat = get_seqmat(seqlist=list(seqlist))
    seqlen = seqmat.shape[1]
    print(seqmat)
    print(seqmat.shape)

    # If sequences are shorter or equal to
    # splitlen, then no splitting required
    if seqlen < splitlen:
        liner.send(' Splitting unnecessary, all sequences synthesizable\n')
        return None

    entvec  = get_entvec(
        seqmat=seqmat)
    print(entvec)
    
    varcont = get_merged_varcont(
        varcont=get_varcont(
            entvec=entvec),
        mergefactor=minhdist // 2)
    print(varcont)

    # Note: Need to check for splittability before processing
    #       as well as when a split region is selected

    # This check is necessary upon ending
    # # If any of the varcont intervals are
    # # smaller than minhdist, then there is
    # # no solution to the present splitting
    # for p,q in varcont:
    #     if q-p < minhdist:
    #         liner.send(
    #             ' Variable ({}, {}) VR, smaller than Dmin\n'.format(
    #                 p,q))
    #         return None

    # If there are two regions are separated
    # by more than splitlen, then there is
    # no solution to the present splitting
    pp, qq = None, None
    for p,q in varcont:
        if not pp is None:
            if p - qq + 2*minhdist > splitlen:
                liner.send(
                    ' Too large separation ({}, {}) and ({}, {}) VR\n'.format(
                        pp, qq, p, q))
                return None
        qq = q
        pp = p

    runmat = get_runmat(
        seqmat=seqmat)
    # print(runmat)

    # Build all Fragments
    split = [] # Store all Splits
    start = 0  # Current  Frag Start
    pcend = -1 # Previous Frag End
    
    # Main Loop
    while True:

        # Is this the last fragment?
        if is_last_fragment(
            start=start,
            splitlen=splitlen,
            seqlen=seqlen):
            break

        # Get Endpoints for Current Fragment
        epq = get_endpoints(
            varcont=varcont,
            start=start,
            pcend=pcend,
            splitlen=splitlen)
        print(epq)
        
        # We don't have any contigs left
        # to build further fragments
        if not epq: break

        # Get the Tm based split
        r, q = get_split(
            seqmat=seqmat,
            runmat=runmat,
            epq=epq,
            mintm=mintm,
            minhdist=minhdist)
        print((r, q))

        # No Tm based split possible
        if r is None:
            break

        break           


    return

    # Build all Fragments
    start, end = 0, splitlen
    while True:

        # Where is the breakpoint given our ending?
        bpstart, bpend = alphavec[end], betavec[end]

        # Adjust ending
        bpend = min(end, bpend)

        # Is the varible region long enough for our Hamming needs?
        if bpend-bpstart < minhdist:
            liner.send(' Variable Region ({}, {}) too small\n'.format(
                bpstart, bpend))

        # If current end is same as previous end
        # no more variable region on right to
        # extend, no viable solution

        print(bpstart, bpend)

        # Select span based on Tm
        # (a full sub-engine)
        tmstart = get_split(
            mintm=mintm,
            runmat=runmat,
            bpstart=bpstart,
            bpend=bpend,
            minspan=minhdist)

        # If span overrides boundary on left
        # there is no viable solution

        # Extend span to meet minhammdist
        # (a full sub-engine)

        # If span overrides boundary on left
        # there is no viable solution

        # Adjust end to build next fragment set

        # Breaking condition execute
        break

# @nb.njit
def get_hdist(
    seqmat,
    i,
    j,
    idx):
    '''
    '''

    # We have store to compare
    if idx:
        return (seqmat[:idx, i:j] != seqmat[idx, i:j]).sum(1).min()

    # Nothing to compare
    return float('inf')

def get_split(seqmat, runmat, epq, mintm, minhdist):
    '''
    Return an integer r (p < r < q) from epq intervals such
    that Tm(runmat[r:q]) > mintm and HD(seqmat) > minhdist.
    Internal use only.
    
    :: seqmat
       type - np.array
       desc - numeric sequence matrix
    :: runmat
       type - np.array
       desc - numeric GC run matrix
    :: epq
       type - deque
       desc - endpoint variable contigs
    :: mintm
       type - float
       desc - melting temperature lower bound
    :: minhdist
       type - integer
       desc - minimum pairwise hamming distance
    '''
    
    # Candidate Results
    r, q = None, None

    # Tm computation loop
    while epq:

        # Fetch Current Interval
        p, q = epq.popleft()
        
        # Span Smaller than minhdist
        if not is_spannable(p, q, minhdist):
            r = None
            continue

        # Maximum r Value
        r = q - minhdist

        # Constraint Match Tracers
        condtm = False
        condhd = False

        # Iterate and Adjust
        idx = 0
        while idx < runmat.shape[0]:

            # Optimize Tm
            if not condtm:

                # Compute Tm for current split
                tm = get_approx_Tm(rv=runmat[idx, :], i=r, j=q)

                print(idx, p, r, q, tm, 'TM', tm >= mintm)

                # Tm was Lower ..
                if tm < mintm:

                    # Minimize r Value
                    r = r - 1

                    # Unresolvable
                    if r < p:
                        # No solution to current split
                        r = None
                        break
                    
                    # Try again ..
                    else:
                        continue
                
                # Tm was OK!
                elif tm >= mintm:
                    condtm = True

            # Optimize Hamming Distance
            if not condhd:

                # Compute Hamming Distance for current split
                hd = get_hdist(seqmat=seqmat, i=r, j=q, idx=idx)

                print(idx, p, r, q, hd, 'HD', hd >= minhdist)

                # HDist was Lower ..
                if hd < minhdist:

                    # Minimize r Value
                    r = r - minhdist + hd
                    
                    # Unresolvable
                    if r < p:
                        # No solution to current split
                        r = None
                        break
                    
                    # Try again ..
                    else:
                        continue

                # HDist was OK!
                elif hd >= minhdist:
                    condhd = True

            # Both Conditions Met!
            if condtm and condhd:
                # Move to next sequence
                idx += 1
                condtm = False # Reset Tm    Tracer
                condhd = False # Reset HDist Tracer
                continue

        # Do we have a solution?
        if not r is None:
            return r, q
    
    # Return Results
    return r, q


def get_Tm_start(mintm, runvec, i, j, minspan):
    maxstart = i
    p = j-minspan
    while p >= maxstart:
        tm = get_approx_Tm(
            rv=runvec,
            i=p,
            j=j)
        print(tm, runvec[p:j], p, j)
        if tm > mintm:
            return p
        p -= 1
    return None


def test1():
    entvec = cx.deque(map(
        int,'00000000000111111111000000011011111111110000000011111111110000'))
    varcont = get_merged_varcont(
        varcont=get_varcont(
            entvec=entvec),
        mergefactor=2)#minhdist // 2)

    print()
    print(varcont)
    print(get_endpoints(varcont, 0, 0, 20))
    print(get_endpoints(varcont, 0, 0, 24))
    print(get_endpoints(varcont, 0, 0, 34))
    print(get_endpoints(varcont, 0, 0, 45))
    print(get_endpoints(varcont, 0, 0, 10))
    print(get_endpoints(varcont, 11, 20, 15))
    print(get_endpoints(varcont, 27, 35, 22))
    print(get_endpoints(varcont, 14, 20, 32))

    entvec = cx.deque(map(
        #               |                     |                             
        #                  |                  |                             
        #                                      |            |
        int,'00000000000111111111111111111111111111111111111110000000000000'))
    varcont = get_merged_varcont(
        varcont=get_varcont(
            entvec=entvec),
        mergefactor=2)#minhdist // 2)
    print()
    print(varcont)
    print(get_endpoints(varcont, 0, 0, 34))
    print(get_endpoints(varcont, 14, 34, 34))
    print(get_endpoints(varcont, 24, 32, 34))
    print(is_last_fragment(24, 38, 62))

    print()
    print(is_splittable([(20, 30), (60, 70)], 30, 0))
    print(is_splittable([(20, 30), (60, 70)], 30, 5))

def get_DNA(l):
    fld = ['A', 'T', 'G', 'C']
    return ''.join(np.random.choice(fld) for _ in range(l))

def main():
    #       Constant                                      Constant
    #       ------------------                  ------------------
    seq1 = 'CCATAGTCAGACGCATCGAGAGAAGGCTGAGAATGAAATCTGCGCATATCGACG'
    seq2 = 'CCATAGTCAGACGCATCGCCGACTCCAATCCTAGACAATCTGCGCATATCGACG'
    seq3 = 'CCATAGTCAGACGCATCGGATCTGATTCCTAGAACTAATCTGCGCATATCGACG'
    seq4 = 'CCATAGTCAGACGCATCGTTCGGCGAGGCGTCACTGAATCTGCGCATATCGACG'
    #       ----------------------++++++++++
    #                             ++++++++++----------------------

    seqlist = [seq1, seq2, seq3, seq4]
    split_engine(
        seqlist=seqlist,
        splitlen=32,
        mintm=5,
        minhdist=10)

    # seq1 = 'TGATTCCTAG'
    # seqlist = [seq1]
    # runmat = get_runmat(
    #     seqmat=get_seqmat(
    #         seqlist=seqlist))
    # print(runmat)
    # print(get_approx_Tm(runmat[0, :], 0, len(seq1)))

    #                        40
    #      ----------------------------------------
    #                          ||||||||||||||||||||
    #                          --------------------
    seq = 'ACGAGAGTAGGCTGAGTGCTAATTTAATCTGCTTTTTAAAAGCGGACTCCCATTAGACGCXX' # 60
    #                            --------------------
    #                            ||||||||||||||||||||
    #                            ----------------------------------------
    #                                                40

    #                        20
    #      --------------------
    seq = 'ACGAGAGTAGGCTGAGTGCTAATTTAATCTGCTTTTTAAAAGCGGACTCCCATTAGACGC' # 60
    #                                              --------------------
    #                                              20

if __name__ == '__main__':
    main()