import time  as tt

import collections as cx

import numpy as np
import numba as nb

import utils as ut

def get_spanlen(minoverlap, minhdist):
    '''
    Return the minimum required split span length.
    Internal use only.

    :: minoverlap
       type - integer
       desc - minimum fragment overlap length
    :: minhdist
       type - integer
       desc - minimum pairwise hamming distance
    '''

    # Automatically fulfills minoverlap
    return max(minoverlap, minhdist)

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

def get_seqmat(seqlist, seqlen, liner):
    '''
    Return the numeric representation
    of seqlist. Internal use only.

    :: seqlist
       type - iterable
       desc - list of sequences to split
    :: seqlen
       type - integer
       desc - length of sequences to split
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Show Segment
    liner.send('\n[Building Sequence Matrix]\n')
    
    # Setup Store
    seqmat = np.zeros(
        (len(seqlist), seqlen),
        dtype=np.float64)

    # Update-keeping
    short = True if seqlen > 10 else False

    # Time-keeping
    t0 = tt.time()

    # Fill Store
    for idx,seq in enumerate(seqlist):
        seqvec = get_seqvec(seq=seq)
        seqmat[idx, :] = seqvec

        if short:
            liner.send(' Storing Vectorized Sequence {}: {}..{}'.format(
                idx, seqvec[:5], seqvec[-5:]))
        else:
            liner.send(' Storing Vectorized Sequence {}: {}'.format(
                idx, seqvec))

    # Final Updates
    liner.send('   Vectorized: {} Sequences\n'.format(idx+1))
    liner.send(' Time Elapsed: {:.2f} sec\n'.format(tt.time()-t0))
    
    # Return Results
    return np.array(seqmat, dtype=np.float64)

def get_entvec(seqmat, seqlen, liner):
    '''
    Return entropy of sequence matrix.
    Internal use only.

    :: seqmat
       type - np.array
       desc - numeric sequence matrix
    :: seqlen
       type - integer
       desc - length of sequences to split
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Show Segment
    liner.send('\n[Building Entropy Vector]\n')

    # Update-keeping
    short = True if seqlen > 10 else False

    # Time-keeping
    t0 = tt.time()

    # Show Update
    liner.send(' Computing Count Vector ...')
    
    # Setup Data Structures
    d = np.zeros(
        (4, seqmat.shape[1]),
        dtype=np.float64)
    p = np.zeros(4)
    
    # Compute Count Frequency
    for idx in range(seqmat.shape[1]):
        
        # Count Symbol Counts at i-th Index
        m = np.unique(
            seqmat[:, idx],
            return_counts=True)[1]
        
        # Update Data Structure and Show Updates
        p[:m.shape[0]] = m # Absorb Local

        liner.send(' Index {} Unique Count: {}'.format(
            idx, p))

        d[:, idx] += p     # Update Global
        p *= 0.            # Reset  Local

    # Normalize Counts
    d = d / d.sum(0)

    # Show Updates
    liner.send('   Count Vector Normalized in {:.2f} sec\n'.format(
        tt.time() - t0))

    # Reset Time-keeping
    t0 = tt.time()

    # Show Updates
    liner.send(' Computing Entropy Vector ...')
    
    # Compute Entropy Vector
    entvec = cx.deque(np.abs(
        (d*(np.log(d, where=d>0.) / np.log(4))).sum(0)))
    
    # Show Updates
    liner.send(' Entropy Vector Computed   in {:.2f} sec\n'.format(
        tt.time() - t0))

    # Return Results
    return entvec

def get_base_varcont(entvec):
    '''
    Return all base variable region span indices.
    Internal use only.

    :: entvec
       type - cx.deque
       desc - positional entropy
    '''

    # Setup Parsing
    varcont = cx.deque()
    start   = None
    end     = None
    idx     = -1

    # Time-keeping
    t0 = tt.time()

    # Parse Contigs
    while entvec:

        # Update index
        idx += 1

        # Extract entropy
        ent = entvec.popleft()

        # Constant Region
        if ent <= 0.25: # Constant Upper Bound

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
                start = idx
            if end is None:
                end = idx
            
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

    # Do we merge?
    if not varcont:
        return varcont

    # Setup Data Structures
    merged   = cx.deque()
    previous = list(varcont.popleft())

    # Time-keeping
    t0 = tt.time()
    
    # Merge Contigs
    while varcont:
        current = varcont.popleft()
        if current[0]-previous[1] <= mergefactor:
            previous[1] = current[1]
        else:
            merged.append(tuple(previous))
            previous = list(current)
    merged.append(tuple(previous))
    
    # Return Result
    return merged

def is_spannable(p, q, spanlen):
    '''
    Determine if a contig span satisfies spanlen.
    Internal use only.

    :: p
       type - integer
       desc - span start
    :: q
       type - integer
       desc - span end
    :: spanlen
       type - integer
       desc - minimum required split span length
    '''

    if q - p < spanlen:
        return False
    return True

def get_filtered_varcont(varcont, spanlen):
    '''
    Filter all variable regions shorter than
    the required spanlen. Internal use only.

    :: varcont
       type - cx.deque
       desc - a deque of tuple with start and
              end coordinates of variable regions
    :: spanlen
       type - integer
       desc - minimum required split span length
    '''

    # Do we filter varcont?
    if not varcont:
        return varcont
    
    # Setup Data Structure
    filtered = cx.deque()

    # Time-keeping
    t0 = tt.time()

    # Filter Contigs
    for p,q in varcont:
        if is_spannable(p, q, spanlen):
            filtered.append((p, q))

    # Return Result
    return filtered

def get_varcont(entvec, minhdist, spanlen, liner):

    # Show Segment
    liner.send('\n[Computing Variable Contigs]\n')

    # Extract varcont
    liner.send(' Extracting Variable Contigs ...')
    t0 = tt.time()
    varcont = get_base_varcont(
        entvec=entvec)
    liner.send(' Extracted Variable Contigs: {} ({:.2f} sec)\n'.format(
        len(varcont), tt.time()-t0))

    # Merge varcont
    liner.send('   Merging Variable Contigs ...')
    t0 = tt.time()
    varcont = get_merged_varcont(
        varcont=varcont,
        mergefactor=minhdist // 2)
    liner.send('    Merged Variable Contigs: {} ({:.2f} sec)\n'.format(
        len(varcont), tt.time()-t0))

    # Filter varcont
    liner.send(' Filtering Variable Contigs ...')
    t0 = tt.time()
    varcont = get_filtered_varcont(
        varcont=varcont,
        spanlen=spanlen)
    liner.send('  Filtered Variable Contigs: {} ({:.2f} sec)\n'.format(
        len(varcont), tt.time()-t0))

    # Return Results
    return varcont

def is_varcont_feasible(varcont, seqlen, splitlen, spanlen, liner):
    '''
    Determine if the variable contigs are
    within splitlen range, otherwise there
    is no solution to problem instance.
    Internal use only.

    :: varcont
       type - cx.deque
       desc - all variable region span indices
    :: seqlen
       type - integer
       desc - length of sequences to split
    :: splitlen
       type - integer
       desc - maximum split length
    :: spanlen
       type - integer
       desc - minimum required split span length
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    ci     = iter(varcont)
    pp, qq = next(ci)

    # Check for the first fragment
    if pp + spanlen > splitlen:
        liner.send(
            ' Verdict: Infeasible (First Contig (Start={}, End={}) too Far from Beginning)\n'.format(
                pp, qq))
        return False

    # Check between contigs
    for p,q in varcont:
        if p - qq + 2*spanlen > splitlen:
            liner.send(
                ' Verdict: Infeasible (Adjacent Contigs (Start={}, End={}) and (Start={}, End={}) Far Apart)\n'.format(
                    pp, qq, p, q))
            return False
        qq = q
        pp = p

    # Check for the last fragment
    if qq - spanlen + splitlen < seqlen:
        liner.send(
            ' Verdict: Infeasible (Last Contig (Start={}, End={}) too far from Ending)\n'.format(
                pp, qq))
        return False

    # No problems found
    return True

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
    :: pcend
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

def is_last_fragment(pcstart, splitlen, seqlen):
    '''
    Determine if the current split is the last
    frament in split. Internal use only.

    :: pcstart
       type - integer
       desc - previous contig interval start
    :: splitlen
       type - integer
       desc - maximum split length
    :: seqlen
       type - integer
       desc - complete oligo length
    '''
    if pcstart+splitlen >= seqlen:
        return True
    return False

def get_split(seqlist, seqmat, epq, mintmelt, minhdist):
    '''
    Return an integer r (p < r < q) from epq intervals such
    that Tm(runmat[r:q]) > mintmelt and HD(seqmat) > minhdist.
    Internal use only.
    
    :: seqlist
       type - list
       desc - list of sequences to split
    :: seqmat
       type - np.array
       desc - numeric sequence matrix
    :: epq
       type - deque
       desc - endpoint variable contigs
    :: mintmelt
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

        # Unresolvable
        if r < p:
            # No solution to current split
            r = None
            break

        # Constraint Match Tracers
        condtm = False
        condhd = False

        # Iterate and Adjust
        idx = 0
        while idx < seqmat.shape[0]:

            # Optimize Tm
            if not condtm:

                # Compute Tm for current split
                tmelt = ut.get_tmelt(
                    seq=seqlist[idx],
                    i=r,
                    j=q)

                print(idx, p, r, q, tmelt, 'TM', tmelt >= mintmelt)

                # Tm was Lower ..
                if tmelt < mintmelt:

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
                elif tmelt >= mintmelt:
                    condtm = True

            # Optimize Hamming Distance
            if not condhd:

                # Compute Hamming Distance for current split
                # hd = get_hdist(seqmat=seqmat, i=r, j=q, idx=idx)
                hdist = ut.get_hdist(
                    store=seqmat,
                    idx=idx,
                    i=r,
                    j=q,
                    direction=0)

                print(idx, p, r, q, hdist, 'HD', hdist >= minhdist)

                # HDist was Lower ..
                if hdist < minhdist:

                    # Minimize r Value
                    r = r - minhdist + hdist
                    
                    # Unresolvable
                    if r < p:
                        # No solution to current split
                        r = None
                        break
                    
                    # Try again ..
                    else:
                        continue

                # HDist was OK!
                elif hdist >= minhdist:
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

def split_engine(
    seqlist,
    splitlen=170,
    mintmelt=50,
    minhdist=10,
    minoverlap=20,
    liner=None):
    '''
    Return a list of split coordinates.

    :: seqlist
       type - list
       desc - list of sequences to split
    :: splitlen
       type - integer
       desc - maximum oligo length for splitting
    :: mintmelt
       type - float
       desc - minimum melting temperature of
              split regions
    :: minhdist
       type - integer
       desc - minimum pairwise hamming distance
              between all regions at a split
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Absorb seqlist
    seqlist = list(seqlist)
    seqlen  = len(seqlist[0])

    # spanlen is the minimum span of a split
    liner.send('\n[Computing Split Region Span Length]\n')
    spanlen = get_spanlen(
        minoverlap=minoverlap,
        minhdist=minhdist)
    liner.send(' Split Span Length: {} bp\n'.format(spanlen))

    # Note to future self: seqlen perhaps a parameter to engine
    # when driver checks for seqlen equality in input (must).

    # Splitting Feasibility Check: Necessity
    liner.send('\n[Checking Splitting Feasibility]\n')
    if seqlen <= splitlen: # Splitting unnecessary
        liner.send(' Verdict: Unnecessary (seqlen <= splitlen)\n')
        return None # No Solution
    else:
        liner.send(
            ' Verdict: At Least {} Fragments per Sequence\n'.format(
                int(np.ceil(seqlen / (splitlen * 1.))))) # Fragment count lowerbound

    # Build seqmat
    seqmat = get_seqmat(
        seqlist=seqlist,
        seqlen=seqlen,
        liner=liner)

    # Build entvec
    entvec = get_entvec(
        seqmat=seqmat,
        seqlen=seqlen,
        liner=liner)

    # Build varcont
    varcont = get_varcont(
        entvec=entvec,
        minhdist=minhdist,
        spanlen=spanlen,
        liner=liner)

    # Splitting Feasibility Check: >= 1 Variable Contigs
    liner.send('\n[Checking Splitting Feasibility]\n')
    if len(varcont) == 0: # No region to split
        liner.send(' Verdict: Infeasible (Sequences in Pool are Similar)\n')
        return None # No Solution
    if not is_varcont_feasible(
        varcont=varcont,
        seqlen=seqlen,
        splitlen=splitlen,
        spanlen=spanlen,
        liner=liner):
        return None # No Solution
    liner.send(' Verdict: Feasible Contigs Found\n')
    
    return

    # Note: Need to check for splittability before processing
    #       as well as when a split region is selected

    # Build all Fragments
    split = [] # Store all Splits
    start = 0  # Current  Frag Start
    pcend = 0  # Previous Frag End

    unsoln = False
    
    # Main Loop
    while True:

        # Get Endpoints for Current Fragment
        epq = get_endpoints(
            varcont=varcont,
            start=start,
            pcend=pcend,
            splitlen=splitlen)
        print(epq)
        
        # We don't have any contigs left
        # to build further fragments
        if not epq:
            unsoln = True
            break

        # Get the Tm and HDist based split
        r, q = get_split(
            seqlist=seqlist,
            seqmat=seqmat,
            epq=epq,
            mintmelt=mintmelt,
            minhdist=minhdist)
        print((r, q))

        # No Tm based split possible
        if r is None:
            unsoln = True
            break

        # Prepare for next split
        else:
            split.append((r, q))
            pcstart = r
            pcend   = q
            start   = q + 1

        # Is this the last fragment?
        if is_last_fragment(
            pcstart=pcstart,
            splitlen=splitlen,
            seqlen=seqlen):
            print('LF!')
            # split.append((pc, seqlen))
            break

    if unsoln:
        return None
    else:
        return split

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

    liner = ut.liner_engine()
    seqlist = [seq1, seq2, seq3, seq4]
    split = split_engine(
        seqlist=seqlist,
        splitlen=32,#26,
        mintmelt=10,#10,
        minhdist=9,#5,
        minoverlap=5,
        liner=liner)
    
    print()
    print(split)

    print(seq1[28:36])

if __name__ == '__main__':
    main()