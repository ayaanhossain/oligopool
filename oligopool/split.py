import time  as tt

import collections as cx

import numpy as np
import numba as nb

import utils as ut

def get_spanlen(
    minoverlap,
    minhdist):
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

def get_seqvec(
    seq):
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

def get_seqmat(
    seqlist,
    seqlen,
    liner):
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

def get_entvec(
    seqmat,
    seqlen,
    liner):
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

def get_base_varcont(
    entvec):
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

def get_merged_varcont(
    varcont,
    mergegap):
    '''
    Return a merged varcont, merging variable regions
    separated by at most mergegap constant bases.
    Internal use only.

    :: varcont
       type - cx.deque
       desc - all variable region span indices
    :: mergegap
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
        if current[0]-previous[1] <= mergegap:
            previous[1] = current[1]
        else:
            merged.append(tuple(previous))
            previous = list(current)
    merged.append(tuple(previous))
    
    # Return Result
    return merged

def is_spannable(
    p,
    q,
    spanlen):
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

def get_filtered_varcont(
    varcont,
    spanlen):
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

def get_varcont(
    entvec,
    minhdist,
    spanlen,
    liner):
    '''
    Return all valid variable contig span regions.
    Internal use only.

    :: entvec
       type - cx.deque
       desc - positional entropy
    :: minhdist
       type - integer
       desc - minimum pairwise hamming distance
    :: spanlen
       type - integer
       desc - minimum required split span length
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Show Segment
    liner.send('\n[Computing Variable Contigs]\n')

    # Extract varcont
    liner.send(' Extracting Variable Contigs ...')
    t0 = tt.time()
    varcont = get_base_varcont(
        entvec=entvec)
    liner.send(
        ' Extracted Variable Contigs: {} ({:.2f} sec)\n'.format(
            len(varcont), tt.time()-t0))

    # Merge varcont
    liner.send('   Merging Variable Contigs ...')
    t0 = tt.time()
    varcont = get_merged_varcont(
        varcont=varcont,
        mergegap=min(3, minhdist // 4))
    liner.send(
        '    Merged Variable Contigs: {} ({:.2f} sec)\n'.format(
            len(varcont), tt.time()-t0))

    # Filter varcont
    liner.send(' Filtering Variable Contigs ...')
    t0 = tt.time()
    varcont = get_filtered_varcont(
        varcont=varcont,
        spanlen=spanlen)
    liner.send(
        '  Filtered Variable Contigs: {} ({:.2f} sec)\n'.format(
            len(varcont), tt.time()-t0))

    # Return Results
    return varcont

def is_varcont_feasible(
    varcont,
    seqlen,
    splitlen,
    spanlen,
    liner):
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

def get_splitend(
    fstart,
    splitlen,
    seqlen):
    '''
    Return the fragment end coordinates.
    Internal use only.

    :: fstart
       type - integer
       desc - current fragment starting point
    :: splitlen
       type - integer
       desc - maximum split length
    :: seqlen
       type - integer
       desc - complete oligo length
    '''

    return min(fstart+splitlen, seqlen)

def get_splitqueue(
    varcont,
    fstart,
    sstart,
    splitlen,
    seqlen):
    '''
    Return all variable contigs splitpoints given
    current oligo starting point and splitlen.
    Internal use only.

    :: varcont
       type - cx.deque
       desc - all variable region span indices
    :: fstart
       type - integer
       desc - current fragment starting point
    :: sstart
       type - integer
       desc - current split starting point
    :: splitlen
       type - integer
       desc - maximum split length
    :: seqlen
       type - integer
       desc - complete oligo length
    '''

    # Do we compute splitpoints?
    if not varcont:
        return varcont

    # Setup Endpoint Queue
    spq = cx.deque()
    end = get_splitend(
        fstart=fstart,
        splitlen=splitlen,
        seqlen=seqlen)

    # Determine Potential Splitpoints from Contigs
    for p,q in varcont:

        # Skip Condition:
        # Contig Ends before Start
        if q <= sstart:
            continue # Skip!

        # Absorb Condition:
        # Contig Starts before End
        if p < end:
            
            # Add Splitpoint
            spq.appendleft((
                max(p, sstart),
                min(q, end)))

            # Contig before End
            if q <  end:
                continue # Next!

            # Contig encompasses End
            if end <= q:
                break # We're done!
        
        # Exhaustion Condition:
        # Contig Starts after End
        if p >= end:
            break # We're done!

    # Return Results
    return spq

def continue_splitting(
    fstart,
    splitlen,
    seqlen):
    '''
    Determine if the current split engenders the last
    frament from splitting. Internal use only.

    :: fstart
       type - integer
       desc - current fragment starting point
    :: splitlen
       type - integer
       desc - maximum split length
    :: seqlen
       type - integer
       desc - complete oligo length
    '''

    # Starting at fstart, we'll
    # reach end of sequence
    if get_splitend(
        fstart=fstart,
        splitlen=splitlen,
        seqlen=seqlen) == seqlen:
        return False # No more splitting required
    
    return True # Splitting may be required

def get_split(
    seqlist,
    seqmat,
    spq,
    mintmelt,
    minhdist,
    spanlen,
    liner):
    '''
    Return an integer r (p < r < q) from spq intervals such
    that Tm(runmat[r:q]) > mintmelt and HD(seqmat) > minhdist.
    Internal use only.
    
    :: seqlist
       type - list
       desc - list of sequences to split
    :: seqmat
       type - np.array
       desc - numeric sequence matrix
    :: spq
       type - deque
       desc - endpoint variable contigs
    :: mintmelt
       type - float
       desc - melting temperature lower bound
    :: minhdist
       type - integer
       desc - minimum pairwise hamming distance
    :: spanlen
       type - integer
       desc - minimum required split span length
    :: liner
       type - coroutine
       desc - dynamic printing
    '''
    
    # Book-keeping
    r, q  = None, None # Split Coordinates
    state = False      # Solution State 

    # Splitting Loop
    while spq:

        # Fetch Current Splitpoints
        p, q = spq.popleft()

        # Show Update
        liner.send('    Attempting Split in Region: (Start={}, End={})\n'.format(
            p, q))
        
        # Is current splitpoint feasigle?
        if not is_spannable(p, q, spanlen):
            liner.send('      Current Split Region of {} bp Infeasible ... Skipping\n'.format(
                q-p))
            state = False # No Solution (Yet)
            continue
        else:
            liner.send('      Current Split Region of {} bp maybe Feasible ... Optimizing\n'.format(
                q-p))

        # Adjusted r Value
        r = q - spanlen # Takes care of spanlen

        # Constraint Match Tracers
        condtm = False
        condhd = False

        # Split Adjustment Loop
        idx = 0 # Current Sequence for Verification
        while idx < seqmat.shape[0]: # Adjust/Verify wrt all Sequences

            # tt.sleep(1)

            # Optimize Tm
            if not condtm:

                # Compute Tm for current split
                tmelt = ut.get_tmelt(
                    seq=seqlist[idx],
                    i=r,
                    j=q)

                # Tm was Lower ..
                if tmelt < mintmelt:

                    # Show Update
                    liner.send(
                        '      Sequence {}: Region (Start={}, End={}) is NOT Tm ({:.2f} C) Optimal'.format(
                            idx, r, q, tmelt))

                    # Minimize r Value
                    r = r - 1
                
                # Tm was OK!
                elif tmelt >= mintmelt:
                    condtm = True

                    # Show Update
                    liner.send(
                        '      Sequence {}: Region (Start={}, End={}) is Tm ({:.2f} C) Optimal'.format(
                            idx, r, q, tmelt))

            # Optimize Hamming Distance (after Tm)
            if condtm and not condhd:

                # Compute Hamming Distance for current split
                hdist = ut.get_hdist(
                    store=seqmat,
                    idx=idx,
                    i=r,
                    j=q,
                    direction=0)

                # HDist was Lower ..
                if hdist < minhdist:

                    # Show Update
                    liner.send(
                        '      Sequence {}: Region (Start={}, End={}) is NOT HDist ({}) Optimal'.format(
                            idx, r, q, hdist))

                    # Minimize r Value
                    r = r - minhdist + hdist

                # HDist was OK!
                elif hdist >= minhdist:
                    condhd = True

                    # Show Update
                    liner.send(
                        '      Sequence {}: Region (Start={}, End={}) is HDist ({}) Optimal'.format(
                            idx, r, q, hdist))

            # Both Conditions Met!
            if condtm and condhd:
                idx   += 1     # Move to next sequence
                condtm = False # Reset Tm    Tracer
                condhd = False # Reset HDist Tracer

                # Show Update
                liner.send(
                    '      Sequence {}: Region (Start={}, End={}) is Tm ({:.2f} C) and HDist ({}) Optimal'.format(
                        idx-1, r, q, tmelt, hdist))

                # Analyze Next Sequence!
                continue

            else:
                # Unresolvable
                if r < p:
                    liner.send('\*      Current Split Region Infeasible for Sequence {} ... Skipping\n'.format(
                            idx))
                    r = None # No solution to current split
                    break    # Try Next Split Region ..
                
                # Try again ..
                else:
                    continue

        # Do we have a solution?
        if not r is None:
            liner.send('\*      Current Split Region Optimized for All {} Sequences\n'.format(
                    idx))
            state  = True  # Solution Found!
            break # We're done!
    
    # Return Results
    if not state:
        return None # No Solution
    else:
        return r,q

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
    liner.send(' Split Region Minimum Span Length: {} bp\n'.format(spanlen))

    # Note to future self: seqlen perhaps a parameter to engine
    # when driver checks for seqlen equality in input (must).

    # Splitting Feasibility Check: Necessity
    liner.send('\n[Checking Splitting Feasibility]\n')
    if seqlen <= splitlen: # Splitting unnecessary
        liner.send(' Verdict: Unnecessary (seqlen <= splitlen)\n')
        return None # No Solution
    else:
        liner.send(
            ' Verdict: Minimum {} Fragments per Sequence\n'.format(
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

    # Note: Need to check for splittability before processing
    #       as well as when a split region is selected

    # Book-keeping
    split  = []    # Split Coordinate Storage
    fstart = 0     # Current Fragment  Start
    sstart = 0     # Current Split     Start
    state  = False # Solution State

    # Show Segment
    liner.send('\n[Computing Split Fragment Coordinates]\n')
    
    # Compute Split Coordinates
    mt0 = tt.time() # Total Split Time-keeping
    while True:

        # Show Update
        liner.send('\n  Now Computing Split for Fragment: {}\n'.format(
            len(split)+1))
        liner.send('    Initial Fragment Coordinates: (Start={}, End={})\n'.format(
            fstart,
            get_splitend(
                fstart=fstart,
                splitlen=splitlen,
                seqlen=seqlen)))

        # Do we need to split any more?
        if not continue_splitting(
            fstart=fstart,
            splitlen=splitlen,
            seqlen=seqlen):

            # Store Final Fragment Coordinates
            split.append((fstart, seqlen))

            # Show Updates
            liner.send('    Split Required? No\n')
            liner.send('    Final Fragment Coordinates: (Start={}, End={})\n'.format(
                    *split[-1]))
            
            # Book-keeping Update
            state = True # Problem Solved!
            break # No more splitting required

        else:

            # Show Updates
            liner.send('    Split Required? Yes\n')

            # Instance Time-keeping
            t0 = tt.time()

            # Get Splitpoints for Current Fragment
            liner.send('    Finding Splittable Regions ...')
            spq = get_splitqueue(
                varcont=varcont,
                fstart=fstart,
                sstart=sstart,
                splitlen=splitlen,
                seqlen=seqlen)
            
            # Did we find split regions?
            if not spq: # No Split Regions Found
                liner.send('    No Splittable Regions Found ... Terminating\n')
                state = False # No Solution
                break
            else:       # Split Regions Found
                liner.send('    Splittable Regions Found: {} (in {:.2f} sec)\n'.format(
                    len(spq),
                    tt.time() - t0))

            # Get the Tm and HDist based split
            rq = get_split(
                seqlist=seqlist,
                seqmat=seqmat,
                spq=spq,
                mintmelt=mintmelt,
                minhdist=minhdist,
                spanlen=spanlen,
                liner=liner)

            # Did we find feasible split regions?
            if rq is None: # No Feasible Split Found
                liner.send('    No Feasible Splits Found ... Terminating\n')
                state = False # No Solution
                break
            
            else:          # Feasigle Split Found
                r,q = rq   # Parse Split
                # Store Current Fragment Coordinates
                split.append((fstart, q))
                
                # Book-keeping Update
                fstart = r # Next Fragment Start Coordinate
                sstart = q # Next Split    Start Coordinate

                # Show Updates
                liner.send('    Split Region Selected: (Start={}, End={}) (in {:.2f} sec)\n'.format(
                    r, q, tt.time()-t0))
                liner.send('    Final Fragment Coordinates: (Start={}, End={})\n'.format(
                    *split[-1]))

    # Final Updates
    liner.send('\n  Time Elapsed: {:.2f} sec\n'.format(
        tt.time()-mt0))

    # Return Results
    if not state:
        return None # No Solution
    return split

def test1():
    liner = ut.liner_engine()
    entvec = cx.deque(map(
        int,'00000000000111111111000000011011111111110000000011111111110000'))
    varcont = get_varcont(
        entvec=entvec,
        minhdist=4,
        spanlen=9,
        liner=liner)

    print()
    print(varcont)
    print(get_splitqueue(varcont, 0, 0, 20))
    print(get_splitqueue(varcont, 0, 0, 24))
    print(get_splitqueue(varcont, 0, 0, 34))
    print(get_splitqueue(varcont, 0, 0, 45))
    print(get_splitqueue(varcont, 0, 0, 10))
    print(get_splitqueue(varcont, 11, 20, 30))
    print(get_splitqueue(varcont, 20, 19, 15))
    print(get_splitqueue(varcont, 35, 34, 22))
    print(get_splitqueue(varcont, 20, 19, 32))

    entvec = cx.deque(map(
        #               |                     |                             
        #                  |                  |                             
        #                                      |             |
        int,'00000000000111111111111111111111111111111111111110000000000000'))
           # 0000000000011111111111111111111111
                                             # 111111111111111
    varcont = get_varcont(
        entvec=entvec,
        minhdist=4,
        spanlen=9,
        liner=liner)
    
    print(varcont)
    print(get_splitqueue(varcont, 0,  0, 34))
    print(get_splitqueue(varcont, 34, 33, 34))
    print(get_splitqueue(varcont, 32, 31, 34))
    print(continue_splitting(24, 38, 62))

def get_DNA(l):
    fld = ['A', 'T', 'G', 'C']
    return ''.join(np.random.choice(fld) for _ in range(l))

def main():
    liner = ut.liner_engine()
    seqlist= []
    for i in range(2500):
        seq = get_DNA(300)
        if (i+1)%1000 == 0:
            print(i, seq[:20], '...')
        seqlist.append(seq[:150] + 'CCATAGTCAGACGCATCGAG' + seq[-150:])
    split = split_engine(
        seqlist=seqlist,
        splitlen=170,
        mintmelt=35,
        minhdist=25,
        minoverlap=25,
        liner=liner)
    print()
    print(split)
    return

    # test1()
    # return

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
        splitlen=27,
        mintmelt=0,
        minhdist=5,
        minoverlap=5,
        liner=liner)
    # split = split_engine(
    #     seqlist=seqlist,
    #     splitlen=32,#26,
    #     mintmelt=10,#10,
    #     minhdist=9,#5,
    #     minoverlap=5,
    #     liner=liner)
    
    print()
    print(split)

    print(seq1[28:36])

if __name__ == '__main__':
    main()