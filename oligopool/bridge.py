import time as tt

import collections as cx
import itertools   as ix

import numpy as np

import utils as ut

def enque_exmotifs(exmotifs, liner):
    '''
    Return all exmotifs sorted by
    length. Internal use only.

    :: exmotifs
       type - list
       desc - list of all motifs
              to be excluded
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
    dq = cx.deque(exmotifs)
    liner.send(' Enqued {} Motifs in {:.2f} sec\n'.format(
        len(exmotifs), tt.time()-t0))

    # Return Results
    return dq

def get_bll(exmotifs, liner):
    '''
    Return the bridge-mer length limit (BLL).
    Internal use only.

    :: exmotifs
       type - deque
       desc - list of all motifs
              to be excluded
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Minimum BLL (default)
    bll = 0

    # BLL is the length of the last element
    if exmotifs:
        maxel = exmotifs.pop()
        bll   = len(maxel) - 1
        exmotifs.append(maxel)

    # Show update
    liner.send(' Bridge Length Limit: {} bp\n'.format(
        bll))

    # Return Results
    return bll

def get_bases(rng):
    '''
    Return a random permutation of bases.
    Internal use only.

    :: rng
       type - np.default_rng
       desc - RNG object
    '''

    bases = ['A', 'T', 'G', 'C']
    rng.shuffle(bases)
    return bases

def partition_exmotifs(lt, ge, bl, liner):
    '''
    Partition the excluded motifs into shorter or
    longer queues, based on bridge-mer length.

    :: lt
       type - deque
       desc - queue of shorter than bl motifs
    :: ge
       type - deque
       desc - queue of longer or equal motifs
    :: bl
       type - integer
       desc - bridge-mer length
    :: liner
       type - coroutine
       desc - dynamic printing
    '''
    
    # Time-keeping
    t0 = tt.time()

    # Partition Loop
    while ge:

        # Get Motif
        motif = ge.popleft()
        
        # Shorter Queue
        if len(motif) < bl:
            liner.send(' {} in Shorter Queue'.format(
                motif))
            lt.append(motif)

        # Longer Queue
        else:
            liner.send(' {} in Longer  Queue'.format(
                motif))
            ge.appendleft(motif)
            break

    # Show Updates
    liner.send(' Shorter Queue Count: {} Motifs\n'.format(
        len(lt)))
    liner.send('  Longer Queue Count: {} Motifs\n'.format(
        len(ge)))
    liner.send('        Time Elapsed: {:.2f} sec\n'.format(
        tt.time()-t0))

    # Return Results
    return lt, ge

def is_lt_embedded(lt, bridge):
    '''
    Determind if candidate contains a smaller
    exmotif. Internal use only.

    :: lt
       type - deque
       desc - queue of shorter than bl motifs
    :: bridge
       type - string
       desc - candidate bridge
    '''
    
    for motif in lt:
        if motif in bridge:
            return True
    return False

def is_ge_embedded(ge, bridge):
    '''
    Determine if bridge contained in a longer
    exmotif. Internal use only.

    :: ge
       type - deque
       desc - queue of longer or equal motifs
    :: bridge
       type - string
       desc - bridge bridge
    '''
    
    for motif in ge:
        if bridge in motif:
            return True
    return False

def get_status(lt, ge, bridge):
    '''
    Return candidate bridge feasibility.
    Internal use only.

    :: lt
       type - deque
       desc - queue of shorter than bl motifs
    :: ge
       type - deque
       desc - queue of longer or equal motifs
    :: bridge
       type - string
       desc - candidate bridge
    '''

    # Is bridge grown in wrong path?
    ltf = is_lt_embedded(lt, bridge)
    if ltf:
        return -1 # Rejection

    # Is bridge potentially in correct path?
    gef = is_ge_embedded(ge, bridge)
    if gef:
        return 0  # Partial Acceptance

    # Candidate is completely accepted
    return 1 # Complete Acceptance

def get_bounded_bridges(lt, ge, bridges, liner):
    '''
    Return all bounded bridge-mers. Internal
    use only.

    :: lt
       type - deque
       desc - queue of shorter than bl motifs
    :: ge
       type - deque
       desc - queue of longer or equal motifs
    :: bridges
       type - list
       desc - list of candidate bridges
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    partial  = [] # Potential for growth
    complete = [] # Completely grown

    # Time-keeping
    t0 = tt.time()

    # Bounding Loop
    for bridge in bridges:

        # Show Update
        liner.send(' Evaluating Bridge {} ...'.format(
            bridge))

        # Is bridge feasible?
        status = get_status(
            lt=lt,
            ge=ge,
            bridge=bridge)

        # Candidate Rejected
        if status == -1:
            liner.send(' Bridge {} Rejected'.format(bridge))

        # Candidate Partially Accepted
        if status == 0:
            liner.send(' Bridge {} Paritally Accepted'.format(bridge))
            partial.append(bridge)
        
        # Candidate Accepted
        if status == 1:
            liner.send(' Bridge {} Accepted'.format(bridge))
            complete.append(bridge)

    # Show Updates
    liner.send('  Partial Candidates: {} Bridge(s) ({:6.2f} %)\n'.format(
        len(partial),  (100. * len(partial)) / len(bridges)))
    liner.send(' Complete Candidates: {} Bridge(s) ({:6.2f} %)\n'.format(
        len(complete), (100. * len(complete)) / len(bridges)))
    liner.send('        Time Elapsed: {:.2f} sec\n'.format(
        tt.time()-t0))

    # Return Results
    return partial, complete

def get_branched_bridges(rng, partial, liner):
    '''
    Return all branched candidate bridges from
    partial solutions. Internal use only.

    :: rng
       type - np.default_rng
       desc - RNG object
    :: partial
       type - list
       desc - list of all partially accepted
              bridge candidates
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    branched = [] # Explorable

    # Time-keeping
    t0 = tt.time()

    # Branching Loop
    for bridge in partial:
        
        # Generate Branches
        for base in get_bases(rng):
            branch = bridge+base    # Next Ppossible Bridge
            branched.append(branch) # Store

            # Show Updates
            liner.send(' Branched Bridge {} -> {}'.format(
                bridge, branch))

    # Show Updates
    liner.send(' Branched Candidates: {} Bridge(s)\n'.format(
        len(branched)))
    liner.send('        Time Elapsed: {:.2f} sec\n'.format(
        tt.time()-t0))
    
    # Return Results
    return branched

def bridge_engine(
    exmotifs,
    liner):
    
    # Preprocessing and BLL Computation
    liner.send('\n[Preprocessing Excluded Motifs]\n')    
    exmotifs = enque_exmotifs(
        exmotifs=exmotifs,
        liner=liner)
    bll = get_bll(
        exmotifs=exmotifs,
        liner=liner)

    # Initialize RNG
    rng = np.random.default_rng()

    # Initialize Bridge Candidates
    bridges = list(get_bases(rng))

    # Parition Exmotifs
    liner.send('\n[Partitioning Exmotifs]\n')
    bl = 1
    lt, exmotifs = partition_exmotifs(
        lt=cx.deque(),
        ge=exmotifs,
        bl=bl,
        liner=liner)

    # Bridge-mer Loop
    while bl <= bll:

        # Bound Candidates
        liner.send('\n[Bounding Level {} Bridges]\n'.format(bl))
        partial, complete = get_bounded_bridges(
            lt=lt,
            ge=exmotifs,
            bridges=bridges,
            liner=liner)

        # We have complete solutions!
        if complete:            
            return complete # Recovered Solutions

        # We cannot branch anymore, RIP!
        if not partial:
            return None # No Solution
        
        # Branch Candidates
        liner.send('\n[Branching Level {} Bridges]\n'.format(bl+1))
        bridges = get_branched_bridges(
            rng=rng,
            partial=partial,
            liner=liner)

        # Partition Exmotifs
        liner.send('\n[Partitioning Exmotifs]\n')
        bl += 1
        lt, exmotifs = partition_exmotifs(
            lt=lt,
            ge=exmotifs,
            bl=bl,
            liner=liner)

    # Nothing Computed
    return None

def main():

    liner = ut.liner_engine()

    exmotifs = ['AAAA', 'GGGG', 'CCCC', 'TTTT', 'GGATCC', 'TCTAGA', 'GAATTC'][::-1]

    exmotifs = []
    exmotifs += [''.join(x) for x in ix.product('ATGC', repeat=4)]
    exmotifs += [''.join(x) for x in ix.product('ATGC', repeat=5)]
    exmotifs += [''.join(x) for x in ix.product('ATGC', repeat=6)]
    np.random.shuffle(exmotifs)
    exmotifs = exmotifs[:256]

    # print(len(exmotifs))

    # exmotifs = ['AA',  'TT', 'GG', 'CC', 'A', 'T', 'G', 'C']

    soln = bridge_engine(exmotifs, liner)
    if soln:
        print(len(soln))
    else:
        print(soln)

if __name__ == '__main__':
    main()