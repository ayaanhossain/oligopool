import time as tt

import collections as cx
import itertools   as ix

import numpy as np
import numba as nb

import utils as ut

def enque_exmotifs(exmotifs):

    # Sort motifs by lengeh
    exmotifs.sort(key=len)
    # print(exmotifs)
    # Return deque
    return cx.deque(exmotifs)

def get_bll(exmotifs):

    bll = 2
    if exmotifs:
        maxel = exmotifs.pop()
        bll   = len(maxel)
        exmotifs.append(maxel)
    return bll

def partition_exmotifs(lt, ge, bl):
    
    while ge:
        motif = ge.popleft()
        if len(motif) < bl:
            lt.append(motif)
        else:
            ge.appendleft(motif)
            break

    return lt, ge

def is_lt_embedded(lt, candidate):
    
    for motif in lt:
        if motif in candidate:
            return True
    return False

def is_ge_embedded(ge, candidate):
    
    for motif in ge:
        if candidate in motif:
            return True
    return False

def get_verdict(lt, ge, candidate):

    ltf = is_lt_embedded(lt, candidate)
    # print(candidate, 'ltf', ltf)

    if ltf:
        return -1 # Rejection I

    gef = is_ge_embedded(ge, candidate)
    # print(candidate, 'gtf', ltf)

    if gef:
        return 0 # Rejection II

    return 1

def get_bases(rng):

    bases = ['A', 'T', 'G', 'C']
    rng.shuffle(bases)
    return bases

def get_branched_candidates(rng, candidates):
    
    expanded = cx.deque()
    for candidate in candidates:
        bases = get_bases(rng)
        for base in bases:
            expanded.append(
                candidate+base)
    return expanded

def get_bounded_candidates(lt, ge, candidates):

    grow   = cx.deque()
    accept = cx.deque()

    for candidate in candidates:

        verdict = get_verdict(
            lt=lt,
            ge=ge,
            candidate=candidate)

        print(candidate, verdict)

        if verdict == 0:
            grow.append(candidate)
        
        if verdict == 1:
            accept.append(candidate)

    return grow, accept    

def bridge_engine(
    exmotifs,
    liner):
    
    exmotifs = enque_exmotifs(exmotifs)
    bll      = get_bll(exmotifs)

    rng = np.random.default_rng()

    candidates = cx.deque(get_bases(rng))
    bl         = 1
    lt, ge = partition_exmotifs(
        lt=cx.deque(),
        ge=exmotifs,
        bl=bl)

    print(lt)
    print(ge)

    while bl < bll:

        print(candidates)

        grow, accept = get_bounded_candidates(lt, ge, candidates)

        if accept:
            print(accept, 'Accepted')
            return accept

        if not grow:
            print('Impossible')
            return None

        candidates = get_branched_candidates(rng, grow)

        bl += 1

        lt, ge = partition_exmotifs(
            lt=lt,
            ge=ge,
            bl=bl)

        print(lt)
        print(ge)

def main():

    liner = ut.liner_engine()

    exmotifs = ['AAAA', 'GGGG', 'CCCC', 'TTTT', 'GGATCC', 'TCTAGA', 'GAATTC'][::-1]

    # exmotifs = []
    # exmotifs += [''.join(x) for x in ix.product('ATGC', repeat=3)][:8]
    # np.random.shuffle(exmotifs)
    # exmotifs = exmotifs[:7]

    print(len(exmotifs))

    # exmotifs += ['AA',  'TT', 'GG', 'CC', 'A', 'T', 'G', 'C']

    soln = bridge_engine(exmotifs, liner)
    if soln:
        print(len(soln))
    else:
        print(soln)

if __name__ == '__main__':
    main()