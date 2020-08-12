import time  as tt

import collections as cx

import numpy as np
import numba as nb

import utils as ut


def get_gc_run_vec(seq):
    '''
    Return GC run vector for melting temp calc.
    Internal use only.

    :: seq
       type - string
       desc - sequence to decompose
    '''

    # Compute run vector
    rv = np.zeros(len(seq)+1, dtype=np.float64)
    i = 0
    while i < len(seq):
        rv[i+1] = rv[i] + (1. if seq[i] in 'GC' else 0.)
        i += 1

    # Return result
    return rv

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
       type - integer
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
    t3 = 0.41 * 100. * ((rv[j] - rv[i]) / (j-i))
    t4 = 600. / (j-i)

    # Calculate correction
    cr = 0.
    if (j-i) < 60:
        if rv[j]   - rv[j-1] > 0 or \
           rv[i+1] - rv[i]   > 0:
            cr += 1.2
        if rv[j-1] - rv[j-2] > 0 or \
           rv[i+2] - rv[i+1] > 0:
            cr += 0.8
        if cr == 0:
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
    
    return np.array(tuple(ord(nt) for nt in seq))

def get_seqmat(seqlist, numseq, seqlen):
    '''
    Return the numeric representation
    of seqlist. Internal use only.

    :: seqlist
       type - iterable
       desc - list of sequences to split
    :: numseq
       type - integer
       desc - total number of elements in
              seqlist
    :: seqlen
       type - integer
       desc - length of each sequence
    '''
    
    seqmat = np.zeros((numseq, seqlen), dtype=np.float64)
    for idx,seq in enumerate(seqlist):
        seqmat[idx, :] = get_seqvec(seq)
    return seqmat

def get_entvec(seqmat):
    '''
    Return entropy of sequence matrix.
    Internal use only.

    :: seqmat
       type - np.array
       desc - numeric sequence matrix
    '''

    d = seqmat / seqmat.sum(0)
    return 1. - ((d*np.log2(d)).sum(0) * -1)

def get_varcont(entvec):
    '''
    Return all variable region indices.
    Internal use only.

    :: entvec
       type - np.array
       desc - positional entropy
    '''

    # Setup Parsing
    varcont = []
    start   = None
    end     = None
    entvec  = cx.deque(entvec)
    i       = -1

    # Parse Contigs
    while entvec:

        # Update index
        i += 1

        # Extract entropy
        ent = entvec.popleft()

        # Constant Region
        if ent == 0.:

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

def get_num_oligos(seqlen, splitlen):
    if seqlen < splitlen:
        return seqlen
    balance = seqlen % splitlen
    if balance > seqlen // 2:
        return seqlen // splitlen + 2
    return seqlen // splitlen + 1

def split_engine(
    seqlist,
    splitlen=170,
    mintm=50,
    minhdist=10):
    seqmat = get_seqmat(
        seqlist=seqlist,
        numseq=len(seqlist),
        seqlen=len(seqlist[0]))
    print(seqmat)
    entvec = get_entvec(seqmat)
    print(entvec)

    print(get_varcont(entvec))

    # Assess if splitting is feasible



def get_DNA(l):
    fld = ['A', 'T', 'G', 'C']
    return ''.join(np.random.choice(fld) for _ in range(l))

def main():
    print(get_num_oligos(90, 30))

    #       Constant                                         Constant
    #       ------------------                  ------------------
    seq1 = 'CCATAGTCAGACGCATCGAGAGTAGGCTGAGAGTGAAATCTGCGCATATCGACG'
    seq2 = 'CCATAGTCAGACGCATCGCCGACTCCAATCCTAGACAATCTGCGCATATCGACG'

    seqlist = [seq1, seq2]
    # split_engine(seqlist)

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