# Code Snippets that are good ideas, but unused


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