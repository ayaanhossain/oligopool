import sys
import os
import shutil  as su
import time    as tt

import uuid    as uu
import atexit  as ae

import collections as cx

import numpy   as np
import numba   as nb
import primer3 as p3


# Global Lookups
complement_table = str.maketrans(
    'ACGTURYSWKMBVDHN',
    'TGCATKYWSRMBDHVN')

dna_alpha = set(
    'ATGC')

dna_space = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'S': {'G', 'C'},
    'W': {'A', 'T'},
    'K': {'G', 'T'},
    'M': {'A', 'C'},
    'B': {'C', 'G', 'T'},
    'V': {'A', 'C', 'G'},
    'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'},
    'N': {'A', 'T', 'G', 'C'}}

typeIIS_dict = {
        'acui'    : 'CTGAAG'  + 'N' * 16,
        'alwi'    : 'GGATC'   + 'N' *  5,
        'bbsi'    : 'GAAGAC'  + 'N' *  6,
        'bcci'    : 'CCATC'   + 'N' *  5,
        'bceai'   : 'ACGGC'   + 'N' * 14,
        'bcivi'   : 'GTATCC'  + 'N' *  6,
        'bcodi'   : 'GTCTC'   + 'N' *  5,
        'bmri'    : 'ACTGGG'  + 'N' *  5,
        'bpuei'   : 'CTTGAG'  + 'N' * 16,
        'bsai'    : 'GGTCTC'  + 'N' *  5,
        'bseri'   : 'GAGGAG'  + 'N' * 10,
        'bsmai'   : 'GTCTC'   + 'N' *  5,
        'bsmbi'   : 'CGTCTC'  + 'N' *  5,
        'bsmfi'   : 'GGGAC'   + 'N' * 14,
        'bsmi'    : 'GAATGC'  + 'N' *  1,
        'bspcni'  : 'CTCAG'   + 'N' *  9,
        'bspqi'   : 'GCTCTTC' + 'N' *  4,
        'bsrdi'   : 'GCAATG'  + 'N' *  2,
        'bsri'    : 'ACTGG'   + 'N' *  1,
        'btgzi'   : 'GCGATG'  + 'N' * 14,
        'btsci'   : 'GGATG'   + 'N' *  2,
        'btsi'    : 'GCAGTG'  + 'N' *  2,
        'btsimuti': 'CAGTG'   + 'N' *  2,
        'eari'    : 'CTCTTC'  + 'N' *  4,
        'ecii'    : 'GGCGGA'  + 'N' * 11,
        'esp3i'   : 'CGTCTC'  + 'N' *  5,
        'faui'    : 'CCCGC'   + 'N' *  6,
        'hgai'    : 'GACGC'   + 'N' * 10,
        'hphi'    : 'GGTGA'   + 'N' *  8,
        'hpyav'   : 'CCTTC'   + 'N' *  6,
        'mlyi'    : 'GAGTC'   + 'N' *  5,
        'mnli'    : 'CCTC'    + 'N' *  7,
        'sapi'    : 'GCTCTTC' + 'N' *  4,
        'sfani'   : 'GCATC'   + 'N' *  9,
    }

# Decorators

def coroutine(func):
    '''
    Decorator used to prime a coroutine.
    Internal use only.

    :: func
       type - coroutine
       desc - a coroutine we want online
    '''

    def online(*args, **kwargs):
        _coroutine = func(*args, **kwargs)
        next(_coroutine)
        return _coroutine
    
    return online

# Printing / Logging

@coroutine
def liner_engine(online=True):
    '''
    Return a coroutine to print stuff
    to stdout and log run information.
    Internal use only.

    :: online
       type - boolean
       desc - if True will print received
              strings to stdout
              (default=True)
    '''
    
    # Book-keeping
    clrlen = 0

    # Online
    try:
        while True:
            
            # Receive String
            printstr = (yield)

            # Line Preservation
            if printstr.startswith('\*'):
                sys.stdout.write('\n')
                clrlen   = 0
                printstr = printstr.lstrip('\*')
            
            # Line String
            if online:
                sys.stdout.write('\r' + ' '*clrlen)
                sys.stdout.write('\r' + printstr)
                clrlen = len(printstr)
                sys.stdout.flush()
    
    # Closure
    except GeneratorExit:
        sys.stdout.write('\n')

# Numeric Functions

def get_trials(prob):
    '''
    Return the number of required trials
    given a probability of success.
    Internal use only.

    :: prob
       type - float
       desc - probability of a rate event
    '''
    growth_rate = -100.0 # Lower is Rarer
    # print(-10 / np.log(1. - prob))
    return 2 * int(np.ceil(
        growth_rate / np.log(1. - prob)))

def get_prob(success, trials):
    '''
    Return the probability of success given
    the number of successes in some finite
    trials. Internal use only.

    :: success
       type - integer
       desc - total number of successes
    :: trials
       type - integer
       desc - total number of attempts made
    '''
    return (success + 1.) / (trials + 2.)

# Oligo Functions

@nb.njit
def get_store_hdist(
    store,
    idx,
    i=0,
    j=None,
    direction=0):
    '''
    Return the minimum pairwise hamming
    distance between store[idx, i:j] and
    store[*: i:j] depending on direction.
    Internal use only.

    :: store
       type - np.array
       desc - numeric sequence array
    :: idx
       type - integer
       desc - location index of sequence
              being compared
    :: i
       type - integer
       desc - starting index of comparison
    :: j
       type - integer
       desc - ending   index of comparison
    :: direction
       type - integer
       desc - direction flag;
              0 = upward   comparison
              1 = downward comparison
              2 = all-pair comparison
    '''

    # Default Result
    hdist = store.shape[1]
    if i >= 0 and not j is None:
        hdist = j-i

    # Something to compare against?
    if idx > 0:

        # Upward / All-Pair Comparison
        if direction == 0 or direction == 2:
            hdist = min(
                hdist,
                (store[:idx, i:j] != store[idx, i:j]).sum(1).min())

    # Something to compare against?
    if idx < store.shape[0] - 1:

        # Downward / All-Pair Comparison
        if direction == 1 or direction == 2:
                hdist = min(
                    hdist,
                    (store[idx+1:, i:j] != store[idx, i:j]).sum(1).min())

    # Return Result
    return hdist

def get_tmelt(
    seq,
    i=0,
    j=None,
    mvc=50.,
    dvc=0.,
    ntc=0.8,
    olc=50.):
    '''
    Return the melting temperature of seq.
    Internal use only.

    Note: primer3-py API docs

    :: seq
       type - string
       desc - DNA string in context
    :: i
       type - integer
       desc - starting index of comparison
              (default=)
    :: j
       type - integer
       desc - ending index of comparison
    :: mvc
       type - float
       desc - monovalent cation conc. (mM)
              (default=50.0 nM)
    :: dvc
       type - float
       desc - divalent cation conc. (mM)
              (default=0.0 nM)
    :: ntc
       type - float
       desc - dNTP conc. (mM)
              (default=0.8 mM)
    :: olc
       type - float
       desc - oligo conc. (nM)
              (default=50 nM)
    '''
    return p3.bindings.calcTm(
        seq=seq[i:j],
        mv_conc=mvc,
        dv_conc=dvc,
        dntp_conc=ntc,
        dna_conc=olc)

# Motif Functions

def prep_exmotifs(exmotifs, packing, liner):
    '''
    Return all exmotifs sorted by
    length. Internal use only.

    :: exmotifs
       type - list
       desc - list of all motifs
              to be excluded
    :: packing
       type - function
       desc - factory function for
              wrapping exmotifs
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
    dq = packing(exmotifs)
    liner.send(' Enqued {} Motifs in {:.2f} sec\n'.format(
        len(exmotifs), tt.time()-t0))

    # Return Results
    return dq

def get_context_len(exmotifs):
    '''
    Get the exmotif context length.
    Internal use only.

    :: exmotifs
       type - iterable / None
       desc - iterable of motifs to exlude,
              None otherwise
    '''

    return len(exmotifs[-1]) - 1

def get_context_type(context):
    # All contexts unique
    if isinstance(context, list):
        return 1

    # All contexts constant
    elif isinstance(context, str):
        return 2
    
    # No context
    elif context == None:
        return 3

    # Unknown context
    else:
        return -1

def get_context_num(context):

    cntxnum = 0
    
    cntxtype = get_context_type(
        context=context)

    # All contexts unique
    if cntxtype   == 1:
        cntxnum = max(
            cntxnum,
            len(context))
    
    # All contexts constant
    elif cntxtype == 2:
        cntxnum = max(
            cntxnum,
            1)

    return cntxnum

def get_context_inference_fn(context):
    '''
    Return selector functions for building
    barcode context. Internal use only.

    :: context
       type - list / string / None
       desc - list of context sequence,
              or a single context sequence,
              or None
    '''

    # Get context type
    cntxtype = get_context_type(context)

    # All contexts unique
    if cntxtype == 1:
        cifn = lambda x: context[x]

    # All contexts constant
    elif cntxtype == 2:
        cifn = lambda x: context
    
    # No context
    elif cntxtype == 3:
        cifn = lambda x: ''

    # Unknown packing
    else:
        raise ValueError(
            'Context Packing Unknown')

    # Return Selector
    return cifn

def get_motif_conflict(
    seq,
    seqlen,
    exmotifs,
    partial=False,
    checkall=False):
    '''
    Determine if the barcode does not contain or is
    contained in one of the exluded motifs (motif
    feasibility). Internal use only.

    :: seq
       type - string
       desc - sequence for conflict checking
    :: seqlen
       type - integer
       desc - length of sequence
    :: exmotifs
       type - iterable / None
       desc - list of motifs to exclude
              in designed barcodes,
              None otherwise
    :: partial
       type - boolean
       desc - if True will not check for
              conflicts for motifs longer
              than seqlen, otherwise run
              all checks
              (default=False)
    :: checkall
       type - boolean
       desc - if True will check for all
              possible conflicts due to
              exmotifs
              (default=False)
    '''

    # Book-keeping
    status = True
    pmotif = set() if checkall else None

    # Do we have anything to exclude?
    if exmotifs:

        # Loop through all motifs
        for motif in exmotifs:

            # Is sequence partial?
            if partial:
                if len(motif) > seqlen:
                    break # No need to check further

            # Embedding and Embedded Conflict
            if (len(motif) <= seqlen and motif in seq) or \
               (len(motif)  > seqlen and seq   in motif):
                
                # We got Conflict!
                status = False
                
                # What is the Conflict?
                if checkall:
                    pmotif.add(motif)
                else:
                    pmotif = motif
                    break
    
    # No Motif Conflct
    return (status, pmotif)

def get_assignment_index(
    seq,
    exmotifs,
    lcifn,
    rcifn,
    cntxlen,
    carr):
    '''
    Determine the context assignment index
    for designed sequence. Internal use only.

    :: seq
       type - string
       desc - designed sequence to be
              assigned
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
    
    # Book-keeping
    i      = 0
    cfails = cx.Counter()

    # Loop through all contexts for assignment
    while i < len(carr):

        # Fetch Context
        aidx = carr.popleft()

        # Fetch left and right contexts
        lcntx = lcifn(aidx)[-cntxlen:]
        rcntx = rcifn(aidx)[:+cntxlen]

        # Build out in-context sequence
        incntxseq = lcntx + seq + rcntx

        # Determine context feasibility
        mcond, motif = get_motif_conflict(
            seq=incntxseq,
            seqlen=len(incntxseq),
            exmotifs=exmotifs,
            partial=False)

        # We have a winner folks!
        if mcond:
            return (True, aidx, None)

        # We lost an assignment, record cause
        else:
            cfails[motif] += 1

        # Update Iteration
        i += 1
        
        # We will try this context again
        carr.append(aidx)

    # Failed Assignment
    return (False, None, cfails)

def is_all_assignable(
    seq,
    exmotifs,
    cifn,
    cntxtype,
    cntxnum,
    cntxlen):
    
    aidx  = 0

    while aidx < cntxnum:

        if cntxtype == 0:
            lcntx = cifn(aidx)[-cntxlen:]
            incntxseq = lcntx + seq
        else:
            rcntx = cifn(aidx)[:+cntxlen]
            incntxseq = seq + rcntx

        # Determine context feasibility
        mcond, motif = get_motif_conflict(
            seq=incntxseq,
            seqlen=len(incntxseq),
            exmotifs=exmotifs,
            partial=True)

        if not mcond:
            mstart = incntxseq.find(motif)

            # Context on Left
            if cntxtype == 0:
                # print((lcntx, seq, 'left'))
                tloc = len(motif) - len(lcntx) + mstart - 1
            # Context on Right
            else:
                # print((seq, rcntx, 'right'))
                tloc = mstart

            return False, tloc

        aidx += 1

    return True, None

# Sequence Analysis

def is_DNA(seq):
    '''
    Determine if seq is a DNA string.
    Internal use only.

    :: seq
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    '''

    if isinstance(seq, str) and \
       set(seq) <= dna_alpha:
        return True
    return False

def get_comp(seq):
    '''
    Return the complement of seq.
    Internal use only.

    :: seq
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    '''
    
    return seq.translate(complement_table)

def get_revcomp(seq):
    '''
    Return the reverse complement of seq.
    Internal use only.

    :: seq
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    '''
    
    return get_comp(seq)[::-1]

def stream_spectrum(seq, k):
    '''
    Stream the k-mer spectrum of seq.
    Internal use only.

    :: seq
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: k
       type - integer
       desc - k-mer length / k-value
    '''
    
    return (seq[i:i+k] for i in range(len(seq)-k+1))

def stream_canon_spectrum(seq, k):
    '''
    Stream the canonical k-mer spectrum
    of seq. Internal use only.

    :: seq
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: k
       type - integer
       desc - k-mer length / k-value
    '''
    
    return map(
        lambda x: min(x, get_revcomp(x)),
        stream_spectrum(seq=seq, k=k))

def get_hdist(seq1, seq2, max_hd=None):
    '''
    Compute the hamming distance between
    seq1 and seq2, within max_hd bound.
    Internal use only.

    :: seq1
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: seq2
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: max_hd
       type - integer / None
       desc - compute hamming distance
              between seq1 and seq2,
              within max_hd bound,
              return None otherwise
              (default=None)
    '''
    
    seq_len = min(len(seq1), len(seq2))
    i = 0
    hdist = 0
    while i < seq_len:
        hdist += seq1[i] != seq2[i]
        i += 1
        if not max_hd is None and \
           hdist > max_hd:
            return None
    return hdist

def get_edist(seq1, seq2, max_ed=None, mode='NW'):
    '''
    Compute the mode edit distance between
    seq1 and seq2, within max_ed bound.
    Internal use only.

    :: seq1
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: seq2
       type - string
       desc - a string in alphabet
              {A, T, G, C}
    :: max_hd
       type - integer / None
       desc - compute edit distance
              between seq1 and seq2,
              within max_ed bound,
              return None otherwise
    :: mode
       type - string
       desc - 'NW' / 'HW' for global
              / infix alignments
    '''
    
    if len(seq2) > len(seq2):
        seq1, seq2 = seq2, seq1
    i = 0
    if max_ed is None:
        max_ed = -1
    ed_align = ed.align(
        seq2,
        seq1,
        mode=mode,
        task='distance',
        k=max_ed)
    edist = ed_align['editDistance']
    if edist == -1:
        return None
    return edist

# Workspace Functions

def get_adjusted_path(path, suffix):
    '''
    Adjust outfile with suffix.
    Internal use only.

    :: path
       type - string
       desc - some path which must
              have given suffix
    :: suffix
       type - string / None
       desc - suffix to concatenate
              with given path
    '''

    if path.endswith('/'):
        path = path.rstrip('/')
    if not suffix is None and \
       not path.endswith(suffix):
        path += str(suffix)
    return path

def get_path_status(
    path,
    suffix=None,
    readable=True,
    writable=False,
    creatable=False):
    '''
    Arbitrary path evaluation.
    Internal use only.

    :: path
       type - string
       desc - path to store or output
              information
    :: suffix
       type - string / None
       desc - required path suffix
              (default=None)
    :: readable
       type - boolean
       desc - if True check is path is
              readable  (file/dir   read)
              (default=True)
    :: writable
       type - boolean
       desc - if True check if path is
              writable  (file/dir  write)
              (default=False)
    :: creatable
       type - boolean
       desc - if True check if path is
              creatable (file/dir create)
              (default=False)
    '''

    # Status:
    # 0  = Path is invalid      input
    # 1  = Path is unreadable   file
    # 2  = Path is unwritable   file
    # 3  = Path is empty        file
    # 4  = Path is non-empty    file
    # 5  = Path is unreadable   directory
    # 6  = Path is unwritable   directory
    # 7  = Path is empty        directory
    # 8  = Path is non-empty    directory
    # 9  = Path is non-existent
    # 10 = Path is creatable
    # 11 = Path is uncreatable
    # X  = Path is existent     special file

    # Is the path valid?
    if not isinstance(path, str) or not path:
        return 0 # only non-empty strings are paths

    # Adjust Path for Checks
    path = get_adjusted_path(
        path=path,
        suffix=suffix)
    
    # Is the path in existence?
    if os.path.exists(path):
        
        # Path points to a file?
        if os.path.isfile(path):

            # Read accessible?
            if  readable and not \
                os.access(path, os.R_OK):
                    return 1 # Unreadable file

            # Is the file empty?
            if not os.stat(path).st_size:

                # Can we write to this file?
                if  writable and not \
                    os.access(path, os.W_OK):
                        return 2 # Unwritable file
                
                return 3 # Empty file
            
            else:
                return 4 # Non-empty file

        # Path points to a directory?
        elif os.path.isdir(path):

            # What's the accessibility like?
            if  readable and not \
                os.access(path, os.R_OK):
                    return 5 # Unreadable directory

            # Directory empty?
            if not any(os.scandir(path)):

                # Can we write to this directory?
                if  writable and not \
                    os.access(path, os.W_OK):
                        return 6 # Unwritable directory

                return 7 # Empty Directory
            
            else:
                return 8 # Non-empty Directory

        # Path points to a special file such as
        # socket, FIFO, device file etc.
        else:
            return 'X' # Special File Identifier

    # Is the path creatable?
    if creatable:

        try:
            # Build parent directory path
            parentpath = os.path.dirname(path) or \
                         os.getcwd()
            # Try building parent path
            os.makedirs(
                name=parentpath,
                exist_ok=True)
        
        except:
            return 11 # Cannot create path
        
        else:
            return 10 # Can create path

    # Path is non-existent
    return 9

def setup_directory(dirpath):
    '''
    Build dirpath. Internal use only.

    :: dirpath
       type - string
       desc - a directory built to
              store information
    '''

    dirpath = get_adjusted_path(
        path=dirpath,
        suffix=None)
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)

def setup_workspace(
    outfile,
    outfile_suffix):
    '''
    Adjust given output file and setup
    output directory, and schedule its
    removal at exit. Internal use only.

    :: outfile
       type - string
       desc - output file
    :: outfile_suffix
       type - string
       desc - output file suffix
    '''

    # Adjust outfile with suffix
    outfile = get_adjusted_path(
        path=outfile,
        suffix=outfile_suffix)

    # Define outdir
    outdir = '.'.join([
        outfile,
        str(uu.uuid4()),
        'dir'])

    # Schedule outdir deletion
    ae.register(
        remove_directory,
        outdir)

    # Setup outdir
    setup_directory(
        dirpath=outdir)

    # Return Workspace
    return outfile, outdir

def remove_file(filepath):
    '''
    Delete file stored at filepath.
    Internal use only.

    :: filepath
       type - string
       desc - a file to be removed
    '''

    filepath = get_adjusted_path(
        path=filepath,
        suffix=None)
    if os.path.exists(filepath):
        os.remove(filepath)

def remove_directory(dirpath):
    '''
    Delete dirpath and its contents.
    Internal use only.

    :: dirpath
       type - string
       desc - a directory (with content)
              to be removed
    '''
    
    dirpath = get_adjusted_path(
        path=dirpath,
        suffix=None)
    if os.path.exists(dirpath):
        su.rmtree(dirpath)

