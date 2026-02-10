import os

import time as tt

import itertools   as ix
import collections as cx
import random      as rn

import numpy    as np
import numba    as nb
import parasail as ps
import edlib    as ed

from . import utils as ut


# Edit Distance Matrix

editmat = ps.matrix_create("ACGTN-", 1, -1)

# Mismatch Density Threshold
# Maximum allowed mismatch rate in overlap region

MAX_MISMATCH_DENSITY = 0.10  # 10% mismatch rate

# Parser and Setup Functions

def get_extracted_overlap_parameters(
    r1file,
    r2file,
    r1type,
    r2type,
    liner):
    '''
    Determine the priority of merging
    functions from FastQ readfiles.
    Internal use only.

    :: r1file
       type - string
       desc - path to FastQ file storing
              R1 reads
    :: r2file
       type - string / None
       desc - path to FastQ file storing
              R2 reads
    :: r1type
       type - integer
       desc - R1 read type identifier
              0 = Forward Reads
              1 = Reverse Reads
    :: r2type
       type - integer / None
       desc - R2 read type identifier
              0 = Forward Reads
              1 = Reverse Reads
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Start Timer
    t0 = tt.time()

    # Setup R1 Streamer
    reader1 = ut.stream_fastq_engine(
        filepath=r1file,
        invert=r1type,
        qualvec=False,
        skipcount=0,
        coreid=1,
        ncores=2)

    # Setup R2 Streamer
    reader2 = ut.stream_fastq_engine(
        filepath=r2file,
        invert=r2type,
        qualvec=False,
        skipcount=0,
        coreid=1,
        ncores=2)

    # Book-keeping
    samplecount = 0
    sampletotal = 0.5 * 10**5
    reachcount  = 0
    inniecount  = 0
    outiecount  = 0
    mergefnhigh = None
    truncread   = False

    # Process Read Samples
    while samplecount < sampletotal:

        # Get Paired Reads
        try:
            r1,_ = next(reader1)
            r2,_ = next(reader2)
        except StopIteration:
            break

        # Exhausted / Truncated File Pair
        if r1 is None and r2 is None:
            break
        if (r1 is None) or (r2 is None):
            truncread = True
            break

        # Update Book-keeping
        samplecount += 1
        reachcount  += 1

        # Show Update
        if reachcount == 10000:
            liner.send(
                ' Reads Analyzed: {:,}'.format(
                    samplecount))
            reachcount = 0

        # Find Innie Overlap Location
        innieend = ps.sg_qb_de_striped_16(
            s1=r1,         # Query
            s2=r2,         # Reference
            open=1,        # Gap Open
            extend=1,      # Gap Extension
            matrix=editmat # +1/-1 Edit Matrix
        ).end_ref + 1

        # Compute Innie Score
        if innieend > 10:
            inniescore = ed.align(
                query=r2[:innieend],
                target=r1,
                mode='HW',
                task='distance',
                k=10)['editDistance']
            if inniescore == -1:
                inniescore = float('inf')
        else:
            inniescore = float('inf')

        # Find Outie Overlap Location
        outieend = ps.sg_qb_de_striped_16(
            s1=r2,         # Query
            s2=r1,         # Reference
            open=1,        # Gap Open
            extend=1,      # Gap Extension
            matrix=editmat # +1/-1 Edit Matrix
        ).end_ref + 1

        # Compute Outie Score
        if outieend > 10:
            outiescore = ed.align(
                query=r1[:outieend],
                target=r2,
                mode='HW',
                task='distance',
                k=10)['editDistance']
            if outiescore == -1:
                outiescore = float('inf')
        else:
            outiescore = float('inf')

        # Which one did the best?
        if   inniescore < outiescore: # Innie Winner
            inniecount += 1
        elif outiescore < inniescore: # Outie Winner
            outiecount += 1
        else:
            if inniescore < float('inf') and \
               outiescore < float('inf') and \
               inniescore == outiescore:
                # No Clear Winner ...
                inniecount += 1
                outiecount += 1

    # Normalize Counts
    totalcount = inniecount + outiecount
    if totalcount > 0:
        inniecount = (100. * inniecount) / totalcount
        outiecount = 100. - inniecount
    else:
        # Default conservatively when no valid overlap evidence is available.
        inniecount = 100.
        outiecount = 0.

    # Show Update
    liner.send(
        ' Reads Analyzed: {:,}\n'.format(
            samplecount))
    liner.send(
        ' Innie Events  : {:6.2f} %\n'.format(
            inniecount))
    liner.send(
        ' Outie Events  : {:6.2f} %\n'.format(
            outiecount))
    if truncread:
        liner.send(
            ' Pair Truncation Detected During Overlap Sampling\n')
    if totalcount == 0:
        liner.send(
            ' Overlap Evidence: Insufficient, Using Default Priority\n')

    # Compute Results
    if inniecount >= outiecount:
        mergefnhigh = get_innie_merged
        mergefnlow  = get_outie_merged
        liner.send(
            ' Priority: Innie Merging\n')
    else:
        mergefnhigh = get_outie_merged
        mergefnlow  = get_innie_merged
        liner.send(
            ' Priority: Outie Merging\n')

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time()-t0))

    # Return Results
    return mergefnhigh, mergefnlow

# Engine Helper Functions

@nb.njit
def get_exact_innie_overlap(r1_bytes, r2_bytes, min_overlap=10):
    '''
    Find exact suffix-prefix overlap between R1 and R2.
    Internal use only.

    :: r1_bytes
       type - np.array (uint8)
       desc - R1 read as byte array
    :: r2_bytes
       type - np.array (uint8)
       desc - R2 read as byte array
    :: min_overlap
       type - integer
       desc - minimum overlap length to consider

    Returns overlap length if found, -1 otherwise.
    '''

    # Compute max overlap
    max_overlap = min(len(r1_bytes), len(r2_bytes))

    # Search from longest to shortest (greedy)
    for olen in range(max_overlap, min_overlap - 1, -1):
        # Check suffix of R1 against prefix of R2
        match = True
        r1_start = len(r1_bytes) - olen
        for i in range(olen):
            if r1_bytes[r1_start + i] != r2_bytes[i]:
                match = False
                break
        if match:
            return olen

    return -1  # No exact match found

@nb.njit
def get_exact_outie_overlap(r1_bytes, r2_bytes, min_overlap=10):
    '''
    Find exact suffix-prefix overlap for outie orientation.
    Internal use only.

    :: r1_bytes
       type - np.array (uint8)
       desc - R1 read as byte array
    :: r2_bytes
       type - np.array (uint8)
       desc - R2 read as byte array
    :: min_overlap
       type - integer
       desc - minimum overlap length to consider

    Returns overlap length if found, -1 otherwise.
    '''

    # Compute max overlap
    max_overlap = min(len(r1_bytes), len(r2_bytes))

    # Search from longest to shortest (greedy)
    for olen in range(max_overlap, min_overlap - 1, -1):
        # Check prefix of R1 against suffix of R2
        match = True
        r2_start = len(r2_bytes) - olen
        for i in range(olen):
            if r1_bytes[i] != r2_bytes[r2_start + i]:
                match = False
                break
        if match:
            return olen

    return -1  # No exact match found

@nb.njit
def qualpass(qvec, threshold):
    '''
    Return Q-Score vector.
    Internal use only.

    :: qvec
       type - np.array
       desc - Q-Score vector
    :: threshold
       type - integer
       desc - quality threshold
    '''

    return np.mean(qvec) >= threshold

def get_concatenated_reads(r1, r2):
    '''
    Return (r1, r2) tuple for storage.
    Internal use only.

    :: r1
       type - string / bytes
       desc - R1 read sequence
    :: r2
       type - string / bytes
       desc - R2 read sequence

    Returns tuple (r1, r2) or (r1, None) if r2 is None.
    '''

    # Note: upstream FastQ streaming may supply bytes for speed.
    # Pack files are expected to store read sequences as strings.
    if type(r1) is bytes:
        r1 = r1.decode('ascii')
    if type(r2) is bytes:
        r2 = r2.decode('ascii')
    return (r1, r2)

@nb.njit
def get_innie_consensus_fast(r1_bytes, r2_bytes, qv1, qv2, olen):
    '''
    Return consensus innie merged read using byte arrays.
    Internal use only.

    :: r1_bytes
       type - np.array (uint8)
       desc - R1 read as byte array
    :: r2_bytes
       type - np.array (uint8)
       desc - R2 read as byte array
    :: qv1
       type - np.array (int8)
       desc - R1 Q-Score vector
    :: qv2
       type - np.array (int8)
       desc - R2 Q-Score vector
    :: olen
       type - integer
       desc - overlap length
    '''

    # Pre-allocate result array
    r1_len = len(r1_bytes)
    r2_len = len(r2_bytes)
    result_len = r1_len + r2_len - olen
    result = np.empty(result_len, dtype=np.uint8)

    # Non-overlapping prefix from R1
    prefix_len = r1_len - olen
    for i in range(prefix_len):
        result[i] = r1_bytes[i]

    # Overlap region - quality-based selection
    r1_start = r1_len - olen
    for k in range(olen):
        if qv1[r1_start + k] < qv2[k]:
            result[prefix_len + k] = r2_bytes[k]
        else:
            result[prefix_len + k] = r1_bytes[r1_start + k]

    # Non-overlapping suffix from R2
    for i in range(olen, r2_len):
        result[r1_len + i - olen] = r2_bytes[i]

    return result

def get_innie_consensus(r1, r2, qv1, qv2, olen):
    '''
    Return consensus innie merged read.
    Internal use only.

    :: r1
       type - string / bytes
       desc - R1 read sequence
    :: r2
       type - string / bytes
       desc - R2 read sequence
    :: qv1
       type - np.array
       desc - R1 Q-Score vector
    :: qv2
       type - np.array
       desc - R2 Q-Score vector
    :: olen
       type - integer
       desc - overlap length
    '''

    # Check if input is bytes
    is_bytes = type(r1) is bytes

    # Convert to byte arrays
    if is_bytes:
        r1_bytes = np.frombuffer(r1, dtype=np.uint8)
        r2_bytes = np.frombuffer(r2, dtype=np.uint8)
    else:
        r1_bytes = np.frombuffer(r1.encode(), dtype=np.uint8)
        r2_bytes = np.frombuffer(r2.encode(), dtype=np.uint8)

    # Compute consensus using fast path
    result_bytes = get_innie_consensus_fast(
        r1_bytes, r2_bytes, qv1, qv2, olen)

    # Return same type as input
    if is_bytes:
        return result_bytes.tobytes()
    else:
        return result_bytes.tobytes().decode('ascii')

def get_innie_merged(r1, r2, qv1, qv2):
    '''
    Return innie-merged r1 and r2.
    Internal use only.

    :: r1
       type - string / bytes
       desc - R1 read sequence
    :: r2
       type - string / bytes
       desc - R2 read sequence
    :: qv1
       type - np.array
       desc - R1 Q-Score vector
    :: qv2
       type - np.array
       desc - R2 Q-Score vector
    '''

    # Check if input is bytes (optimized path)
    is_bytes = type(r1) is bytes

    # Fast-path: Check for exact overlap first
    if is_bytes:
        r1_bytes = np.frombuffer(r1, dtype=np.uint8)
        r2_bytes = np.frombuffer(r2, dtype=np.uint8)
    else:
        r1_bytes = np.frombuffer(r1.encode(), dtype=np.uint8)
        r2_bytes = np.frombuffer(r2.encode(), dtype=np.uint8)

    exact_olen = get_exact_innie_overlap(
        r1_bytes=r1_bytes,
        r2_bytes=r2_bytes,
        min_overlap=10)

    # Exact match found - direct concatenation
    if exact_olen > 0:
        merged = r1 + r2[exact_olen:]
        return merged, 0  # Score 0 = perfect match

    # Slow-path: Use alignment for reads with mismatches
    # Note: edlib expects strings; keep byte-optimized consensus paths but
    # decode once for alignment when upstream streaming uses bytes.
    if is_bytes:
        r1_ed = r1.decode('ascii')
        r2_ed = r2.decode('ascii')
    else:
        r1_ed = r1
        r2_ed = r2

    # Adaptive seed length based on read length
    seed_len = max(10, min(20, len(r1) // 10))

    # Adaptive k based on estimated overlap and mismatch density
    est_overlap = int(0.4 * min(len(r1), len(r2)))
    k_seed = max(4, int(est_overlap * MAX_MISMATCH_DENSITY))

    aln = ed.align(
        query=r1_ed[-seed_len:],
        target=r2_ed,
        mode='HW',
        task='distance',
        k=k_seed)

    # Bad Alignment?
    if aln['editDistance'] == -1:
        # Too many mutations,
        # Skip this Read
        aln.clear() # Memory Management
        del aln     # Memory Management
        return None, float('inf')

    # Extract End Locations
    end_locs = aln['locations']

    # Clear Memory
    aln.clear() # Memory Management
    del aln     # Memory Management

    # Too Many Overlap Locations?
    if len(end_locs) > 10:
        # Likely Lots of Mutation,
        # Skip this Read
        end_locs.clear() # Memory Management
        del end_locs     # Memory Management
        return None, None

    # Find Overlap End Coordinate
    end = end_locs[0][-1] + 1

    # Too Small of an Overlap?
    if end < 10:
        # Not High Confidence,
        # Skip this Read
        end_locs.clear() # Memory Management
        del end_locs     # Memory Management
        return None, None

    # Define Query Slice
    qslice = r2_ed[:end]

    # Density-based k for validation alignment
    k_validate = max(10, int(end * MAX_MISMATCH_DENSITY))

    # Infix Align r1 and r2
    aln = ed.align(
        query=qslice,
        target=r1_ed,
        mode='HW',
        task='path',
        k=k_validate)

    # Bad Alignment?
    if aln['editDistance'] == -1:
        # Too many mutations,
        # Skip this Read
        del qslice  # Memory Management
        aln.clear() # Memory Management
        del aln     # Memory Management
        return None, float('inf')

    # Do we have Indels?
    if 'D' in aln['cigar'] or \
       'I' in aln['cigar']:
        # We don't care about Indels,
        # Skip this Read
        del qslice  # Memory Management
        aln.clear() # Memory Management
        del aln     # Memory Management
        return None, float('inf')

    # Compute Merged Read
    score = aln['editDistance']
    if score == 0:
        merged = r1 + r2[end:]
    else:
        merged = get_innie_consensus(
            r1=r1,
            r2=r2,
            qv1=qv1,
            qv2=qv2,
            olen=end)

    # Memory Management
    del qslice
    aln.clear()
    del aln

    # Return Results
    return merged, score

@nb.njit
def get_outie_consensus_fast(r1_bytes, r2_bytes, qv1, qv2, olen):
    '''
    Return consensus outie merged read using byte arrays.
    Internal use only.

    :: r1_bytes
       type - np.array (uint8)
       desc - R1 read as byte array
    :: r2_bytes
       type - np.array (uint8)
       desc - R2 read as byte array
    :: qv1
       type - np.array (int8)
       desc - R1 Q-Score vector
    :: qv2
       type - np.array (int8)
       desc - R2 Q-Score vector
    :: olen
       type - integer
       desc - overlap length
    '''

    # Pre-allocate result array (overlap region only)
    result = np.empty(olen, dtype=np.uint8)

    # Overlap region - quality-based selection
    r2_start = len(r2_bytes) - olen
    for k in range(olen):
        if qv1[k] < qv2[r2_start + k]:
            result[k] = r2_bytes[r2_start + k]
        else:
            result[k] = r1_bytes[k]

    return result

def get_outie_consensus(r1, r2, qv1, qv2, olen):
    '''
    Return consensus outie merged read.
    Internal use only.

    :: r1
       type - string / bytes
       desc - R1 read sequence
    :: r2
       type - string / bytes
       desc - R2 read sequence
    :: qv1
       type - np.array
       desc - R1 Q-Score vector
    :: qv2
       type - np.array
       desc - R2 Q-Score vector
    :: olen
       type - integer
       desc - overlap length
    '''

    # Check if input is bytes
    is_bytes = type(r1) is bytes

    # Convert to byte arrays
    if is_bytes:
        r1_bytes = np.frombuffer(r1, dtype=np.uint8)
        r2_bytes = np.frombuffer(r2, dtype=np.uint8)
    else:
        r1_bytes = np.frombuffer(r1.encode(), dtype=np.uint8)
        r2_bytes = np.frombuffer(r2.encode(), dtype=np.uint8)

    # Compute consensus using fast path
    result_bytes = get_outie_consensus_fast(
        r1_bytes, r2_bytes, qv1, qv2, olen)

    # Return same type as input
    if is_bytes:
        return result_bytes.tobytes()
    else:
        return result_bytes.tobytes().decode('ascii')

def get_outie_merged(r1, r2, qv1, qv2):
    '''
    Return outie-merged r1 and r2.
    Internal use only.

    :: r1
       type - string / bytes
       desc - R1 read sequence
    :: r2
       type - string / bytes
       desc - R2 read sequence
    :: qv1
       type - np.array
       desc - R1 Q-Score vector
    :: qv2
       type - np.array
       desc - R2 Q-Score vector
    '''

    # Check if input is bytes (optimized path)
    is_bytes = type(r1) is bytes

    # Fast-path: Check for exact overlap first
    if is_bytes:
        r1_bytes = np.frombuffer(r1, dtype=np.uint8)
        r2_bytes = np.frombuffer(r2, dtype=np.uint8)
    else:
        r1_bytes = np.frombuffer(r1.encode(), dtype=np.uint8)
        r2_bytes = np.frombuffer(r2.encode(), dtype=np.uint8)

    exact_olen = get_exact_outie_overlap(
        r1_bytes=r1_bytes,
        r2_bytes=r2_bytes,
        min_overlap=10)

    # Exact match found - return overlap region
    if exact_olen > 0:
        merged = r1[:exact_olen]
        return merged, 0  # Score 0 = perfect match

    # Slow-path: Use alignment for reads with mismatches
    # Note: edlib expects strings; keep byte-optimized consensus paths but
    # decode once for alignment when upstream streaming uses bytes.
    if is_bytes:
        r1_ed = r1.decode('ascii')
        r2_ed = r2.decode('ascii')
    else:
        r1_ed = r1
        r2_ed = r2

    # Adaptive seed length based on read length
    seed_len = max(10, min(20, len(r2) // 10))

    # Adaptive k based on estimated overlap and mismatch density
    est_overlap = int(0.4 * min(len(r1), len(r2)))
    k_seed = max(4, int(est_overlap * MAX_MISMATCH_DENSITY))

    aln = ed.align(
        query=r2_ed[-seed_len:],
        target=r1_ed,
        mode='HW',
        task='distance',
        k=k_seed)

    # Bad Alignment?
    if aln['editDistance'] == -1:
        # Too many mutations,
        # Skip this Read
        aln.clear() # Memory Management
        del aln     # Memory Management
        return None, float('inf')

    # Extract End Locations
    end_locs = aln['locations']

    # Clear Memory
    aln.clear() # Memory Management
    del aln     # Memory Management

    # Too Many Overlap Locations?
    if len(end_locs) > 10:
        # Likely Lots of Mutation,
        # Skip this Read
        end_locs.clear() # Memory Management
        del end_locs     # Memory Management
        return None, None

    # Find Overlap End Coordinate
    end = end_locs[0][-1] + 1

    # Too Small of an Overlap?
    if end < 10:
        # Not High Confidence,
        # Skip this Read
        end_locs.clear() # Memory Management
        del end_locs     # Memory Management
        return None, None

    # Define Query Slice (edlib expects strings)
    qslice_ed = r1_ed[:end]

    # Density-based k for validation alignment
    k_validate = max(10, int(end * MAX_MISMATCH_DENSITY))

    # Infix Align r1 and r2
    aln = ed.align(
        query=qslice_ed,
        target=r2_ed,
        mode='HW',
        task='path',
        k=k_validate)

    # Bad Alignment?
    if aln['editDistance'] == -1:
        # Too many mutations,
        # Skip this Read
        del qslice_ed  # Memory Management
        aln.clear() # Memory Management
        del aln     # Memory Management
        return None, float('inf')

    # Do we have Indels?
    if 'D' in aln['cigar'] or \
       'I' in aln['cigar']:
        # We don't care about Indels,
        # Skip this Read
        del qslice_ed  # Memory Management
        aln.clear() # Memory Management
        del aln     # Memory Management
        return None, float('inf')

    # Compute Merged Read
    score = aln['editDistance']
    if score == 0:
        merged = r1[:end]
    else:
        merged = get_outie_consensus(
            r1=r1,
            r2=r2,
            qv1=qv1,
            qv2=qv2,
            olen=end)

    # Memory Management
    del qslice_ed
    aln.clear()
    del aln

    # Return Results
    return merged, score

def get_merged_reads(r1, r2, qv1, qv2,
    mergefnhigh,
    mergefnlow):
    '''
    Return assembled read from r1 and r2 as tuple.
    Internal use only.

    :: r1
       type - string / bytes
       desc - R1 read sequence
    :: r2
       type - string / bytes
       desc - R2 read sequence
    :: qv1
       type - np.array
       desc - R1 Q-Score vector
    :: qv2
       type - np.array
       desc - R2 Q-Score vector
    :: mergefnhigh
       type - function
       desc - high priority merging
              function
    :: mergefnlow
       type - function
       desc - low priority merging
              function

    Returns (merged_read, None) on success, None on failure.
    '''

    # Do we not need to assemble?
    if r2 is None:
        # Nothing to Merge
        if type(r1) is bytes:
            r1 = r1.decode('ascii')
        return (r1, None)

    # High Priority Assembly
    (high_assembled,
    high_score) = mergefnhigh(
        r1=r1, r2=r2, qv1=qv1, qv2=qv2)

    # Is High Priority Enough?
    if high_score == 0:
        if type(high_assembled) is bytes:
            high_assembled = high_assembled.decode('ascii')
        return (high_assembled, None)

    # Low Priority Assembly
    (low_assembled,
    low_score) = mergefnlow(
        r1=r1, r2=r2, qv1=qv1, qv2=qv2)

    # Only High Priority Assembly was Successful
    if not high_assembled is None and \
       low_assembled is None:
       if type(high_assembled) is bytes:
           high_assembled = high_assembled.decode('ascii')
       return (high_assembled, None)

    # Only Low Priority Assembly was Successful
    if not low_assembled is None and \
       high_assembled is None:
       if type(low_assembled) is bytes:
           low_assembled = low_assembled.decode('ascii')
       return (low_assembled, None)

    # Nothing Assembled?
    if high_assembled is None and \
       low_assembled is None:
       return None

    # Both Equally Assembled?
    if high_score == low_score:
        return None

    # Return Assembly with Best Score
    if high_score < low_score:
        if type(high_assembled) is bytes:
            high_assembled = high_assembled.decode('ascii')
        return (high_assembled, None)
    if type(low_assembled) is bytes:
        low_assembled = low_assembled.decode('ascii')
    return (low_assembled, None)

def get_skipcount(
    previousreads,
    coreid,
    ncores):
    '''
    Return the number of reads to be
    skipped for current core.
    Internal use only.

    :: previousreads
       type - integer
       desc - total number of read pairs scanned
              previously from r1file and r2file
    :: coreid
       type - integer
       desc - current core integer id
    :: ncores
       type - integer
       desc - total number of readers
              concurrently initiated
    '''

    # Book-keeping
    n = coreid
    m = ncores
    p = previousreads
    d = ncores - coreid - 1

    # Compute and Return Results
    if p <= 0:
        return 0
    return ((p-1) * m) + (n + 1) + d

def stream_processed_fastq(
    r1file,
    r2file,
    r1type,
    r2type,
    r1length,
    r2length,
    r1qual,
    r2qual,
    packtype,
    assemblyparams,
    r1truncfile,
    r2truncfile,
    previousreads,
    scannedreads,
    ambiguousreads,
    shortreads,
    survivedreads,
    verbiagetarget,
    coreid,
    ncores,
    memlimit,
    launchtime,
    restart,
    shutdown,
    liner):
    '''
    Scan, process and stream reads from
    FastQ files. Internal use only.

    :: r1file
       type - string
       desc - path to FastQ file storing
              R1 reads
    :: r2file
       type - string / None
       desc - path to FastQ file storing
              R2 reads
    :: r1type
       type - integer
       desc - R1 read type identifier
              0 = Forward Reads
              1 = Reverse Reads
    :: r2type
       type - integer / None
       desc - R2 read type identifier
              0 = Forward Reads
              1 = Reverse Reads
    :: r1length
       type - integer
       desc - minimum required R1 read
              length for packing
    :: r2length
       type - integer / None
       desc - minimum required R2 read
              length for packing
    :: r1qual
       type - integer
       desc - minimum required R1 read
              Q-Score for packing
    :: r2qual
       type - integer / None
       desc - minimum required R2 read
              Q-Score for packing
    :: packtype
       type - integer
       desc - packing operation identifier
              0 = concatenated / joined reads
              1 = assembled / merged reads
    :: assemblyparams
       type - dict / None
       desc - dictionary storing read overlap
              merging parameters
    :: r1truncfile
       type - Event
       desc - multiprocessing Event triggered
              on truncated r1file
    :: r2truncfile
       type - Event
       desc - multiprocessing Event triggered
              on truncated r2file
    :: filecomplete
       type - mp.Event
       desc - multiprocessing Event triggered
              when r1file and r2file fully
              scanned
    :: insofarreads
       type - integer
       desc - total number of read pairs scanned
              so far from r1file and r2file
    :: previousreads
       type - integer
       desc - total number of read pairs scanned
              previously from r1file and r2file
    :: scannedreads
       type - SafeCounter
       desc - total number of existing read
              pairs in r1file and r2file
    :: ambiguousreads
       type - SafeCounter
       desc - total number of rejected read
              pairs in r1file and r2file due
              to ambiguous bases in the reads
    :: shortreads
       type - SafeCounter
       desc - total number of rejected read
              pairs in r1file and r2file due
              to shorter than expected lengths
    :: survivedreads
       type - SafeCounter
       desc - total number of filtered read
              pairs in r1file and r2file
    :: verbiagetarget
       type - integer
       desc - total number of reads to scan per
              batch before updating scannedreads
    :: coreid
       type - integer
       desc - current core integer id
    :: ncores
       type - integer
       desc - total number of readers
              concurrently initiated
    :: memlimit
       type - nu.Real
       desc - total amount of memory
              allowed per core
    :: launchtime
       type - time
       desc - initial launch timestamp
    :: restart
       type - mp.Event
       desc - multiprocessing Event
              triggered to restart
    :: shutdown
       type - mp.Event
       desc - multiprocessing Event
              triggered to shutdown
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Start Timer
    t0 = launchtime

    # Book-keeping Variables
    verbiagereach    = 0
    c_previousreads  = previousreads.value()
    c_scannedreads   = 0
    c_ambiguousreads = 0
    c_shortreads     = 0
    c_survivedreads  = 0
    clen             = ut.get_printlen(value=ncores)
    truncable        = not r2file is None
    batchreach       = 0
    batchsize        = 499979 # Largest Prime < 500k
    skipcount        = get_skipcount(
        previousreads=c_previousreads,
        coreid=coreid,
        ncores=ncores)

    # Setup R1 Streamer (use bytes for faster processing)
    reader1 = ut.stream_fastq_engine(
        filepath=r1file,
        invert=r1type,
        qualvec=r1qual > 0,
        readbytes=True,
        skipcount=skipcount,
        coreid=coreid,
        ncores=ncores,
        fastqid=1,
        liner=liner)

    # Setup R2 Streamer (use bytes for faster processing)
    if not r2file is None:
        reader2 = ut.stream_fastq_engine(
            filepath=r2file,
            invert=r2type,
            qualvec=r2qual > 0,
            readbytes=True,
            skipcount=skipcount,
            coreid=coreid,
            ncores=ncores,
            fastqid=2,
            liner=liner)
    else:
        reader2 = ix.repeat((None, None))

    # Setup Merging Functions
    if packtype:
        mergefnhigh = assemblyparams['mergefnhigh']
        mergefnlow  = assemblyparams['mergefnlow']

    # Read Scanning Loop
    while True:

        # Concurrent Reading
        r1,qv1 = next(reader1)
        r2,qv2 = next(reader2)

        # # Simulate Early Termination
        # if c_scannedreads > 5 * 10**6:
        #     shutdown.set()
        #     break

        # Exhausted File Pair
        if r1 is None and \
           r2 is None:
            shutdown.set()
            break

        # Can Read Files be Truncated?
        if truncable:

            # Truncated R1 File
            if       r1 is None and \
                 not r2 is None:
                 r1truncfile.set()
                 shutdown.set()
                 break

            # Truncated R2 File
            elif     r2 is None and \
                 not r1 is None:
                 r2truncfile.set()
                 shutdown.set()
                 break

        # Apply Ambiguity, Length, and Quality Filters
        filterpass = True

        # Filter I on R1 (works with both bytes and strings)
        if   r1.count(ut.N_BYTE) > 10:
            c_ambiguousreads += 1 # Update Ambiguous Read Count based on R1
            filterpass = False    # This Read will be Skipped
        elif len(r1) < r1length:
            c_shortreads     += 1 # Update Short Read Count based on R1
            filterpass = False    # This Read will be Skipped
        elif len(set(r1)) < 4:
            c_ambiguousreads += 1 # Update Ambiguous Read Count based on R1
            filterpass = False    # This Read will be Skipped

        # Filter I on R2 (works with both bytes and strings)
        if filterpass and (not r2 is None):
            if   r2.count(ut.N_BYTE) > 10:
                c_ambiguousreads += 1 # Update Ambiguous Read Count based on R2
                filterpass = False    # This Read will be Skipped
            elif len(r2) < r2length:
                c_shortreads     += 1 # Update Short Read Count based on R2
                filterpass = False    # This Read will be Skipped
            elif len(set(r2)) < 4:
                c_ambiguousreads += 1 # Update Ambiguous Read Count based on R2
                filterpass = False    # This Read will be Skipped

        # Filter II
        if filterpass:

            # Filter II on R1
            if  (r1qual > 0) and \
                (not qualpass(
                    qvec=qv1,
                    threshold=r1qual)):
                c_ambiguousreads += 1 # Update Ambiguous Read Count based on R1
                filterpass = False    # This Read will be Skipped

            # Filter II on R2
            elif (not r2 is None):
                if  (r2qual > 0) and \
                    (not qualpass(
                        qvec=qv2,
                        threshold=r2qual)):
                        c_ambiguousreads += 1 # Update Ambiguous Read Count based on R2
                        filterpass = False    # This Read will be Skipped

        # Read(-Pair) Operation
        if filterpass:

            # Assemble R1 and R2
            if packtype:
                read = get_merged_reads(
                    r1=r1, r2=r2, qv1=qv1, qv2=qv2,
                    mergefnhigh=mergefnhigh,
                    mergefnlow=mergefnlow)

            # Join R1 and R2
            else:
                read = get_concatenated_reads(
                    r1=r1, r2=r2)

            # Do we stream this Read?
            if not read is None:
                # High Quality Read
                c_survivedreads  += 1 # Update Survived Read Count
                yield read            # This Read will be Stored
                del   read            # Memory Management
            else:
                c_ambiguousreads += 1 # Update Ambiguous Read Count

        # Update Book-keeping
        c_scannedreads += 1
        verbiagereach  += 1
        batchreach     += 1

        # Time to show updates?
        if verbiagereach >= verbiagetarget:

            # Show Updates
            liner.send(
                ' Core {:{},d}: Scanned {:,d} Reads in {:.2f} sec'.format(
                    coreid,
                    clen,
                    c_previousreads + c_scannedreads,
                    tt.time()-t0))

            # Update Book-keeping
            verbiagereach = 0

        # Did we complete a batch?
        if batchreach == batchsize:
            # Need to Restart?
            if ut.needs_restart(
                memlimit=memlimit):
                restart.set() # Enable Restart
                break # Release, your Memory Real Estate, biatch!
            else:
                batchreach = 0 # Reset Counter

        # # Simulate Asynchornous Restarts
        # if batchreach == ((coreid+1) * 10000):
        #     restart.set()
        #     break

    # Final Book-keeping Update
    previousreads.increment(incr=c_scannedreads)
    scannedreads.increment(incr=c_scannedreads)
    ambiguousreads.increment(incr=c_ambiguousreads)
    shortreads.increment(incr=c_shortreads)
    survivedreads.increment(incr=c_survivedreads)

    # Show Final Updates?
    if survivedreads.value() == 0:
        # Nothing Survived, Important
        # to Show Reads were Scanned!
        # Show Updates
        liner.send(
            ' Core {:{},d}: Scanned {:,d} Reads in {:.2f} sec\n'.format(
                coreid,
                clen,
                c_previousreads + c_scannedreads,
                tt.time()-t0))

    # Shutdown!
    if shutdown.is_set():
        # Show Updates
        liner.send(' Core {:{},d}: Shutting Down\n'.format(
            coreid,
            clen))
    # Restart, We Must!
    elif restart.is_set():
        # Show Updates
        liner.send(' Core {:{},d}: Restarting ...\n'.format(
            coreid,
            clen))

def pack_aggregator(
    metaqueue,
    packqueue,
    packedreads,
    packsbuilt,
    liner):
    '''
    Aggregate all meta read packs referenced
    in metaqueue by merging common reads in
    meta read packs, and save them as usual
    read packs and their paths in packqueue.

    : metaqueue
       type - SafeQueue
       desc - queue storing meta read pack
              file paths once saved into
              packdir
    :: packqueue
       type - SafeQueue
       desc - queue storing read pack file
              paths once saved into packdir
    :: packedreads
       type - SafeCounter
       desc - total number of uniquely
              packed reads in batches
    :: packsbuilt
       type - SafeCounter
       desc - total number of read packs
              built from r1file and r2file
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping
    aggcount = 0

    # Do we have Meta Packs to aggregate?
    if metaqueue.empty():
        return None

    # Fetch First Meta Pack
    mpath = metaqueue.get()

    # Meta Pack hasn't Materialized?
    while not os.path.isfile(mpath):
        tt.sleep(0)
        continue

    # Load First Meta Pack
    mpack = ut.loadmeta(
        mfile=mpath)
    aggcount += 1

    # Merge Remaining Meta Packs
    while not metaqueue.empty():

        # Release Control
        tt.sleep(0)

        # Fetch Another Meta Pack
        cpath = metaqueue.get()

        # Meta Pack hasn't Materialized?
        while not os.path.isfile(cpath):
            tt.sleep(0)
            continue

        # Load Meta Pack
        cpack = ut.loadmeta(
            mfile=cpath)
        aggcount += 1

        # Show Update
        cpackid = ut.removestarfix(
            string=cpath.split('/')[-1],
            fix='.meta',
            loc=1)
        liner.send(
            ' Aggregating {} in Progress'.format(
                cpackid))

        # Reduced Meta Pack
        npack = cx.Counter()
        npath = ut.removestarfix(
            string=cpath,
            fix='.meta',
            loc=1)

        # Merge Fetched Meta Pack
        while cpack:

            # Pop Reads and their Count
            reads,count = cpack.popitem()

            # Reads Common to First Meta Pack
            if reads in mpack:

                # Add Count to First Meta Pack
                mpack[reads] += count

                # Subtract from Packed Reads
                packedreads.decrement()

            # Uncommon Read
            else:

                # Store in Reduced Meta Pack
                npack[reads] = count

        # Unable to fully reduce Fetched Meta Pack?
        if npack:

            # Dump Reduced Meta Pack
            ut.savepack(
                pobj=npack,
                filepath=npath)

            # Enqueue Pack Path
            packqueue.put(npath)

        # We fully reduced Fetched Meta Pack
        else:

            # One less Read Pack to Count
            packsbuilt.decrement()

    # Show Update
    mpackid = ut.removestarfix(
        string=mpath.split('/')[-1],
        fix='.meta',
        loc=1)
    liner.send(
        ' Aggregating {} in Progress'.format(
            mpackid))

    # First Meta Pack Path
    mpath = ut.removestarfix(
        string=mpath,
        fix='.meta',
        loc=1)

    # Dump First Meta Pack
    ut.savepack(
        pobj=mpack,
        filepath=mpath)

    # Clear Line
    liner.send(' ')

    # Enqueue Pack Path and Done!
    packqueue.put(mpath)
    packqueue.put(None)

    # Final Update
    liner.send(
        ' Meta Aggregated: {:,} Read Pack(s)\n'.format(aggcount))

def pack_engine(
    r1file,
    r2file,
    r1type,
    r2type,
    r1length,
    r2length,
    r1qual,
    r2qual,
    packdir,
    metaqueue,
    packqueue,
    packtype,
    assemblyparams,
    packsize,
    r1truncfile,
    r2truncfile,
    previousreads,
    scannedreads,
    ambiguousreads,
    shortreads,
    survivedreads,
    packedreads,
    packsbuilt,
    coreid,
    batchid,
    ncores,
    nactive,
    memlimit,
    launchtime,
    restart,
    shutdown,
    liner):
    '''
    Process r1file and r2file and pack read
    pairs by frequency. Internal use only.

    Scan, process and stream reads from
    FastQ files. Internal use only.

    :: r1file
       type - string
       desc - path to FastQ file storing
              R1 reads
    :: r2file
       type - string / None
       desc - path to FastQ file storing
              R2 reads
    :: r1type
       type - integer
       desc - R1 read type identifier
              0 = Forward Reads
              1 = Reverse Reads
    :: r2type
       type - integer / None
       desc - R2 read type identifier
              0 = Forward Reads
              1 = Reverse Reads
    :: r1length
       type - integer
       desc - minimum required R1 read
              length for packing
    :: r2length
       type - integer / None
       desc - minimum required R2 read
              length for packing
    :: r1qual
       type - integer
       desc - minimum required R1 read
              Q-Score for packing
    :: r2qual
       type - integer / None
       desc - minimum required R2 read
              Q-Score for packing
    :: packdir
       type - string
       desc - path to directory temporarily
              storing packed reads
    :: metaqueue
       type - SafeQueue
       desc - queue storing meta read pack
              file paths once saved into
              packdir
    :: packqueue
       type - SafeQueue
       desc - queue storing read pack file
              paths once saved into packdir
    :: packtype
       type - integer
       desc - packing operation identifier
              0 = concatenated / joined reads
              1 = assembled / merged reads
    :: assemblyparams
       type - dict / None
       desc - dictionary storing read overlap
              merging parameters
    :: packsize
       type - integer
       desc - maximum number of reads (in
              millions) stored per pack
    :: r1truncfile
       type - mp.Event
       desc - multiprocessing Event triggered
              on truncated r1file
    :: r2truncfile
       type - mp.Event
       desc - multiprocessing Event triggered
              on truncated r2file
    :: previousreads
       type - SafeCounter
       desc - total number of read pairs scanned
              previously from r1file and r2file
              by current core
    :: scannedreads
       type - SafeCounter
       desc - total number of read pairs scanned
              from r1file and r2file
    :: ambiguousreads
       type - SafeCounter
       desc - total number of rejected read
              pairs in r1file and r2file due
              to ambiguous bases in the reads
    :: shortreads
       type - SafeCounter
       desc - total number of rejected read
              pairs in r1file and r2file due
              to shorter than expected lengths
    :: survivedreads
       type - SafeCounter
       desc - total number of filtered read
              pairs in r1file and r2file
    :: packedreads
       type - SafeCounter
       desc - total number of uniquely
              packed reads in batches
    :: packsbuilt
       type - SafeCounter
       desc - total number of read packs
              built from r1file and r2file
    :: coreid
       type - integer
       desc - current core integer id
    :: batchid
       type - integer
       desc - current batch integer id
    :: ncores
       type - integer
       desc - total number of packers
              concurrently initiated
    :: nactive
       type - SafeCounter
       desc - total number of packers
              concurrently active
    :: memlimit
       type - nu.Real
       desc - total amount of memory
              allowed per core
    :: launchtime
       type - time
       desc - initial launch timestamp
    :: restart
       type - mp.Event
       desc - multiprocessing event when
              process needs to restart
    :: shutdown
       type - mp.Event
       desc - multiprocessing event when
              process needs to shutdown
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Book-keeping Variables
    t0 = tt.time()
    truncated = False

    c_packsbuilt  = 0
    c_packedreads = 0

    min_dump_reach  = 0
    min_dump_target = packsize // 2

    max_dump_reach  = 0
    max_dump_target = packsize

    # Setup Verbiage Variables
    if packtype:
        verbiagefactor = 0.2
    else:
        verbiagefactor = 1
    verbiagetarget = rn.randint(
        *map(round, (min_dump_target * 0.080 * verbiagefactor,
                     min_dump_target * 0.120 * verbiagefactor)))

    clen = ut.get_printlen(value=ncores)
    plen = ut.get_printlen(value=packsize)

    # Current and Meta Pack Storage
    cpack = cx.Counter()
    mpack = {}
    mpath = None

    # Read Packing Loop
    for reads in stream_processed_fastq(
        r1file=r1file,
        r2file=r2file,
        r1type=r1type,
        r2type=r2type,
        r1length=r1length,
        r2length=r2length,
        r1qual=r1qual,
        r2qual=r2qual,
        packtype=packtype,
        assemblyparams=assemblyparams,
        r1truncfile=r1truncfile,
        r2truncfile=r2truncfile,
        previousreads=previousreads,
        scannedreads=scannedreads,
        ambiguousreads=ambiguousreads,
        shortreads=shortreads,
        survivedreads=survivedreads,
        verbiagetarget=verbiagetarget,
        coreid=coreid,
        ncores=ncores,
        memlimit=memlimit,
        launchtime=launchtime,
        restart=restart,
        shutdown=shutdown,
        liner=liner):

        # Truncated Files?
        if  r1truncfile.is_set() or \
            r2truncfile.is_set():
            truncated = True
            shutdown.set()
            break

        # Pack Read
        if reads in mpack:
            mpack[reads] += 1
        else:
            cpack[reads] += 1

        # Update Book-keeping
        min_dump_reach  = len(cpack)
        max_dump_reach += 1

        # Time to dump a pack?
        # Dump packs per packsize / 2 (e.g. 1.5M) unique
        # entries or packsize (3.0M) total reads processed,
        # whichever is reached earlier
        if  max_dump_reach >= max_dump_target or \
            min_dump_reach >= min_dump_target:

            # Show Updates
            liner.send(
                ' Core {:{},d}: Built Pack {}.{}.{} w/ {:{},d} Reads in {:.2f} sec\n'.format(
                    coreid,
                    clen,
                    coreid,
                    c_packsbuilt,
                    batchid,
                    len(cpack),
                    plen,
                    tt.time()-t0))

            # Define Pack Path
            cpath = '{}/{}.{}.{}.pack'.format(
                packdir,
                coreid,
                c_packsbuilt,
                batchid)

            # Update Book-keeping
            t0 = tt.time()
            c_packedreads  += len(cpack)
            c_packsbuilt   += 1
            max_dump_reach  = 0
            min_dump_reach  = 0

            # Dump Current Pack
            if mpack:

                # Persist Current Pack to Disk
                ut.savepack(
                    pobj=cpack,
                    filepath=cpath)

                # Enqueue Current Pack Path
                packqueue.put(cpath)

                # Release Control
                tt.sleep(0)

                # Cleanup Current Pack
                cpack.clear()
                del cpack
                cpack = cx.Counter()
                ut.free_mem()

            # Setup Meta Pack
            else:
                mpack = cpack
                mpath = cpath
                cpack = cx.Counter()
                ut.free_mem()

    # Final Dumping
    if not truncated:

        # Dump Final Pack?
        if cpack:

            # Show Updates
            liner.send(
                ' Core {:{},d}: Built Pack {}.{}.{} w/ {:{},d} Reads in {:05.2f} sec\n'.format(
                    coreid,
                    clen,
                    coreid,
                    c_packsbuilt,
                    batchid,
                    len(cpack),
                    plen,
                    tt.time()-t0))

            # Dump Read Pack Appropriately

            # Meta pack
            if c_packsbuilt == 0:

                # Define Pack Path
                cpath = '{}/{}.{}.{}.pack.meta'.format(
                    packdir,
                    coreid,
                    c_packsbuilt,
                    batchid,)

                # Dump Meta Pack
                ut.savemeta(
                    pobj=cpack,
                    filepath=cpath)

                # Enqueue Meta Pack Path
                metaqueue.put(cpath)

                # Release Control
                tt.sleep(0)

            # Non-Meta Pack
            else:

                # Define Pack Path
                cpath = '{}/{}.{}.{}.pack'.format(
                    packdir,
                    coreid,
                    c_packsbuilt,
                    batchid)

                # Dump Pack
                ut.savepack(
                    pobj=cpack,
                    filepath=cpath)

                # Enqueue Pack Path
                packqueue.put(cpath)

                # Release Control
                tt.sleep(0)

            # Update Book-keeping
            c_packedreads += len(cpack)
            c_packsbuilt  += 1

            # Cleanup Final Pack
            cpack.clear()
            del cpack
            ut.free_mem()

        # Dump Meta Pack?
        if mpack:

            # Define Pack Path
            mpath = '{}.meta'.format(
                mpath)

            # Dump Meta Pack
            ut.savemeta(
                pobj=mpack,
                filepath=mpath)

            # Enqueue Meta Pack Path
            metaqueue.put(mpath)

            # Release Control
            tt.sleep(0)

            # Cleanup Meta Pack
            mpack.clear()
            del mpack
            ut.free_mem()

        # Update Read Packing Book-keeping
        packedreads.increment(incr=c_packedreads)
        packsbuilt.increment(incr=c_packsbuilt)

    # Packing Completed
    nactive.decrement()
    if shutdown.is_set():
        packqueue.put(None)

    # Release Control
    tt.sleep(0)
