import time as tt

import itertools   as ix
import collections as cx

import numpy   as np
import numba   as nb
import edlib   as ed
import pyfastx as pf

import utils as ut


# Engine Helper Functions

@nb.njit
def qual2qvec(qualstring):
    '''
    Return Q-Score vector.
    Internal use only.

    :: qualstring
       type - string
       desc - quality string
    '''

    qvec = np.zeros(len(qualstring))
    for i in range(len(qualstring)):
        qvec[i] = ord(qualstring[i]) - 33.
    return qvec

@nb.njit
def qualpass_merge(qvec, threshold):
    '''
    Return Q-Score vector.
    Internal use only.

    :: qvec
       type - np.array
       desc - Q-Score vector
    '''

    return np.mean(qvec) >= threshold

@nb.njit
def qualpass_concat(qualstring, threshold):
    '''
    Determine if qualstring mean is above
    threshold. Internal use only.

    :: qualstring
       type - string
       desc - quality string
    :: threshold
       type - integer
       desc - quality threshold
    '''

    tot = 0.
    for qc in qualstring:
        tot += ord(qc) - 33.
    return (tot / len(qualstring)) >= threshold

def get_concatenated_reads(r1, r2):
    '''
    Return r1+r2 joined with '-' delimiters.
    Internal use only.

    :: r1
       type - string
       desc - R1 read sequence
    :: r2
       type - string
       desc - R2 read sequence
    '''

    return r1 if r2 is None else r1 + '-----' + r2

@nb.njit
def get_innie_consensus(r1, r2, qv1, qv2, olen):
    '''
    Return consensus innie merged read.
    Internal use only.

    :: r1
       type - string
       desc - R1 read sequence
    :: r2
       type - string
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

    # Book-keeping
    i = len(r1)-olen
    j = 0

    # Split Merged Parts
    ml = nb.typed.List(r1)
    ml.extend(r2[olen:])

    # Compute Merged Region based
    # on Q-Score Values
    while i < len(r1):
        if qv1[i] < qv2[j]:
            ml[i] = r2[j]
        i += 1
        j += 1

    # Return Results
    return ''.join(ml)

def get_innie_merged(r1, r2, qv1, qv2):
    '''
    Return innie-merged r1 and r2.
    Internal use only.

    :: r1
       type - string
       desc - R1 read sequence
    :: r2
       type - string
       desc - R2 read sequence
    :: qv1
       type - np.array
       desc - R1 Q-Score vector
    :: qv2
       type - np.array
       desc - R2 Q-Score vector
    '''

    # Find Overlap Locations
    end_locs = ed.align(
        query=r1,
        target=r2,
        mode='SHW',
        task='distance')['locations']

    # Too Many Overlap Locations?
    if len(end_locs) > 10:
        # Likely Lots of Mutation,
        # Skip this Read
        return None, None

    # Find Overlap End Coordinate
    end = end_locs[-1][-1] + 1

    # Too Small of an Overlap?
    if end < 10:
        # Not High Confidence,
        # Skip this Read
        return None, None

    # Infix Align r1 and r2
    aln = ed.align(
        query=r2[:end],
        target=r1,
        mode='HW',
        task='path',
        k=10)

    # Bad Alignment?
    if aln['editDistance'] == -1:
        # Too many mutations,
        # Skip this Read
        return None, float('inf')

    # Do we have Indels?
    cigar = set(aln['cigar'])
    if 'D' in cigar or \
       'I' in cigar:
        # We don't care about Indels,
        # Skip this Read
        return None, float('inf')

    # Compute Merged Read
    score = aln['editDistance']
    if aln['editDistance'] == 0:
        merged = r1 + r2[end:]
    else:
        merged = get_innie_consensus(
            r1=r1,
            r2=r2,
            qv1=qv1,
            qv2=qv2,
            olen=end)

    # Return Results
    return merged, score

@nb.njit
def get_outie_consensus(r1, r2, qv1, qv2, olen):
    '''
    Return consensus outie merged read.
    Internal use only.

    :: r1
       type - string
       desc - R1 read sequence
    :: r2
       type - string
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

    # Book-keeping
    i = 0
    j = len(r2)-olen

    # Split Merged Parts
    ml = nb.typed.List(r1[:olen])

    # Compute Merged Region based
    # on Q-Score Values
    while i < olen:
        if qv1[i] < qv2[j]:
            ml[i] = r2[j]
        i += 1
        j += 1

    # Return Results
    return ''.join(ml)

def get_outie_merged(r1, r2, qv1, qv2):
    '''
    Return outie-merged r1 and r2.
    Internal use only.

    :: r1
       type - string
       desc - R1 read sequence
    :: r2
       type - string
       desc - R2 read sequence
    :: qv1
       type - np.array
       desc - R1 Q-Score vector
    :: qv2
       type - np.array
       desc - R2 Q-Score vector
    '''

    # Find Overlap Locations
    end_locs = ed.align(
        query=r2,
        target=r1,
        mode='SHW',
        task='distance')['locations']

    # Too Many Overlap Locations?
    if len(end_locs) > 10:
        # Likely Lots of Mutation,
        # Skip this Read
        return None, None

    # Find Overlap End Coordinate
    end = end_locs[-1][-1] + 1

    # Too Small of an Overlap?
    if end < 10:
        # Not High Confidence,
        # Skip this Read
        return None, None

    # Infix Align r1 and r2
    aln = ed.align(
        query=r1[:end],
        target=r2,
        mode='HW',
        task='path',
        k=10)

    # Bad Alignment?
    if aln['editDistance'] == -1:
        # Too many mutations,
        # Skip this Read
        return None, float('inf')

    # Do we have Indels?
    cigar = set(aln['cigar'])
    if 'D' in cigar or \
       'I' in cigar:
        # We don't care about Indels,
        # Skip this Read
        return None, float('inf')

    # Compute Merged Read
    score = aln['editDistance']
    if aln['editDistance'] == 0:
        merged = r1[:end]
    else:
        merged = get_outie_consensus(
            r1=r1,
            r2=r2,
            qv1=qv1,
            qv2=qv2,
            olen=end)

    # Return Results
    return merged, score

def get_merged_reads(r1, r2, qv1, qv2):
    '''
    Return assembled read from r1 and r2.
    Internal use only.

    :: r1
       type - string
       desc - R1 read sequence
    :: r2
       type - string
       desc - R2 read sequence
    :: qv1
       type - np.array
       desc - R1 Q-Score vector
    :: qv2
       type - np.array
       desc - R2 Q-Score vector
    '''

    # Is r2 None?
    if r2 is None:
        # Nothing to Merge
        return r1

    # Innie Assembly
    (innie_aseembled,
    innie_score) = get_innie_merged(
        r1=r1, r2=r2, qv1=qv1, qv2=qv2)

    # Is Innie Enough?
    if innie_score == 0:
        return innie_aseembled

    # Outie Assembly
    (outie_aseembled,
    outie_score) = get_outie_merged(
        r1=r1, r2=r2, qv1=qv1, qv2=qv2)

    # Only Innie Assembly was Successful
    if not innie_aseembled is None and \
       outie_aseembled is None:
       return innie_aseembled

    # Only Outie Assembly was Successful
    if not outie_aseembled is None and \
       innie_aseembled is None:
       return outie_aseembled

    # Nothing Assembled?
    if innie_aseembled is None and \
       outie_aseembled is None:
       return None

    # Both Equally Assembled?
    if innie_score == outie_score:
        return None

    # Return Assembly with Best Score
    if innie_score < outie_score:
        return innie_aseembled
    return outie_aseembled

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
    r1truncfile,
    r2truncfile,
    scannedreads,
    ambiguousreads,
    shortreads,
    survivedreads,
    verbagetarget,
    coreid,
    ncores,
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
    :: r1truncfile
       type - Event
       desc - multiprocessing Event triggered
              on truncated r1file
    :: r2truncfile
       type - Event
       desc - multiprocessing Event triggered
              on truncated r2file
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
    :: verbagetarget
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
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Start Timer
    t0 = tt.time()

    # Book-keeping Variables
    verbagereach     = 0
    c_scannedreads   = 0
    c_ambiguousreads = 0
    c_shortreads     = 0
    c_survivedreads  = 0
    clen             = len(str(ncores))
    truncable        = not r2file is None

    # Setup R1 Streamer
    reader1 = ut.stream_fastq_engine(
        filepath=r1file,
        coreid=coreid,
        ncores=ncores)

    # Setup R2 Streamer
    if not r2file is None:
        reader2 = ut.stream_fastq_engine(
            filepath=r2file,
            coreid=coreid,
            ncores=ncores)
    else:
        reader2 = ix.repeat((None, None))

    # Read Scanning Loop
    while True:

        # Concurrent Reading
        r1,q1 = next(reader1)
        r2,q2 = next(reader2)

        # Exhausted File Pair
        if r1 is None and \
           r2 is None:
            break

        # Can Read Files be Truncated?
        if truncable:

            # Truncated R1 File
            if       r1 is None and \
                 not r2 is None:
                 r1truncfile.set()
                 break

            # Truncated R2 File
            elif     r2 is None and \
                 not r1 is None:
                 r2truncfile.set()
                 break

        # Correct Orientation for Reads
        if r1type:
            r1 = ut.get_revcomp(r1)
        if r2type:
            r2 = ut.get_revcomp(r2)

        # Apply Ambiguity, Length, and Quality Filters
        filterpass = True

        # Filter I on R1
        if   'N' in r1:
            c_ambiguousreads += 1 # Update Ambiguous Read Count based on R1
            filterpass = False    # This Read will be Skipped
        elif len(r1) < r1length:
            c_shortreads     += 1 # Update Short Read Count based on R1
            filterpass = False    # This Read will be Skipped
        elif len(set(r1)) < 4:
            c_ambiguousreads += 1 # Update Ambiguous Read Count based on R1
            filterpass = False    # This Read will be Skipped

        # Filter I on R2
        if filterpass and (not r2 is None):
            if   'N' in r2:
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

            # Merging Operation Quality Check
            if packtype:

                # Get Quality Vectors for Merging
                qv1 = qual2qvec(qualstring=q1)
                # R2 Quality Vector
                if r2 is None:
                    qv2 = None
                else:
                    qv2 = qual2qvec(qualstring=q2)

                # Filter II on R1
                if   not qualpass_merge(
                        qvec=qv1,
                        threshold=r1qual):
                    c_ambiguousreads += 1 # Update Ambiguous Read Count based on R1
                    filterpass = False    # This Read will be Skipped

                # Filter II on R2
                elif (not r2 is None) and \
                     (not qualpass_merge(
                        qvec=qv2,
                        threshold=r2qual)):
                    c_ambiguousreads += 1 # Update Ambiguous Read Count based on R2
                    filterpass = False    # This Read will be Skipped

            # Concatenation Operation Quality Check
            else:

                # Filter II on R1
                if   not qualpass_concat(
                        qualstring=q1,
                        threshold=r1qual):
                    c_ambiguousreads += 1 # Update Ambiguous Read Count based on R1
                    filterpass = False    # This Read will be Skipped

                # Filter II on R2
                elif (not r2 is None) and \
                     (not qualpass_concat(
                        qualstring=q2,
                        threshold=r2qual)):
                    c_ambiguousreads += 1 # Update Ambiguous Read Count based on R2
                    filterpass = False    # This Read will be Skipped

        # Read(-Pair) Operation
        if filterpass:

            # Assemble R1 and R2
            if packtype:
                read = get_merged_reads(
                    r1=r1, r2=r2, qv1=qv1, qv2=qv2)

            # Join R1 and R2
            else:
                read = get_concatenated_reads(
                    r1=r1, r2=r2)

            # Do we stream this Read?
            if not read is None:
                # High Quality Read
                c_survivedreads  += 1 # Update Survived Read Count
                yield read            # This Read will be Stored
            else:
                c_ambiguousreads += 1 # Update Ambiguous Read Count

        # Update Book-keeping
        c_scannedreads += 1
        verbagereach   += 1

        # Time to show updates?
        if verbagereach >= verbagetarget:

            # Show Updates
            liner.send(
                ' Core {:{},d}: Scanned {:,d} Reads in {:.2f} sec'.format(
                    coreid,
                    clen,
                    c_scannedreads,
                    tt.time()-t0))

            # Update Book-keeping
            verbagereach = 0

    # Final Book-keeping Update
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
                c_scannedreads,
                tt.time()-t0))

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
       type - SimpleQueue
       desc - queue storing meta read pack
              file paths once saved into
              packdir
    :: packqueue
       type - SimpleQueue
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

    # Do we have Meta Packs to aggregate?
    if metaqueue.empty():
        return None

    # Fetch First Meta Pack
    mpath = metaqueue.get()
    mpack = ut.loadmeta(
        filepath=mpath)

    # Merge Remaining Meta Packs
    while not metaqueue.empty():

        # Fetch Another Meta Pack
        cpath = metaqueue.get()
        cpack = ut.loadmeta(
            filepath=cpath)

        # Show Update
        cpackid = cpath.split('/')[-1].rstrip(
            '.meta')
        liner.send(
            ' Aggregating {} in Progress'.format(
                cpackid))

        # Reduced Meta Pack
        npack = cx.Counter()
        npath = cpath.rstrip('.meta')

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

    # First Meta Pack Path
    mpath = mpath.rstrip('.meta')

    # Dump First Meta Pack
    ut.savepack(
        pobj=mpack,
        filepath=mpath)

    # Enqueue Pack Path
    packqueue.put(mpath)
