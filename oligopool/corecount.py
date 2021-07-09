import time as tt

import collections as cx
import zipfile     as zf

import numpy  as np
import numba  as nb
import edlib  as ed

import utils as ut


# Parser and Setup Functions

def pack_loader(
    packfile,
    packqueue,
    enqueuecomplete,
    ncores,
    liner):
    '''
    Enqueue all read pack names from
    packfile. Internal use only.

    :: packfile
       type - string
       desc - filename storing read packs
    :: packqueue
       type - SafeQueue
       desc - queue storing read pack file
              paths
    :: enqueuecomplete
       type - mp.Event
       desc - multiprocessing Event set
              when all pack names stored
              in pack queue
    :: ncores
       type - integer
       desc - total number of counting
              processes engaged
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Start Timer
    t0 = tt.time()

    # Show Update
    liner.send(' Loading Read Packs ...')

    # Open Archive
    archive = ut.get_archive(
        arcfile=packfile)

    # Enque Read Packs
    packqueue.multiput(filter(
        lambda arcentry: arcentry.endswith('.pack'),
        ('0.0.0.pack',)))
        # archive.namelist()))

    # Show Final Updates
    liner.send(
        ' Pack Queue  : {:,} Read Pack(s) Loaded\n'.format(
            len(packqueue)))
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - t0))

    # Insert Exhaustion Tokens
    packqueue.multiput((None for _ in range(ncores)))

    # Unblock Queue
    enqueuecomplete.set()

# Engine Helper Functions

def exoneration_procedure(
    exoread,
    exofreq,
    phiXkval,
    phiXspec,
    cctrs):
    '''
    Determine exoneration for given
    read and update core counters.
    Internal use only.

    TBD
    '''

    # Adjust Read
    if '-' in exoread:
        exoread = exoread.replace('-', '')

    # Exoneration Flag
    exonerated = False

    # PhiX Match?
    for kmer in ut.stream_spectrum(
        seq=exoread,
        k=phiXkval):

        if kmer in phiXspec:
            cctrs['phiXreads'] += exofreq
            exonerated = True
            break

    # Nucleotide Composition
    acount = exoread.count('A')
    tcount = exoread.count('T')
    gcount = exoread.count('G')
    ccount = exoread.count('C')

    # Dinucleotide Match?
    if not exonerated:
        # Composition Largely Dinucleotide?
        dinukethresh = 0.75 * len(exoread)
        if (ccount + tcount >= dinukethresh) or \
           (acount + ccount >= dinukethresh) or \
           (ccount + gcount >= dinukethresh) or \
           (acount + gcount >= dinukethresh) or \
           (acount + tcount >= dinukethresh) or \
           (tcount + gcount >= dinukethresh):
            cctrs['lowcomplexreads'] += exofreq
            exonerated = True

    # Mononucleotide Match?
    if not exonerated:
        # Composition Large Mononucleotide?
        mononukethresh = 0.50 * len(exoread)
        if ccount >= mononukethresh or \
           acount >= mononukethresh or \
           gcount >= mononukethresh or \
           tcount >= mononukethresh:
            cctrs['lowcomplexreads'] += exofreq
            exonerated = True

    # Trinucleotide Match?
    if not exonerated:
        # Composition Largely Dinucleotide?
        trinukethresh = 1.00 * len(exoread)
        if (ccount + acount + tcount >= trinukethresh) or \
           (ccount + acount + gcount >= trinukethresh) or \
           (acount + tcount + gcount >= trinukethresh) or \
           (tcount + ccount + gcount >= trinukethresh):
            cctrs['lowcomplexreads'] += exofreq
            exonerated = True

    # Uncategorized Discard
    if not exonerated:
        cctrs['incalcreads'] += exofreq

def get_anchored_read(
    read,
    metamap,
    cctrs):
    '''
    TBD
    '''

    # Too Short of a Read!
    if len(read) < (len(metamap['anchor']) - metamap['anchortval']):
        return read, False

    # Quick Anchor!
    qcount = 0
    qtype  = None
    if metamap['anchor'] in read:
        # cctrs['directcount'] += 1
        qcount += 1
        qtype   = 0
        return read, True
    if metamap['revanchor'] in read:
        # cctrs['inversioncount'] += 1
        qcount += 1
        qtype   = 1

    # Anchor Ambiguous?
    if qcount > 1:
        return read, False

    # Anchor Resolved
    if qcount == 1:
        if qtype:
            return ut.get_revcomp(
                seq=read), True
        else:
            return read, True

    # Define Anchor Score
    alnfscore = float('inf')

    # Compute Anchor Alignment
    alnf = ed.align(
        query=metamap['anchor'],
        target=read,
        mode='HW',
        task='distance',
        k=metamap['anchortval'])

    # Update Anchor Score
    alnfd = alnf['editDistance']
    if alnfd > -1:
        alnfscore = alnfd

    # Anchor Perfect Fit
    if alnfscore == 0:
        # cctrs['directcount'] += 1
        return read, True

    # Define Reverse Score
    alnrscore = float('inf')

    # Compute Reverse Alignment
    alnr = ed.align(
        query=metamap['revanchor'],
        target=read,
        mode='HW',
        task='distance',
        k=metamap['anchortval'])

    # Update Reverse Score
    alnrd = alnr['editDistance']
    if alnrd > -1:
        alnrscore = alnrd

    # Reverse Perfect Fit
    if alnrscore == 0:
        # cctrs['inversioncount'] += 1
        return ut.get_revcomp(
            seq=read), True

    # Return Results
    if alnfscore == alnrscore:
        return read, False
    elif alnfscore < alnrscore:
        # cctrs['directcount'] += 1
        return read, True
    else:
        # cctrs['inversioncount'] += 1
        return ut.get_revcomp(
            seq=read), True

def get_trimmed_read(
    read,
    const,
    constype,
    constval):
    '''
    TBD
    '''

    # Nothing to Trim!
    if const is None:
        return read, True
    if len(read) == 0:
        return read, False

    # Too Short of a Read!
    if len(read) < (len(const) - constval):
        return read, False

    # Quick Find!
    if constype <= 0:
        start = read.rfind(const)
    else:
        start = read.find(const)

    # Match Found!
    if start > -1:

        # Compute Extraction Coordinates
        if constype <= 0:
            start = start + len(const)
            stop  = len(read)
        else:
            stop  = start
            start = 0

    # No Direct Match ..
    else:

        # Compute Constant Alignment
        aln = ed.align(
            query=const,
            target=read,
            mode='HW',
            task='locations',
            k=constval)

        # Constant Absent
        if aln['editDistance'] == -1:
            return read, False

        # Extract Locations
        locs = aln['locations']

        # Compute Extraction Coordinates
        if constype <= 0:
            start = locs[-1][-1]+1
            stop  = len(read)
        else:
            stop  = locs[+0][+0]+0
            start = 0

    # Compute and Return Trimmed Read
    extracted_seq = read[start:stop]
    return extracted_seq, True

def get_associate_match(
    read,
    associate,
    associatetval):
    '''
    TBD
    '''

    # Query-Reference Adjustment
    if len(read) < len(associate):
        query  = read
        target = associate
    else:
        query  = associate
        target = read

    # Quick Match!
    if query in target:
        return True

    # Compute Associate Alignment
    aln = ed.align(
        query=query,
        target=target,
        mode='HW',
        task='distance',
        k=associatetval)

    # Return Results
    if aln['editDistance'] == -1:
        return False
    return True

def count_aggregator(
    prepfile,
    countqueue,
    countfile,
    prodcount):
    '''
    Aggregate all count matrices stored in
    countqueue in to a pandas DataFrame and
    optionally write it out as countfile
    CSV file. Internal use only.

    :: prepfile
       type - string
       desc - path to compressed zipfile
              storing prepared structures
              and data models for amplicon
              and sample data
    :: countqueue
       type - SimpleQueue
       desc - queue storing count matrices
              aggregating one or more read
              packs processed by each core
    :: countfile
       type - string / None
       desc - path to CSV file storing
              amplicon read counts of
              prepared samples; if None
              then no CSV file is written
              and only an in-memory pandas
              DataFrame of read counts is
              returned to the user
    :: prodcount
       type - integer
       desc - total number of producers
              scheduled to add to countqueue
    '''

    # Define Aggregate Matrix and DataFrame
    agmatrix = None
    agdf     = None

    # Open prepfile
    prepfile = zf.ZipFile(
        file=prepfile)

    # Load Meta Map and Sample Features
    metamap = ut.loaddict(
        archive=prepfile,
        dfile='meta.map')
    samfeats = ut.loaddict(
        archive=prepfile,
        dfile='sam.feats')

    # Close prepfile
    prepfile.close()

    # Waiting on Count Matrices
    while countqueue.empty():
        continue

    # Count Matrix Aggregation Loop
    while prodcount:

        # Fetch Count Matrix / Token
        cqtoken = countqueue.get()

        # Exhaustion Token
        if cqtoken is None:
            prodcount -= 1
            continue

        # Update Aggregate Matrix
        if not agmatrix is None:
            agmatrix += cqtoken
        else:
            agmatrix  = cqtoken

    # Empty Aggregate Matrix?
    if np.sum(agmatrix) == 0.:
        agmatrix = None

    # Do we have Aggregate Matrix?
    if not agmatrix is None:

        # Build Aggregate DataFrame
        agdf = pd.DataFrame(
            data=agmatrix.astype(np.int64),
            columns=metamap['headers'],
            index=metamap['SampleID'])
        agdf.index.name = 'SampleID'

        # Build Meta Feature DataFrame
        ftdf = pd.DataFrame.from_dict(
            data=samfeats,
            orient='index')
        ftdf.index.name = 'SampleID'

        # Merge DataFrames
        agdf = agdf.join(
            ftdf, on='SampleID')

        # Dump Aggregate DataFrame
        if not countfile is None:
            agdf.to_csv(
                path_or_buf=countfile,
                sep=',')

    # Return Aggregate DataFrame
    return agdf
