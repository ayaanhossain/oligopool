import time    as tt

import zipfile as zf

import edlib   as ed
import leveldb as lv

import utils   as ut


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
    packqueue.multiput(map(lambda x: x.removesuffix('.pack'),
        filter(lambda arcentry: arcentry.endswith('.pack'),
            archive.namelist())))

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

    :: exoread
       type - string
       desc - read to be exonerated
    :: exofreq
       type - integer
       desc - absolute frequency of
              the read in given pack
    :: phiXkval
       type - integer
       desc - k-value for phiX matching
    :: phiXspec
       type - set
       desc - set of all possible phiX
              k-mers (includes revcomp)
    :: cctrs
       type - dict
       desc - dictionary of core counters
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
    metamap):
    '''
    Return orientation corrected reads,
    in case the reads are flipped.
    Internal use only.

    :: read
       type - string
       desc - read to be anchored
    :: metamap
       type - dict
       desc - ditionary containing index
              meta information
    '''

    # Too Short of a Read!
    if len(read) < (len(metamap['anchor']) - metamap['anchortval']):
        return read, False

    # Quick Anchor!
    qcount = 0
    qtype  = None
    if metamap['anchor'] in read:
        qcount += 1
        qtype   = 0
    if metamap['revanchor'] in read:
        qcount += 1
        qtype   = 1

    # Anchor Ambiguous?
    if qcount > 1:
        return read, False

    # Resolve Anchor
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
        return ut.get_revcomp(
            seq=read), True

    # Return Results
    if alnfscore == alnrscore:
        return read, False
    elif alnfscore < alnrscore:
        return read, True
    return ut.get_revcomp(
        seq=read), True

def get_trimmed_read(
    read,
    const,
    constype,
    constval):
    '''
    Return read after trimming
    constant sequence.
    Internal use only.

    :: read
       type - string
       desc - read to be trimmed
    :: const
       type - string / None
       desc - constant to trim
    :: constype
       type - integer
       desc - constant type identifier
              0 = prefix constant
              1 = suffix constant
    :: constval
       type - integer / None
       desc - constant t-value to be
              tolerated for matching
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
    return read[start:stop], True

def get_barcode_index(
    barcoderead,
    metamap,
    model):
    '''
    Determine and return barcode identifier
    index in given anchored read.
    Internal use only.

    :: barcoderead
       type - string
       desc - read containing barcode
              information
    :: metamap
       type - dict
       desc - ditionary containing index
              meta information
    :: model
       type - Scry
       desc - Scry model for classifying
              barcode
    '''

    # Book-keeping
    trimpolicy   = metamap['trimpolicy']
    pxtrimstatus = True
    sxtrimstatus = True
    trimstatus   = True

    # Prefix Trim
    if trimpolicy <= 1.5:
        # Prefix Trim Read
        barcoderead, pxtrimstatus = get_trimmed_read(
            read=barcoderead,
            const=metamap['barcodeprefix'],
            constype=0,
            constval=metamap['bpxtval'])
        # Trim was Unsuccessful
        if (trimpolicy == 1) and (pxtrimstatus is False):
            return None

    # Suffix Trim
    if trimpolicy >= 1.5:
        # Suffix Trim Read
        barcoderead, sxtrimstatus = get_trimmed_read(
            read=barcoderead,
            const=metamap['barcodesuffix'],
            constype=1,
            constval=metamap['bsxtval'])
        # Trim was Unsuccessful
        if (trimpolicy == 2) and (sxtrimstatus is False):
            return None

    # Policy Trim
    if trimpolicy == 1:
        barcoderead = barcoderead[:+(metamap['barcodelen'] + metamap['barcodetval']) ]
    if trimpolicy == 2:
        barcoderead = barcoderead[ -(metamap['barcodelen'] + metamap['barcodetval']):]
    if trimpolicy == 1.5:
        trimstatus = pxtrimstatus and sxtrimstatus
        if not trimstatus:
            return None

    # Compute Barcode Index
    return model.predict(
        x=barcoderead)[0]

def is_associate_match(
    associate,
    associateread,
    associatetval):
    '''
    Determine if associate exists
    in given anchored read.
    Internal use only.

    :: associate
       type - string
       desc - full associate sequence
              to match
    :: associateread
       type - string
       desc - read containing associate
              information
    :: associatetval
       type - integer
       desc - associate t-value to be
              tolerated for matching
    '''

    # Query-Reference Adjustment
    if len(associateread) <= len(associate):
        query  = associateread
        target = associate
    else:
        query  = associate
        target = associateread

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

def get_associate_match(
    associateread,
    associateerrors,
    associatedict,
    index,
    metamap):
    '''
    Process associate read and determine if
    it contains required associate.
    Internal use only.

    :: associateread
       type - string
       desc - read containing associate
              information
    :: associaterrors
       type - integer
       desc - maximum number of mismatches
              between reference and read
              associate
    :: associatedict
       type - dict
       desc - dictionary containing associate
              information from index
    :: index
       type - integer
       desc - associate index identifier
    :: metamap
       type - dict
       desc - ditionary containing index
              meta information
    '''

    # Trim Associate Prefix
    associateread, trimstatus = get_trimmed_read(
        read=associateread,
        const=metamap['associateprefix'],
        constype=0,
        constval=metamap['apxtval'])

    # Associate Prefix Absent
    if not trimstatus:
        return False, 0

    # Trim Associate Suffix
    associateread, trimstatus = get_trimmed_read(
        read=associateread,
        const=metamap['associatesuffix'],
        constype=1,
        constval=metamap['asxtval'])

    # Associate Suffix Absent
    if not trimstatus:
        return False, 0

    # Match Associate
    associate, associatetval = associatedict[index]
    if (associateerrors > -1) and \
       (associateerrors < associatetval):
        associatetval = associateerrors

    # Return Results
    return (is_associate_match(
        associate=associate,
        associateread=associateread,
        associatetval=associatetval),
        1)

def count_aggregator(
    countqueue,
    countdir,
    prodcount,
    prodactive,
    liner):
    '''
    TBD
    '''

    # Define CountDB file
    countdb = '{}/countdb'.format(
        countdir)

    # Open CountDB Instance
    countdb = lv.LevelDB(
        filename=countdb,
        create_if_missing=True,
        error_if_exists=True)

    # Waiting on Count Matrices
    while countqueue.empty():
        tt.sleep(0)
        continue

    # Count Dicionary Aggregation Loop
    while prodcount:

        # Fetch Count Path / Token
        cqtoken = countqueue.get()

        # Exhaustion Token
        if cqtoken is None:
            prodcount -= 1
            continue

        # Build Count Path
        cpath = '{}/{}'.format(
            countdir,
            cqtoken)

        # Load Count List
        fname = cqtoken
        countlist = ut.loadcount(
            cfile=cpath)

        # Remove Count List
        ut.remove_file(
            filepath=cpath)

        # Show Updates
        if prodactive and \
           prodactive.value() == 0:
            liner.send(
                ' Aggregating: {}'.format(
                    fname))

        # Create New Batch
        batch = lv.WriteBatch()

        # Update Batch
        while countlist:
            k,v  = countlist.pop()
            strk = str(k).encode()
            try:
                w = eval(countdb.Get(strk))
            except:
                w = 0
            strv = str(w+v).encode()
            batch.Put(strk, strv)

        # Update CountDB
        countdb.Write(
            batch,
            sync=True)

        # Release Control
        tt.sleep(0)

def write_count(
    indexfiles,
    countdir,
    countfile,
    liner):
    '''
    TBD
    '''

    # Book-keeping
    IDdicts    = []
    indexnames = []
    t0 = tt.time()

    # Show Update
    liner.send(' Loading IDs ...')

    # Load IDdicts
    for indexfile in indexfiles:

        # Update Index Names
        indexnames.append(indexfile.split(
            '/')[-1].removesuffix(
                '.oligopool.index'))

        # Open indexfile
        indexfile = zf.ZipFile(
            file=indexfile)

        # Update IDdict
        IDdicts.append(ut.loaddict(
            archive=indexfile,
            dfile='ID.map'))

        # Close indexfile
        indexfile.close()

    # Show Update
    liner.send(' Writing Count Matrix ...')

    # Define CountDB file
    countdb = '{}/countdb'.format(
        countdir)

    # Open CountDB Instance
    countdb = lv.LevelDB(
        filename=countdb,
        create_if_missing=False,
        error_if_exists=False)

    # Count Matrix Loop
    with open(countfile, 'w') as outfile:

        # Write Header
        outfile.write(','.join('{}.ID'.format(
            idxname) for idxname in indexnames) + ',Counts\n')

        # Book-keeping
        entrycount = 0
        batchsize  = 10**4
        batchreach = 0

        # Loop through CountDB Entries
        for indextuple, count in countdb.RangeIter(
            key_from=None, key_to=None):

            # Update Book-keeping
            entrycount += 1
            batchreach += 1

            # Build Entry
            rows = []
            for IDx, index in enumerate(eval(indextuple)):
                if index == '-':
                    rows.append('-')
                else:
                    rows.append(IDdicts[IDx][index])
            rows.append(count.decode('ascii'))

            # Write Entry to Count Matrix
            outfile.write(','.join(rows) + '\n')

            # Show Update
            if batchreach == batchsize:
                liner.send(
                    ' Rows Written: {:,} Unique ID / Combination(s)'.format(
                        entrycount))
                batchreach = 0

    # Final Updates
    liner.send(
        ' Rows Written: {:,} Unique ID / Combination(s)\n'.format(
            entrycount))
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - t0))
