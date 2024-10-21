import os

import time as tt

import zipfile as zf
import collections as cx
import random as rn

import numpy  as np
import pandas as pd
import edlib  as ed

import plyvel

from . import utils as ut
from . import phiX  as px


# Parser and Setup Functions

def get_random_DNA(length):
    '''
    Return a random DNA sequence
    of required length.
    Internal use only.

    :: length
       type - integer
       desc - length of required
              DNA sequence
    '''

    return ''.join(np.random.choice(
        ('A', 'T', 'G', 'C')) for _ in range(length))

def get_meta_info(packfile):
    '''
    Compute the number of reads to
    stream for callback parsing, and
    extract the meta read pack.

    :: packfile
       type - string
       desc - path to packfile storing
              read packs
    '''

    # Open Pack File
    packfile = ut.get_archive(
        arcfile=packfile)

    # List all packs
    allpacknames = [pkname for pkname in packfile.namelist() \
        if pkname != 'packing.stat']

    # Load Header Pack
    packzero = ut.loadpack(
        archive=packfile,
            pfile=allpacknames[0])

    # Close Packfile
    packfile.close()

    # Determine Stream Count
    count = min(1000, len(packzero))

    # Return Results
    return (packzero,
        count)

def stream_random_reads(count):
    '''
    Stream several random reads.
    Internal use only.

    :: count
       type - integer
       desc - total number of reads
              to be streamed
    '''

    for _ in range(count):
        yield get_random_DNA(
            length=np.random.randint(
                10, 1000))

def stream_packed_reads(packzero, count):
    '''
    Stream several experiment reads.
    Internal use only.

    :: packfile
       type - string
       desc - meta read pack containing
              most frequent reads
    :: count
       type - integer
       desc - total number of reads
              to be streamed
    '''

    # Build Streaming Index Vector
    indexvec = np.arange(len(packzero))
    np.random.shuffle(indexvec)
    indexvec = indexvec[:count]

    # Read Streaming Loop
    for readidx in indexvec:
        yield packzero[readidx][0]

def stream_IDs(indexfiles, count):
    '''
    Stream random ID combination.
    Internal use only.

    :: indexfiles
       type - tuple
       desc - tuple of file paths
              to applied indices
    :: count
       type - integer
       desc - total number of reads
              to be streamed
    '''

    # Open Index Files
    indexfiles = list(
        ut.get_archive(idxfile) for idxfile in indexfiles)

    # Load IDdics
    IDdicts = [ut.loaddict(
        archive=indexfile,
        dfile='ID.map') for indexfile in indexfiles]

    # Close Index Files
    for indexfile in indexfiles:
        indexfile.close()

    # Stream Random ID Combination
    for _ in range(count):

        # Build ID Tuple
        IDtuple = []
        for IDdict in IDdicts:
            if len(indexfiles) == 1:
                entry = IDdict[
                    np.random.randint(len(IDdict))]
            else:
                toss = np.random.randint(10)
                if toss == 0:
                    entry = '-'
                else:
                    entry = IDdict[
                        np.random.randint(len(IDdict))]
            IDtuple.append(entry)

        # Stream ID Tuple
        yield tuple(IDtuple)

def get_parsed_callback(
    indexfiles,
    callback,
    packfile,
    ncores,
    liner):
    '''
    Determine if given callback
    function is valid.
    Internal use only.

    :: indexfiles
       type - tuple
       desc - tuple of file paths
              to applied indices
    :: callback
       type - function
       desc - callback function
              specified
    :: packfile
       type - string
       desc - path to packfile storing
              read packs
    :: ncores
       type - integer
       desc - total number of cores to
              be used in counting
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Start Timer
    t0 = tt.time()

    # Show Update
    liner.send(' Loading Read Generators ...')

    # Get Meta Information
    (packzero,
    count) = get_meta_info(
        packfile=packfile)

    # Setup Streamers
    randstreamer = stream_random_reads(
        count=count)
    packstreamer = stream_packed_reads(
        packzero=packzero,
        count=count)
    streamers  = [randstreamer, packstreamer]
    IDstreamer = stream_IDs(
        indexfiles=indexfiles,
        count=count)

    # Parsing Loop
    pass_count   = 0
    failedinputs = []
    for readcount in range(count):

        # Get a Read!
        read = next(streamers[np.random.randint(2)])

        # Execute Callback
        try:
            callback_input = {
                  'read': read,
                    'ID': next(IDstreamer),
                 'count': np.random.randint(1, 10001),
                'coreid': np.random.randint(ncores)}
            result = callback(**callback_input)
            if not isinstance(result,bool):
                raise
        except:
            failedinputs.append(callback_input)
            iter_valid  = False
        else:
            iter_valid  = True
            pass_count += 1

        # Show Update
        liner.send(
            ' Callback Evaluation: {:,} / 1,000 {}'.format(
                readcount,
                ('Rejected', 'Accepted')[iter_valid]))

    # Close Generators
    randstreamer.close()
    packstreamer.close()
    IDstreamer.close()

    # Final Update
    plen = ut.get_printlen(
        value=count)
    liner.send(
        ' Callback Passed: {:{},} Input(s)\n'.format(
            pass_count,
            plen))
    liner.send(
        ' Callback Failed: {:{},} Input(s)\n'.format(
            count - pass_count,
            plen))
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - t0))

    # Show Final Verdict
    parsestatus = len(failedinputs) == 0
    if not parsestatus:
        liner.send(
            ' Verdict: Counting Infeasible due to Callback Function\n')
    else:
        liner.send(
            ' Verdict: Counting Possibly Feasible\n')

    # Cleanup
    del packzero

    # Return results
    return (parsestatus,
        failedinputs)

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
    packqueue.multiput(map(
        lambda x: ut.removestarfix(
            string=x,
            fix='.pack',
            loc=1),
        filter(
            lambda arcentry: arcentry.endswith('.pack'),
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
    if '-----' in exoread:
        exoread = exoread.replace('-----', '')

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

    # Gap Trim
    if metamap['barcodegapped']:
        # Localise Gap Lengths
        pregap, postgap = (metamap['barcodepregap'],
            metamap['barcodepostgap'])
        # Pregap Trim
        if pregap:
            barcoderead = barcoderead[ +pregap:]
        # Postgap Trim
        if postgap:
            barcoderead = barcoderead[:-postgap]
        # Nothing Remains after Gap Trim
        if not barcoderead:
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

def get_failed_inputs(
    packqueue,
    countdir,
    liner):
    '''
    Empty packqueue and return failed inputs
    for callback. Internal use only.

    :: packqueue
       type - SafeQueue
       desc - queue storing read pack file
              paths
    :: countdir
       type - string
       desc - filepath to counting workspace
              directory
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Show Update
    liner.send(
        ' Emptying Pack Queue ...')

    # Consume Pack Queue
    for item in packqueue.multiget():
        pass

    # Show Update
    liner.send(
        ' Extracting Failed Input(s) ...')

    # Aggregate Callback Dumps
    failedinputs = []
    for entry in os.listdir(countdir):
        if entry.endswith('.callback.dump'):
            failedinputs.append(
                ut.loaddump(dfile='{}/{}'.format(
                    countdir,
                    entry)))

    # Return Results
    return failedinputs

def count_aggregator(
    countqueue,
    countdir,
    prodcount,
    prodactive,
    assoc,
    liner):
    '''
    Aggregate count dictionaries in
    countqueue into a count database.
    Internal use only.

    :: countqueue
       type - SafeQueue
       desc - count queue storing
              count dictionaries
    :: countdir
       type - string
       desc - path to work space
              count directory
    :: prodcount
       type - integer
       desc - total number of producers
              scheduled to add to
              countqueue
    :: prodactive
       type - SafeCounter
       desc - total number of producers
              actively adding to
              countqueue
    :: assoc
       type - boolean
       desc - if True, we are aggregating
              association counts of tuples
    :: liner
       type - coroutine
       desc - dynamic printing and
              logging
    '''

    # Define CountDB file
    countdb = '{}/countdb'.format(
        countdir)

    # Open CountDB Instance
    countdb = plyvel.DB(
        countdb,
        create_if_missing=True,
        error_if_exists=True)

    # Waiting on Count Matrices
    while countqueue.empty():
        tt.sleep(0)
        continue

    # Count Dictionary Aggregation Loop
    while prodcount:

        # Fetch Count Path / Token
        cqtoken = countqueue.get()

        # Exhaustion Token
        if cqtoken is None:
            prodcount -= 1
            continue

        # Build Count Path
        cpath = '{}/{}'.format(countdir, cqtoken)

        # Load Count List
        fname = cqtoken
        countlist = ut.loadcount(cfile=cpath)

        # Remove Count List
        ut.remove_file(filepath=cpath)

        # Show Updates
        if prodactive and \
           prodactive.value() == 0:
            liner.send(' Aggregating: {}'.format(fname))

        # Create New Batch
        batch = countdb.write_batch(transaction=True)

        # Update Batch
        while countlist:
            k,v  = countlist.pop()
            strk = str(k).encode()
            w = countdb.get(strk)
            match assoc:
                case True:
                    if w is None:
                        w = [0, 0]
                    else:
                        w = eval(w.decode('ascii'))
                    v[0] += w[0]
                    v[1] += w[1]
                    strv = str(v).encode()
                case False:
                    if w is None:
                        w = 0
                    else:
                        w = int(w.decode('ascii'))
                    strv = str(w+v).encode()
            batch.put(strk, strv)

        # Update CountDB
        batch.write()

        # Release Control
        tt.sleep(0)

def write_count(
    indexfiles,
    countdir,
    countfile,
    assoc,
    liner):
    '''
    Write entries in count database
    to a final count file / matrix.
    Internal use only.

    :: indexfiles
       type - tuple
       desc - tuple of index file
              paths
    :: countdir
       type - string
       desc - path to workspace
              count directory
    :: countfile
       type - string
       desc - path to CSV file storing
              read counts for discovered
              ID combinations
    :: assoc
       type - boolean
       desc - if True, we are writing
              association counts
    :: liner
       type - coroutine
       desc - dynamic printing
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
        indexnames.append(ut.removestarfix(
            string=indexfile.split('/')[-1],
            fix='.oligopool.index',
            loc=1))

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
    countdb = plyvel.DB(
        countdb,
        create_if_missing=False,
        error_if_exists=False)

    # Count Matrix Loop
    with open(countfile, 'w') as outfile:

        # Write Header
        match assoc:
            case True:
                outfile.write(','.join('{}.ID'.format(
                    idxname) for idxname in indexnames) + ',BarcodeCounts,AssociationCounts\n')
            case False:
                outfile.write(','.join('{}.ID'.format(
                    idxname) for idxname in indexnames) + ',Counts\n')

        # Book-keeping
        entrycount = 0
        batchsize  = 10**4
        batchreach = 0

        # Loop through CountDB Entries
        for indextuple,count in countdb:

            # Update Book-keeping
            entrycount += 1
            batchreach += 1

            # Build Entry
            rows = []
            for IDx,index in enumerate(eval(indextuple)):
                if index == '-':
                    rows.append('-')
                else:
                    rows.append(IDdicts[IDx][int(index)])
            match assoc:
                case True:
                    rows.append(','.join(str(x) for x in eval(count.decode('ascii'))))
                case False:
                    rows.append(count.decode('ascii'))

            # Write Entry to Count Matrix
            outfile.write(','.join(map(str,rows)) + '\n')

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

    # Read Back the CSV
    liner.send('|* Reading Count Matrix ...')
    df = pd.read_csv(countfile)

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - t0))

    # Return DataFrame
    return df

def acount_engine(
    indexfile,
    packfile,
    packqueue,
    countdir,
    countqueue,
    maptype,
    barcodeerrors,
    associateerrors,
    previousreads,
    analyzedreads,
    phiXreads,
    lowcomplexreads,
    misassocreads,
    falsereads,
    incalcreads,
    experimentreads,
    callback,
    callbackerror,
    coreid,
    batchid,
    ncores,
    nactive,
    memlimit,
    restart,
    shutdown,
    launchtime,
    liner):
    '''
    Association count read packs stored in packfile
    and enqueue the final read count matrix in to
    countqueue. Internal use only.

    :: indexfile
       type - string
       desc - archive of indexed objects storing barcode
              and associate data and models
    :: packfile
       type - string
       desc - path to compressed zipfile
              storing read packs
    :: packqueue
       type - SimpleQueue
       desc - queue storing path to read pack
              files stored in packfile
    :: countdir
       type - string
       desc - filepath to counting workspace
              directory
    :: countqueue
       type - SimpleQueue
       desc - queue storing count matrices
              aggregating one or more read
              packs processed by each core
    :: maptype
       type - integer
       desc - mapping type identifier
              0 = fast mapping
              1 = exact mapping
    :: barcodeerrors
       type - Real / None
       desc - maximum number of mutations
              tolerated in barcodes before
              discarding reads from counting
    :: associateerrors
       type - Real / None
       desc - maximum number of mutations
              tolerated in associates before
              discarding reads from counting
    :: previousreads
       type - SafeCounter
       desc - total number of reads processed
              previously
    :: analyzedreads
       type - SafeCounter
       desc - total number of reads processed
              so far during counting
    :: phiXreads
       type - SafeCounter
       desc - total number of reads processed
              attributed to PhiX contamination
    :: lowcomplexreads
       type - SafeCounter
       desc - total number of reads processed
              attributed to low complexity products
    :: misassocreads
       type - SafeCounter
       desc - total number of reads processed
              attributed to mis-association
    :: falsereads
       type - SafeCounter
       desc - total number of reads processed
              that were rejected by callback function
    :: incalcreads
       type - SafeCounter
       desc - total number of invalid reads or ones
              not usable for quantification
    :: experimentreads
       type - SafeCounter
       desc - total number of valid reads used for
              quantification
    :: callback
       type - function
       desc - callback function to invoke during
              counting for concurrent processing
    :: callbackerror
       type - mp.Event
       desc - If True, signals error in callback
              execution
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
    :: restart
       type - mp.Event
       desc - multiprocessing event when
              process needs to restart
    :: shutdown
       type - mp.Event
       desc - multiprocessing event when
              process needs to shutdown
    :: launchtime
       type - time
       desc - initial launch timestamp
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Open indexfile
    indexfile = ut.get_archive(
        arcfile=indexfile)

    # Load Barcode Model
    model = ut.loadmodel(
        archive=indexfile,
        mfile='barcode.model')

    # Prime Barcode Model
    model.prime(
        t=barcodeerrors,
        mode=maptype)

    # Load Maps
    IDdict = ut.loaddict(
        archive=indexfile,
        dfile='ID.map')
    metamap = ut.loaddict(
        archive=indexfile,
        dfile='meta.map')
    associatedict = ut.loaddict(
        archive=indexfile,
        dfile='associate.map')

    # Close indexfile
    indexfile.close()

    # Define PhiX Spectrum
    phiXkval = 30
    phiXspec = px.get_phiX_spectrum(
        k=phiXkval)

    # Open packfile
    packfile = ut.get_archive(
        arcfile=packfile)

    # Load packing.stat
    packstat = ut.loaddict(
        archive=packfile,
        dfile='packing.stat')

    # Book-keeping Variables
    cctrs = {
        'previousreads'  : previousreads.value(),
        'analyzedreads'  : 0,
        'phiXreads'      : 0,
        'lowcomplexreads': 0,
        'misassocreads'  : 0,
        'falsereads'     : 0,
        'incalcreads'    : 0,
        'experimentreads': 0}
    callbackabort = False

    # Compute Printing Lengths
    clen = ut.get_printlen(value=ncores)
    slen = plen = ut.get_printlen(
        value=packstat['survived_reads'])

    # Read Count Storage
    countdict = cx.defaultdict(lambda: [0, 0])
    countpath = '{}/{}.{}.count'.format(
        countdir, coreid, batchid)

    # Pack Counting Loop
    while not packqueue.empty():

        # Callback Error Somewhere?
        if not callback is None:
            if callbackerror.is_set():
                callbackabort = True
                shutdown.set()
                break

        # Fetch Pack Name / Token
        packname = packqueue.get()

        # Exhaustion Token
        if packname is None:
            shutdown.set()
            break

        # Fetch Read Pack
        cpack = ut.loadpack(
            archive=packfile,
            pfile='{}.pack'.format(
                packname))

        # Book-keeping Variables
        exoread   = None
        exofreq   = None
        readcount = 0

        verbagereach  = 0
        verbagetarget = rn.randint(
            *map(round, (len(cpack) * 0.080,
                         len(cpack) * 0.120)))

        # Start Timer
        t0 = tt.time()

        # Read Counting Loop
        while True:

            # Callback Error Somewhere?
            if not callback is None:
                if callbackerror.is_set():
                    callbackabort = True
                    break

            # Exoneration Block
            if not exoread is None:

                # Run Exoneration Procedure
                exoneration_procedure(
                    exoread=exoread,
                    exofreq=exofreq,
                    phiXkval=phiXkval,
                    phiXspec=phiXspec,
                    cctrs=cctrs)

                # Clear for Next Exoneration
                exoread = None
                exofreq = None

            # Time to Show Update?
            if verbagereach >= verbagetarget:

                # Show Update
                liner.send(
                    ' Core {:{},d}: Analyzed {:{},d} Reads in {:.2f} sec'.format(
                        coreid,
                        clen,
                        cctrs['analyzedreads'] + cctrs['previousreads'],
                        slen,
                        tt.time()-launchtime))

                # Update Book-keeping
                verbagereach = 0

            # Continue processing pack?
            if not cpack:
                break

            # Fetch Read and Frequency
            (read,
            freq) = cpack.pop()

            # Update Book-keeping
            cctrs['analyzedreads'] += freq
            verbagereach += 1
            readcount    += 1

            # Anchor Read
            (anchoredread,
            anchorstatus) = get_anchored_read(
                read=read,
                metamap=metamap)

            # Anchor Absent
            if not anchorstatus:
                exoread = read
                exofreq = freq
                continue

            # Setup Read References
            barcoderead   = anchoredread
            associateread = anchoredread

            # Compute Barcode Index
            index = get_barcode_index(
                barcoderead=barcoderead,
                metamap=metamap,
                model=model)

            # Barcode Absent
            if index is None:
                cctrs['incalcreads'] += freq
                continue

            # Compute Associate Match
            associatematch, basalmatch = get_associate_match(
                associateread=associateread,
                associateerrors=associateerrors,
                associatedict=associatedict,
                index=index,
                metamap=metamap)

            # Associate Absent / Incorrect
            if not associatematch:

                # Associate Constants Missing
                if not basalmatch:
                    cctrs['incalcreads'] += freq
                    continue
                # Associate Mismatches with Reference
                else:
                    cctrs['misassocreads'] += freq

            # Compute Callback Evaluation
            if associatematch and (not callback is None):
                try:
                    # Execute Callback
                    evaluation = callback(
                        read=anchoredread,
                        ID=(IDdict[index],),
                        count=freq,
                        coreid=coreid)

                    # Uh oh ... Non-boolean Output
                    if not ((evaluation is True) or \
                            (evaluation is False)):
                        raise
                except:
                    # We have a failed evaluation!
                    callbackerror.set()
                    callbackabort = True

                    # Dump input for Stats
                    ut.savedump(dobj={
                          'read': anchoredread,
                            'ID': (index,),
                         'count': freq,
                        'coreid': coreid},
                        filepath='{}/{}.callback.dump'.format(
                            countdir,
                            coreid))

                    # Break out!
                    break

                else:
                    # Callback Evaluation is False
                    if not evaluation:
                        cctrs['falsereads'] += freq
                        continue

            # Tally Read Counts
            if associatematch:
                cctrs['experimentreads'] += freq
                countdict[((index),)][1] += freq # Association Count
            countdict[((index),)][0] += freq     # Barcode Count

        # Show Final Updates
        liner.send(
            ' Core {:{},d}: Counted Pack {} w/ {:{},d} Reads in {:05.2f} sec\n'.format(
                coreid,
                clen,
                packname,
                readcount,
                plen,
                tt.time()-t0))

        # Free Memory
        ut.free_mem()

        # Release Control
        tt.sleep(0)

        # Did we Abort due to Callback?
        if callbackabort:
            shutdown.set()
            break

        # Need to Restart?
        if ut.needs_restart(
            memlimit=memlimit):
            restart.set() # Enable Restart
            break # Release, your Memory Real Estate

    # Pack Queue is Empty
    else:
        shutdown.set()


    # Close packfile
    packfile.close()

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

    # We didn't Abort right?
    if not callbackabort:

        # Do we have counts?
        if countdict:

            # Save Count Dictionary
            ut.savecount(
                cobj=countdict,
                filepath=countpath)

            # Release Control
            tt.sleep(0)

            # Queue Count Dicionary Path
            countqueue.put(
                countpath.split('/')[-1])

            # Release Control
            tt.sleep(0)

    # Update Read Counting Book-keeping
    previousreads.increment(incr=cctrs['analyzedreads'])
    analyzedreads.increment(incr=cctrs['analyzedreads'])
    phiXreads.increment(incr=cctrs['phiXreads'])
    lowcomplexreads.increment(incr=cctrs['lowcomplexreads'])
    misassocreads.increment(incr=cctrs['misassocreads'])
    falsereads.increment(incr=cctrs['falsereads'])
    incalcreads.increment(incr=cctrs['incalcreads'])
    experimentreads.increment(incr=cctrs['experimentreads'])

    # Counting Completed
    nactive.decrement()
    if shutdown.is_set():
        countqueue.put(None)

    # Release Control
    tt.sleep(0)

def xcount_engine(
    indexfiles,
    packfile,
    packqueue,
    countdir,
    countqueue,
    maptype,
    barcodeerrors,
    previousreads,
    analyzedreads,
    phiXreads,
    lowcomplexreads,
    falsereads,
    incalcreads,
    experimentreads,
    callback,
    callbackerror,
    coreid,
    batchid,
    ncores,
    nactive,
    memlimit,
    restart,
    shutdown,
    launchtime,
    liner):
    '''
    Combinatorial count barcodes from multiple
    indexes in packed reads and enque the count
    entries to countqueue. Internal use only.

    :: indexfiles
       type - string / list
       desc - one or more archive of indexed objects
              storing barcode data and models
    :: packfile
       type - string
       desc - path to compressed zipfile
              storing read packs
    :: packqueue
       type - SimpleQueue
       desc - queue storing path to read pack
              files stored in packfile
    :: countdir
       type - string
       desc - filepath to counting workspace
              directory
    :: countqueue
       type - SimpleQueue
       desc - queue storing count matrices
              aggregating one or more read
              packs processed by each core
    :: maptype
       type - integer
       desc - mapping type identifier
              0 = fast mapping
              1 = exact mapping
    :: barcodeerrors
       type - Real / None
       desc - maximum number of mutations
              tolerated in barcodes before
              discarding reads from counting
    :: associateerrors
       type - Real / None
       desc - maximum number of mutations
              tolerated in associates before
              discarding reads from counting
    :: previousreads
       type - SafeCounter
       desc - total number of reads processed
              previously
    :: analyzedreads
       type - SafeCounter
       desc - total number of reads processed
              so far during counting
    :: phiXreads
       type - SafeCounter
       desc - total number of reads processed
              attributed to PhiX contamination
    :: lowcomplexreads
       type - SafeCounter
       desc - total number of reads processed
              attributed to low complexity products
    :: misassocreads
       type - SafeCounter
       desc - total number of reads processed
              attributed to mis-association
    :: falsereads
       type - SafeCounter
       desc - total number of reads processed
              that were rejected by callback function
    :: incalcreads
       type - SafeCounter
       desc - total number of invalid reads or ones
              not usable for quantification
    :: experimentreads
       type - SafeCounter
       desc - total number of valid reads used for
              quantification
    :: callback
       type - function
       desc - callback function to invoke during
              counting for concurrent processing
    :: callbackerror
       type - mp.Event
       desc - If True, signals error in callback
              execution
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
    :: restart
       type - mp.Event
       desc - multiprocessing event when
              process needs to restart
    :: shutdown
       type - mp.Event
       desc - multiprocessing event when
              process needs to shutdown
    :: launchtime
       type - time
       desc - initial launch timestamp
    :: liner
       type - coroutine
       desc - dynamic printing
    '''

    # Aggregate Index Files
    models   = []
    metamaps = []
    numindex = len(indexfiles)
    for indexfile in indexfiles:

        # Open indexfile
        indexfile = ut.get_archive(
            arcfile=indexfile)

        # Load Barcode Model
        model = ut.loadmodel(
            archive=indexfile,
            mfile='barcode.model')

        # Prime Barcode Model
        model.prime(
            t=barcodeerrors,
            mode=maptype)

        # Load Maps
        IDdict = ut.loaddict(
            archive=indexfile,
            dfile='ID.map')
        metamap = ut.loaddict(
            archive=indexfile,
            dfile='meta.map')

        # Close indexfile
        indexfile.close()

        # Append Index Objects
        models.append(model)
        metamaps.append(metamap)

    # Define PhiX Spectrum
    phiXkval = 30
    phiXspec = px.get_phiX_spectrum(
        k=phiXkval)

    # Open packfile
    packfile = ut.get_archive(
        arcfile=packfile)

    # Load packing.stat
    packstat = ut.loaddict(
        archive=packfile,
        dfile='packing.stat')

    # Book-keeping Variables
    cctrs = {
        'previousreads'  : previousreads.value(),
        'analyzedreads'  : 0,
        'phiXreads'      : 0,
        'lowcomplexreads': 0,
        'falsereads'     : 0,
        'incalcreads'    : 0,
        'experimentreads': 0}
    callbackabort = False

    # Compute Printing Lengths
    clen = ut.get_printlen(value=ncores)
    slen = plen = ut.get_printlen(
        value=packstat['survived_reads'])

    # Read Count Storage
    countdict = cx.Counter()
    countpath = '{}/{}.{}.count'.format(
        countdir, coreid, batchid)

    # Pack Counting Loop
    while not packqueue.empty():

        # Callback Error Somewhere?
        if not callback is None:
            if callbackerror.is_set():
                callbackabort = True
                shutdown.set()
                break

        # Fetch Pack Name / Token
        packname = packqueue.get()

        # Exhaustion Token
        if packname is None:
            shutdown.set()
            break

        # Fetch Read Pack
        cpack = ut.loadpack(
            archive=packfile,
            pfile='{}.pack'.format(
                packname))

        # Book-keeping Variables
        exoread   = None
        exofreq   = None
        readcount = 0

        verbagereach  = 0
        verbagetarget = rn.randint(
            *map(round, (len(cpack) * 0.080,
                         len(cpack) * 0.120)))

        # Start Timer
        t0 = tt.time()

        # Read Counting Loop
        while True:

            # Callback Error Somewhere?
            if not callback is None:
                if callbackerror.is_set():
                    callbackabort = True
                    break

            # Exoneration Block
            if not exoread is None:

                # Run Exoneration Procedure
                exoneration_procedure(
                    exoread=exoread,
                    exofreq=exofreq,
                    phiXkval=phiXkval,
                    phiXspec=phiXspec,
                    cctrs=cctrs)

                # Clear for Next Exoneration
                exoread = None
                exofreq = None

            # Time to Show Update?
            if verbagereach >= verbagetarget:

                # Show Update
                liner.send(
                    ' Core {:{},d}: Analyzed {:{},d} Reads in {:.2f} sec'.format(
                        coreid,
                        clen,
                        cctrs['analyzedreads'] + cctrs['previousreads'],
                        slen,
                        tt.time()-launchtime))

                # Update Book-keeping
                verbagereach = 0

            # Continue processing pack?
            if not cpack:
                break

            # Fetch Read and Frequency
            (read,
            freq) = cpack.pop()

            # Update Book-keeping
            cctrs['analyzedreads'] += freq
            verbagereach += 1
            readcount    += 1

            # Barcode Mapping Loop
            indextuple = []
            partialanc = False
            partialmap = False
            for idx in range(numindex):

                # Anchor Read
                (anchoredread,
                anchorstatus) = get_anchored_read(
                    read=read,
                    metamap=metamaps[idx])

                # Anchor Absent
                if not anchorstatus:
                    indextuple.append('-')
                    continue

                # Anchoring Successful
                partialanc = True

                # Compute Barcode Index
                index = get_barcode_index(
                    barcoderead=anchoredread,
                    metamap=metamaps[idx],
                    model=models[idx])

                # Barcode Absent
                if index is None:
                    indextuple.append('-')
                    continue

                # Found a Barcode!
                indextuple.append(
                    index)

                # Mapping Successful
                partialmap = True

            # All Anchors Absent
            if not partialanc:
                exoread = read
                exofreq = freq
                continue

            # All Barcodes Absent
            if not partialmap:
                cctrs['incalcreads'] += freq
                continue

            # Convert Tuples to Indexes
            indextuple = tuple(indextuple)

            # Compute Callback Evaluation
            if not callback is None:
                try:
                    # Execute Callback
                    evaluation = callback(
                        read=anchoredread,
                        ID=tuple(IDdict[it] if not it == '-' else None for it in indextuple),
                        count=freq,
                        coreid=coreid)

                    # Uh oh ... Non-boolean Output
                    if not ((evaluation is True) or \
                            (evaluation is False)):
                        raise
                except:
                    # We have a failed evaluation!
                    callbackerror.set()
                    callbackabort = True

                    # Dump input for Stats
                    ut.savedump(dobj={
                          'read': anchoredread,
                            'ID': indextuple,
                         'count': freq,
                        'coreid': coreid},
                        filepath='{}/{}.callback.dump'.format(
                            countdir,
                            coreid))

                    # Break out!
                    break

                else:
                    # Callback Evaluation is False
                    if not evaluation:
                        cctrs['falsereads'] += freq
                        continue

            # All Components Valid
            countdict[indextuple]    += freq
            cctrs['experimentreads'] += freq

        # Show Final Updates
        liner.send(
            ' Core {:{},d}: Counted Pack {} w/ {:{},d} Reads in {:05.2f} sec\n'.format(
                coreid,
                clen,
                packname,
                readcount,
                plen,
                tt.time()-t0))

        # Free Memory
        ut.free_mem()

        # Release Control
        tt.sleep(0)

        # Did we Abort due to Callback?
        if callbackabort:
            shutdown.set()
            break

        # Need to Restart?
        if ut.needs_restart(
            memlimit=memlimit):
            restart.set() # Enable Restart
            break # Release, your Memory Real Estate

    # Pack Queue is Empty
    else:
        shutdown.set()


    # Close packfile
    packfile.close()

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

    # We didn't Abort right?
    if not callbackabort:

        # Do we have counts?
        if countdict:

            # Save Count Dictionary
            ut.savecount(
                cobj=countdict,
                filepath=countpath)

            # Release Control
            tt.sleep(0)

            # Queue Count Dicionary Path
            countqueue.put(
                countpath.split('/')[-1])

            # Release Control
            tt.sleep(0)

    # Update Read Counting Book-keeping
    previousreads.increment(incr=cctrs['analyzedreads'])
    analyzedreads.increment(incr=cctrs['analyzedreads'])
    phiXreads.increment(incr=cctrs['phiXreads'])
    lowcomplexreads.increment(incr=cctrs['lowcomplexreads'])
    falsereads.increment(incr=cctrs['falsereads'])
    incalcreads.increment(incr=cctrs['incalcreads'])
    experimentreads.increment(incr=cctrs['experimentreads'])

    # Counting Completed
    nactive.decrement()
    if shutdown.is_set():
        countqueue.put(None)

    # Release Control
    tt.sleep(0)
