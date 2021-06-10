import time as tt
import gc

import collections as cx
import random      as rn
import atexit      as ae

import ctypes                       as ct
import multiprocessing              as mp
import multiprocessing.sharedctypes as mpsct

import utils    as ut
import valparse as vp
import corepack as cp


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
    packsize,
    r1truncfile,
    r2truncfile,
    scannedreads,
    ambiguousreads,
    shortreads,
    survivedreads,
    packedreads,
    packsbuilt,
    coreid,
    ncores,
    nactive,
    liner):
    '''
    TBD
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

    verbagetarget   = rn.randint(
        *map(round, (min_dump_target * 0.080,
                     min_dump_target * 0.120)))

    clen = len(str(ncores))
    plen = len(str(packsize)) + \
        int(ut.safelog10(packsize) / 3)

    # Current and Meta Pack Storage
    cpack = cx.Counter()
    mpack = {}
    mpath = None

    # Read Packing Loop
    for reads in cp.stream_processed_fastq(
        r1file=r1file,
        r2file=r2file,
        r1type=r1type,
        r2type=r2type,
        r1length=r1length,
        r2length=r2length,
        r1qual=r1qual,
        r2qual=r2qual,
        packtype=packtype,
        r1truncfile=r1truncfile,
        r2truncfile=r2truncfile,
        scannedreads=scannedreads,
        ambiguousreads=ambiguousreads,
        shortreads=shortreads,
        survivedreads=survivedreads,
        verbagetarget=verbagetarget,
        coreid=coreid,
        ncores=ncores,
        liner=liner):

        # Truncated Files?
        if  r1truncfile.is_set() or \
            r2truncfile.is_set():
            truncated = True
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
                ' Core {:{},d}: Built Pack {}.{} w/ {:{},d} Reads in {:.2f} sec\n'.format(
                    coreid,
                    clen,
                    coreid,
                    c_packsbuilt,
                    len(cpack),
                    plen,
                    tt.time()-t0))

            # Define Pack Path
            cpath = '{}/{}.{}.pack'.format(
                packdir,
                coreid,
                c_packsbuilt)

            # Dump Current Pack or Setup Meta Pack
            if mpack:
                ut.savepack(
                    pobj=cpack,
                    filepath=cpath)
                packqueue.put(cpath)
            else:
                mpack = cpack
                mpath = cpath

            # Update Book-keeping
            t0 = tt.time()
            c_packedreads  += len(cpack)
            c_packsbuilt   += 1
            max_dump_reach  = 0
            min_dump_reach  = 0

            # Setup Next Pack
            del cpack
            cpack = cx.Counter()
            gc.collect()

    # Final Dumping
    if not truncated:

        # Dump Final Pack?
        if cpack:

            # Show Updates
            liner.send(
                ' Core {:{},d}: Built Pack {}.{} w/ {:{},d} Reads in {:05.2f} sec\n'.format(
                    coreid,
                    clen,
                    coreid,
                    c_packsbuilt,
                    len(cpack),
                    plen,
                    tt.time()-t0))

            # Dump Read Pack Appropriately

            if c_packsbuilt == 0: # Meta pack

                # Define Pack Path
                cpath = '{}/{}.{}.pack.meta'.format(
                    packdir,
                    coreid,
                    c_packsbuilt)

                # Dump Meta Pack
                ut.savemeta(
                    pobj=cpack,
                    filepath=cpath)

                # Enqueue Meta Pack Path
                metaqueue.put(cpath)

            else:                 # Non-Meta Pack

                # Define Pack Path
                cpath = '{}/{}.{}.pack'.format(
                    packdir,
                    coreid,
                    c_packsbuilt)

                # Dump Pack
                ut.savepack(
                    pobj=cpack,
                    filepath=cpath)

                # Enqueue Pack Path
                packqueue.put(cpath)

            # Update Book-keeping
            c_packedreads += len(cpack)
            c_packsbuilt   = max(1, c_packsbuilt)

            # Cleanup Final Pack
            del cpack

        # Dump Meta Pack?
        if mpack:

            # # Dump Meta Pack
            # ut.savepack(
            #     pobj=mpack,
            #     filepath=mpath)
            # packqueue.put(mpath)

            # Define Pack Path
            mpath = '{}.meta'.format(
                mpath)

            # Dump Meta Pack
            ut.savemeta(
                pobj=mpack,
                filepath=mpath)

            # Enqueue Meta Pack Path
            metaqueue.put(mpath)

            # Update Book-keeping
            c_packsbuilt += 1

            # Cleanup Meta Pack
            del mpack

    # Update Read Packing Book-keeping
    packedreads.increment(incr=c_packedreads)
    packsbuilt.increment(incr=c_packsbuilt)
    nactive.decrement()

    # Packing Completed
    packqueue.put(None)

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
            ' Reducing {} in Progress'.format(
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

def pack(
    r1file,
    r1type,
    packfile,
    r1length=1,
    r1qual=30,
    r2file=None,
    r2type=None,
    r2length=None,
    r2qual=None,
    packtype=0,
    packsize=3.0,
    ncores=0,
    verbose=True):
    '''
    TBD
    '''

    # Start Liner
    liner = ut.liner_engine(online=verbose)

    # Packing Verbage Print
    liner.send('\n[Oligopool Calculator: Analysis Mode - Pack]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # Full r1file Validation
    r1valid = vp.get_readfile_validity(
        readfile=r1file,
        readfile_field='   R1 File   ',
        paired_readfile=None,
        liner=liner)

    # Full r1type Validation
    t1valid = vp.get_categorical_validity(
        category=r1type,
        category_field='   R1 Type   ',
        category_pre_desc=' R1 has ',
        category_post_desc='',
        category_dict={
            0: 'Forward Reads 5\' ---F--> 3\'',
            1: 'Reverse Reads 3\' <--R--- 5\''},
        liner=liner)

    # Full packfile Validation
    packfile_valid = vp.get_outfile_validity(
        outfile=packfile,
        outfile_suffix='.oligopool.pack',
        outfile_field=' Pack File   ',
        liner=liner)

    # Adjust packfile Suffix
    packfile = ut.get_adjusted_path(
        path=packfile,
        suffix='.oligopool.pack')

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full r1length Validation
    l1valid = vp.get_numeric_validity(
        numeric=r1length,
        numeric_field='   R1 Length ',
        numeric_pre_desc=' Use Reads of Length ',
        numeric_post_desc=' bp or Longer',
        minval=1,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full r1qual Validation
    q1valid = vp.get_numeric_validity(
        numeric=r1qual,
        numeric_field='   R1 Quality',
        numeric_pre_desc=' Use Reads w/ Mean Q-Score of ',
        numeric_post_desc=' or Higher',
        minval=30,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full r2file Validation
    r2valid = vp.get_optional_readfile_validity(
        readfile=r2file,
        readfile_field='   R2 File   ',
        paired_readfile=r1file,
        liner=liner)

    # Full r2type Validation
    if r2valid and (not r2file is None):
        validfn = vp.get_categorical_validity
    else:
        validfn = vp.get_optional_categorical_validity
    t2valid = validfn(
        category=r2type,
        category_field='   R2 Type   ',
        category_pre_desc=' R2 has ',
        category_post_desc='',
        category_dict={
            0: 'Forward Reads 5\' ---F--> 3\'',
            1: 'Reverse Reads 3\' <--R--- 5\''},
        liner=liner)

    # Full r2length Validation
    if r2valid and (not r2file is None):
        validfn = vp.get_numeric_validity
    else:
        validfn = vp.get_optional_numeric_validity
    l2valid = validfn(
        numeric=r2length,
        numeric_field='   R2 Length ',
        numeric_pre_desc=' Use Reads of Length ',
        numeric_post_desc=' bp or Longer',
        minval=1,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full r2qual Validation
    if r2valid and (not r2file is None):
        validfn = vp.get_numeric_validity
    else:
        validfn = vp.get_optional_numeric_validity
    q2valid = validfn(
        numeric=r2qual,
        numeric_field='   R2 Quality',
        numeric_pre_desc=' Use Reads w/ Mean Q-Score of ',
        numeric_post_desc=' or Higher',
        minval=30,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full packtype Validation
    packtype_valid = vp.get_categorical_validity(
        category=packtype,
        category_field=' Pack Type   ',
        category_pre_desc=' ',
        category_post_desc='',
        category_dict={
            0: 'Store Concatenated / Joined Reads',
            1: 'Store Assembled / Merged Reads'},
        liner=liner)

    # Full packsize Validation
    packsize_valid = vp.get_numeric_validity(
        numeric=packsize,
        numeric_field=' Pack Size   ',
        numeric_pre_desc=' Store up to ',
        numeric_post_desc=' Million Reads per Pack',
        minval=0.10,
        maxval=10.0,
        precheck=False,
        liner=liner)

    # Full num_core Parsing and Validation
    (ncores,
    ncores_valid) = vp.get_parsed_core_info(
        ncores=ncores,
        core_field='  Num Cores  ',
        default=None if packtype == 0 else float('inf'),
        liner=liner)

    # First Pass Validation
    if not all([
        r1valid,
        t1valid,
        packfile_valid,
        l1valid,
        q1valid,
        r2valid,
        t2valid,
        l2valid,
        q2valid,
        packtype_valid,
        packsize_valid,
        ncores_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Setup Warning Dictionary
    warns = {}

    # Setup Workspace
    (packfile,
    packdir) = ut.setup_workspace(
        outfile=packfile,
        outfile_suffix='.oligopool.pack')

    # Schedule packfile deletion
    pkdeletion = ae.register(
        ut.remove_file,
        packfile)

    # Expand packsize
    packsize = int(packsize * (10.**6))

    # Read Pack Queues
    metaqueue = mp.SimpleQueue()
    packqueue = mp.SimpleQueue()

    # Truncated Read Files Event
    r1truncfile = mp.Event()
    r2truncfile = mp.Event()

    # Read Packing Book-keeping
    nactive        = ut.SafeCounter(initval=ncores)
    scannedreads   = ut.SafeCounter()
    ambiguousreads = ut.SafeCounter()
    shortreads     = ut.SafeCounter()
    survivedreads  = ut.SafeCounter()
    packedreads    = ut.SafeCounter()
    packsbuilt     = ut.SafeCounter()

    # Launching Read Packing
    liner.send('\n[Step 1: Computing Read Packs]\n')

    # Define Packing Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 1,
        'stepname': 'computing-read-packs',
        'vars'    : {
            'r1truncated': False,
            'r2truncated': False,
               'packsize': packsize,
              'packcount': int(packsbuilt.value()),
                'scanned': int(scannedreads.value()),
               'survived': int(survivedreads.value()),
                 'packed': int(packedreads.value())},
        'warns'   : warns}

    # Packer Process Store
    readpackers = []

    # Fire off Single-Core Read Packer
    if ncores == 1:

        pack_engine(
            r1file=r1file,
            r2file=r2file,
            r1type=r1type,
            r2type=r2type,
            r1length=r1length,
            r2length=r2length,
            r1qual=r1qual,
            r2qual=r2qual,
            packdir=packdir,
            metaqueue=metaqueue,
            packqueue=packqueue,
            packtype=packtype,
            packsize=packsize,
            r1truncfile=r1truncfile,
            r2truncfile=r2truncfile,
            scannedreads=scannedreads,
            ambiguousreads=ambiguousreads,
            shortreads=shortreads,
            survivedreads=survivedreads,
            packedreads=packedreads,
            packsbuilt=packsbuilt,
            coreid=0,
            ncores=ncores,
            nactive=nactive,
            liner=liner)

    # Fire off Multi-Core Read Packers
    else:

        # Queue All Read Packers
        coreid = 0
        while coreid < ncores:

            # Define Packer
            readpacker = mp.Process(
                target=pack_engine,
                args=(r1file,
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
                    packsize,
                    r1truncfile,
                    r2truncfile,
                    scannedreads,
                    ambiguousreads,
                    shortreads,
                    survivedreads,
                    packedreads,
                    packsbuilt,
                    coreid,
                    ncores,
                    nactive,
                    liner,))

            # Start Packer
            readpacker.start()

            # Update Book-keeping
            readpackers.append(readpacker)
            coreid += 1

    # Archive Non-Meta Read Packs
    ut.archive(
        objqueue=packqueue,
        arcfile=packfile,
        prodcount=ncores,
        prodactive=nactive,
        liner=liner)

    # Join All Read Packers
    for readpacker in readpackers:
        readpacker.join()

    # Handle Truncated Read Files
    if r1truncfile.is_set():
        liner.send('\n')
        ut.remove_file(
            filepath=packfile)
        raise RuntimeError(
            'R1 File is Truncated or Incompatible with R2 File.')

    if r2truncfile.is_set():
        liner.send('\n')
        ut.remove_file(
            filepath=packfile)
        raise RuntimeError(
            'R2 File is Truncated or Incompatible with R1 File.')

    # Aggregate Meta Read Packs
    pack_aggregator(
        metaqueue=metaqueue,
        packqueue=packqueue,
        packedreads=packedreads,
        packsbuilt=packsbuilt,
        liner=liner)

    # Define Archiver
    archiver = ut.archive_engine(
        arcfile=packfile,
        mode='a')

    # Archive Meta Read Packs
    while not packqueue.empty():
        xpath = packqueue.get()
        archiver.send((xpath, liner))

    # Packing Status
    if packsbuilt.value() > 0:
        packstatus = (True,  'Successful')
    else:
        packstatus = (False, 'Failed')

    # Read Packing Stats
    liner.send('\n[Packing Stats]\n')

    plen = ut.get_printlen(
        value=scannedreads.value())

    liner.send(
        '   Packing  Status: {}\n'.format(
            packstatus[1]))
    liner.send(
        '   Scanned   Reads: {:{},d}\n'.format(
            scannedreads.value(),
            plen))
    liner.send(
        ' Ambiguous   Reads: {:{},d} ({:6.2f} %)\n'.format(
            ambiguousreads.value(),
            plen,
            ut.safediv(
                A=(100. * ambiguousreads.value()),
                B=scannedreads.value())))
    liner.send(
        '     Short   Reads: {:{},d} ({:6.2f} %)\n'.format(
            shortreads.value(),
            plen,
            ut.safediv(
                A=(100. * shortreads.value()),
                B=scannedreads.value())))
    liner.send(
        '  Survived   Reads: {:{},d} ({:6.2f} %)\n'.format(
             survivedreads.value(),
            plen,
            ut.safediv(
                A=(100. * survivedreads.value()),
                B=scannedreads.value())))
    liner.send(
        '    Packed   Reads: {:{},d}\n'.format(
            packedreads.value(),
            plen))
    liner.send(
        '   Packing   Ratio: {:.2f} to 1\n'.format(
            ut.safediv(
                A=(1. * survivedreads.value()),
                B=packedreads.value())))
    liner.send(
        '   Packing   Order: {:.2f} %\n'.format(
            (100. * (
                1. - ut.safediv(
                    A=packedreads.value(),
                    B=survivedreads.value())))))
    liner.send(
        '     Packs   Built: {} Read Packs\n'.format(
            packsbuilt.value()))
    liner.send(
        '      Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - t0))

    # Archive Packing Stats
    packstat = {
        'packtype' : packtype,
        'packsize' : packsize,
        'packcount': int(packsbuilt.value()),
        'scanned'  : int(scannedreads.value()),
        'survived' : int(survivedreads.value()),
        'packed'   : int(packedreads.value())
    }

    statpath = '{}/packing.stat'.format(
        packdir)

    ut.savedict(
        dobj=packstat,
        filepath=statpath)

    archiver.send((statpath, None))

    # Close Archiver
    archiver.close()

    # # Remove Workspace
    ut.remove_directory(
        dirpath=packdir)

    # Unschedule packfile deletion
    if packstatus[0]:
        ae.unregister(pkdeletion)

    # Close Liner
    liner.close()