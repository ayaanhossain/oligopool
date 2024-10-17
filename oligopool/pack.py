import time as tt

import atexit as ae

import multiprocessing as mp

from .base import utils as ut
from .base import validation_parsing as vp
from .base import core_pack as cp


def pack(
    r1_fastq_file:str,
    r1_read_type:int,
    pack_type:int,
    pack_file:str,
    minimum_r1_read_length:int=1,
    minimum_r1_read_quality:int=20,
    r2_fastq_file:str|None=None,
    r2_read_type:int|None=None,
    minimum_r2_read_length:int|None=None,
    minimum_r2_read_quality:int|None=None,
    pack_size:float=3.0,
    core_count:int=0,
    memory_limit:float=0.0,
    verbose:bool=True) -> dict:
    '''
    Pack NGS reads with indexed barcodes and optional associated variants. Processes, filters, and
    optionally merges reads before packing. The resulting pack file represents a single characterization
    experiment, outcome, or replicate for subsequent index-based counting.

    Required Parameters:
        - `r1_fastq_file` (`str`): Path to R1 FastQ file (may be gzipped).
        - `r1_read_type` (`int`): Orientation of R1 reads (0 for forward, 1 for reverse).
        - `pack_type` (`int`): Packing storage type (0 for concatenated, 1 for merged)
        - `pack_file` (`str`): Filename for output pack file.

    Optional Parameters:
        - `minimum_r1_read_length` (`int`): Minimum R1 read length (≥ 1).
        - `minimum_r1_read_quality` (`int`): Minimum average R1 quality (default: 20).
        - `r2_fastq_file` (`str`): Path to R2 FastQ file (default: `None`).
        - `r2_read_type` (`int`): Orientation of R2 reads (0 for forward, 1 for reverse) (default: `None`).
        - `minimum_r2_read_length` (`int`): Minimum R2 read length (default: `None`).
        - `minimum_r2_read_quality` (`int`): Minimum average R2 quality (default: `None`).
        - `pack_size` (`float`): Million unique reads per pack (default: 3.0, range: 0.1 to 5.0).
        - `core_count` (`int`): CPU cores to use (0: auto-infer, default: 0).
        - `memory_limit` (`float`): GB of memory per core (0: auto-infer, default: 0)
        - `verbose` (`bool`): If `True`, logs updates to stdout (default: `True`).

    Returns:
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - For single-end reads, use R1 arguments only; set R2 arguments to None.
        - Read quality is average Phred score.
        - Both reads must pass criteria for paired-end acceptance.
        - Concatenated storage (`pack_type` = 0) is IO bound; use 2 cores.
        - Pack size balances memory usage and computation speed.
    '''

    # Alias Arguments
    r1file   = r1_fastq_file
    r2file   = r2_fastq_file
    r1type   = r1_read_type
    r2type   = r2_read_type
    r1length = minimum_r1_read_length
    r2length = minimum_r2_read_length
    r1qual   = minimum_r1_read_quality
    r2qual   = minimum_r2_read_quality
    packtype = pack_type
    packsize = pack_size
    packfile = pack_file
    ncores   = core_count
    memlimit = memory_limit
    verbose  = verbose

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
        minval=0,
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
        minval=0,
        maxval=float('inf'),
        precheck=False,
        liner=liner)

    # Full packsize Validation
    packsize_valid = vp.get_numeric_validity(
        numeric=packsize,
        numeric_field=' Pack Size   ',
        numeric_pre_desc=' Store up to ',
        numeric_post_desc=' Million Reads per Pack',
        minval=0.10,
        maxval=5.00,
        precheck=False,
        liner=liner)

    # Full num_core Parsing and Validation
    (ncores,
    ncores_valid) = vp.get_parsed_core_info(
        ncores=ncores,
        core_field='  Num Cores  ',
        default=None if not packtype else mp.cpu_count() // 3,
        offset=2,
        liner=liner)

    # Full num_core Parsing and Validation
    (memlimit,
    memlimit_valid) = vp.get_parsed_memory_info(
        memlimit=memlimit,
        memlimit_field='  Mem Limit  ',
        ncores=ncores,
        ncores_valid=ncores_valid,
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
        ncores_valid,
        memlimit_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Adjust Numeric Parameters
    if r2file is None:
        packtype = 0

    # Setup Warning Dictionary
    warns = {}

    # Some Assembly Required?
    if packtype:

        # Compute Assembly Parameters
        liner.send('\n[Step 1: Extracting Assembly Parameters]\n')

        # Extract Overlap Parameters
        (mergefnhigh,
        mergefnlow) = cp.get_extracted_overlap_parameters(
            r1file=r1file,
            r2file=r2file,
            r1type=r1type,
            r2type=r2type,
            liner=liner)

        # Store Parameters
        assemblyparams = {
            'mergefnhigh': mergefnhigh,
             'mergefnlow': mergefnlow}

    # Vanilla Storage
    else:
        # Nothing to see here ...
        assemblyparams = None

    # Free Memory
    ut.free_mem()

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
    metaqueue = ut.SafeQueue()
    packqueue = ut.SafeQueue()

    # Read File Processing Book-keeping
    r1truncfile  = mp.Event()
    r2truncfile  = mp.Event()
    restarts  = [
        mp.Event() for _ in range(ncores)]
    shutdowns = [
        mp.Event() for _ in range(ncores)]

    # Read Packing Book-keeping
    nactive        = ut.SafeCounter(initval=ncores)
    scannedreads   = ut.SafeCounter()
    ambiguousreads = ut.SafeCounter()
    shortreads     = ut.SafeCounter()
    survivedreads  = ut.SafeCounter()
    packedreads    = ut.SafeCounter()
    packsbuilt     = ut.SafeCounter()
    batchids       = [0] * ncores
    previousreads  = [
        ut.SafeCounter() for _ in range(ncores)]

    # Launching Read Packing
    liner.send('\n[Step 2: Computing Read Packs]\n')

    # Define Packing Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 2,
        'step_name': 'computing-read-packs',
        'vars'    : {
             'meta_agg_count': 0,
               'r1_truncated': False,
               'r2_truncated': False,
                  'pack_size': packsize,
                 'pack_count': int(packsbuilt.value()),
              'scanned_reads': int(scannedreads.value()),
            'ambiguous_reads': int(ambiguousreads.value()),
                'short_reads': int(shortreads.value()),
             'survived_reads': int(survivedreads.value()),
               'packed_reads': int(packedreads.value())},
        'warns'   : warns}

    # Engine Timer
    et = tt.time()

    # Define Archiver (Non-Meta Read Packs)
    archiver = mp.Process(
        target=ut.archive,
        args=(packqueue,
            packfile,
            'x',
            ncores,
            nactive,
            liner,))

    # Fire-off Archiver
    archiver.start()

    # Define Packer Process Store
    readpackers = []

    # Fire-off Initial Read Packers
    coreid = 0
    clen = ut.get_printlen(value=ncores)
    while coreid < ncores:

        # Define Packer
        readpacker = mp.Process(
            target=cp.pack_engine,
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
                assemblyparams,
                packsize,
                r1truncfile,
                r2truncfile,
                previousreads[coreid],
                scannedreads,
                ambiguousreads,
                shortreads,
                survivedreads,
                packedreads,
                packsbuilt,
                coreid,
                batchids[coreid],
                ncores,
                nactive,
                memlimit,
                et,
                restarts[coreid],
                shutdowns[coreid],
                liner,))

        # Show Update
        liner.send(
            ' Core {:{},d}: Starting Up\n'.format(
                coreid,
                clen))

        # Start Packer
        readpacker.start()

        # Update Book-keeping
        readpackers.append(readpacker)
        coreid += 1

    # Packer Management
    coreid = 0
    activepackers = ncores
    while activepackers:

        # Had Packer Finished?
        if readpackers[coreid] is None:
            pass

        # Has Packer Shutdown?
        elif shutdowns[coreid].is_set():
            # Cleanup
            readpackers[coreid].join()
            readpackers[coreid].close()
            # Update
            readpackers[coreid] = None
            activepackers -= 1
            ut.free_mem()
            # Reset
            restarts[coreid].clear()
            shutdowns[coreid].clear()

        # Must Packer Restart?
        elif restarts[coreid].is_set():
            # Cleanup
            readpackers[coreid].join()
            readpackers[coreid].close()
            # Update
            readpackers[coreid] = None
            batchids[coreid] += 1
            ut.free_mem()
            # Reset
            restarts[coreid].clear()
            shutdowns[coreid].clear()
            readpackers[coreid] = mp.Process(
                target=cp.pack_engine,
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
                    assemblyparams,
                    packsize,
                    r1truncfile,
                    r2truncfile,
                    previousreads[coreid],
                    scannedreads,
                    ambiguousreads,
                    shortreads,
                    survivedreads,
                    packedreads,
                    packsbuilt,
                    coreid,
                    batchids[coreid],
                    ncores,
                    nactive,
                    memlimit,
                    et,
                    restarts[coreid],
                    shutdowns[coreid],
                    liner,))
            readpackers[coreid].start()

        # Next Iteration
        coreid = (coreid + 1) % ncores
        tt.sleep(0)

    # Join Archiver
    archiver.join()
    archiver.close()

    # Free Memory
    ut.free_mem()

    # Handle Truncated Read Files
    if r1truncfile.is_set():
        liner.send(
            ' R1 File Truncated or Incompatible with R2 File\n')
        ut.remove_file(
            filepath=packfile)
        stats['vars']['r1_truncated'] = True
    elif r2truncfile.is_set():
        liner.send(
            ' R2 File Truncated or Incompatible with R1 File\n')
        ut.remove_file(
            filepath=packfile)
        stats['vars']['r2_truncated'] = True

    # Packing Successful
    else:
        stats['status'] = True
        stats['basis']  = 'solved'

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - et))

    # Did we succeed?
    if stats['status']:

        # Launching Read Packing
        liner.send('\n[Step 3: Aggregating Meta Packs]\n')

        # Update Stats
        stats['step'] = 3
        stats['step_name'] = 'aggregating-meta-packs'
        stats['vars']['meta_agg_count'] = len(metaqueue)

        # Aggregation Timer
        at = tt.time()

        # Aggregate Meta Read Packs
        # Define Meta Pack Aggregator
        # on a Separate Core
        aggregator = mp.Process(
            target=cp.pack_aggregator,
            args=(metaqueue,
                packqueue,
                packedreads,
                packsbuilt,
                liner,))

        # Start Aggregator
        aggregator.start()

        # Update Producer Count
        nactive.increment()

        # Archive Meta Read Packs
        ut.archive(
            objqueue=packqueue,
            arcfile=packfile,
            mode='a',
            prodcount=1,
            prodactive=nactive,
            liner=liner)

        # Update Producer Count
        nactive.decrement()

        # Join and Close Aggregator
        aggregator.join()
        aggregator.close()

        # Show Time Elapsed
        liner.send(
            ' Time Elapsed: {:.2f} sec\n'.format(
                tt.time() - at))

    # Update Stats
    stats['vars']['pack_count'] = int(packsbuilt.value())
    stats['vars']['scanned_reads']   = int(scannedreads.value())
    stats['vars']['ambiguous_reads'] = int(ambiguousreads.value())
    stats['vars']['short_reads']     = int(shortreads.value())
    stats['vars']['survived_reads']  = int(survivedreads.value())
    stats['vars']['packed_reads']    = int(packedreads.value())

    # Packing Status
    if stats['status']:
        packstatus = 'Successful'
    else:
        packstatus = 'Failed'

    # Read Packing Stats
    liner.send('\n[Packing Statistics]\n')

    plen = ut.get_printlen(
        value=scannedreads.value())

    liner.send(
        '   Packing Status: {}\n'.format(
            packstatus))
    liner.send(
        '     R1 Truncated: {}\n'.format(
            ['No', 'Yes'][r1truncfile.is_set()]))
    liner.send(
        '     R2 Truncated: {}\n'.format(
            ['No', 'Yes'][r2truncfile.is_set()]))
    liner.send(
        '   Scanned Reads : {:{},d}\n'.format(
            scannedreads.value(),
            plen))
    liner.send(
        ' Ambiguous Reads : {:{},d} ({:6.2f} %)\n'.format(
            ambiguousreads.value(),
            plen,
            ut.safediv(
                A=(100. * ambiguousreads.value()),
                B=scannedreads.value())))
    liner.send(
        '     Short Reads : {:{},d} ({:6.2f} %)\n'.format(
            shortreads.value(),
            plen,
            ut.safediv(
                A=(100. * shortreads.value()),
                B=scannedreads.value())))
    liner.send(
        '  Survived Reads : {:{},d} ({:6.2f} %)\n'.format(
             survivedreads.value(),
            plen,
            ut.safediv(
                A=(100. * survivedreads.value()),
                B=scannedreads.value())))
    liner.send(
        '    Packed Reads : {:{},d}\n'.format(
            packedreads.value(),
            plen))
    liner.send(
        '   Packing Ratio : {:.2f} to 1\n'.format(
            ut.safediv(
                A=(1. * survivedreads.value()),
                B=packedreads.value())))
    liner.send(
        '   Packing Order : {:.2f} %\n'.format(
            (100. * (
                1. - ut.safediv(
                    A=packedreads.value(),
                    B=survivedreads.value())))))
    liner.send(
        '     Packs Built : {} Read Packs\n'.format(
            packsbuilt.value()))
    liner.send(
        '  Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - t0))

    # Did we succeed?
    if stats['status']:

        # Archive Packing Stats
        packstat = {
            'pack_type' : packtype,
            'pack_size' : packsize,
            'pack_count': int(packsbuilt.value()),
            'scanned_reads'  : int(scannedreads.value()),
            'survived_reads' : int(survivedreads.value()),
            'packed_reads'   : int(packedreads.value())}

        # Prepare Archiving Sentinels
        statpath = '{}/packing.stat'.format(
            packdir)
        packqueue.put(statpath)
        packqueue.put(None)

        # Dump Pack Stat
        ut.savedict(
            dobj=packstat,
            filepath=statpath)

        # Update Producer Count
        nactive.increment()

        # Archive Pack Stat
        ut.archive(
            objqueue=packqueue,
            arcfile=packfile,
            mode='a',
            prodcount=1,
            prodactive=nactive,
            liner=liner)

        # Update Producer Count
        nactive.decrement()

        # Close Queues
        metaqueue.close()
        packqueue.close()

    # Remove Workspace
    ut.remove_directory(
        dirpath=packdir)

    # Unschedule packfile deletion
    if packstatus == 'Successful':
        ae.unregister(pkdeletion)

    # Close Liner
    liner.close()

    # Return Statistics
    return stats
