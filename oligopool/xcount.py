import time as tt

import atexit as ae
import multiprocess as mp

import pandas as pd

from .base import utils as ut
from .base import validation_parsing as vp
from .base import core_count as cc

from typing import Callable, Tuple


def xcount(
    index_files:str,
    pack_file:str,
    count_file:str,
    mapping_type:int|str=0,
    barcode_errors:int=-1,
    callback:Callable[[str, str|None, Tuple, int, int], bool]|None=None,
    core_count:int=0,
    memory_limit:float=0.0,
    failed_reads_file:str|None=None,
    failed_reads_sample_size:int=1000,
    verbose:bool=True) -> Tuple[pd.DataFrame, dict]:
    '''
    Count barcodes (single or combinatorial) using one or more indices; no associate verification.
    Writes a count matrix to disk and returns it as a DataFrame (callbacks available via Python API).

    Required Parameters:
        - `index_files` (`str` / `list`): A single (or a list of) index filename(s).
        - `pack_file` (`str`): Pack file path.
        - `count_file` (`str`): Output count matrix filename.

    Optional Parameters:
        - `mapping_type` (`int` / `str`): Barcode classification mode (default: 0). See Notes.
        - `barcode_errors` (`int`): Maximum errors in barcodes (-1: auto-infer, default: -1).
        - `callback` (`callable`): Custom read processing function (default: `None`).
        - `core_count` (`int`): CPU cores to use (0: auto-infer, default: 0).
        - `memory_limit` (`float`): GB of memory per core (0: auto-infer, default: 0)
        - `failed_reads_file` (`str` / `None`): Output CSV path for failed read samples (None: disabled, default: `None`).
        - `failed_reads_sample_size` (`int`): Maximum samples per failure category (default: 1000).
        - `verbose` (`bool`): If `True`, logs updates to stdout (default: `True`).

    Returns:
        - A pandas DataFrame of barcode combination counts.
        - A dictionary of stats from the last step in pipeline.

    Notes:
        - Partial and missing combinations are included in counts.
        - Reads are retained if at least one barcode maps; missing barcodes are represented as
            `'-'` (gaps) in the output combination.
        - `mapping_type`:
            0 or 'fast' for fast, 1 or 'sensitive' for sensitive
            (aliases: 'quick', 'near-exact', 'sens', 'accurate', 'slow').
        - CLI note: callback functions are not currently supported via the `op`/`oligopool` CLI
            (the CLI always runs with `callback=None`); use the Python API to supply callbacks.
        - Callback signature: `callback(r1, r2, ID, count, coreid) -> bool` (return `True` to accept the read);
            `r2` is `None` for merged/single-end reads.
        - Associate information in indexes is ignored.
        - Barcodes can be isolated or be sub-barcodes of a larger combinatorial assembly.
        - If an anchor appears multiple times, the best-scoring one is used; ties with multiple barcodes rejected.
        - Failed reads sampling collects representative samples from each failure category for diagnostics.
            Categories: phix_match, low_complexity, anchor_missing, barcode_absent, barcode_ambiguous,
            callback_false, incalculable.
    '''

    # Alias Arguments
    indexfiles = index_files
    packfile   = pack_file
    countfile  = count_file
    maptype    = mapping_type
    barcodeerrors = barcode_errors
    callback   = callback
    ncores     = core_count
    memlimit   = memory_limit
    failedreadsfile = failed_reads_file
    failedreadssamplesize = failed_reads_sample_size
    verbose    = verbose

    # Start Liner
    liner = ut.liner_engine(online=verbose)

    # Counting Verbiage Print
    liner.send('\n[Oligopool Calculator: Analysis Mode - Combinatorial Count]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # Full indexfiles Validation
    indexfiles_valid = vp.get_indexfiles_validity(
        indexfiles=indexfiles,
        indexfiles_field='     Index File(s)',
        associated=False,
        liner=liner)

    # Normalize indexfiles for Combinatorial + Regular Usecases
    if isinstance(indexfiles,str):
        indexfiles = [indexfiles]

    # Adjust indexfile Suffix
    if indexfiles_valid:
        indexfiles = [ut.get_adjusted_path(
            path=indexfile,
            suffix='.oligopool.index') for indexfile in indexfiles]

    # Full packfile Validation
    (packfile_valid,
    packcount) = vp.get_parsed_packfile(
        packfile=packfile,
        packfile_field='      Pack File   ',
        liner=liner)
    input_rows = packcount if packfile_valid else 0

    # Adjust packfile Suffix
    if packfile_valid:
        packfile = ut.get_adjusted_path(
            path=packfile,
            suffix='.oligopool.pack')

    # Full countfile Validation
    countfile_valid = vp.get_outfile_validity(
        outfile=countfile,
        outfile_suffix='.oligopool.xcount.csv',
        outfile_field='     Count File   ',
        liner=liner)

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full maptype Validation
    (maptype,
    maptype_valid) = vp.get_typed_categorical_validity(
        category=maptype,
        category_field='   Mapping Type   ',
        category_pre_desc=' ',
        category_post_desc=' Classification',
        type_name='mapping_type',
        liner=liner)

    # Full barcodeerrors Validation
    (barcodeerrors,
    barcodeerrors_valid) = vp.get_errors_validity(
        errors=barcodeerrors,
        errors_field='   Barcode Errors ',
        errors_pre_desc=' At most ',
        errors_post_desc=' Mutations per Barcode',
        errors_base='B',
        indexfiles_valid=indexfiles_valid,
        indexfiles=indexfiles,
        liner=liner)

    # Full callback Validation
    callback_valid = vp.get_callback_validity(
        callback=callback,
        callback_field='  Callback Method ',
        liner=liner)

    # Full num_core Parsing and Validation
    (ncores,
    ncores_valid) = vp.get_parsed_core_info(
        ncores=ncores,
        core_field='       Num Cores  ',
        default=packcount,
        offset=2,
        liner=liner)

    # Full num_core Parsing and Validation
    (memlimit,
    memlimit_valid) = vp.get_parsed_memory_info(
        memlimit=memlimit,
        memlimit_field='       Mem Limit  ',
        ncores=ncores,
        ncores_valid=ncores_valid,
        liner=liner)

    # Full failedreadsfile Validation
    failedreadsfile_valid = vp.get_optional_outfile_validity(
        outfile=failedreadsfile,
        outfile_suffix='.oligopool.xcount.failed_reads.csv',
        outfile_field='    Failed Reads  ',
        liner=liner)

    # Adjust failedreadsfile Suffix
    if failedreadsfile_valid and failedreadsfile is not None:
        failedreadsfile = ut.get_adjusted_path(
            path=failedreadsfile,
            suffix='.oligopool.xcount.failed_reads.csv')

    # Full failedreadssamplesize Validation
    failedreadssamplesize_valid = True
    if failedreadsfile is not None:
        failedreadssamplesize_valid = vp.get_sample_size_validity(
            sample_size=failedreadssamplesize,
            sample_size_field='    Sample Size   ',
            liner=liner)
    else:
        liner.send('    Sample Size   : N/A (Disabled)\n')

    # First Pass Validation
    if not all([
        indexfiles_valid,
        packfile_valid,
        countfile_valid,
        maptype_valid,
        barcodeerrors_valid,
        callback_valid,
        ncores_valid,
        memlimit_valid,
        failedreadsfile_valid,
        failedreadssamplesize_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Counting Book-keeping
    stats = None
    outdf = None
    warns = {}

    # Parse Callback Function
    if not callback is None:
        liner.send('\n[Step 1: Parsing Callback Method]\n')

        # Parse callback
        (parsestatus,
        failedinputs) = cc.get_parsed_callback(
            indexfiles=indexfiles,
            callback=callback,
            packfile=packfile,
            ncores=ncores,
            liner=liner)

        # callback infeasible
        if not parsestatus:

            # Prepare stats
            stats = {
                'status'  : False,
                'basis'   : 'infeasible',
                'step'    : 1,
                'step_name': 'parsing-callback-method',
                'vars'    : {
                    'failed_inputs': failedinputs},
                'warns'   : warns}
            stats = ut.stamp_stats(
                stats=stats,
                module='xcount',
                input_rows=input_rows,
                output_rows=0)

            # Return results
            liner.close()
            return (outdf, stats)

    # Enqueue Read Packs
    liner.send('\n[Step 2: Enqueing Read Packs]\n')

    # Define Queing Sentinels
    packqueue = ut.SafeQueue()
    enqueuecomplete = mp.Event()

    # Define Pack Enquer
    packenqueuer = mp.Process(
        target=cc.pack_loader,
        args=(packfile,
            packqueue,
            enqueuecomplete,
            ncores,
            liner,))

    # Start Enquer
    packenqueuer.start()

    # Setup Workspace
    (countfile,
    countdir) = ut.setup_workspace(
        outfile=countfile,
        outfile_suffix='.oligopool.xcount.csv')

    # Schedule countfile deletion
    ctdeletion = ae.register(
        ut.remove_file,
        countfile)

    # Adjust Errors
    barcodeerrors = round(barcodeerrors)

    # Define countqueue
    countqueue = mp.SimpleQueue()

    # Define failuresamplequeue (if enabled)
    failuresamplequeue = None
    if failedreadsfile is not None:
        failuresamplequeue = ut.SafeQueue()

    # Pack File Processing Book-keeping
    callbackerror = mp.Event()
    restarts  = [
        mp.Event() for _ in range(ncores)]
    shutdowns = [
        mp.Event() for _ in range(ncores)]

    # Read Counting Book-keeping
    nactive         = ut.SafeCounter(initval=ncores)
    analyzedreads   = ut.SafeCounter()
    phiXreads       = ut.SafeCounter()
    lowcomplexreads = ut.SafeCounter()
    falsereads      = ut.SafeCounter()
    incalcreads     = ut.SafeCounter()
    experimentreads = ut.SafeCounter()
    batchids        = [0] * ncores
    previousreads   = [
        ut.SafeCounter() for _ in range(ncores)]

    # Wait on Enqueuing
    enqueuecomplete.wait()

    # Launching Read Counting
    liner.send('\n[Step 3: Counting Read Packs]\n')

    # Define Counting Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 3,
        'step_name': 'counting-read-packs',
        'vars'    : {
            'callback_error': False,
             'failed_inputs': None,
                  'analyzed_reads': int(analyzedreads.value()),
                      'phiX_reads': int(phiXreads.value()),
               'low_complex_reads': int(lowcomplexreads.value()),
            'callback_false_reads': int(falsereads.value()),
                    'incalc_reads': int(incalcreads.value()),
                'experiment_reads': int(experimentreads.value())},
        'warns'   : warns}

    # Engine Timer
    et = tt.time()

    # Define Aggregator
    assoc = False
    aggregator = mp.Process(
        target=cc.count_aggregator,
        args=(countqueue,
            countdir,
            ncores,
            nactive,
            assoc,
            liner,))

    # Fire-off Aggregator
    aggregator.start()

    # Define Counter Process Store
    readcounters = []

    # Fire-off Initial Read Counters
    coreid = 0
    clen = ut.get_printlen(value=ncores)
    while coreid < ncores:

        # Define Counter
        readcounter = mp.Process(
            target=cc.xcount_engine,
            args=(indexfiles,
                packfile,
                packqueue,
                countdir,
                countqueue,
                maptype,
                barcodeerrors,
                previousreads[coreid],
                analyzedreads,
                phiXreads,
                lowcomplexreads,
                falsereads,
                incalcreads,
                experimentreads,
                callback,
                callbackerror,
                coreid,
                batchids[coreid],
                ncores,
                nactive,
                memlimit,
                restarts[coreid],
                shutdowns[coreid],
                et,
                liner,
                failuresamplequeue,
                failedreadssamplesize,))

        # Show Update
        liner.send(
            ' Core {:{},d}: Starting Up\n'.format(
                coreid,
                clen))

        # Start Counter
        readcounter.start()

        # Update Book-keeping
        readcounters.append(readcounter)
        coreid += 1

    # Counter Management
    coreid = 0
    activecounters = ncores
    while activecounters:

        # Had Counter Finished?
        if readcounters[coreid] is None:
            pass

        # Has Counter Shutdown?
        elif shutdowns[coreid].is_set():
            # Cleanup
            readcounters[coreid].join()
            readcounters[coreid].close()
            # Update
            readcounters[coreid] = None
            activecounters -= 1
            ut.free_mem()
            # Reset
            restarts[coreid].clear()
            shutdowns[coreid].clear()

        # Must Counter Restart?
        elif restarts[coreid].is_set():
            # Cleanup
            readcounters[coreid].join()
            readcounters[coreid].close()
            # Update
            readcounters[coreid] = None
            batchids[coreid] += 1
            ut.free_mem()
            # Reset
            restarts[coreid].clear()
            shutdowns[coreid].clear()
            readcounters[coreid] = mp.Process(
                target=cc.xcount_engine,
                args=(indexfiles,
                    packfile,
                    packqueue,
                    countdir,
                    countqueue,
                    maptype,
                    barcodeerrors,
                    previousreads[coreid],
                    analyzedreads,
                    phiXreads,
                    lowcomplexreads,
                    falsereads,
                    incalcreads,
                    experimentreads,
                    callback,
                    callbackerror,
                    coreid,
                    batchids[coreid],
                    ncores,
                    nactive,
                    memlimit,
                    restarts[coreid],
                    shutdowns[coreid],
                    et,
                    liner,
                    failuresamplequeue,
                    failedreadssamplesize,))
            readcounters[coreid].start()

        # Next Iteration
        coreid = (coreid + 1) % ncores
        tt.sleep(0)

    # Join Aggregator
    aggregator.join()
    aggregator.close()

    # Free Memory
    ut.free_mem()

    # Handle Callback Error
    if callbackerror.is_set():
        failedinputs = cc.get_failed_inputs(
            packqueue=packqueue,
            countdir=countdir,
            liner=liner)
        liner.send(
            ' Callback Function Erroneous\n')
        ut.remove_file(
            filepath=countfile)
        stats['vars']['callback_error'] = True
        stats['vars']['failed_inputs']  = failedinputs

    # Handle Unmapped Reads
    elif experimentreads.value() <= 0:
        liner.send(
            ' No Reads Mapped Successfully\n')

    # Counting Successful
    else:
        stats['status'] = True
        stats['basis']  = 'solved'

    # Join Enqueuer
    packenqueuer.join()
    packenqueuer.close()

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - et))

    # Did we succeed?
    if stats['status']:

        # Launching Count Matrix Writing
        liner.send('\n[Step 4: Writing Count Matrix]\n')

        # Write Count Matrix
        outdf = cc.write_count(
            indexfiles=indexfiles,
            countdir=countdir,
            countfile=countfile,
            assoc=False,
            liner=liner)

        # Update Stats
        stats['step'] = 4
        stats['step_name'] = 'writing-count-matrix'

    # Aggregate failed reads samples if enabled (Step 5)
    if failuresamplequeue is not None:
        cc.aggregate_failure_samples(
            failuresamplequeue=failuresamplequeue,
            countdir=countdir,
            failed_reads_file=failedreadsfile,
            failed_reads_sample_size=failedreadssamplesize,
            prodcount=ncores,
            liner=liner)
        stats['step'] = 5
        stats['step_name'] = 'writing-failed-reads'

    # Update Stats
    stats['vars']['analyzed_reads']       = int(analyzedreads.value())
    stats['vars']['phiX_reads']           = int(phiXreads.value())
    stats['vars']['low_complex_reads']    = int(lowcomplexreads.value())
    stats['vars']['callback_false_reads'] = int(falsereads.value())
    stats['vars']['incalc_reads']         = int(incalcreads.value())
    stats['vars']['experiment_reads']     = int(experimentreads.value())

    # Counting Status
    if stats['status']:
        countstatus = 'Successful'
    else:
        countstatus = 'Failed'

    # Read Counting Stats
    liner.send('\n[Combinatorial Counting Stats]\n')

    plen = ut.get_printlen(
        value=analyzedreads.value())

    liner.send(
        '       Counting Status: {}\n'.format(
            countstatus))
    liner.send(
        '       Analyzed Reads : {:{},d}\n'.format(
            analyzedreads.value(),
            plen))
    liner.send(
        '           PhiX Reads : {:{},d} ({:6.2f} %)\n'.format(
            phiXreads.value(),
            plen,
            ut.safediv(
                A=100. * phiXreads.value(),
                B=analyzedreads.value())))
    liner.send(
        ' Low-Complexity Reads : {:{},d} ({:6.2f} %)\n'.format(
            lowcomplexreads.value(),
            plen,
            ut.safediv(
                A=100. * lowcomplexreads.value(),
                B=analyzedreads.value())))
    liner.send(
        ' Callback-False Reads : {:{},d} ({:6.2f} %)\n'.format(
            falsereads.value(),
            plen,
            ut.safediv(
                A=100. * falsereads.value(),
                B=analyzedreads.value())))
    liner.send(
        '   Incalculable Reads : {:{},d} ({:6.2f} %)\n'.format(
            incalcreads.value(),
            plen,
            ut.safediv(
                A=100. * incalcreads.value(),
                B=analyzedreads.value())))
    liner.send(
        '     Experiment Reads : {:{},d} ({:6.2f} %)\n'.format(
            experimentreads.value(),
            plen,
            ut.safediv(
                A=(100. * experimentreads.value()),
                B=analyzedreads.value())))
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - t0))

    # Remove Workspace
    ut.remove_directory(
        dirpath=countdir)

    # Unschedule countfile deletion
    if countstatus == 'Successful':
        ae.unregister(ctdeletion)

    # Close Liner
    liner.close()

    # Return Counts and Statistics
    stats = ut.stamp_stats(
        stats=stats,
        module='xcount',
        input_rows=input_rows,
        output_rows=len(outdf.index) if outdf is not None else 0)
    return (outdf, stats)
