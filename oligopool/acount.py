import time as tt

import collections as cx
import random      as rn
import atexit      as ae

import multiprocessing as mp

import utils     as ut
import valparse  as vp
import phiX      as px
import corecount as cc


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
    incalcreads,
    experimentreads,
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
    TBD
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
        'incalcreads'    : 0,
        'experimentreads': 0}

    # Compute Printing Lengths
    clen = ut.get_printlen(value=ncores)
    slen = plen = ut.get_printlen(
        value=packstat['survived'])

    # Read Count Storage
    countdict = cx.Counter()
    countpath = '{}/{}.{}.count'.format(
        countdir, coreid, batchid)

    # Initial Update
    if packqueue.empty():
        liner.send(
            ' Core {:{},d}: Shutting Down\n'.format(
                coreid,
                clen))
        shutdown.set()

    # Pack Counting Loop
    while not packqueue.empty():

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

            # Exoneration Block
            if not exoread is None:

                # Run Exoneration Procedure
                cc.exoneration_procedure(
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
            anchorstatus) = cc.get_anchored_read(
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
            index = cc.get_barcode_index(
                barcoderead=barcoderead,
                metamap=metamap,
                model=model)

            # Barcode Absent
            if index is None:
                cctrs['incalcreads'] += freq
                continue

            # Compute Associate Match
            associatematch, basalmatch = cc.get_associate_match(
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
                # Associate Mismatches Reference
                else:
                    cctrs['misassocreads'] += freq
                continue

            # All Components Valid
            countdict[(index,)] += freq
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

        # Need to Restart?
        if ut.needs_restart(
            memlimit=memlimit):
            restart.set() # Enable Restart
            break # Release, your Memory Real Estate, biatch!

    # Close packfile
    packfile.close()

    # Restart, We Must!
    if restart.is_set():
        # Show Updates
        liner.send(' Core {:{},d}: Restarting ...\n'.format(
            coreid,
            clen))
    # Shutdown!
    elif shutdown.is_set():
        # Show Updates
        liner.send(' Core {:{},d}: Shutting Down\n'.format(
            coreid,
            clen))

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
    incalcreads.increment(incr=cctrs['incalcreads'])
    experimentreads.increment(incr=cctrs['experimentreads'])

    # Counting Completed
    nactive.decrement()
    if shutdown.is_set():
        countqueue.put(None)

    # Release Control
    tt.sleep(0)

def acount(
    indexfile,
    packfile,
    countfile,
    maptype=0,
    barcodeerrors=-1,
    associateerrors=-1,
    ncores=0,
    memlimit=0,
    verbose=True):
    '''
    TBD
    '''

    # Start Liner
    liner = ut.liner_engine(online=verbose)

    # Packing Verbage Print
    liner.send('\n[Oligopool Calculator: Analysis Mode - Associate Count]\n')

    # Required Argument Parsing
    liner.send('\n Required Arguments\n')

    # Full indexfile Validation
    indexfile_valid = vp.get_indexfile_validity(
        indexfile=indexfile,
        indexfile_field='     Index File  ',
        associated=True,
        liner=liner)

    # Adjust indexfile Suffix
    if indexfile_valid:
        indexfile = ut.get_adjusted_path(
            path=indexfile,
            suffix='.oligopool.index')

    # Full packfile Validation
    (packfile_valid,
    packcount) = vp.get_parsed_packfile(
        packfile=packfile,
        packfile_field='      Pack File  ',
        liner=liner)

    # Adjust packfile Suffix
    if packfile_valid:
        packfile = ut.get_adjusted_path(
            path=packfile,
            suffix='.oligopool.pack')

    # Full countfile Validation
    countfile_valid = vp.get_outfile_validity(
        outfile=countfile,
        outfile_suffix='.oligopool.acount.csv',
        outfile_field='     Count File  ',
        liner=liner)

    # Adjust countfile Suffix
    if not countfile is None:
        countfile = ut.get_adjusted_path(
            path=countfile,
            suffix='.oligopool.acount.csv')

    # Optional Argument Parsing
    liner.send('\n Optional Arguments\n')

    # Full maptype Validation
    maptype_valid = vp.get_categorical_validity(
        category=maptype,
        category_field='   Mapping Type  ',
        category_pre_desc=' ',
        category_post_desc=' Classification',
        category_dict={
            0: 'Fast / Near-Exact',
            1: 'Slow / Sensitive'},
        liner=liner)

    # Full barcodeerrors Validation
    (barcodeerrors,
    barcodeerrors_valid) = vp.get_errors_validity(
        errors=barcodeerrors,
        errors_field='   Barcode Errors',
        errors_pre_desc=' At most ',
        errors_post_desc=' Mutations per Barcode',
        errors_base='B',
        indexfile_valid=indexfile_valid,
        indexfile=indexfile,
        liner=liner)

    # Full associateerrors Validation
    (associateerrors,
    associateerrors_valid) = vp.get_errors_validity(
        errors=associateerrors,
        errors_field=' Associate Errors',
        errors_pre_desc=' At most ',
        errors_post_desc=' Mutations per Associate',
        errors_base='A',
        indexfile_valid=indexfile_valid,
        indexfile=indexfile,
        liner=liner)

    # Full num_core Parsing and Validation
    (ncores,
    ncores_valid) = vp.get_parsed_core_info(
        ncores=ncores,
        core_field='       Num Cores ',
        default=packcount,
        offset=2,
        liner=liner)

    # Full num_core Parsing and Validation
    (memlimit,
    memlimit_valid) = vp.get_parsed_memory_info(
        memlimit=memlimit,
        memlimit_field='       Mem Limit ',
        ncores=ncores,
        ncores_valid=ncores_valid,
        liner=liner)

    # First Pass Validation
    if not all([
        indexfile_valid,
        packfile_valid,
        countfile_valid,
        maptype_valid,
        barcodeerrors_valid,
        associateerrors_valid,
        ncores_valid,
        memlimit_valid]):
        liner.send('\n')
        raise RuntimeError(
            'Invalid Argument Input(s).')

    # Start Timer
    t0 = tt.time()

    # Compute Assembly Parameters
    liner.send('\n[Step 1: Enqueing Read Packs]\n')

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
        outfile_suffix='.oligopool.acount.csv')

    # Schedule countfile deletion
    ctdeletion = ae.register(
        ut.remove_file,
        countfile)

    # Adjust Errors
    barcodeerrors   = round(barcodeerrors)
    associateerrors = round(associateerrors)

    # Define countqueue
    countqueue = mp.SimpleQueue()

    # Pack File Processing Book-keeping
    restarts  = [
        mp.Event() for _ in range(ncores)]
    shutdowns = [
        mp.Event() for _ in range(ncores)]

    # Read Counting Book-keeping
    nactive         = ut.SafeCounter(initval=ncores)
    analyzedreads   = ut.SafeCounter()
    phiXreads       = ut.SafeCounter()
    lowcomplexreads = ut.SafeCounter()
    misassocreads   = ut.SafeCounter()
    incalcreads     = ut.SafeCounter()
    experimentreads = ut.SafeCounter()
    batchids        = [0] * ncores
    previousreads   = [
        ut.SafeCounter() for _ in range(ncores)]

    # Wait on Enqueuing
    enqueuecomplete.wait()

    # Launching Read Packing
    liner.send('\n[Step 2: Computing Read Packs]\n')

    # Define Packing Stats
    stats = {
        'status'  : False,
        'basis'   : 'unsolved',
        'step'    : 2,
        'stepname': 'computing-read-packs',
        'vars'    : {
              'analyzed': int(analyzedreads.value()),
                  'phiX': int(phiXreads.value()),
            'lowcomplex': int(lowcomplexreads.value()),
              'misassoc': int(misassocreads.value()),
                'incalc': int(incalcreads.value()),
            'experiment': int(experimentreads.value())},
        'warns'   : {}}

    # Engine Timer
    et = tt.time()

    # Define Aggregator
    aggregator = mp.Process(
        target=cc.count_aggregator,
        args=(countqueue,
            countdir,
            ncores,
            nactive,
            liner,))

    # Fire-off Aggregator
    aggregator.start()

    # Define Counter Process Store
    readcounters = []

    # Fire-off Initial Read Counters
    coreid = 0
    clen = ut.get_printlen(value=ncores)
    while coreid < ncores:

        # Define Packer
        readcounter = mp.Process(
            target=acount_engine,
            args=(indexfile,
                packfile,
                packqueue,
                countdir,
                countqueue,
                maptype,
                barcodeerrors,
                associateerrors,
                previousreads[coreid],
                analyzedreads,
                phiXreads,
                lowcomplexreads,
                misassocreads,
                incalcreads,
                experimentreads,
                coreid,
                batchids[coreid],
                ncores,
                nactive,
                memlimit,
                restarts[coreid],
                shutdowns[coreid],
                et,
                liner,))

        # Show Update
        liner.send(
            ' Core {:{},d}: Starting Up\n'.format(
                coreid,
                clen))

        # Start Packer
        readcounter.start()

        # Update Book-keeping
        readcounters.append(readcounter)
        coreid += 1

    # Packer Management
    coreid = 0
    activepackers = ncores
    while activepackers:

        # Had Packer Finished?
        if readcounters[coreid] is None:
            pass

        # Has Packer Shutdown?
        elif shutdowns[coreid].is_set():
            # Cleanup
            readcounters[coreid].join()
            readcounters[coreid].close()
            # Update
            readcounters[coreid] = None
            activepackers -= 1
            ut.free_mem()
            # Reset
            restarts[coreid].clear()
            shutdowns[coreid].clear()

        # Must Packer Restart?
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
                target=acount_engine,
                args=(indexfile,
                    packfile,
                    packqueue,
                    countdir,
                    countqueue,
                    maptype,
                    barcodeerrors,
                    associateerrors,
                    previousreads[coreid],
                    analyzedreads,
                    phiXreads,
                    lowcomplexreads,
                    misassocreads,
                    incalcreads,
                    experimentreads,
                    coreid,
                    batchids[coreid],
                    ncores,
                    nactive,
                    memlimit,
                    restarts[coreid],
                    shutdowns[coreid],
                    et,
                    liner,))
            readcounters[coreid].start()

        # Next Iteration
        coreid = (coreid + 1) % ncores
        tt.sleep(0)

    # Join Enqueuer
    packenqueuer.join()

    # Join Aggregator
    aggregator.join()
    aggregator.close()

    # Free Memory
    ut.free_mem()

    # Compute Success Status
    if experimentreads.value() > 0:
        stats['status'] = True
        stats['basis']  = 'solved'

    # Show Time Elapsed
    liner.send(
        ' Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - et))

    # Did we succeed?
    if stats['status']:

        # Launching Read Packing
        liner.send('\n[Step 3: Writing Count Matrix]\n')

        # Aggregate Meta Read Packs
        cc.write_count(
            indexfiles=(indexfile,),
            countdir=countdir,
            countfile=countfile,
            liner=liner)

    # Update Stats
    stats['vars']['analyzed']   = int(analyzedreads.value())
    stats['vars']['phiX']       = int(phiXreads.value())
    stats['vars']['lowcomplex'] = int(lowcomplexreads.value())
    stats['vars']['misassoc']   = int(misassocreads.value())
    stats['vars']['incalc']     = int(incalcreads.value())
    stats['vars']['experiment'] = int(experimentreads.value())

    # Packing Status
    if stats['status']:
        countstatus = 'Successful'
    else:
        countstatus = 'Failed'

    # Read Counting Stats
    liner.send('\n[Associate Counting Stats]\n')

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
        ' Mis-Associated Reads : {:{},d} ({:6.2f} %)\n'.format(
            misassocreads.value(),
            plen,
            ut.safediv(
                A=100. * misassocreads.value(),
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

    # # Return Statistics
    return stats
