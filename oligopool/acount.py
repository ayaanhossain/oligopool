import os
import time as tt

import collections as cx
import zipfile     as zf
import random      as rn
import atexit      as ae

import multiprocessing as mp

import psutil as pu

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
    incalcreads,
    experimentreads,
    coreid,
    batchid,
    ncores,
    memlimit,
    launchtime,
    restart,
    shutdown,
    liner):
    '''
    TBD
    '''

    # Open indexfile
    indexfile = zf.ZipFile(
        file=indexfile)

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
    packfile = zf.ZipFile(
        file=packfile)

    # Load packing.stat
    packstat = ut.loaddict(
        archive=packfile,
        dfile='packing.stat')

    # Book-keeping Variables
    cctrs = {
        'analyzedreads'  : previousreads.value(),
        'phiXreads'      : 0,
        'lowcomplexreads': 0,
        'incalcreads'    : 0,
        'experimentreads': 0}

    # Compute Printing Lengths
    clen = len(str(ncores))
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

        # Fetch Pack Path / Token
        cpath = packqueue.get()

        # Exhaustion Token
        if cpath is None:
            shutdown.set()
            break

        # Fetch Read Pack
        cpack = ut.loadpack(
            archive=packfile,
            pfile=cpath)

        # Book-keeping Variables
        exoread = None
        exofreq = None

        verbagereach  = 0
        verbagetarget = rn.randint(
            *map(round, (len(cpack) * 0.080,
                         len(cpack) * 0.120)))

        readcount = 0
        packname  = cpath.rstrip('.pack')

        # Start Timer
        t0 = tt.time()

        print(len(cpack))

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
                        cctrs['analyzedreads'],
                        slen,
                        tt.time()-t0))

                # Update Book-keeping
                verbagereach = 0

            # Continue processing pack?
            if not cpack:
                break

            # Fetch Read and Frequency
            (read,
            freq) = cpack.pop()

            # read = ut.get_revcomp(read)

            # print('--------------')

            # print(read)

            # Update Book-keeping
            cctrs['analyzedreads'] += freq
            verbagereach += 1
            readcount    += 1

            # Anchor Read
            anchoredread, anchorstatus = cc.get_anchored_read(
                read=read,
                metamap=metamap,
                cctrs=cctrs)

            # Read References
            barcoderead   = anchoredread
            associateread = anchoredread

            # print('Barcode Part')
            # print(barcoderead)
            # print(anchor)

            # Anchor Absent
            if not anchorstatus:
                exoread = barcoderead
                exofreq = freq
                continue

            # Trim Barcode Prefix
            barcoderead, trimstatus = cc.get_trimmed_read(
                read=barcoderead,
                const=metamap['barcodeprefix'],
                constype=0,
                constval=metamap['bpxtval'])

            # print(barcoderead)
            # print(barcodeprefix)

            # Barcode Prefix Absent
            if not trimstatus:
                cctrs['incalcreads'] += freq
                continue

            # Trim Barcode Suffix
            barcoderead, trimstatus = cc.get_trimmed_read(
                read=barcoderead,
                const=metamap['barcodesuffix'],
                constype=1,
                constval=metamap['bsxtval'])

            # print(barcoderead)
            # print(barcodesuffix)

            # Barcode Suffix Absent
            if not trimstatus:
                cctrs['incalcreads'] += freq
                continue

            # Classify Barcode
            idx,score = model.predict(
                x=barcoderead)

            # print(idx,score)

            # Barcode Absent
            if idx is None:
                cctrs['incalcreads'] += freq
                continue

            # print('Associate Part')
            # print(associateread)

            # Trim Associate Prefix
            associateread, trimstatus = cc.get_trimmed_read(
                read=associateread,
                const=metamap['associateprefix'],
                constype=0,
                constval=metamap['apxtval'])

            # print(associateread)
            # print(associateprefix)

            # Associate Prefix Absent
            if not trimstatus:
                cctrs['incalcreads'] += freq
                continue

            # Trim Associate Suffix
            associateread, trimstatus = cc.get_trimmed_read(
                read=associateread,
                const=metamap['associatesuffix'],
                constype=1,
                constval=metamap['asxtval'])

            # print(associateread)
            # print(associatesuffix)

            # Associate Suffix Absent
            if not trimstatus:
                cctrs['incalcreads'] += freq
                continue

            # Match Associate
            associate, associatetval = associatedict[idx]
            if associateerrors < associatetval:
                associatetval = associateerrors
            associatematch = cc.get_associate_match(
                read=associateread,
                associate=associate,
                associatetval=associatetval)

            # print(associatematch)

            # Associate Absent
            if not associatematch:
                cctrs['incalcreads'] += freq
                continue

            # All Components Valid
            countdict[idx]    += freq
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

    # Save Count Dictionary
    ut.savecount(
        cobj=countdict,
        filepath=countpath)

    # countdict = ut.loadcount(countpath)

    print(len(countdict))

    # Release Control
    tt.sleep(0)

    # Update Read Counting Book-keeping
    previousreads.increment(incr=cctrs['analyzedreads'])
    analyzedreads.increment(incr=cctrs['analyzedreads'])
    phiXreads.increment(incr=cctrs['phiXreads'])
    lowcomplexreads.increment(incr=cctrs['lowcomplexreads'])
    incalcreads.increment(incr=cctrs['incalcreads'])
    experimentreads.increment(incr=cctrs['experimentreads'])

    # Counting Completed
    countqueue.put(countpath)
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
    barcodeerrors_valid = vp.get_errors_validity(
        errors=barcodeerrors,
        errors_field='   Barcode Errors',
        errors_pre_desc=' At most ',
        errors_post_desc=' Mutations per Barcode',
        errors_base='B',
        indexfile_valid=indexfile_valid,
        indexfile=indexfile,
        liner=liner)

    # Full associateerrors Validation
    associateerrors_valid = vp.get_errors_validity(
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
        offset=3,
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
        'vars'    : {},
        'warns'   : {}}

    # Engine Timer
    et = tt.time()

    # Execute Counting
    acount_engine(
        indexfile,
        packfile,
        packqueue,
        countdir,
        countqueue,
        maptype,
        barcodeerrors,
        associateerrors,
        previousreads[0],
        analyzedreads,
        phiXreads,
        lowcomplexreads,
        incalcreads,
        experimentreads,
        0,
        batchids[0],
        ncores,
        memlimit,
        et,
        restarts[0],
        shutdowns[0],
        liner)

    pot = 1
    while pot:
        item = countqueue.get()
        print(item, pot)
        if item == None:
            pot -= 1

    # Join Enqueuer
    packenqueuer.join()

    print(analyzedreads.value())
    print(phiXreads.value())
    print(lowcomplexreads.value())
    print(incalcreads.value())
    print(experimentreads.value())

    # Remove Workspace
    ut.remove_directory(
        dirpath=countdir)

    # # Unschedule packfile deletion
    # if countstatus == 'Successful':
    #     ae.unregister(ctdeletion)

    # Close Liner
    liner.close()

    # # Return Statistics
    # return stats

    return

    # Launching Read Counting
    liner.send('\n[Counting Read Packs]\n')

    # Start Timer
    t0 = tt.time()

    # Counter Process Store
    readcounters = []

    # Fire off Single-Core Read Counter
    if ncores == 1:
        count_engine(
            prepfile=prepfile,
            packfile=packfile,
            packqueue=packqueue,
            countqueue=countqueue,
            amperrors=amperrors,
            barerrors=barerrors,
            coreid=0,
            ncores=ncores,
            analyzedreads=analyzedreads,
            diagnosticreads=diagnosticreads,
            phiXreads=phiXreads,
            mononukereads=mononukereads,
            dinukereads=dinukereads,
            trinukereads=trinukereads,
            incalcreads=incalcreads,
            liner=liner)

    # Fire off Multi-Core Read Counters
    else:

        # Queue All Read Counters
        coreid = 0
        while coreid < ncores:

            # Define Counter
            readcounter = mp.Process(
                target=count_engine,
                args=(prepfile,
                    packfile,
                    packqueue,
                    countqueue,
                    amperrors,
                    barerrors,
                    coreid,
                    ncores,
                    analyzedreads,
                    diagnosticreads,
                    phiXreads,
                    mononukereads,
                    dinukereads,
                    trinukereads,
                    incalcreads,
                    liner,))

            # Start Counter
            readcounter.start()

            # Update Book-keeping
            readcounters.append(readcounter)
            coreid += 1

    # Aggregate All Count Matrices
    agdf = count_aggregator(
        prepfile=prepfile,
        countqueue=countqueue,
        countfile=countfile,
        prodcount=ncores)

    # Join All Read Counters
    for readcounter in readcounters:
        readcounter.join()

    # Counting Status
    if diagnosticreads.value() > 0:
        countstatus = (True,  'Successful')
    else:
        countstatus = (False, 'Failed')

    # Read Counting Stats
    liner.send('\n[Counting Stats]\n')

    plen = ut.get_printlen(
        value=analyzedreads.value())

    liner.send(
        '       Counting  Status: {}\n'.format(
            countstatus[1]))
    liner.send(
        '       Analyzed   Reads: {:{},d}\n'.format(
            analyzedreads.value(),
            plen))
    liner.send(
        '           PhiX   Reads: {:{},d} ({:6.2f} %)\n'.format(
            phiXreads.value(),
            plen,
            ut.safediv(
                A=100. * phiXreads.value(),
                B=analyzedreads.value())))
    liner.send(
        ' Mononucleotide   Reads: {:{},d} ({:6.2f} %)\n'.format(
            mononukereads.value(),
            plen,
            ut.safediv(
                A=100. * mononukereads.value(),
                B=analyzedreads.value())))
    liner.send(
        '   Dinucleotide   Reads: {:{},d} ({:6.2f} %)\n'.format(
            dinukereads.value(),
            plen,
            ut.safediv(
                A=100. * dinukereads.value(),
                B=analyzedreads.value())))
    liner.send(
        '  Trinucleotide   Reads: {:{},d} ({:6.2f} %)\n'.format(
            trinukereads.value(),
            plen,
            ut.safediv(
                A=100. * trinukereads.value(),
                B=analyzedreads.value())))
    liner.send(
        '   Incalculable   Reads: {:{},d} ({:6.2f} %)\n'.format(
            incalcreads.value(),
            plen,
            ut.safediv(
                A=100. * incalcreads.value(),
                B=analyzedreads.value())))
    liner.send(
        '     Diagnostic   Reads: {:{},d} ({:6.2f} %)\n'.format(
            diagnosticreads.value(),
            plen,
            ut.safediv(
                A=(100. * diagnosticreads.value()),
                B=analyzedreads.value())))
    liner.send(
        '           Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - t0))

    # Unschedule countfile deletion
    if countstatus[0]:
        ae.unregister(ctdeletion)

    # Close Liner
    liner.close()

    # Return Aggregate DataFrame
    return agdf